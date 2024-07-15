import gzip
import os
import re
from threading import Thread
from pathlib import Path
from typing import Any, IO, NamedTuple
from common.config import StratoMod, RefsetKey, VCFFile
from common.io import bgzip_file, spawn_stream, check_processes

NORMCMD = ["bcftools", "norm", "--multiallelics", "-"]


class Env(NamedTuple):
    refsetkey: RefsetKey
    vcf: VCFFile
    config: StratoMod


def fix_DV_refcall(filter_col: bytes, sample_col: bytes) -> bytes:
    return (
        sample_col.replace(b"./.", b"0/1").replace(b"0/0", b"0/1")
        if filter_col == b"RefCall"
        else sample_col
    )


def strip_format_fields(
    fields: set[str],
    format_col: bytes,
    sample_col: bytes,
) -> tuple[bytes, bytes]:
    f, s = zip(
        *[
            (f, s)
            for f, s in zip(format_col.split(b":"), sample_col.split(b":"))
            if f.decode() not in fields
        ]
    )
    return (b":".join(f), b":".join(s))


def filter_file(env: Env, fi: IO[bytes], fo: IO[bytes], log: IO[str]) -> None:
    chr_prefix = env.vcf.chr_prefix
    cs = env.config.refsetkey_to_chr_indices(env.refsetkey)
    chr_mapper = {c.chr_name_full(chr_prefix): c.value for c in cs}

    for ln in fi:
        # for commented lines, either print it as-is or do this weird dance if
        # it has a contig parameter (these must be changed to reflect the
        # integer chromosome numbers we are about to assign as some tools will
        # complain if there is a mismatch)
        if ln.startswith(b"#"):
            # If a contig line, change the ID field (and hopefully leave
            # everything else as-is).
            if cmatch := re.match("##contig=<(.*)>", ln.decode()):
                contig_matches = [
                    re.match("(.*)=(.*)", c) for c in cmatch[1].split(",")
                ]
                if any([x for x in contig_matches if x is None]):
                    log.write("Invalid contig line: %s" % ln.decode())
                    exit(1)
                contig_map: dict[str, str] = {
                    x[1]: x[2] for x in contig_matches if x is not None
                }
                if "ID" not in contig_map:
                    log.write("ID field not in contig line: %s" % ln.decode())
                    exit(1)
                try:
                    contig_map["ID"] = chr_mapper[contig_map["ID"]]
                    new_contig = ",".join([f"{k}={v}" for k, v in contig_map.items()])
                    fo.write((f"##contig=<{new_contig}>\n").encode())
                except KeyError:
                    pass
                except Exception as e:
                    log.write(str(e))
            else:
                fo.write(ln)
        else:
            ls = ln.rstrip().split(b"\t")[:10]
            # CHROM = 0
            # POS = 1
            # ID = 2
            # REF = 3
            # ALT = 4
            # QUAL = 5
            # FILTER = 6
            # INFO = 7
            # FORMAT = 8
            # SAMPLE = 9
            try:
                ls[0] = str(chr_mapper[ls[0].decode()]).encode()
                # optionally fix the GT field for filtered variants
                if env.vcf.corrections.fix_refcall_gt:
                    ls[9] = fix_DV_refcall(ls[6], ls[9])
                # optionally remove some fields in FORMAT/SAMPLE
                if len(env.vcf.corrections.strip_format_fields) > 0:
                    ls[8], ls[9] = strip_format_fields(
                        env.vcf.corrections.strip_format_fields,
                        ls[8],
                        ls[9],
                    )
                fo.write(b"\t".join(ls) + b"\n")
            except KeyError:
                pass


def main(smk: Any, config: StratoMod) -> None:
    env = Env(
        refsetkey=RefsetKey(smk.wildcards["refset_key"]),
        vcf=smk.params.vcf,
        config=config,
    )
    out = Path(smk.output[0])

    with gzip.open(smk.input[0], "rb") as i, open(smk.log["filtered"], "w") as lo:
        r, w = os.pipe()
        r0 = os.fdopen(r, "rb")
        w0 = os.fdopen(w, "wb")
        try:

            def go() -> None:
                try:
                    filter_file(env, i, w0, lo)
                except Exception:
                    r0.close()
                finally:
                    w0.close()

            t1 = Thread(target=go)
            t1.start()

            if env.vcf.split_biallelics:
                p1, o1 = spawn_stream(NORMCMD, r0)  # NORM!!!!!!
                bgzip_file(o1, out)
                check_processes([p1], Path(smk.log["split"]))
            else:
                bgzip_file(r0, out)
        finally:
            r0.close()
            w0.close()


main(snakemake, snakemake.config)  # type: ignore
