import re
import os
from threading import Thread
from pathlib import Path
from typing import Any, IO
import subprocess as sp
import common.config as cfg
from common.io import spawn_stream, check_processes, bgzip_file
from common.functional import DesignError


def stream_fasta(ipath: str, chr_names: list[str]) -> tuple[sp.Popen[bytes], IO[bytes]]:
    return spawn_stream(["samtools", "faidx", ipath, *chr_names])


def stream_sdf(ipath: str, chr_names: list[str]) -> tuple[sp.Popen[bytes], IO[bytes]]:
    cmd = [
        *["rtg", "sdf2fasta", "--no-gzip", "--line-length=70"],
        *["--input", ipath],
        *["--output", "-"],
        *["--names", *chr_names],
    ]
    return spawn_stream(cmd)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    cs = sconf.refsetkey_to_chr_indices(rsk)
    prefix = sconf.refsetkey_to_ref(rsk).sdf.chr_prefix

    chr_mapper = {c.chr_name_full(prefix): c.value for c in cs}
    chr_names = [*chr_mapper]

    # Read from a fasta or sdf depending on what we were given; in either
    # case, read only the chromosomes we want in sorted order and return a
    # fasta text stream
    def choose_input(i: Any) -> tuple[sp.Popen[bytes], IO[bytes]]:
        try:
            return stream_fasta(i.fasta[0], chr_names)
        except AttributeError:
            try:
                return stream_sdf(i.sdf[0], chr_names)
            except AttributeError:
                assert False, "unknown input key, this should not happen"

    with open(smk.log["index"], "w") as lo:
        p, o = choose_input(smk.input)

        r, w = os.pipe()
        r1 = os.fdopen(r, "rb")
        w1 = os.fdopen(w, "wb")

        def go() -> None:
            # Stream the fasta and replace the chromosome names in the header with
            # its integer index
            try:
                for i in o:
                    if i.startswith(b">"):
                        m = re.match(">([^ \n]+)", i.decode())
                        if m is None:
                            lo.write("could get chrom name from FASTA header")
                            exit(1)
                        try:
                            w1.write(f">{chr_mapper[m[1]]}\n".encode())
                        except KeyError:
                            raise DesignError(f"could not convert '{m[1]}' to index")
                    else:
                        w1.write(i)
            except Exception:
                r1.close()
            finally:
                w1.close()

        try:
            t = Thread(target=go)
            t.start()
            bgzip_file(r1, Path(smk.output[0]))
        finally:
            r1.close()
            w1.close()

        check_processes([p], Path(smk.log["convert"]))


main(snakemake, snakemake.config)  # type: ignore
