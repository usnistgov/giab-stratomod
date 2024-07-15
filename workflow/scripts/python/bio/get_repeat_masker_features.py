import os
import gzip
from threading import Thread
from pathlib import Path
from typing import Any
import common.config as cfg
from common.io import spawn_stream, check_processes
from common.functional import maybe

# The repeat masker database is documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema


def main(smk: Any, config: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    rk = config.refsetkey_to_refkey(rsk)
    src = config.references[rk].feature_data.repeat_masker

    def merge_and_write_group(clsname: str, famname: str | None = None) -> None:
        col = config.feature_definitions.repeat_masker.fmt_name(src, clsname, famname)
        c = clsname.encode()
        f = maybe(None, lambda s: s.encode(), famname)
        i = smk.input[0]
        o = smk.output[0]
        with gzip.open(i, "r") as ih, gzip.open(o, "w") as oh:
            # write header
            oh.write(("\t".join([*cfg.BED_COLS, col]) + "\n").encode())

            try:
                r, w = os.pipe()
                r0 = os.fdopen(r, "rb")
                w0 = os.fdopen(w, "wb")

                # 1. Filter stream for class and/or family (if present). Assume
                # stream has columns: chrom, start, end, class, family. Only
                # write the first 3 columns.
                def filter_base() -> None:
                    try:
                        for x in ih:
                            s = x.rstrip().split(b"\t")
                            if c == s[3] and (f is None or f == s[4]):
                                w0.write(b"\t".join(s[:3]) + b"\n")
                    except Exception:
                        r0.close()
                    finally:
                        w0.close()

                t1 = Thread(target=filter_base)
                t1.start()

                p1, o1 = spawn_stream(["mergeBed", "-i", "stdin"], r0)

                for x in o1:
                    s = x.rstrip().split(b"\t")
                    start = int(s[1])
                    end = int(s[2])
                    newline = [*s, str(end - start).encode()]
                    oh.write(b"\t".join(newline) + b"\n")

                check_processes([p1], Path(smk.log[0]))
            finally:
                r0.close()
                w0.close()

    cls = smk.wildcards.rmsk_class

    try:
        fam = smk.wildcards.rmsk_family
        merge_and_write_group(cls, fam)
    except AttributeError:
        merge_and_write_group(cls)


main(snakemake, snakemake.config)  # type: ignore
