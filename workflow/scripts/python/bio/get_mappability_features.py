import gzip
from pathlib import Path
import pandas as pd
import common.config as cfg
from typing import Any
from common.tsv import write_tsv
from common.bed import read_bed, bed_to_stream
from common.io import spawn_stream, check_processes


def main(smk: Any, config: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    cs = config.refsetkey_to_chr_indices(rsk)
    mapconf = config.refsetkey_to_ref(rsk).feature_data.mappability
    mapmeta = config.feature_definitions.mappability

    def read_map_bed(p: Path, ps: cfg.BedFileParams, col: str) -> pd.DataFrame:
        df = read_bed(p, ps, {}, cs)
        df[col] = 1
        return df

    # first, read, sort/filter, and write the high mappability file
    high = read_map_bed(smk.input["high"][0], mapconf.high.params, mapmeta.high)
    highout = smk.output["high"]
    write_tsv(highout, high)

    # once high is fully processed and on-disk, read/sort/filter the low file,
    # then subtract off the high mappability file (since the latter is a subset
    # of the former) using a stream
    cmd = ["subtractBed", "-a", "stdin", "-b", highout, "-g", smk.input["genome"][0]]
    with gzip.open(smk.output["low"], "wt") as lowout:
        low = read_map_bed(smk.input["low"][0], mapconf.low.params, mapmeta.low)

        # write header first
        lowout.write("\t".join(low.columns) + "\n")

        with bed_to_stream(low) as s:
            p, o = spawn_stream(cmd, s)
            for x in o:
                lowout.write(x.decode())
            check_processes([p], Path(smk.log[0]))


main(snakemake, snakemake.config)  # type: ignore
