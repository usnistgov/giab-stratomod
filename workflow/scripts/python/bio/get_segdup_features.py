import gzip
import pandas as pd
from pathlib import Path
from typing import Any
import common.config as cfg
from common.bed import read_bed, merge_and_apply_stats, bed_to_stream
from common.io import check_processes

# This database is documented here:
# http://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg38&hgta_doSchemaTable=genomicSuperDups


def read_segdups(
    smk: Any,
    config: cfg.StratoMod,
    path: Path,
    fconf: cfg.SegDupsGroup,
) -> pd.DataFrame:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    rk = config.refsetkey_to_refkey(rsk)
    s = config.references[rk].feature_data.segdups
    ocs = s.other_cols
    feature_cols = {
        ocs.align_L: str(fconf.fmt_col(lambda x: x.alignL)[0]),
        ocs.frac_match_indel: str(fconf.fmt_col(lambda x: x.fracMatchIndel)[0]),
    }
    cs = config.refsetkey_to_chr_indices(rsk)
    return read_bed(path, s.params, feature_cols, cs)


def main(smk: Any, config: cfg.StratoMod) -> None:
    fconf = config.feature_definitions.segdups
    df = read_segdups(smk, config, smk.input[0], fconf)
    with bed_to_stream(df) as s:
        p, o, header = merge_and_apply_stats(fconf, s, df.columns.tolist()[3:])
        with gzip.open(smk.output[0], "wt") as oh:
            oh.write("\t".join(header) + "\n")
            for x in o:
                oh.write(x.decode())
        check_processes([p], smk.log[0])


main(snakemake, snakemake.config)  # type: ignore
