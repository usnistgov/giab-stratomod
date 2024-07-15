import gzip
from pathlib import Path
import pandas as pd
from typing import Any
import common.config as cfg
from common.bed import read_bed, merge_and_apply_stats, bed_to_stream
from common.io import check_processes, spawn_stream

# Input dataframe documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema

SLOP = 5


def read_tandem_repeats(
    smk: Any,
    path: Path,
    fconf: cfg.TandemRepeatGroup,
    sconf: cfg.StratoMod,
) -> pd.DataFrame:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    rk = sconf.refsetkey_to_refkey(rsk)
    ss = sconf.references[rk].feature_data.tandem_repeats
    ocs = ss.other_cols
    fmt_col = fconf.fmt_col
    perc_a_col = str(fconf.A[0])
    perc_t_col = str(fconf.T[0])
    perc_c_col = str(fconf.C[0])
    perc_g_col = str(fconf.G[0])
    unit_size_col = fmt_col(lambda x: x.period)[0]
    feature_cols = {
        ocs.period: unit_size_col,
        ocs.copy_num: fmt_col(lambda x: x.copyNum)[0],
        ocs.per_match: fmt_col(lambda x: x.perMatch)[0],
        ocs.per_indel: fmt_col(lambda x: x.perIndel)[0],
        ocs.score: fmt_col(lambda x: x.score)[0],
        ocs.per_A: perc_a_col,
        ocs.per_C: perc_c_col,
        ocs.per_G: perc_g_col,
        ocs.per_T: perc_t_col,
    }
    cs = sconf.refsetkey_to_chr_indices(rsk)
    df = read_bed(path, ss.params, feature_cols, cs)
    base_groups = [
        (fconf.AT[0], perc_a_col, perc_t_col),
        (fconf.AG[0], perc_a_col, perc_g_col),
        (fconf.CT[0], perc_c_col, perc_t_col),
        (fconf.GC[0], perc_c_col, perc_g_col),
    ]
    for double, single1, single2 in base_groups:
        df[double] = df[single1] + df[single2]
    # Filter out all TRs that have period == 1, since those by definition are
    # homopolymers. NOTE, there is a difference between period and consensusSize
    # in this database; however, it turns out that at least for GRCh38 that the
    # sets of TRs where either == 1 are identical, so just use period here
    # since I can easily refer to it.
    # logger.info("Removing TRs with unitsize == 1")
    return df[df[unit_size_col] > 1]


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    i = smk.input
    fconf = sconf.feature_definitions.tandem_repeats
    len_col = fconf.length[0]
    slopcmd = ["slopBed", "-i", "stdin", "-b", str(SLOP), "-g", i.genome[0]]
    df = read_tandem_repeats(smk, Path(i.src[0]), fconf, sconf)
    with bed_to_stream(df) as s:
        p0, o0, header = merge_and_apply_stats(fconf, s, df.columns.tolist()[3:])
        p1, o1 = spawn_stream(slopcmd, o0)
        with gzip.open(smk.output[0], "wt") as oh:
            oh.write("\t".join([*header, len_col]) + "\n")
            for x in o1:
                y = x.decode().rstrip().split("\t")
                start = int(y[1])
                end = int(y[2])
                newline = [*y, str(end - start - SLOP * 2)]
                oh.write("\t".join(newline) + "\n")
        check_processes([p0, p1], smk.log[0])


main(snakemake, snakemake.config)  # type: ignore
