from typing import Any
import common.config as cfg
from common.bed import read_bed


def main(smk: Any, config: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    rk = config.refsetkey_to_refkey(rsk)
    src = config.references[rk].feature_data.repeat_masker
    cs = config.refsetkey_to_chr_indices(rsk)

    cols = {11: "_class", 12: "_family"}
    df = read_bed(smk.input[0], src.params, cols, cs)
    df.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
