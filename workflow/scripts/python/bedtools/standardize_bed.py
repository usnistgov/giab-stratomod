import gzip
import io
import re
from typing import Callable, Any
from common.config import StratoMod, ChrIndex, RefsetKey

# I could use pandas for all this, but vcfeval will complain if I strip out the
# headers (which pandas will do). Imperative loop it is...


def filter_file(smk: Any, config: StratoMod, fi: io.TextIOWrapper) -> None:
    chr_prefix = smk.params.chr_prefix
    chr_indices = config.refsetkey_to_chr_indices(
        RefsetKey(smk.wildcards["refset_key"])
    )
    fs = tuple(["#", *[f"{i.chr_name}\t" for i in chr_indices]])

    def make_sub(i: ChrIndex) -> Callable[[str], str]:
        pat = re.compile(f"^{i.chr_name}")
        return lambda s: pat.sub(str(i.value), s)

    subX = make_sub(ChrIndex.CHRX)
    subY = make_sub(ChrIndex.CHRY)

    def tolines(f: io.TextIOWrapper) -> None:
        f.writelines(
            (
                subY(subX(y))
                for x in fi
                if (y := x.removeprefix(chr_prefix)).startswith(fs)
            ),
        )

    if str(smk.output[0]).endswith(".gz"):
        with gzip.open(smk.output[0], "wt") as fo:
            tolines(fo)
    else:
        with open(smk.output[0], "wt") as fo:
            tolines(fo)


def main(smk: Any, config: StratoMod) -> None:
    if str(smk.input[0]).endswith(".gz"):
        with gzip.open(smk.input[0], "rt") as f:
            filter_file(smk, config, f)
    else:
        with open(smk.input[0], "rt") as f:
            filter_file(smk, config, f)


main(snakemake, snakemake.config)  # type: ignore
