import gzip
import common.config as cfg


def main(opath: str, regions: list[cfg.BedRegion]) -> None:
    with gzip.open(opath, "wt") as f:
        for r in sorted(regions):
            f.write(r.fmt() + "\n")


main(snakemake.output[0], snakemake.params.regions)  # type: ignore
