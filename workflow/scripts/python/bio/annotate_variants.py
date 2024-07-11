import gzip
from functools import reduce
from pathlib import Path
from typing import Any, Generator, NamedTuple, IO
from common.io import setup_logging
from common.functional import maybe
import common.config as cfg

logger = setup_logging(snakemake.log[0])  # type: ignore


Stream = Generator[list[str], None, None]


class Line(NamedTuple):
    chrom: int
    start: int
    end: int
    rest: list[str]


def line_to_strs(x: Line) -> list[str]:
    return [str(x.chrom), str(x.start), str(x.end)] + x.rest


def strs_to_line(xs: list[str]) -> Line:
    return Line(int(xs[0]), int(xs[1]), int(xs[2]), xs[3:])


def split_line(s: str) -> list[str]:
    return s.rstrip().split("\t")


def read_variants(p: Path) -> Stream:
    with gzip.open(p, "rt") as f:
        for x in f:
            yield split_line(x)


def left_outer_intersect(left: Stream, p: Path) -> Stream:
    """Perform a left-outer join of 'left' with contents in 'p'.

    Assume that 'left' and 'p' both represent bed-like files with
    bed coordinates in the first three columns, and also have headers.
    Also assume both are sorted.

    Join left and right if both have intervals that overlap. This is
    roughly equivalent to "intersectBed -a 'left' -b 'p' -loj" (except
    it can't deal with headers, hence this function)
    """
    rline_buffer: list[Line] = []

    def left_line(s: Stream) -> Line | None:
        return maybe(None, strs_to_line, next(s, None))

    def right_line(s: IO[str]) -> Line | None:
        try:
            return rline_buffer.pop(0)
        except IndexError:
            return maybe(
                None, lambda x: strs_to_line(x.rstrip().split("\t")), next(s, None)
            )

    with gzip.open(p, "rt") as f:
        # TODO this will cryptically fail if the header happens to be blank
        lheader = next(left)
        rheader = next(f).rstrip().split("\t")[3:]
        yield lheader + rheader

        blank = [""] * len(rheader)
        lline = left_line(left)
        rline = right_line(f)
        while lline is not None and rline is not None:
            # If right has greater chromosome, there is nothing in the right
            # side to join so return left line with blanks
            if lline.chrom < rline.chrom:
                yield line_to_strs(lline) + blank
                lline = left_line(left)

            # Left chrom >= right chrom
            else:
                # Throw away all the right side until chroms are equal and the
                # end of the right interval is greater than the left start.
                while rline is not None and not (
                    lline.chrom == rline.chrom and lline.start < rline.end
                ):
                    rline = right_line(f)

                # Bail if we happen to exhaust the right side completely
                if rline is None:
                    break

                # If the left end is greater than the right start, then they
                # overlap. In this case return the left with the right appended.
                if lline.end > rline.start:
                    first_rline = rline
                    tmp: list[Line] = []
                    yield line_to_strs(lline) + rline.rest
                    # Continue to loop through right side for this left interval
                    # since multiple on the right may overlap.
                    while (rline := right_line(f)) is not None:
                        tmp.append(rline)
                        if lline.end > rline.start:
                            yield line_to_strs(lline) + rline.rest
                        else:
                            break
                    # Since the next left interval might start at the exact same
                    # position as the current left interval, save all the right
                    # intervals through which we iterated here by adding them
                    # to a temp buffer which will be read before the next line
                    # in the file.
                    rline_buffer = tmp + rline_buffer
                    # Start at the exact same right interval next time
                    rline = first_rline
                    lline = left_line(left)
                # Otherwise, the right interval is beyond the left,
                # in which case return left with blanks and get a new left
                # interval.
                else:
                    yield line_to_strs(lline) + blank
                    lline = left_line(left)

    # if we run out of lines on the right, loop through the rest of the left
    # and return blanks
    while lline is not None:
        yield line_to_strs(lline) + blank
        lline = left_line(left)


def write_annotations(xs: Stream, p: Path) -> None:
    with gzip.open(p, "wt") as f:
        header = next(xs)
        f.write("\t".join([cfg.VAR_IDX, *header]) + "\n")
        for i, x in enumerate(xs):
            f.write("\t".join([str(i), *x]) + "\n")


def main(smk: Any, config: cfg.StratoMod) -> None:
    fs = [Path(p) for p in smk.input.features]
    vcf = Path(smk.input.variants[0])
    logger.info("Adding annotations to %s\n", vcf)
    write_annotations(
        reduce(left_outer_intersect, fs, read_variants(vcf)),
        Path(smk.output[0]),
    )


main(snakemake, snakemake.config)  # type: ignore
