import contextlib
import os
import subprocess as sp
from threading import Thread
from pathlib import Path
from typing import TypeVar, IO, Generator
from itertools import product
import pandas as pd
import common.config as cfg
from common.io import spawn_stream
from common.functional import DesignError


X = TypeVar("X")


def read_bed(
    path: Path,
    b: cfg.BedFileParams = cfg.BedFileParams(),
    more: dict[int, str] = {},
    chr_indices: list[cfg.ChrIndex] = [],
) -> pd.DataFrame:
    """Read a bed file as a pandas dataframe.

    Return a dataframe where the first three columns are numbered 0, 1, 2 and
    typed str, int, int (first is str regardless of how the chr names are
    formated). Columns from 'more' are appended to the end of the dataframe
    in the order given starting from 3.
    """
    bedcols = {**b.bed_cols.indexed, **more}
    df = pd.read_table(
        path,
        header=None,
        usecols=[*bedcols],
        sep=b.sep,
        comment="#",
        skiprows=b.skip_lines,
        # satisfy type checker :/
        dtype={k: v for k, v in b.bed_cols.typed.items()},
    )
    df.columns = pd.Index(bedcols.values())
    if len(chr_indices) > 0:
        f = cfg.ChrFilter(b.chr_prefix, chr_indices)
        return filter_sort_bed(f, df)
    else:
        return df


def filter_sort_bed(cfilt: cfg.ChrFilter, df: pd.DataFrame, n: int = 3) -> pd.DataFrame:
    from_map = {i.chr_name_full(cfilt.prefix): i.value for i in cfilt.indices}
    return filter_sort_bed_inner(from_map, df, n)


def filter_sort_bed_inner(
    from_map: dict[str, int],
    df: pd.DataFrame,
    n: int = 3,
) -> pd.DataFrame:
    chr_col = df.columns.tolist()[0]
    df[chr_col] = df[chr_col].map(from_map)
    df = df.dropna(subset=[chr_col]).astype({chr_col: int})
    return sort_bed_numerically(df, n)


def sort_bed_numerically(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Sort a bed file encoded by a dataframe.

    Assumes the first three columns correspond to coordinates, and that all are
    integer typed. Use 'n = 2' to sort only by chr/start, and 'n=1' to sort only
    by chr.

    """
    cols = df.columns.tolist()
    bycols = [cols[i] for i in range(0, n)]
    return df.sort_values(by=bycols, axis=0, ignore_index=True)


def merge_and_apply_stats(
    fconf: cfg.MergedFeatureGroup[X],
    bed_stream: IO[bytes],
    stat_cols: list[str],
) -> tuple[sp.Popen[bytes], IO[bytes], list[str]]:
    res = [
        (i + 4, m.value, fconf.fmt_merged_feature(s, m))
        for (i, s), m in product(enumerate(stat_cols), fconf.operations)
    ]
    cols = [r[0] for r in res]
    opts = [r[1] for r in res]
    headers = [r[2] for r in res]

    # just use one column for count since all columns will produce the same
    # number
    full_opts = ",".join(["count", *opts])
    # 4 since we only need to count the first of the extra columns
    full_cols = ",".join(map(str, [4, *cols]))
    full_headers: list[str] = [*cfg.BED_COLS, fconf.count_feature[0], *headers]

    cmd = ["mergeBed", "-i", "stdin", "-c", full_cols, "-o", full_opts]
    p, o = spawn_stream(cmd, bed_stream)
    return p, o, full_headers


@contextlib.contextmanager
def bed_to_stream(df: pd.DataFrame) -> Generator[IO[bytes], None, None]:
    """Transform bed-like dataframe into a bytestream.

    This is useful for using a pandas dataframe in a bedtools subprocess.
    """
    r, w = os.pipe()
    _r = os.fdopen(r, "rb")
    _w = os.fdopen(w, "w")

    # NOTE: python threads can't pass exceptions between themselves and the
    # calling thread, so need to hack something to signal when something bad
    # happens. Do this by closing read side of the thread (which will also make
    # any process that depends on the pipe die since it can't be opened after
    # we close it), and then testing if the read end is closed after the context
    # block.
    def read_df() -> None:
        try:
            df.to_csv(_w, sep="\t", header=False, index=False)
            # write_bed_stream(_w, df)
        except Exception:
            _r.close()
        finally:
            _w.close()

    t = Thread(target=read_df)
    t.start()

    try:
        yield _r
        if _r.closed:
            raise DesignError("bed stream thread exited unexpectedly")
    finally:
        _r.close()
        _w.close()
