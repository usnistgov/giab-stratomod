import gzip
import os
from tempfile import mkdtemp
import subprocess as sp
import hashlib
from typing import Callable, TextIO, TypeVar, Generator, IO
from pathlib import Path
from logging import Logger
from contextlib import contextmanager

X = TypeVar("X")


@contextmanager
def temp_fifo() -> Generator[Path, None, None]:
    """Context Manager for creating named pipes with temporary names."""
    tmpdir = mkdtemp()
    filename = Path(tmpdir) / "fifo"
    os.mkfifo(filename)
    try:
        yield filename
    finally:
        os.unlink(filename)
        os.rmdir(tmpdir)


def with_gzip_maybe(f: Callable[[TextIO, TextIO], X], i: str, o: str) -> X:
    hi = (
        gzip.open(i, "rt", encoding="latin1")
        if i.endswith(".gz")
        else open(i, "rt", encoding="latin1")
    )
    ho = gzip.open(o, "wt") if o.endswith(".gz") else open(o, "wt")
    with hi as fi, ho as fo:
        return f(fi, fo)


def get_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


def get_md5_dir(path: Path) -> str:
    """Get the "MD5" of a directory.

    This is equivalent to the following bash command:

    find <dir> -type f -exec md5sum {} \\; | \\
        grep -v snakemake_timestamp | \\
        LC_COLLATE=C sort -k 2 | \\
        cut -d' ' -f1 | \\
        tr -d '\n' | \\
        md5sum

    Will only get the md5 of the first layer of files (eg not recursive). This
    is appropriate for SDF files which have this property.
    """
    h = hashlib.md5()
    ps = sorted(
        [
            p
            for p in path.iterdir()
            if p.is_file() and "snakemake_timestamp" not in p.name
        ]
    )
    for p in ps:
        h.update(get_md5(p).encode())
    return h.hexdigest()


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with gzip.open(p, "r") as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False


# set up basic logger that prints to both console and a file (the log directive
# from snakemake) and captures warnings so those don't go unnoticed
def setup_logging(path: str, console: bool = False) -> Logger:
    import logging

    logging.basicConfig(filename=path, level=logging.INFO)
    logging.captureWarnings(True)
    logger = logging.getLogger()
    if console:
        logger.addHandler(logging.StreamHandler())
    return logger


def spawn_stream(
    cmd: list[str],
    i: IO[bytes] | int | None = None,
) -> tuple[sp.Popen[bytes], IO[bytes]]:
    p = sp.Popen(cmd, stdin=i, stdout=sp.PIPE, stderr=sp.PIPE)
    # ASSUME since we typed the inputs so that the stdin/stdout can only take
    # file descriptors or file streams, the return for each will never be
    # none
    if p.stdout is None:
        raise
    return p, p.stdout


def check_processes(
    ps: list[sp.Popen[bytes] | sp.CompletedProcess[bytes]], log: Path
) -> None:
    some_error = False
    # TODO make parent directory if it doesn't exist? probably won't be necessary
    with open(log, "w") as lf:
        for p in ps:
            if isinstance(p, sp.CompletedProcess):
                err = p.stderr
            else:
                _, err = p.communicate()  # break deadlocks if there are any
            if p.returncode != 0:
                some_error = True

                args = p.args
                if isinstance(args, list):
                    cmd = " ".join(args)
                elif isinstance(args, bytes):
                    cmd = args.decode()
                else:
                    cmd = str(args)
                if err:
                    lf.write(f"{cmd}: {err.decode()}\n")
                else:
                    lf.write(f"{cmd}: return code {p.returncode}\n")

    if some_error:
        exit(1)
