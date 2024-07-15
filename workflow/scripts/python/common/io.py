import gzip
import subprocess as sp
import hashlib
from typing import Callable, TextIO, TypeVar, IO
from pathlib import Path
from logging import Logger
from common.functional import DesignError

X = TypeVar("X")


def is_gzip_stream(i: IO[bytes]) -> bool:
    return i.read(2) == b"\x1f\x8b"


def is_bgzip_stream(i: IO[bytes]) -> bool:
    return i.read(4) == b"\x1f\x8b\x08\x04"


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with open(p, "rb") as f:
        return is_gzip_stream(f)


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        return is_bgzip_stream(f)


def bgzip(i: IO[bytes], o: IO[bytes]) -> sp.CompletedProcess[bytes]:
    """Stream bgzip to endpoint.

    NOTE: this will block since this is almost always going to be the
    final step in a pipeline.
    """
    return sp.run(["bgzip", "-c"], stdin=i, stdout=o)


def bgzip_file(
    i: IO[bytes],
    p: Path,
    parents: bool = True,
) -> sp.CompletedProcess[bytes]:
    # make parent directory as necessary since snakemake won't make the parent
    # in the case of checkpoints where the output files aren't fully known at
    # parsetime
    p.parent.mkdir(parents=parents, exist_ok=True)
    with open(p, "wb") as f:
        return bgzip(i, f)


def gunzip(i: Path) -> tuple[sp.Popen[bytes], IO[bytes]]:
    """Stream bgzip to endpoint.

    NOTE: this will block since this is almost always going to be the
    final step in a pipeline.
    """
    p = sp.Popen(["gunzip", "-c", i], stdout=sp.PIPE)
    if p.stdout is None:
        raise DesignError()
    return (p, p.stdout)


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
