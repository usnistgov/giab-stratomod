import os
import threading as tr
import gzip
from pathlib import Path
import subprocess as sp
import common.config as cfg
from typing import Any, IO


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


# logger = setup_logging(snakemake.log[0])  # type: ignore

# temporary columns used for dataframe processing

SLOP = 1


def merge_base(
    config: cfg.StratoMod,
    i: Path,
    o: Path,
    genome: Path,
    base: cfg.Base,
    log: Path,
) -> None:
    baseval = f"unit={base.value}".encode()
    hgroup = config.feature_definitions.homopolymers
    length_col = hgroup.fmt_name(base, lambda x: x.len)
    frac_col = hgroup.fmt_name(base, lambda x: x.imp_frac)
    cols = [x.encode() for x in [*cfg.BED_COLS, length_col, frac_col]]

    mergecmd = ["mergeBed", "-i", "stdin", "-d", "1", "-c", "4", "-o", "sum"]
    slopcmd = ["slopBed", "-i", "stdin", "-b", str(SLOP), "-g", str(genome)]

    # ASSUME input is already sorted
    with gzip.open(i, "r") as ih, gzip.open(o, "w") as oh:
        # write header to output first
        oh.write(b"\t".join(cols) + b"\n")

        r, w = os.pipe()
        r0 = os.fdopen(r, "rb")
        w0 = os.fdopen(w, "wb")

        try:
            # 1. Filter stream for lines that match 'base' and calculate length.
            # Resulting stream will have 4 columns (chrom, start, end, length)
            def filter_base() -> None:
                try:
                    for x in ih:
                        s = x.rstrip().split(b"\t")
                        if baseval == s[3]:
                            start = int(s[1])
                            end = int(s[2])
                            newline = [*s[:3], str(end - start).encode()]
                            w0.write(b"\t".join(newline) + b"\n")
                except Exception:
                    r0.close()
                finally:
                    w0.close()

            t1 = tr.Thread(target=filter_base)
            t1.start()

            # 2. Merge this stream using mergeBed and sum the 4th column
            p1, r1 = spawn_stream(mergecmd, r0)

            # 3. Add slop to stream
            p2, r2 = spawn_stream(slopcmd, r1)

            # 4. Compute real length of (possibly imperfect) homopolymer and
            # imperfect fraction. Write the result. Assume input has 4 columns
            # (chrom, start, end, sum_of_merged_length_before_slop)
            for x in r2:
                s = x.rstrip().split(b"\t")
                start = int(s[1])
                end = int(s[2])
                perfect_length = int(s[3])
                real_length = end - start - SLOP * 2
                frac = 1 - perfect_length / real_length
                newline = [*s[:3], str(real_length).encode(), str(frac).encode()]
                oh.write(b"\t".join(newline) + b"\n")

            check_processes([p1, p2], log)
        finally:
            r0.close()
            w0.close()


def main(smk: Any, config: cfg.StratoMod) -> None:
    merge_base(
        config,
        Path(smk.input["bed"][0]),
        Path(smk.output[0]),
        Path(smk.input["genome"][0]),
        cfg.Base(smk.wildcards["base"]),
        Path(smk.log[0]),
    )


main(snakemake, snakemake.config)  # type: ignore
