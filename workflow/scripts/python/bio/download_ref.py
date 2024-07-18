from pathlib import Path
import subprocess as sp
from typing import Any
from typing_extensions import assert_never
import common.io as io
import common.config as cfg


def main(smk: Any, params: Any) -> None:
    src: cfg.FileSrc = params.src
    opath = Path(smk.output[0])
    is_fasta = smk.params.is_fasta
    log = smk.log[0]

    if isinstance(src, cfg.LocalSrc):
        # ASSUME this is in the format we indicate (TODO be more paranoid)
        opath.symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HTTPSrc):
        url = str(src.url)
        if is_fasta:
            with open(opath, "wb") as f:
                if io.curl_test(url, io.is_bgzip_stream, log):
                    io.curl(url, f, log)
                elif io.curl_test(url, io.is_gzip_stream, log):
                    p1, o1 = io.spawn_stream(io.curl_cmd(url))
                    p2, o2 = io.gunzip_stream(o1)
                    p3 = io.bgzip(o2, f)
                    o1.close()
                    o2.close()
                    io.check_processes([p1, p2, p3], log)
                else:
                    io.curl_gzip(url, f, log, True)

        else:
            tarcmd = [
                *["bsdtar", "-xf", "-"],
                *["--directory", str(opath)],
                "--strip-component=1",
            ]

            opath.mkdir(parents=True)

            _p1, _o1 = io.spawn_stream(io.curl_cmd(url))
            _p2 = sp.run(tarcmd, stdin=_p1.stdout)
            io.check_processes([_p1, _p2], log)

    else:
        assert_never(src)

    if (
        src.md5 is not None
        and (actual := io.get_md5(opath, True) if is_fasta else io.get_md5_dir(opath))
        != src.md5
    ):
        with open(log, "a") as f:
            f.write(f"md5s don't match; wanted {src.md5}, actual {actual}")
            exit(1)


main(snakemake, snakemake.params)  # type: ignore
