from pathlib import Path
from typing_extensions import assert_never
import common.config as cfg
import common.io as io


def main(opath: Path, src: cfg.FileSrc | None, log: Path) -> None:
    if isinstance(src, cfg.LocalSrc):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        Path(opath).symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HTTPSrc):
        url = str(src.url)
        with open(opath, "wb") as f:
            # TODO add bigbed eventually
            if io.curl_test(url, io.is_gzip_stream, log):
                io.curl(url, f, log)
            else:
                io.curl_gzip(url, f, log, False)

    elif src is None:
        assert False, "file src is null; this should not happen"
    else:
        assert_never(src)

    if src.md5 is not None and src.md5 != (actual := io.get_md5(opath, True)):
        with open(log, "a") as f:
            f.write(f"md5s don't match; wanted {src.md5}, actual {actual}")
            exit(1)


main(Path(snakemake.output[0]), snakemake.params.src, snakemake.log[0])  # type: ignore
