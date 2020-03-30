"""
Gather gzipped JSON per-ZMW metrics from CCS, output as zipped JSON
"""

from zipfile import ZipFile
import logging
import gzip
import json
import re
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)
__version__ = "0.1"


def gather_chunks(chunks, output_file):
    ccs_zmws = []
    for file_name in chunks:
        with gzip.open(file_name, mode="rt") as gz_in:
            d = json.loads(gz_in.read())
            ccs_zmws.extend(d["zmws"])
    # gzip is more C++ friendly, but zip is more SMRTLink-friendly
    with ZipFile(output_file, mode="w") as zip_out:
        base_name = op.basename(re.sub(".zip$", "", output_file))
        zip_out.writestr(base_name, json.dumps({"zmws": ccs_zmws}))
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("merged", help="Name of merged json.zip")
    p.add_argument("chunks", nargs="+", help="Chunk outputs")
    return p


def run_args(args):
    return gather_chunks(args.chunks, args.merged)


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
