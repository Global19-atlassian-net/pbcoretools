"""
Gather gzipped JSON per-ZMW metrics from CCS
"""

from collections import namedtuple
import logging
import gzip
import json
from json import JSONEncoder
import re
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)
__version__ = "0.1"

# XXX these are the only fields we need for pbreports; storing the individual
# dicts as namedtuples conserves memory
FIELDS = "insert_size polymerase_length predicted_accuracy num_full_passes status"

ZmwInfo = namedtuple("ZmwInfo", FIELDS.split())

def _to_zmw_info_hook(d):
    if "insert_size" in d:
        return ZmwInfo(d["insert_size"],
                       d["polymerase_length"],
                       d["predicted_accuracy"],
                       d["num_full_passes"],
                       d["status"])
    return d


class ZmwInfoEncoder(JSONEncoder):
    def default(self, o):  # pylint: disable=method-hidden
        return {getattr(o, name) for name in FIELDS.split()}


def gather_chunks(chunks, output_file):
    ccs_zmws = []
    for file_name in chunks:
        with gzip.open(file_name, mode="rt") as gz_in:
            d = json.loads(gz_in.read(), object_hook=_to_zmw_info_hook)
            ccs_zmws.extend(d["zmws"])
    with gzip.open(output_file, mode="wt") as gz_out:
        gz_out.write(json.dumps({"zmws": ccs_zmws}, cls=ZmwInfoEncoder))
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("merged", help="Name of merged json.gz")
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
