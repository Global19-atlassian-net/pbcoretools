"""
Dataset filtering tool that incorporates downsampling; replaces old pbsmrtpipe
task pbcoretools.tasks.filter_dataset.
"""

import logging
import re
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcore.util.statistics import phred_qv_as_accuracy

from pbcoretools.filters import run_filter_dataset

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    return run_filter_dataset(
        in_file=args.dataset,
        out_file=args.xml_out,
        read_length=args.min_read_length,
        other_filters=args.filters,
        downsample_factor=args.downsample,
        min_rq=args.min_rq)


def _phred_qv_as_accuracy(x):
    """
    Wrapper for phred_qv_as_accuracy that allows values of -1, which is
    ignored (no filter added).  This is useful for CCS applications where we
    want to allow every read through, even the ones with RQ = -1.
    """
    x = int(x)
    if x == -1:
        return -1
    else:
        return phred_qv_as_accuracy(x)


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("dataset", help="Dataset XML input file")
    p.add_argument("xml_out", help="Output dataset XML file")
    p.add_argument("filters", help="Filters as string")
    p.add_argument("--downsample", action="store", type=int, default=0,
                   help="Downsampling factor")
    p.add_argument("--min-read-length", action="store", type=int, default=0,
                   help="Minimum read length")
    p.add_argument("--min-rq", action="store", type=float, default=None,
                   help="Minimum read score/quality (range 0.0-1.0)")
    p.add_argument("--min-qv",
                   dest="min_rq",
                   type=_phred_qv_as_accuracy,
                   help="Alternative to --min-rq, as integer on Phred scale (0-60)")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
