#!/usr/bin/env python

"""
Generate a combined ZIP file for LAA output.
"""

import logging
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbcoretools.file_utils import make_combined_laa_zip

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    return make_combined_laa_zip(
        fastq_file=args.fastq_in,
        summary_csv=args.summary_csv,
        input_subreads=args.subreads_in,
        output_file_name=args.zip_out)


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("fastq_in", help="Input FASTQ file")
    p.add_argument("summary_csv", help="Input Summary CSV file")
    p.add_argument("subreads_in", help="Input SubreadSet XML")
    p.add_argument("zip_out", help="Output ZIP file")
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
