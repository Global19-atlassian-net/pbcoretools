#!/usr/bin/env python

"""
Split LAA outputs by barcode.
"""

import logging
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbcoretools.file_utils import split_laa_fastq_archived

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    # XXX a bit of a hack to support unique file names for the FASTQ tarballs
    return max(split_laa_fastq_archived(args.consensus_in,
                                        args.consensus_out,
                                        args.subreads_in),
               split_laa_fastq_archived(args.chimeras_in,
                                        args.chimeras_out,
                                        args.subreads_in))


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("consensus_in", help="Input consensus FASTQ file")
    p.add_argument("chimeras_in", help="Input chimeras FASTQ file")
    p.add_argument("subreads_in", help="Input SubreadSet XML")
    p.add_argument("consensus_out", help="Output consensus ZIP file")
    p.add_argument("chimeras_out", help="Output chimeras ZIP file")
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
