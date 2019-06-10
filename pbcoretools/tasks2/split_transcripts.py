
"""
Divide transcripts into high- and low-quality datasets
"""

import logging
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbcoretools.filters import split_transcripts, Constants

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    return split_transcripts(
        transcripts=args.transcripts,
        hq_file=args.hq_file,
        lq_file=args.lq_file,
        cutoff=args.cutoff)


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("transcripts", help="TranscriptSet XML file")
    p.add_argument("--lq-file", action="store", default="lq.transcriptset.xml",
                   help="Name of output low-quality TranscriptSet")
    p.add_argument("--hq-file", action="store", default="hq.transcriptset.xml",
                   help="Name of output high-quality TranscriptSet")
    p.add_argument("--cutoff", action="store", type=float,
                   default=Constants.TRANSCRIPT_QV_CUTOFF,
                   help="Minimum read quality (0 to 1.0) of high-quality " +
                        "transcripts")
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
