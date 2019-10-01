
"""
Command line wrapper for running bam2fasta
"""

import logging
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pacbio_args_runner

from pbcoretools.bam2fastx import run_bam_to_fasta
from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


def get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("subreads", help="Input SubreadSet XML")
    p.add_argument("fasta_out", help="Exported FASTA")
    return p


def run_args(args):
    return run_bam_to_fasta(args.subreads, args.fasta_out)


def main(argv=sys.argv):
    return pacbio_args_runner(argv[1:],
                              get_parser(),
                              run_args,
                              log,
                              setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
