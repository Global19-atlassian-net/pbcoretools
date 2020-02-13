
"""
bam2fasta zip file export with tmp dir support.  This module also contains
functionality shared with related tasks.
"""

import functools
import logging
import sys

from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log

from pbcoretools.bam2fastx import run_bam_to_fasta
from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


class Constants:
    FORMAT_NAME = "fasta"
    TOOL_NAME = "bam2fasta"
    READ_TYPE = "subreads"
    DESC = __doc__


def get_parser_impl(constants):
    fmt_name = constants.FORMAT_NAME
    p = get_base_parser(constants.DESC)
    p.add_argument("bam", help="Input PacBio dataset XML or BAM file")
    p.add_argument("{f}_out".format(f=fmt_name),
                   help="Exported {f} as ZIP archive".format(f=fmt_name.upper()))
    return p


def run_args_impl(f, args):
    out = args.fasta_out if hasattr(args, 'fasta_out') else args.fastq_out
    return f(args.bam, out)


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv[1:],
        get_parser_impl(Constants),
        functools.partial(run_args_impl, run_bam_to_fasta),
        log,
        setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
