
"""
bam2fastq zip file export
"""

import functools
import logging
import sys

from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log

from pbcoretools.bam2fastx import run_bam_to_fastq
from pbcoretools.tasks.bam2fasta_archive import get_parser_impl, run_args_impl

log = logging.getLogger(__name__)


class Constants:
    FORMAT_NAME = "fastq"
    TOOL_NAME = "bam2fastq"
    READ_TYPE = "subreads"
    DESC = __doc__


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv[1:],
        get_parser_impl(Constants),
        functools.partial(run_args_impl, run_bam_to_fastq),
        log,
        setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
