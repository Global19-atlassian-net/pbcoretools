
"""
bam2fasta CCS zip file export with tmp dir support
"""

import functools
import logging
import sys

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.bam2fastx import run_bam_to_fasta
from pbcoretools.tasks.bam2fasta_archive import (get_parser_impl, run_args_impl, run_rtc_impl)

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fasta_ccs"
    VERSION = "0.3.1"
    DRIVER = "python -m pbcoretools.tasks.bam2fasta_ccs --resolved-tool-contract"
    FILE_TYPE = FileTypes.DS_CCS
    FORMAT_NAME = "fasta"
    TOOL_NAME = "bam2fasta"
    READ_TYPE = "ccs"


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser_impl(Constants),
                           functools.partial(run_args_impl, run_bam_to_fasta),
                           functools.partial(run_rtc_impl, run_bam_to_fasta),
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
