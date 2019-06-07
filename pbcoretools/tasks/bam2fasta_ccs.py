
"""
Tool contract wrapper for running bam2fasta on CCS reads, with tmp_dir support
and barcoded sample name annotation (requires original SubreadSet input).
"""

import functools
import tempfile
import logging
import gzip
import re
import os.path as op
import os
import sys

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcommand.engine import run_cmd
from pbcommand.utils import walker

from pbcoretools.bam2fastx import run_bam_to_fasta
from pbcoretools.tasks.bam2fasta_archive import get_parser_impl

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fasta_ccs"
    VERSION = "0.3.1"
    DRIVER = "python -m pbcoretools.tasks.bam2fasta_ccs --resolved-tool-contract"
    FILE_TYPE = FileTypes.DS_CCS
    FORMAT_NAME = "fasta"
    TOOL_NAME = "bam2fasta"
    READ_TYPE = "ccs"


def get_parser():
    p = get_parser_impl(Constants)
    p.tool_contract_parser.add_input_file_type(
        FileTypes.DS_SUBREADS,
        "subreads",
        "Raw Subreads",
        "Input SubreadSet XML")
    return p


def run_args(args):
    return run_bam_to_fasta(args.ccs, args.fasta_out)


def run_rtc(rtc):
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0],
                            tmp_dir=rtc.task.tmpdir_resources[0].path,
                            subreads_in=rtc.task.input_files[1])


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           run_args,
                           run_rtc,
                           log,
                           setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
