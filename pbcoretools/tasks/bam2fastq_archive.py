
"""
bam2fastq zip file export with tmp dir support
"""

import logging
import sys

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.tasks.converters import run_bam_to_fastq

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fastq_archive"
    VERSION = "0.3.0"
    DRIVER = "python -m pbcoretools.tasks.bam2fastq_archive --resolved-tool-contract"

def get_parser():
    p = get_pbparser(Constants.TOOL_ID,
                     Constants.VERSION,
                     "bam2fastq export to ZIP",
                     __doc__,
                     Constants.DRIVER,
                     is_distributed=True,
                     resource_types=(ResourceTypes.TMP_DIR,))
    p.add_input_file_type(FileTypes.DS_SUBREADS, "subreads",
                          "Input Subreads",
                          "Input SubreadSet XML")
    p.add_output_file_type(FileTypes.ZIP,
                           "fastq_out",
                           "FASTQ file(s)",
                           description="Exported FASTQ as ZIP archive",
                           default_name="subreads_fastq")
    return p


def run_args(args):
    return run_bam_to_fastq(args.subreads, args.fastq_out)


def run_rtc(rtc):
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0],
                            tmp_dir=rtc.task.tmpdir_resources[0].path)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           run_args,
                           run_rtc,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
