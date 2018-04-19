
"""
bam2fasta export with tmp dir support
"""

import logging
import sys

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.tasks.converters import run_bam_to_fasta

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fasta_transcripts"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.bam2fasta_transcripts --resolved-tool-contract"


def get_parser():
    p = get_pbparser(Constants.TOOL_ID,
                     Constants.VERSION,
                     "bam2fasta TranscriptSet export",
                     __doc__,
                     Constants.DRIVER,
                     is_distributed=True,
                     resource_types=(ResourceTypes.TMP_DIR,))
    p.add_input_file_type(FileTypes.DS_TRANSCRIPT, "hq_transcripts",
                          "HQ Transcripts",
                          "High-Quality TranscriptSet XML")
    p.add_input_file_type(FileTypes.DS_TRANSCRIPT, "lq_transcripts",
                          "LQ Transcripts",
                          "Low-Quality TranscriptSet XML")
    p.add_output_file_type(FileTypes.FASTA,
                           "hq_fasta",
                           "High-Quality Transcripts",
                           description="Exported FASTA containing high-quality transcripts",
                           default_name="hq_transcripts")
    p.add_output_file_type(FileTypes.FASTA,
                           "lq_fasta",
                           "Low-Quality Transcripts",
                           description="Exported FASTA containing low-quality transcripts",
                           default_name="lq_transcripts")
    return p


def run_args(args):
    return max(
        run_bam_to_fasta(args.hq_transcripts, args.hq_fasta),
        run_bam_to_fasta(args.lq_transcripts, args.lq_fasta))


def run_rtc(rtc):
    return max(
        run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0],
                         tmp_dir=rtc.task.tmpdir_resources[0].path),
        run_bam_to_fasta(rtc.task.input_files[1], rtc.task.output_files[1],
                         tmp_dir=rtc.task.tmpdir_resources[0].path))



def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           run_args,
                           run_rtc,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
