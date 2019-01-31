
"""
bam2fastq export for TranscriptSet
"""

import logging
import sys

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.bam2fastx import run_bam_to_fastq
from pbcoretools.tasks.bam2fasta_transcripts import get_prefixes

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fastq_transcripts"
    VERSION = "0.1.2"
    DRIVER = "python -m pbcoretools.tasks.bam2fastq_transcripts --resolved-tool-contract"


def get_parser():
    p = get_pbparser(Constants.TOOL_ID,
                     Constants.VERSION,
                     "bam2fastq TranscriptSet export",
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
    p.add_input_file_type(FileTypes.DS_SUBREADS, "subreads",
                          "Input Subreads",
                          "SubreadSet used to generate transcripts")
    p.add_output_file_type(FileTypes.FASTQ,
                           "hq_fastq",
                           "High-Quality Transcripts",
                           description="Exported FASTQ containing high-quality transcripts",
                           default_name="hq_transcripts")
    p.add_output_file_type(FileTypes.FASTQ,
                           "lq_fastq",
                           "Low-Quality Transcripts",
                           description="Exported FASTQ containing low-quality transcripts",
                           default_name="lq_transcripts")
    return p


def run_args(args):
    hq_prefix, lq_prefix = get_prefixes(args.subreads)
    return max(
        run_bam_to_fastq(args.hq_transcripts, args.hq_fastq,
                         seqid_prefix=hq_prefix),
        run_bam_to_fastq(args.lq_transcripts, args.lq_fastq,
                         seqid_prefix=lq_prefix))


def run_rtc(rtc):
    hq_prefix, lq_prefix = get_prefixes(rtc.task.input_files[2])
    return max(
        run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0],
                         tmp_dir=rtc.task.tmpdir_resources[0].path,
                         seqid_prefix=hq_prefix),
        run_bam_to_fastq(rtc.task.input_files[1], rtc.task.output_files[1],
                         tmp_dir=rtc.task.tmpdir_resources[0].path,
                         seqid_prefix=lq_prefix))


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           run_args,
                           run_rtc,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
