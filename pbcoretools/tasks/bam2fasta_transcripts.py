
"""
bam2fasta export with tmp dir support
"""

import logging
import sys
import re

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import SubreadSet

from pbcoretools.bam2fastx import run_bam_to_fasta
from pbcoretools.file_utils import get_bio_sample_name

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fasta_transcripts"
    VERSION = "0.1.1"
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
    p.add_input_file_type(FileTypes.DS_SUBREADS, "subreads",
                          "Input Subreads",
                          "SubreadSet used to generate transcripts")
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


# Regular expression pattern of sample strings: must be a string
# of length >= 1, the leading character must be a letter or number,
# the remaining characters must be in [a-zA-Z0-9\_\-]
SAMPLE_CHARSET_RE_STR = '[a-zA-Z0-9\-\_]'
SAMPLE_CHARSET_RE = re.compile(r"{}".format(SAMPLE_CHARSET_RE_STR))

def sanitize_sample(sample):
    """Simple method to sanitize sample to match sample pattern
    ...doctest:
        >>> def a_generator(): return 'a'
        >>> sanitize_sample('1' * 20) # no length limit
        '11111111111111111111'
        >>> sanitize_sample('-123')
        '-123'
        >>> sanitize_sample('123 !&?') # all invalid characters go to '_'
        '123____'
    """
    if len(sample) == 0:
        raise ValueError('Sample must not be an empty string')
    sanitized_sample = ''
    for c in sample:
        if not SAMPLE_CHARSET_RE.search(c):
            c = '_'
        sanitized_sample +=  c
    return sanitized_sample


def get_sanitized_bio_sample_name(subreads):
    """Return sanitized biosample name"""
    sample = get_bio_sample_name(subreads)
    ssample = sanitize_sample(sample)
    log.warning("Sanitize biosample name from {!r} to {!r}".format(sample, ssample))
    return ssample


def get_prefixes(subreads_file):
    with SubreadSet(subreads_file) as subreads:
        seqid_prefix = get_sanitized_bio_sample_name(subreads)
        return ("{}_HQ_".format(seqid_prefix), "{}_LQ_".format(seqid_prefix))


def run_args(args):
    hq_prefix, lq_prefix = get_prefixes(args.subreads)
    return max(
        run_bam_to_fasta(args.hq_transcripts, args.hq_fasta,
                         seqid_prefix=hq_prefix),
        run_bam_to_fasta(args.lq_transcripts, args.lq_fasta,
                         seqid_prefix=lq_prefix))


def run_rtc(rtc):
    hq_prefix, lq_prefix = get_prefixes(rtc.task.input_files[2])
    return max(
        run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0],
                         tmp_dir=rtc.task.tmpdir_resources[0].path,
                         seqid_prefix=hq_prefix),
        run_bam_to_fasta(rtc.task.input_files[1], rtc.task.output_files[1],
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
