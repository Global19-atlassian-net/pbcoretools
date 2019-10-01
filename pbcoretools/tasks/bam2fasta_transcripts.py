"""
bam2fasta export with tmp dir support
"""

import logging
import sys

from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log

from pbcoretools.bam2fastx import run_bam_to_fasta
from pbcoretools.file_utils import get_prefixes
from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


def get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("hq_transcripts", help="High-Quality TranscriptSet XML")
    p.add_argument("lq_transcripts", help="Low-Quality TranscriptSet XML")
    p.add_argument("subreads", help="SubreadSet used to generate transcripts")
    p.add_argument("hq_fasta",
                   help="Exported FASTA containing high-quality transcripts")
    p.add_argument("lq_fasta",
                   help="Exported FASTA containing low-quality transcripts")
    return p


def run_args(args):
    hq_prefix, lq_prefix = get_prefixes(args.subreads)
    return max(
        run_bam_to_fasta(args.hq_transcripts, args.hq_fasta,
                         seqid_prefix=hq_prefix),
        run_bam_to_fasta(args.lq_transcripts, args.lq_fasta,
                         seqid_prefix=lq_prefix))


def main(argv=sys.argv):
    return pacbio_args_runner(argv[1:],
                              get_parser(),
                              run_args,
                              log,
                              setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
