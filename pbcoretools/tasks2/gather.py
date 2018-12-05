
"""
Simple gather tool for non-dataset file types
"""

import logging
import sys

from pbcommand.cli import (pacbio_args_runner,
    get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import gather_gff, gather_vcf, gather_csv, gather_fasta_contigset, gather_fastq_contigset

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    if not args.output_file.endswith(args.extension):
        raise ValueError("Output file name {f} does not have expected extension {e}".format(f=args.output_file, e=args.extension))
    args.func(args.chunked_files, args.output_file)
    return 0


def _get_parser():
    modes = [
        ("gff", gather_gff),
        ("vcf", gather_vcf),
        ("csv", gather_csv),
        ("fasta", gather_fasta_contigset),
        ("fastq", gather_fastq_contigset)
    ]
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    sp = p.add_subparsers()
    for mode, runner in modes:
        p_sub = sp.add_parser(mode, help="Gather {m} files".format(m=mode))
        p_sub.set_defaults(func=runner, extension="." + mode)
        p_sub.add_argument("output_file", help="Gathered output file")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
