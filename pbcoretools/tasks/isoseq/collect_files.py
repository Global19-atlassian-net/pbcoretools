"""
Utility to create a datastore for a set of per-sample Iso-Seq Analysis outputs.
"""

import logging
import uuid
import os.path as op
import sys

from pbcommand.models import FileTypes, DataStore, DataStoreFile
from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log
from pbcore.io import openDataFile

from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)

FILE_IDS_AND_NAMES = [
    ("flnc_bam", FileTypes.BAM, "Full-Length Non-Concatemer Reads"),
    ("flnc_report", FileTypes.CSV, "Full-Length Non-Concatemer Report"),
    ("cluster_report_csv", FileTypes.CSV, "Cluster Report"),
    ("barcode_overview_report", FileTypes.CSV, "Isoform Counts by Barcode"),
    ("hq_aln_bam", FileTypes.BAM, "Mapped High-Quality Isoforms"),
    ("hq_aln_bai", FileTypes.BAMBAI, "Mapped High-Quality Isoforms (BAM Index)"),
    ("hq_fasta", FileTypes.FASTA, "High-Quality Isoforms"),
    ("lq_fasta", FileTypes.FASTA, "Low-Quality Isoforms"),
    ("collapse_fasta", FileTypes.FASTA, "Collapsed Filtered Isoforms"),
    ("collapse_gff", FileTypes.GFF, "Collapsed Filtered Isoform GFF"),
    ("collapse_group", FileTypes.TXT, "Collapsed Filtered Isoform Groups"),
    ("collapse_abundance", FileTypes.TXT, "Collapsed Filtered Isoform Counts"),
    ("collapse_readstat", FileTypes.TXT, "Full-Length Non-Concatemer Read Assignments"),
]


def to_datastore_file(file_name, file_id, file_type, label):
    return DataStoreFile(uuid.uuid4(),
                         file_id,
                         file_type.file_type_id,
                         op.abspath(file_name),
                         name=label,
                         description=label)


def run_args(args):
    sample_name = None
    if not args.single_sample and not args.all_samples:
        bam = openDataFile(args.samples_file)
        sample_name = bam.readGroupTable[0].SampleName
        log.info("Sample name is {}".format(sample_name))
    elif args.all_samples:
        sample_name = "All Samples"
    files = []
    for file_id, file_type, label in FILE_IDS_AND_NAMES:
        file_path = getattr(args, file_id)
        if file_path is None:
            log.info("Skipping {}".format(file_id))
            continue
        assert file_path is not None and op.exists(file_path)
        if sample_name:
            label += " ({})".format(sample_name)
        files.append(to_datastore_file(file_path, file_id, file_type, label))
    DataStore(files).write_json(args.datastore)
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("samples_file", help="BAM file with sample name")
    for file_id, _, label in FILE_IDS_AND_NAMES:
        p.add_argument("--{}".format(file_id.replace("_", "-")),
                       default=None,
                       help=label)
    p.add_argument("--single-sample",
                   action="store_true",
                   default=False,
                   help="Indicates whether the analysis has only one sample")
    p.add_argument("--all-samples",
                   action="store_true",
                   default=False,
                   help="Indicates whether the outputs are for all samples in a multi-sample experiment")
    p.add_argument("--datastore",
                   default="isoseq_sample.datastore.json",
                   help="Output JSON file name")
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
