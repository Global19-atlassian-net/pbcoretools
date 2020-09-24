"""
Memory allocation task for lima in SMRT Link demultiplexing workflows
"""

import logging
import math
import os.path as op
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pacbio_args_runner
from pbcore.io import BarcodeSet, openDataSet

from pbcoretools.utils import get_base_parser


log = logging.getLogger(__name__)
DEFAULT_MEM_GB = 2
BAM_COMPRESSION_FACTOR = 5


def estimate_lima_memory(barcodes, dataset, symmetric):
    """
    Returns estimated system memory for lima, with a minimum of 2GB (which is
    wildly excessive).  Memory is increased only if there are more than 500
    unique barcode pairs being identified.
    """
    n_barcodes = len(BarcodeSet(barcodes, strict=True))
    if symmetric and n_barcodes < 500:
        log.info("symmetric with {n} barcodes, will use default memory".format(
                 n=n_barcodes))
        return DEFAULT_MEM_GB
    else:
        with openDataSet(dataset, skipCounts=True) as reads:
            bioSampleBarcodes = []
            for coll in reads.metadata.collections:
                for bioSample in coll.wellSample.bioSamples:
                    bioSampleBarcodes.extend(bioSample.DNABarcodes)
            n_bc_pairs = len(bioSampleBarcodes)
            if n_bc_pairs == 0:
                log.warn("No biosamples defined, assuming all barcodes used")
                n_bc_pairs = n_barcodes
                if not symmetric:
                    n_bc_pairs *= n_bc_pairs
            if 0 < n_bc_pairs < 500:
                log.info(
                    "only {n} sample barcodes, will use default memory".format(n=n_bc_pairs))
                return DEFAULT_MEM_GB
            else:
                bam_files = [er.bam for er in reads.externalResources]
                bam_size_bytes = sum([op.getsize(f) for f in bam_files])
                bam_size_gb = int(math.ceil(bam_size_bytes / 1024**3))
                log.info(
                    "guessing memory from total BAM size ({m} GB)".format(m=bam_size_gb))
                return DEFAULT_MEM_GB + (bam_size_gb * BAM_COMPRESSION_FACTOR)


def run_args(args):
    mem_gb = estimate_lima_memory(args.barcodes,
                                  args.dataset,
                                  args.symmetric_barcodes)
    log.info("Final memory: {}GB".format(mem_gb))
    with open("lima_mem_gb.txt", mode="wt") as txt_out:
        txt_out.write(str(mem_gb))
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("barcodes", help="BarcodeSet XML or FASTA")
    p.add_argument(
        "dataset", help="Reads (Subreads or CCS) PacBio dataset XML")
    p.add_argument("--symmetric",
                   action="store_true",
                   dest="symmetric_barcodes",
                   default=None,
                   help="Symmetric barcode pairs")
    p.add_argument("--asymmetric",
                   action="store_false",
                   dest="symmetric_barcodes",
                   help="Asymmetric barcode pairs")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(argv[1:],
                              _get_parser(),
                              run_args,
                              log,
                              setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
