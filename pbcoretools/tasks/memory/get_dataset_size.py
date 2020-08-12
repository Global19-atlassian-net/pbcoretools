"""
Dump dataset dimensions to text files for use in Cromwell workflow memory
management
"""

from collections import namedtuple
import logging
import math
import os.path as op
import sys

from pbcore.io import openDataSet
from pbcore.io.align.PacBioBamIndex import get_index_size_bytes
from pbcommand.utils import setup_log
from pbcommand.cli import pacbio_args_runner

from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


ATTRS = ["numRecords", "totalLengthMb", "indexSizeGb", "numResources", "numFilters"]
DsInfo = namedtuple("DsInfo", ATTRS)


def get_dataset_size(dataset, get_index_size=False, skip_counts=True):
    ds = openDataSet(dataset, skipCounts=skip_counts)
    length_mb = int(math.ceil(ds.totalLength / 1024**2))
    index_size = 1
    if get_index_size:
        for ext_res in ds.externalResources:
            if ext_res.pbi:
                index_size += get_index_size_bytes(ext_res.pbi)
        # FIXME this is excessive, please reduce it
        index_size = 1 + math.ceil(2 * index_size / 1024**3)
    n_filt = len(ds.filters)
    n_res = len(ds.externalResources)
    return DsInfo(ds.numRecords, length_mb, index_size, n_res, n_filt)


def run_args(args):
    def logf(p): return log.info("Wrote %s", op.abspath(p))
    m = get_dataset_size(args.dataset, args.get_index_size, args.skip_counts)
    ofns = ["numrecords", "totallength", "indexsize", "numresources", "numfilters"]
    for ofn, attr in zip(ofns, ATTRS):
        ofn = ofn + ".txt"
        with open(ofn, "wt") as f:
            log.info("{a}={n}".format(a=attr, n=getattr(m, attr)))
            f.write(str(getattr(m, attr)))
        logf(ofn)
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("dataset", help="PacBio dataset XML")
    p.add_argument("--skip-counts",
                   action="store_true",
                   default=True,
                   help="Don't load database indices")
    p.add_argument("--no-skip-counts",
                   action="store_false",
                   dest="skip_counts",
                   help="Load database indices and recalculate size")
    p.add_argument("--get-index-size",
                   action="store_true",
                   help="Compute file index size")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(argv[1:],
                              _get_parser(),
                              run_args,
                              log,
                              setup_log)


if __name__ == '__main__':
    sys.exit(main())
