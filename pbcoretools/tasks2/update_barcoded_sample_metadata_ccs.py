
"""
Update barcoded sample metadata ccs
"""

import logging
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
    get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbcoretools.file_utils import update_barcoded_sample_metadata

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    ccs = args.ccs
    lima_datastore = args.lima_datastore
    barcodes = args.barcodes
    out_json = args.out_json

    datastore = update_barcoded_sample_metadata(
        base_dir=op.dirname(lima_datastore),
        datastore_file=lima_datastore,
        input_reads=ccs,
        barcode_set=barcodes,
        isoseq_mode=True,
        use_barcode_uuids=False)
    datastore.write_json(out_json)
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("ccs", help="ConsensusReadSet")
    p.add_argument("lima_datastore", help="Datastore json generated by lima to demultiplex ccs.")
    p.add_argument("barcodes", help="BarcodeSet lima used to demultiplex ccs")
    p.add_argument("out_json", help="Output datastore json")
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
