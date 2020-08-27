"""
Create the final output dataset XML for the Trim PCR Adapters workflow.
"""

import logging
import os.path
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcommand.models import DataStore
from pbcore.io import ConsensusReadSet

from pbcoretools.file_utils import sanitize_dataset_tags
from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


def run_args(args):
    dstore = DataStore.load_from_json(os.path.realpath(args.datastore))
    ds_in = ConsensusReadSet(args.ccs_in, trustCounts=True)
    ds_out = ConsensusReadSet(*([f.path for f in dstore.files.values()]),
                              trustCounts=True)
    sanitize_dataset_tags(ds_out, remove_hidden=True)
    ds_out.name = ds_in.name.replace(" (filtered)", "") + " (trimmed)"
    ds_out.subdatasets = []
    ds_out.write("trimmed.consensusreadset.xml")
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("datastore", help="Datastore of trimmed reads")
    p.add_argument("ccs_in", help="Input (untrimmed) CCS reads")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    main(sys.argv)
