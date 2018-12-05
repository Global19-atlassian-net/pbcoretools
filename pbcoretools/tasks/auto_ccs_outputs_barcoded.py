"""
Generate single-file FASTQ outputs from a datastore of demultiplexed
ConsensusReadSets.  Will run individual exports in parallel.
"""

import itertools
import logging
import uuid
import sys
import os.path as op
import re

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser, DataStoreFile, DataStore
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log, pool_map

from pbcoretools.tasks.auto_ccs_outputs import run_ccs_bam_fastq_exports

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.auto_ccs_outputs_barcoded"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.auto_ccs_outputs_barcoded --resolved-tool-contract"
    MAX_NPROC = 8 # just a guess


def _get_parser():
    p = get_pbparser(Constants.TOOL_ID,
                     Constants.VERSION,
                     "Generate primary demultiplexed CCS outputs",
                     __doc__,
                     Constants.DRIVER,
                     is_distributed=True,
                     nproc=Constants.MAX_NPROC)
                     #resource_types=(ResourceTypes.TMP_DIR,))
    p.add_input_file_type(FileTypes.DATASTORE, "datastore_in",
                          "DataStore JSON",
                          "DataStore JSON of ConsensusReadSet files")
    p.add_output_file_type(FileTypes.DATASTORE,
                           "datastore_out",
                           "DataStore JSON",
                           description="DataStore JSON of FASTQ files",
                           default_name="ccs_outputs")
    return p


def __run_ccs_bam_fastq_exports(args):
    return run_ccs_bam_fastq_exports(*args, is_barcoded=True)


def _run_auto_ccs_outputs_barcoded(datastore_in, datastore_out, nproc=Constants.MAX_NPROC):
    base_dir = op.dirname(datastore_out)
    files = DataStore.load_from_json(datastore_in).files.values()
    ccs_files = []
    for ds_file in files:
        if ds_file.file_type_id == FileTypes.DS_CCS.file_type_id and ds_file.file_id == "barcoding.tasks.lima-0":
            ccs_files.append(ds_file.path)
            log.info("Exporting %s", ds_file.path)
    log.info("Exporting %d CCS datasets", len(ccs_files))
    args = [(f, base_dir) for f in ccs_files]
    output_files = itertools.chain.from_iterable(
        pool_map(__run_ccs_bam_fastq_exports, args, nproc))
    DataStore(output_files).write_json(datastore_out)
    return 0


def _run_args(args):
    return _run_auto_ccs_outputs_barcoded(args.datastore_in, args.datastore_out)


def _run_rtc(rtc):
    return _run_auto_ccs_outputs_barcoded(rtc.task.input_files[0],
                                          rtc.task.output_files[0],
                                          nproc=rtc.task.nproc)


def _main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           _get_parser(),
                           _run_args,
                           _run_rtc,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(_main(sys.argv))
