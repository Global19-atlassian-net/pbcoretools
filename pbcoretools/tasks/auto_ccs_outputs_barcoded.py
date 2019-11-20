"""
Generate single-file FASTQ outputs from a datastore of demultiplexed
ConsensusReadSets.  Will run individual exports in parallel.
"""

from zipfile import ZipFile
import itertools
import functools
import logging
import uuid
import sys
import os.path as op
import re

from pbcommand.models import FileTypes, DataStoreFile, DataStore
from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log, pool_map

from pbcoretools.tasks.auto_ccs_outputs import run_ccs_bam_fastq_exports
from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


class Constants:
    MAX_NPROC = 8  # just a guess


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("datastore_in",
                   help="DataStore JSON of ConsensusReadSet files")
    p.add_argument("datastore_out", help="DataStore JSON of FASTQ files")
    p.add_argument("--nproc", type=int, default=1,
                   help="Number of processors to use")
    return p


def __run_ccs_bam_fastq_exports(args):
    return run_ccs_bam_fastq_exports(*args, is_barcoded=True)


def __create_zipped_fastx(file_type_id, source_id, ds_files, output_file):
    fastx_files = [f.path for f in ds_files if f.file_type_id == file_type_id]
    with ZipFile(output_file, "w", allowZip64=True) as zip_out:
        for file_name in fastx_files:
            with open(file_name, "r") as fastx_in:
                zip_out.writestr(op.basename(file_name), fastx_in.read())
    file_type_label = file_type_id.split(".")[-1].upper()
    return DataStoreFile(uuid.uuid4(),
                         source_id,
                         FileTypes.ZIP.file_type_id,
                         op.abspath(output_file),
                         name="All Barcodes ({l})".format(l=file_type_label))


_create_zipped_fasta = functools.partial(__create_zipped_fastx,
                                         FileTypes.FASTA.file_type_id,
                                         "pbcoretools.bc_fasta_zip")

_create_zipped_fastq = functools.partial(__create_zipped_fastx,
                                         FileTypes.FASTQ.file_type_id,
                                         "pbcoretools.bc_fastq_zip")


def _run_auto_ccs_outputs_barcoded(datastore_in, datastore_out, nproc=Constants.MAX_NPROC):
    base_dir = op.dirname(datastore_out)
    files = DataStore.load_from_json(datastore_in).files.values()
    ccs_files = []
    for ds_file in files:
        # FIXME use a better file_id
        if ds_file.file_type_id == FileTypes.DS_CCS.file_type_id and ds_file.file_id == "barcoding.tasks.lima-0":
            ccs_files.append(ds_file.path)
            log.info("Exporting %s", ds_file.path)
    log.info("Exporting %d CCS datasets", len(ccs_files))
    args = [(f, base_dir) for f in ccs_files]
    output_files = list(itertools.chain.from_iterable(
        pool_map(__run_ccs_bam_fastq_exports, args, nproc)))
    output_files.extend([
        _create_zipped_fastq(output_files, "all_barcodes.fastq.zip"),
        _create_zipped_fasta(output_files, "all_barcodes.fasta.zip")
    ])
    DataStore(output_files).write_json(datastore_out)
    return 0


def _run_args(args):
    return _run_auto_ccs_outputs_barcoded(args.datastore_in,
                                          args.datastore_out,
                                          nproc=args.nproc)


def _main(argv=sys.argv):
    return pacbio_args_runner(argv[1:],
                              _get_parser(),
                              _run_args,
                              log,
                              setup_log)


if __name__ == "__main__":
    sys.exit(_main(sys.argv))
