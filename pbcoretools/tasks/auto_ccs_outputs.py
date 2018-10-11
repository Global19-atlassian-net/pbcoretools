"""
Generate single-file CCS BAM and FASTQ outputs from a ConsensusReadSet.
"""

import logging
import uuid
import sys
import os.path as op
import re

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser, DataStoreFile, DataStore
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import ConsensusReadSet

from pbcoretools.bam2fastx import run_bam_to_fastq

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.auto_ccs_outputs"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.auto_ccs_outputs --resolved-tool-contract"
    FASTQ_EXT = ".Q20.fastq"
    BAM_EXT = ".ccs.bam"
    FASTQ_ID = TOOL_ID + "-out-1"
    BAM_ID = TOOL_ID + "-out-2"


def _get_parser():
    p = get_pbparser(Constants.TOOL_ID,
                     Constants.VERSION,
                     "Generate primary CCS outputs",
                     __doc__,
                     Constants.DRIVER,
                     is_distributed=True)
                     #resource_types=(ResourceTypes.TMP_DIR,))
    p.add_input_file_type(FileTypes.DS_CCS, "ccs_dataset",
                          "ConsensusReadSet XML",
                          "ConsensusReadSet XML")
    p.add_output_file_type(FileTypes.DATASTORE,
                           "datastore_out",
                           "DataStore JSON",
                           description="DataStore JSON",
                           default_name="ccs_outputs")
    return p


def run_ccs_bam_fastq_exports(ccs_dataset_file, base_dir, is_barcoded=False):
    """
    Take a ConsensusReadSet and write BAM/FASTQ files to the output
    directory.  If this is a demultiplexed dataset, it is assumed to have
    a single BAM file within a dataset that is already imported in SMRT Link.
    """
    datastore_files = []
    with ConsensusReadSet(ccs_dataset_file, strict=True) as ccs:
        bam_file_name = None
        if is_barcoded:
            barcodes = list(set(zip(ccs.index.bcForward, ccs.index.bcReverse)))
            assert len(barcodes) == 1, "Multiple barcodes found: {b}".format(
                b=", ".join([str(b) for b in barcodes]))
            assert len(ccs.externalResources) == 1
            ccs_bam = ccs.externalResources[0].bam
            log.info("Found a single barcoded BAM %s", ccs_bam)
            file_prefix = re.sub(".ccs.bam$", "", op.basename(ccs_bam))
        else:
            movies = sorted(list({rg.MovieName for rg in ccs.readGroupTable}))
            if len(movies) > 1:
                log.warn("Multiple movies found: %s", movies)
                log.warn("This pipeline was not designed to run on more than")
                log.warn("one instrument movie.  Results not guaranteed.")
            file_prefix = "_".join(movies)
            bam_file_name = op.join(base_dir, file_prefix + Constants.BAM_EXT)
            ccs.consolidate(bam_file_name)
            datastore_files.append(
                DataStoreFile(uuid.uuid4(),
                              Constants.BAM_ID,
                              FileTypes.BAM_CCS.file_type_id,
                              bam_file_name,
                              name=op.basename(bam_file_name), # XXX is this right?
                              description="CCS BAM file"))
        fastq_file = op.join(base_dir, file_prefix + Constants.FASTQ_EXT)
        run_bam_to_fastq(ccs_dataset_file, fastq_file)
        datastore_files.append(
            DataStoreFile(uuid.uuid4(),
                          Constants.FASTQ_ID,
                          FileTypes.FASTQ.file_type_id,
                          fastq_file,
                          name=op.basename(fastq_file),
                          description="Q20 Reads"))
    return datastore_files


def _run_auto_ccs_outputs(ccs_dataset, datastore_out):
    base_dir = op.dirname(datastore_out)
    DataStore(run_ccs_bam_fastq_exports(ccs_dataset,
                                        base_dir)).write_json(datastore_out)
    return 0


def _run_args(args):
    return _run_auto_ccs_outputs(args.ccs_dataset, args.datastore_out)


def _run_rtc(rtc):
    return _run_auto_ccs_outputs(rtc.task.input_files[0],
                                 rtc.task.output_files[0])


def _main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           _get_parser(),
                           _run_args,
                           _run_rtc,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(_main(sys.argv))
