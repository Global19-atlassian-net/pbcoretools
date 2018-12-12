
# FIXME this should probably live somewhere more general, e.g. pbdataset?

"""
Consolidate AlignmentSet .bam files
"""

import tempfile
import logging
import os.path as op
import sys

from pbcommand.models import get_pbparser, FileTypes, ResourceTypes, DataStore, DataStoreFile
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import openDataSet, AlignmentSet, ConsensusAlignmentSet, TranscriptAlignmentSet


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.consolidate_alignments"
    VERSION = "0.3.0"
    DRIVER = "python -m pbcoretools.tasks.consolidate_alignments --resolved-tool-contract "
    CONSOLIDATE_ID = "pbcoretools.task_options.consolidate_aligned_bam"
    N_FILES_ID = "pbcoretools.task_options.consolidate_n_files"
    BAI_FILE_TYPES = {
        FileTypes.BAMBAI.file_type_id,
        FileTypes.I_BAI.file_type_id
    }


def get_parser(tool_id=Constants.TOOL_ID,
               file_type=FileTypes.DS_ALIGN,
               driver_exe=Constants.DRIVER,
               version=Constants.VERSION,
               description=__doc__):
    ds_type = file_type.file_type_id.split(".")[-1]
    p = get_pbparser(tool_id,
                     version,
                     "{t} consolidate".format(t=ds_type),
                     description,
                     driver_exe,
                     is_distributed=True,
                     resource_types=(ResourceTypes.TMP_DIR,))

    p.add_input_file_type(file_type,
                          "align_in",
                          "Input {t}".format(t=ds_type),
                          "Gathered {t} to consolidate".format(t=ds_type))
    p.add_output_file_type(file_type,
                           "ds_out",
                           "Alignments",
                           description="Alignment results dataset",
                           default_name="combined")
    p.add_output_file_type(FileTypes.DATASTORE,
                           "datastore",
                           "JSON Datastore",
                           description="Datastore containing BAM resource",
                           default_name="resources")
    p.add_boolean(Constants.CONSOLIDATE_ID, "consolidate",
                  default=False,
                  name="Consolidate .bam",
                  description="Merge chunked/gathered .bam files")
    p.add_int(Constants.N_FILES_ID, "consolidate_n_files",
              default=1,
              name="Number of .bam files",
              description="Number of .bam files to create in consolidate mode")
    return p


def bam_of_dataset(dataset_fn):
    return op.splitext(dataset_fn)[0] + ".bam"


def get_reads_name(ds_in):
    if isinstance(ds_in, TranscriptAlignmentSet):
        return 'Aligned transcripts'
    if isinstance(ds_in, ConsensusAlignmentSet):
        return 'Aligned consensus reads'
    return 'Aligned reads'


def run_consolidate(dataset_file, output_file, datastore_file,
                    consolidate, n_files, task_id=Constants.TOOL_ID,
                    consolidate_f=lambda ds: ds.consolidate):
    datastore_files = []
    with openDataSet(dataset_file) as ds_in:
        if consolidate:
            if len(ds_in.toExternalFiles()) <= 0:
                raise ValueError("DataSet {} must contain one or more files!".format(dataset_file))
            new_resource_file = bam_of_dataset(output_file)
            consolidate_f(ds_in)(new_resource_file, numFiles=n_files, useTmp=False)
            # always display the BAM/BAI if consolidation is enabled
            # XXX there is no uniqueness constraint on the sourceId, but this
            # seems sloppy nonetheless - unfortunately I don't know how else to
            # make view rule whitelisting work
            reads_name = get_reads_name(ds_in)
            for ext_res in ds_in.externalResources:
                if ext_res.resourceId.endswith(".bam"):
                    ds_file = DataStoreFile(
                        ext_res.uniqueId,
                        task_id + "-out-2",
                        ext_res.metaType,
                        ext_res.bam,
                        name=reads_name,
                        description=reads_name)
                    datastore_files.append(ds_file)
                    # Prevent duplicated index files being added to datastore, since consolidated
                    # dataset may contain multiple indices pointing to the same physical file
                    added_resources = set()
                    for index in ext_res.indices:
                        if (index.metaType in Constants.BAI_FILE_TYPES and
                            index.resourceId not in added_resources):
                            added_resources.add(index.resourceId)
                            ds_file = DataStoreFile(
                                index.uniqueId,
                                task_id + "-out-3",
                                index.metaType,
                                index.resourceId,
                                name="Index of {}".format(reads_name.lower()),
                                description="Index of {}".format(reads_name.lower()))
                            datastore_files.append(ds_file)
        ds_in.newUuid()
        ds_in.write(output_file)
    datastore = DataStore(datastore_files)
    datastore.write_json(datastore_file)
    return 0


def args_runner(args, task_id=Constants.TOOL_ID):
    return run_consolidate(
        dataset_file=args.align_in,
        output_file=args.ds_out,
        datastore_file=args.datastore,
        consolidate=args.consolidate,
        n_files=args.consolidate_n_files,
        task_id=task_id)


def rtc_runner(rtc, task_id=Constants.TOOL_ID):
    tempfile.tempdir = rtc.task.tmpdir_resources[0].path
    return run_consolidate(
        dataset_file=rtc.task.input_files[0],
        output_file=rtc.task.output_files[0],
        datastore_file=rtc.task.output_files[1],
        consolidate=rtc.task.options[Constants.CONSOLIDATE_ID],
        n_files=rtc.task.options[Constants.N_FILES_ID],
        task_id=task_id)


def main(argv=sys.argv):
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger()
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
