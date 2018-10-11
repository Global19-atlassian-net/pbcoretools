
"""
Tool contract wrappers for barcoding-related utilities.
"""

import subprocess
import functools
import logging
import shutil
import os.path as op
import os
import sys

from pbcore.io import (SubreadSet, ConsensusReadSet, openDataSet)
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.utils import get_dataset_metadata
from pbcommand.models import FileTypes, SymbolTypes, OutputFileType, DataStore, DataStoreFile

from pbcoretools.file_utils import (update_barcoded_sample_metadata,
                                    iterate_datastore_read_set_files)


log = logging.getLogger(__name__)


class Constants(object):
    TOOL_NAMESPACE = 'pbcoretools'
    DRIVER_BASE = "python -m pbcoretools.tasks.barcoding "


registry = registry_builder(Constants.TOOL_NAMESPACE, Constants.DRIVER_BASE)


# internal only
@registry("datastore_to_subreads", "0.2.1",
          FileTypes.DATASTORE,
          FileTypes.DS_SUBREADS,
          is_distributed=False,
          nproc=1)
def run_datastore_to_subreads(rtc):
    datasets = list(iterate_datastore_read_set_files(rtc.task.input_files[0]))
    if len(datasets) > 0:
        with SubreadSet(*[f.path for f in datasets], strict=True, skipCounts=True) as ds:
            ds.newUuid()
            ds.write(rtc.task.output_files[0])
    else:
        raise ValueError("Expected one or more SubreadSets in datastore")
    return 0


# internal only
@registry("datastore_to_ccs", "0.1.1",
          FileTypes.DATASTORE,
          FileTypes.DS_CCS,
          is_distributed=False,
          nproc=1)
def run_datastore_to_ccs(rtc):
    datasets = list(iterate_datastore_read_set_files(rtc.task.input_files[0]))
    if len(datasets) > 0:
        with ConsensusReadSet(*[f.path for f in datasets], strict=True, skipCounts=True) as ds:
            ds.newUuid()
            ds.write(rtc.task.output_files[0])
    else:
        raise ValueError("Expected one or more ConsensusReadSets in datastore")
    return 0


@registry("update_barcoded_sample_metadata", "0.5.0",
          (FileTypes.JSON, FileTypes.DS_SUBREADS, FileTypes.DS_BARCODE),
          FileTypes.DATASTORE,
          is_distributed=True,
          nproc=1,
          options={"use_barcode_uuids": True})
def _run_update_barcoded_sample_metadata(rtc):
    use_barcode_uuids = rtc.task.options["pbcoretools.task_options.use_barcode_uuids"]
    base_dir = op.dirname(rtc.task.output_files[0])
    datastore = update_barcoded_sample_metadata(
        base_dir=op.dirname(rtc.task.output_files[0]),
        datastore_file=rtc.task.input_files[0],
        input_reads=rtc.task.input_files[1],
        barcode_set=rtc.task.input_files[2],
        isoseq_mode=False,
        use_barcode_uuids=use_barcode_uuids)
    datastore.write_json(rtc.task.output_files[0])
    return 0


@registry("update_barcoded_sample_metadata_ccs", "0.1.2",
          (FileTypes.JSON, FileTypes.DS_CCS, FileTypes.DS_BARCODE),
          FileTypes.DATASTORE,
          is_distributed=False,
          nproc=1)
def _run_update_barcoded_sample_metadata(rtc):
    base_dir = op.dirname(rtc.task.output_files[0])
    datastore = update_barcoded_sample_metadata(
        base_dir=op.dirname(rtc.task.output_files[0]),
        datastore_file=rtc.task.input_files[0],
        input_reads=rtc.task.input_files[1],
        barcode_set=rtc.task.input_files[2],
        isoseq_mode=True,
        use_barcode_uuids=False)
    datastore.write_json(rtc.task.output_files[0])
    return 0


def _reparent_dataset(input_file, dataset_name, output_file):
    with openDataSet(input_file, strict=True, skipCounts=True) as ds_in:
        if len(ds_in.metadata.provenance) > 0:
            log.warn("Removing existing provenance record: %s",
                     ds_in.metadata.provenance)
            ds_in.metadata.provenance = None
        ds_in.name = dataset_name
        ds_in.newUuid(random=True)
        ds_in.write(output_file)
    return 0


def _reparent_dataset_rtc(rtc):
    NAME_OPT_ID = "pbcoretools.task_options.new_dataset_name"
    if rtc.task.options[NAME_OPT_ID].strip() == "":
        raise ValueError("New dataset name is required")
    return _reparent_dataset(rtc.task.input_files[0],
                             rtc.task.options[NAME_OPT_ID],
                             rtc.task.output_files[0])


ds_name_opt = QuickOpt("", "Name of Output Data Set",
                       "Name of new demultiplexed data set as it appears in " +
                       "SMRT Link")


@registry("reparent_subreads", "0.1.3",
          FileTypes.DS_SUBREADS,
          FileTypes.DS_SUBREADS,
          is_distributed=False,
          nproc=1,
          options={"new_dataset_name": ds_name_opt})
def _run_reparent_subreads(rtc):
    return _reparent_dataset_rtc(rtc)


@registry("reparent_ccs", "0.1.0",
          FileTypes.DS_CCS,
          FileTypes.DS_CCS,
          is_distributed=False,
          nproc=1,
          options={"new_dataset_name": ds_name_opt})
def _run_reparent_subreads(rtc):
    return _reparent_dataset_rtc(rtc)


def _ds_to_datastore(dataset_file, datastore_file,
                     source_id="pbcoretools.tasks.barcoding-out-0"):
    dsmd = get_dataset_metadata(dataset_file)
    ds_file = DataStoreFile(dsmd.uuid, source_id, dsmd.metatype, dataset_file)
    ds_out = DataStore([ds_file])
    ds_out.write_json(datastore_file)
    return 0


@registry("subreads_to_datastore", "0.1.2",
          FileTypes.DS_SUBREADS,
          FileTypes.JSON,
          is_distributed=False,
          nproc=1)
def _run_subreads_to_datastore(rtc):
    return _ds_to_datastore(rtc.task.input_files[0],
                            rtc.task.output_files[0],
                            source_id=rtc.task.task_id + "-out-0")


@registry("ccs_to_datastore", "0.1.2",
          FileTypes.DS_CCS,
          FileTypes.JSON,
          is_distributed=False,
          nproc=1)
def _run_ccs_to_datastore(rtc):
    return _ds_to_datastore(rtc.task.input_files[0],
                            rtc.task.output_files[0],
                            source_id=rtc.task.task_id + "-out-0")


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
