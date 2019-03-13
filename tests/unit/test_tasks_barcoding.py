
"""
Tests for tool contract wrappers in pbcoretools.tasks.barcoding.
"""

import tempfile
import unittest
import logging
import uuid
import os.path as op
import os
import sys

from pbcore.io import SubreadSet, openDataSet, ConsensusReadSet
import pbcommand.testkit
from pbcommand.utils import which
from pbcommand.models.common import DataStore, DataStoreFile, FileTypes

import pbtestdata

from pbcoretools import pbvalidate

from base import get_temp_file
from test_file_utils import (validate_barcoded_datastore_files,
                             split_barcoded_dataset,
                             make_mock_laa_inputs,
                             make_fastq_inputs)

log = logging.getLogger(__name__)


class TestDataStoreToSubreads(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_subreads"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".datastore.json").name]

    @classmethod
    def setUpClass(cls):
        subreads = pbtestdata.get_file("subreads-sequel")
        files = [
            DataStoreFile(uuid.uuid4(), "barcoding.tasks.lima-out-0",
                          FileTypes.DS_SUBREADS.file_type_id, subreads)
        ]
        ds = DataStore(files)
        ds.write_json(cls.INPUT_FILES[0])


class TestDataStoreToCCS(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_ccs"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".datastore.json").name]

    @classmethod
    def setUpClass(cls):
        subreads = pbtestdata.get_file("ccs-barcoded")
        files = [
            DataStoreFile(uuid.uuid4(), "barcoding.tasks.lima-out-0",
                          FileTypes.DS_CCS.file_type_id, subreads)
        ]
        ds = DataStore(files)
        ds.write_json(cls.INPUT_FILES[0])


class TestUpdateBarcodedSampleMetadata(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.update_barcoded_sample_metadata"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".datastore.json").name,
        pbtestdata.get_file("barcoded-subreadset"),
        pbtestdata.get_file("barcodeset")
    ]
    MAX_NPROC = 3
    RESOLVED_NPROC = 3

    @classmethod
    def setUpClass(cls):
        ds = split_barcoded_dataset(cls.INPUT_FILES[1])
        ds.write_json(cls.INPUT_FILES[0])

    def run_after(self, rtc, output_dir):
        datastore = DataStore.load_from_json(rtc.task.output_files[0])
        validate_barcoded_datastore_files(self, self.INPUT_FILES[1], datastore)


class TestUpdateBarcodedSampleMetadataNoUuid(TestUpdateBarcodedSampleMetadata):
    TASK_OPTIONS = {
        "pbcoretools.task_options.use_barcode_uuids": False
    }

    def run_after(self, rtc, output_dir):
        datastore = DataStore.load_from_json(rtc.task.output_files[0])
        validate_barcoded_datastore_files(self, self.INPUT_FILES[1], datastore,
                                          use_barcode_uuids=False)


class TestUpdateBarcodedSampleMetadataCCS(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.update_barcoded_sample_metadata_ccs"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".datastore.json").name,
        pbtestdata.get_file("ccs-barcoded"),
        pbtestdata.get_file("barcodeset")
    ]
    MAX_NPROC = 3
    RESOLVED_NPROC = 3

    @classmethod
    def setUpClass(cls):
        ds = split_barcoded_dataset(cls.INPUT_FILES[1], ".consensusreadset.xml")
        ds.write_json(cls.INPUT_FILES[0])


class TestReparentSubreads(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.reparent_subreads"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [pbtestdata.get_file("subreads-sequel")]
    TASK_OPTIONS = {"pbcoretools.task_options.new_dataset_name": "My Data"}

    def run_after(self, rtc, output_dir):
        with SubreadSet(rtc.task.output_files[0]) as ds_out:
            self.assertEqual(ds_out.name, "My Data")


class TestReparentCCS(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.reparent_ccs"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [pbtestdata.get_file("ccs-barcoded")]
    TASK_OPTIONS = {"pbcoretools.task_options.new_dataset_name": "My Data"}

    def run_after(self, rtc, output_dir):
        with ConsensusReadSet(rtc.task.output_files[0]) as ds_out:
            self.assertEqual(ds_out.name, "My Data")


class TestSubreadsToDataStore(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.subreads_to_datastore"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [pbtestdata.get_file("subreads-sequel")]


class TestCCSToDataStore(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.ccs_to_datastore"
    DRIVER_EMIT = "python -m pbcoretools.tasks.barcoding emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [pbtestdata.get_file("rsii-ccs")]
