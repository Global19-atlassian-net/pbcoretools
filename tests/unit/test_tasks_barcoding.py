
"""
Tests for barcoding-related tasks used in SMRT Link applications.
"""

import subprocess
import tempfile
import unittest
import logging
import uuid
import os.path as op
import os
import sys

from pbcore.io import openDataSet
from pbcommand.models.common import DataStore
from pbcommand.testkit import PbIntegrationBase

import pbtestdata

from base import get_temp_file
from test_file_utils import (validate_barcoded_datastore_files,
                             split_barcoded_dataset,
                             make_mock_laa_inputs,
                             make_fastq_inputs)

log = logging.getLogger(__name__)


class TestUpdateBarcodedSampleMetadata(PbIntegrationBase):

    def _to_args(self, ds_in, extension=".subreadset.xml", use_barcode_uuids=True):
        barcodes = pbtestdata.get_file("barcodeset")
        ds = split_barcoded_dataset(ds_in, extension)
        ds.write_json("input.datastore.json")
        args = [
            "python", "-m",
            "pbcoretools.tasks2.update_barcoded_sample_metadata",
            ds_in,
            "input.datastore.json",
            barcodes,
            "output.datastore.json"
        ]
        if use_barcode_uuids:
            args.append("--use-barcode-uuids")
        return args

    def _run_update_barcoded_sample_metadata(self, use_barcode_uuids):
        ds_in = pbtestdata.get_file("barcoded-subreadset")
        args = self._to_args(ds_in, use_barcode_uuids=use_barcode_uuids)
        subprocess.check_call(args)
        datastore = DataStore.load_from_json("output.datastore.json")
        validate_barcoded_datastore_files(self, ds_in, datastore,
                                          use_barcode_uuids=use_barcode_uuids)

    def test_update_barcoded_sample_metadata_with_uuids(self):
        self._run_update_barcoded_sample_metadata(use_barcode_uuids=True)

    def test_update_barcoded_sample_metadata_default(self):
        self._run_update_barcoded_sample_metadata(use_barcode_uuids=False)

    def test_update_barcoded_sample_metadata_ccs(self):
        args = self._to_args(pbtestdata.get_file("ccs-barcoded"),
                             ".consensusreadset.xml")
        subprocess.check_call(args)


class TestReparent(PbIntegrationBase):
    DATASET_NAME = "My Data {u}".format(u=uuid.uuid4())

    def _run_and_check_output(self, args):
        subprocess.check_call(args)
        with openDataSet(args[-1], strict=True) as ds_out:
            self.assertEqual(ds_out.name, self.DATASET_NAME)
            with openDataSet(args[-3], strict=True) as ds_in:
                self.assertNotEqual(ds_out.uuid, ds_in.uuid)

    def _to_args(self, file_name, output_file_name):
        return [
            "python", "-m", "pbcoretools.tasks2.reparent_dataset",
            file_name,
            self.DATASET_NAME,
            output_file_name
        ]

    def test_reparent_subreads(self):
        args = self._to_args(pbtestdata.get_file("subreads-sequel"),
                             "new_parent.subreadset.xml")
        self._run_and_check_output(args)

    def test_reparent_ccs(self):
        args = self._to_args(pbtestdata.get_file("ccs-barcoded"),
                             "new_parent.consensusreadset.xml")
        self._run_and_check_output(args)
