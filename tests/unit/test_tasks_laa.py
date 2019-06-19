
from zipfile import ZipFile
import tempfile
import logging
import os.path as op
import os
import sys

from pbcore.io import FastqRecord

import pbtestdata

from base import IntegrationBase
from test_file_utils import (make_mock_laa_inputs,
                             make_fastq_inputs)

log = logging.getLogger(__name__)

SUBREADS_IN = pbtestdata.get_file("barcoded-subreadset")


class TestSplitLAATask(IntegrationBase):
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".fastq").name,
        tempfile.NamedTemporaryFile(suffix=".fastq").name,
        SUBREADS_IN
    ]

    @classmethod
    def setUpClass(cls):
        make_fastq_inputs(ofn=cls.INPUT_FILES[0])
        chimera_records = [
            FastqRecord("Barcode3--3_Cluster0_Phase0_NumReads10",
                        "AAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCC",
                        qualityString="&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"),
            FastqRecord("Barcode4--4_Cluster0_Phase0_NumReads11",
                        "AAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTT",
                        qualityString="&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        ]
        make_fastq_inputs(chimera_records, cls.INPUT_FILES[1])

    def test_split_laa_fastq(self):
        consensus_zip = tempfile.NamedTemporaryFile(suffix=".zip").name
        chimeras_zip = tempfile.NamedTemporaryFile(suffix=".zip").name
        args = ["python", "-m", "pbcoretools.tasks2.split_laa_fastq",
                self.INPUT_FILES[0], self.INPUT_FILES[1], self.INPUT_FILES[2],
                consensus_zip, chimeras_zip]
        self._check_call(args)
        with ZipFile(consensus_zip, "r") as zip_out:
            files = zip_out.namelist()
            suffixes = sorted([".".join(of.split('.')[1:]) for of in files])
            self.assertEqual(suffixes, ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq'])


class TestCombinedLAAZip(IntegrationBase):
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".fastq").name,
        tempfile.NamedTemporaryFile(suffix=".csv").name,
        SUBREADS_IN
    ]

    @classmethod
    def setUpClass(cls):
        make_mock_laa_inputs(cls.INPUT_FILES[0], cls.INPUT_FILES[1])

    def test_make_combined_laa_zip(self):
        combined_zip = tempfile.NamedTemporaryFile(suffix=".zip").name
        args = ["python", "-m", "pbcoretools.tasks2.make_combined_laa_zip",
                self.INPUT_FILES[0], self.INPUT_FILES[1], self.INPUT_FILES[2],
                combined_zip]
        self._check_call(args)
        with ZipFile(combined_zip, "r") as zip_out:
            files = zip_out.namelist()
            suffixes = sorted([".".join(of.split('.')[1:]) for of in files])
            self.assertEqual(suffixes, ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq', "csv", "csv"])
