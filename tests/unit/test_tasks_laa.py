
from zipfile import ZipFile
import tempfile
import logging
import os.path as op
import os
import sys

from pbcore.io import FastqRecord
from pbcommand.testkit import PbTestApp

import pbtestdata

from test_file_utils import (make_mock_laa_inputs,
                             make_fastq_inputs)

log = logging.getLogger(__name__)

SUBREADS_IN = pbtestdata.get_file("barcoded-subreadset")


class TestSplitLAATask(PbTestApp):
    TASK_ID = "pbcoretools.tasks.split_laa_fastq"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.laa emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.laa run-rtc '
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

    def run_after(self, rtc, output_dir):
        with ZipFile(rtc.task.output_files[0], "r") as zip_out:
            files = zip_out.namelist()
            suffixes = sorted([".".join(of.split('.')[1:]) for of in files])
            self.assertEqual(suffixes, ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq'])


class TestCombinedLAAZip(PbTestApp):
    TASK_ID = "pbcoretools.tasks.make_combined_laa_zip"
    DRIVER_EMIT = "python -m pbcoretools.tasks.laa emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.laa run-rtc '
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".fastq").name,
        tempfile.NamedTemporaryFile(suffix=".csv").name,
        SUBREADS_IN
    ]

    @classmethod
    def setUpClass(cls):
        make_mock_laa_inputs(cls.INPUT_FILES[0], cls.INPUT_FILES[1])


    def run_after(self, rtc, output_dir):
        with ZipFile(rtc.task.output_files[0], "r") as zip_out:
            files = zip_out.namelist()
            suffixes = sorted([".".join(of.split('.')[1:]) for of in files])
            self.assertEqual(suffixes, ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq', "csv", "csv"])
