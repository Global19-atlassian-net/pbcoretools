
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


class TestSplitLAATask(PbTestApp):
    TASK_ID = "pbcoretools.tasks.split_laa_fastq"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.laa emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.laa run-rtc '
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".fastq").name,
        tempfile.NamedTemporaryFile(suffix=".fastq").name
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


class TestCombinedLAAZip(PbTestApp):
    SUBREADS_IN = pbtestdata.get_file("barcoded-subreadset")
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
