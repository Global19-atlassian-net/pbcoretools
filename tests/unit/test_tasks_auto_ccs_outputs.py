
import unittest
import tempfile
import logging
import os.path as op

from pbcommand.testkit import PbTestApp
from pbcommand.models import DataStore

from pbcoretools.tasks.auto_ccs_outputs import run_ccs_bam_fastq_exports

from base import TESTDATA, skip_if_no_testdata

log = logging.getLogger(__name__)

@skip_if_no_testdata
class TestAutoCCSOutputs(PbTestApp):
    DRIVER_BASE = "python -m pbcoretools.tasks.auto_ccs_outputs"
    INPUT_FILES = [
        op.join(TESTDATA, "auto_ccs_outputs/m54006_180707_211919.consensusreadset.xml")
    ]

    def run_after(self, rtc, output_dir):
        ds = DataStore.load_from_json(rtc.task.output_files[0])
        file_names = sorted([op.basename(f.path) for f in ds.files.values()])
        self.assertEqual(file_names, [
            "m54006_180707_211919.Q20.fastq",
            "m54006_180707_211919.ccs.bam"
        ])


@skip_if_no_testdata
class TestAutoCCSBarcodedOutputs(PbTestApp):
    DRIVER_BASE = "python -m pbcoretools.tasks.auto_ccs_outputs_barcoded"
    INPUT_FILES = [
        op.join(TESTDATA, "auto_ccs_outputs_barcoded/file.datastore.json")
    ]

    def run_after(self, rtc, output_dir):
        ds = DataStore.load_from_json(rtc.task.output_files[0])
        self.assertEqual(len(ds.files), 1)
