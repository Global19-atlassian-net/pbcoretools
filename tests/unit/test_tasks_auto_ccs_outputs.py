
import unittest
import tempfile
import logging
import os.path as op

from pbcommand.testkit import PbTestApp
from pbcommand.models import DataStore

from pbcoretools.tasks.auto_ccs_outputs import run_ccs_bam_fastq_exports

log = logging.getLogger(__name__)

TESTDATA = "/pbi/dept/secondary/siv/testdata/pbcoretools-unittest/data"


@unittest.skipUnless(op.isdir(TESTDATA), "Testdata not found")
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


@unittest.skip("TODO")
class TestAutoCCSBarcodedOutputs(PbTestApp):
    DRIVER_BASE = "python -m pbcoretools.tasks.auto_ccs_outputs_barcoded"
    INPUT_FILES = [
        op.join(TESTDATA, "?")
    ]

    def run_after(self, rtc, output_dir):
        ds = DataStore.load_from_json(rtc.task.output_files[0])
        self.assertEqual(len(ds.files), _N_FILES_TODO)
