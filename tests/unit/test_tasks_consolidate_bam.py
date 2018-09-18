import subprocess
import tempfile
import unittest
import os.path
import sys

import pbcommand.testkit
from pbcommand.models import DataStore
from pbcore.io import AlignmentSet, ConsensusAlignmentSet, openDataSet

import pbtestdata


HAVE_PBMERGE = False
try:
    with tempfile.TemporaryFile() as O, \
         tempfile.TemporaryFile() as E:
        assert subprocess.call(["pbmerge", "--help"], stdout=O, stderr=E) == 0
except Exception as e:
    sys.stderr.write(str(e)+"\n")
    sys.stderr.write("pbmerge missing, skipping test\n")
else:
    HAVE_PBMERGE = True


@unittest.skipUnless(HAVE_PBMERGE, "pbmerge not installed")
class TestConsolidateBam(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m pbcoretools.tasks.consolidate_alignments"
    INPUT_FILES = [pbtestdata.get_file("aligned-ds-2")]
    TASK_OPTIONS = {
        "pbcoretools.task_options.consolidate_aligned_bam": True,
    }

    def run_after(self, rtc, output_dir):
        with openDataSet(rtc.task.output_files[0]) as f:
            f.assertIndexed()
            self.assertEqual(len(f.toExternalFiles()), 1)
            # test for bug 33778
            qnames = set()
            for rec in f:
                qnames.add(rec.qName)
            self.assertEqual(len(qnames), len(f))
        ds = DataStore.load_from_json(rtc.task.output_files[1])
        self.assertEqual(len(ds.files), 2)


@unittest.skipUnless(HAVE_PBMERGE, "pbmerge not installed")
class TestConsolidateBamDisabled(TestConsolidateBam):
    TASK_OPTIONS = {
        "pbcoretools.task_options.consolidate_aligned_bam": False,
    }

    def run_after(self, rtc, output_dir):
        with AlignmentSet(rtc.task.output_files[0]) as f:
            self.assertEqual(len(f.toExternalFiles()), 2)


@unittest.skipUnless(HAVE_PBMERGE, "pbmerge not installed")
class TestConsolidateBamCCS(TestConsolidateBam):
    DRIVER_BASE = "python -m pbcoretools.tasks.consolidate_alignments_ccs"
    INPUT_FILES = [pbtestdata.get_file("rsii-ccs-aligned")]
