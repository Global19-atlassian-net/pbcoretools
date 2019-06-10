"""
Test application of the 'dataset' tool to BAM consolidation.
"""

import subprocess
import tempfile
import unittest
import os.path
import sys

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
class TestConsolidateBam(unittest.TestCase):
    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp_dir = tempfile.mkdtemp()
        os.chdir(self._tmp_dir)

    def tearDown(self):
        os.chdir(self._cwd)

    def _run_and_check_outputs(self, args):
        subprocess.check_call(args)
        self.assertTrue(os.path.isfile("mapped.bam"))
        self.assertTrue(os.path.isfile("mapped.bam.bai"))
        self.assertTrue(os.path.isfile("mapped.bam.pbi"))
        with openDataSet(args[-1]) as f:
            f.assertIndexed()
            self.assertEqual(len(f.toExternalFiles()), 1)
            # test for bug 33778
            qnames = set()
            for rec in f:
                qnames.add(rec.qName)
            self.assertEqual(len(qnames), len(f))

    def test_consolidate_split_alignments(self):
        path = pbtestdata.get_file("aligned-ds-2")
        args = ["dataset", "consolidate", path, "mapped.bam", "mapped.alignmentset.xml"]
        self._run_and_check_outputs(args)

    def test_consolidate_ccs_single_file(self):
        path = pbtestdata.get_file("rsii-ccs-aligned")
        args = ["dataset", "consolidate", path, "mapped.bam", "mapped.consensusalignmentset.xml"]
        self._run_and_check_outputs(args)
