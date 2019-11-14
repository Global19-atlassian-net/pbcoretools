"""
Test application of the 'dataset' tool to BAM consolidation.
"""

import subprocess
import tempfile
import unittest
import os.path as op
import os
import sys

from pbcommand.models import DataStore
from pbcore.io import AlignmentSet, ConsensusAlignmentSet, openDataSet

from pbcoretools.tasks import auto_consolidate

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

    SPLIT_SUBREADS = pbtestdata.get_file("aligned-ds-2")

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp_dir = tempfile.mkdtemp()
        os.chdir(self._tmp_dir)

    def tearDown(self):
        os.chdir(self._cwd)

    @property
    def output_bam(self):
        return op.join(self._tmp_dir, "mapped.bam")

    def _run_and_check_outputs(self, args):
        subprocess.check_call(args)
        self._check_outputs(args[-1])

    def _check_outputs(self, dataset_file):
        self.assertTrue(op.isfile(self.output_bam))
        self.assertTrue(op.isfile(self.output_bam + ".bai"))
        self.assertTrue(op.isfile(self.output_bam + ".pbi"))
        with openDataSet(dataset_file) as f:
            f.assertIndexed()
            self.assertEqual(len(f.toExternalFiles()), 1)
            # test for bug 33778
            qnames = set()
            for rec in f:
                qnames.add(rec.qName)
            self.assertEqual(len(qnames), len(f))

    def test_consolidate_split_alignments(self):
        path = pbtestdata.get_file("aligned-ds-2")
        args = ["dataset", "consolidate", path,
                "mapped.bam", "mapped.alignmentset.xml"]
        self._run_and_check_outputs(args)

    def test_consolidate_ccs_single_file(self):
        path = pbtestdata.get_file("rsii-ccs-aligned")
        args = ["dataset", "consolidate", path,
                "mapped.bam", "mapped.consensusalignmentset.xml"]
        self._run_and_check_outputs(args)

    def _run_auto(self, args):
        base_args = ["auto_consolidate.py"]
        return auto_consolidate.main(base_args + args)

    def _check_datastore(self, file_name):
        ds = DataStore.load_from_json(file_name)
        files = sorted([f.source_id for f in ds.files.values()])
        self.assertEqual(files, ["mapped_bam", "mapped_bam_bai"])

    def test_auto_consolidate_split(self):
        args = [self.SPLIT_SUBREADS, self.output_bam]
        self._run_auto(args)
        xml_file = op.splitext(self.output_bam)[0] + ".alignmentset.xml"
        self._check_outputs(xml_file)
        datastore_file = op.splitext(self.output_bam)[0] + ".datastore.json"
        self._check_datastore(datastore_file)

    def test_auto_consolidate_ccs(self):
        args = [pbtestdata.get_file("rsii-ccs-aligned"), self.output_bam]
        self._run_auto(args)
        xml_file = op.splitext(self.output_bam)[
            0] + ".consensusalignmentset.xml"
        self._check_outputs(xml_file)
        datastore_file = op.splitext(self.output_bam)[0] + ".datastore.json"
        self._check_datastore(datastore_file)

    def test_auto_consolidate_exceeds_cutoff(self):
        args = [self.SPLIT_SUBREADS, self.output_bam, "--max-size", "0"]
        self.assertEqual(self._run_auto(args), 0)
        # this should not have written any files
        self.assertFalse(op.isfile(self.output_bam))
        self.assertFalse(op.isfile(self.output_bam + ".bai"))
        datastore_file = op.splitext(self.output_bam)[0] + ".datastore.json"
        self.assertFalse(op.isfile(datastore_file))
        # now force it
        args = [self.SPLIT_SUBREADS, self.output_bam,
                "--max-size", "0", "--force"]
        self._run_auto(args)
        xml_file = op.splitext(self.output_bam)[0] + ".alignmentset.xml"
        self._check_outputs(xml_file)
        self._check_datastore(datastore_file)
