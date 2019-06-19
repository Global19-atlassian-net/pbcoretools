
import unittest
import tempfile
import logging
import os.path as op
import sys

from pbcore.io import FastaReader
from pbcommand.models import DataStore
from pbcommand.utils import which

from pbcoretools.tasks.auto_ccs_outputs import run_ccs_bam_fastq_exports

import pbtestdata
from base import TESTDATA, skip_if_no_testdata, IntegrationBase

HAVE_PBMERGE = which("pbmerge")
skip_unless_pbmerge = unittest.skipUnless(HAVE_PBMERGE, "Missing pbmerge")

log = logging.getLogger(__name__)

@skip_unless_pbmerge
@skip_if_no_testdata
class TestAutoCCSOutputs(IntegrationBase):
    INPUT_FILE = op.join(TESTDATA, "auto_ccs_outputs/m54006_180707_211919.consensusreadset.xml")

    def setUp(self):
        # FIXME workaround for 'nose' conflict with how we run external cmds
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        IntegrationBase.setUp(self)

    def _check_datastore_files(self, files):
        file_names = sorted([op.basename(f.path) for f in files])
        self.assertEqual(file_names, [
            "m54006_180707_211919.Q20.fasta",
            "m54006_180707_211919.Q20.fastq",
            "m54006_180707_211919.ccs.bam"
        ])

    def _to_args(self, input_file):
        return [
            "python", "-m",
            "pbcoretools.tasks.auto_ccs_outputs",
            input_file,
            "output.datastore.json"
        ]

    def test_auto_ccs_outputs(self):
        args = self._to_args(self.INPUT_FILE)
        self._check_call(args)
        ds = DataStore.load_from_json("output.datastore.json")
        self._check_datastore_files(ds.files.values())

    def test_run_ccs_bam_fastq_exports(self):
        tmp_dir = tempfile.mkdtemp()
        files = run_ccs_bam_fastq_exports(self.INPUT_FILE, tmp_dir)
        self._check_datastore_files(files)

    def test_export_sub_q20(self):
        ds_file = pbtestdata.get_file("rsii-ccs-multi-cell")
        tmp_dir = tempfile.mkdtemp()
        files = run_ccs_bam_fastq_exports(ds_file, tmp_dir)
        self.assertEqual(len(files), 5)
        bam_file = op.basename(files[0].path)
        self.assertEqual(bam_file, "m150404_101626_42267_c100807920800000001823174110291514_s1_p0_m150803_002149_42161_c100745121910000001823165807071563_s1_p0.ccs.bam")
        with FastaReader(files[2].path) as fasta_q20:
            records = [rec.id for rec in fasta_q20]
            self.assertEqual(records, ["m150404_101626_42267_c100807920800000001823174110291514_s1_p0/480/ccs"])
        self.assertEqual(files[4].description, "Q0 Reads")
        with FastaReader(files[4].path) as fasta_lq:
            records = [rec.id for rec in fasta_lq]
            self.assertEqual(records, ["m150404_101626_42267_c100807920800000001823174110291514_s1_p0/480/ccs", "m150803_002149_42161_c100745121910000001823165807071563_s1_p0/137/ccs"])

    def test_auto_ccs_outputs_barcoded(self):
        input_file = op.join(TESTDATA, "auto_ccs_outputs_barcoded/file.datastore.json")
        args = [
            "python", "-m",
            "pbcoretools.tasks.auto_ccs_outputs_barcoded",
            input_file,
            "output.datastore.json"
        ]
        self._check_call(args)
        ds = DataStore.load_from_json("output.datastore.json")
        self.assertEqual(len(ds.files), 2)
