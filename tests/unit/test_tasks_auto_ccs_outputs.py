import tempfile
import logging
import json
import os.path as op
import os
import pytest
import sys

from pbcore.io import FastaReader
from pbcommand.models import DataStore
from pbcommand.utils import which
from pbcommand.testkit import PbIntegrationBase

from pbcoretools.tasks.auto_ccs_outputs import run_ccs_bam_fastq_exports

import pbtestdata
from base import TESTDATA

log = logging.getLogger(__name__)


@pytest.mark.constools
@pytest.mark.internal_data
class TestAutoCCSOutputs(PbIntegrationBase):
    INPUT_FILE = op.join(
        TESTDATA, "auto_ccs_outputs/m54006_180707_211919.consensusreadset.xml")
    OUTPUT_FILES = [
        "m54006_180707_211919.ccs.bam",
        "m54006_180707_211919.hifi.fasta.zip",
        "m54006_180707_211919.hifi.fastq.zip"
    ]

    def setup_method(self, method):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        PbIntegrationBase.setup_method(self, method)

    def _check_datastore_files_exist(self, file_name):
        file_name = op.abspath(file_name)
        # guard against relative paths
        tmp_dir = tempfile.mkdtemp()
        os.chdir(tmp_dir)
        with open(file_name, "r") as json_in:
            d = json.loads(json_in.read())
            files = d["files"]
            assert all([op.isfile(f["path"]) for f in files])

    def _check_datastore_files(self, files, expected_files):
        file_names = sorted([op.basename(f.path) for f in files])
        assert file_names == expected_files

    def _check_all_datastore_files(self, files):
        self._check_datastore_files(files, self.OUTPUT_FILES)

    def _to_args(self, input_file, mode):
        return [
            "python", "-m",
            "pbcoretools.tasks.auto_ccs_outputs",
            mode,
            input_file,
            "output.datastore.json"
        ]

    def test_run_ccs_bam_fastq_exports(self):
        tmp_dir = tempfile.mkdtemp()
        files = run_ccs_bam_fastq_exports(self.INPUT_FILE, tmp_dir)
        self._check_all_datastore_files(files)

    def test_export_sub_q20(self):
        ds_file = pbtestdata.get_file("rsii-ccs-multi-cell")
        tmp_dir = tempfile.mkdtemp()
        files = run_ccs_bam_fastq_exports(ds_file, tmp_dir, no_zip=True,
                                          include_lowq=True)
        assert len(files) == 5
        bam_file = op.basename(files[0].path)
        assert bam_file == "multiple_movies.ccs.bam"

        def _get_file(id_):
            for ds_file in files:
                if ds_file.source_id == id_:
                    return ds_file.path
        with FastaReader(_get_file("ccs_fasta_out")) as fasta_q20:
            records = [rec.id for rec in fasta_q20]
            assert records == [
                             "m150404_101626_42267_c100807920800000001823174110291514_s1_p0/480/ccs"]
        assert files[4].description == "Non-High Fidelity Reads (FASTQ)"
        with FastaReader(_get_file("ccs_fasta_lq_out")) as fasta_lq:
            records = [rec.id for rec in fasta_lq]
            assert records == ["m150803_002149_42161_c100745121910000001823165807071563_s1_p0/137/ccs"]
        # now with a lower cutoff
        tmp_dir = tempfile.mkdtemp()
        files = run_ccs_bam_fastq_exports(ds_file, tmp_dir, min_rq=0)
        assert len(files) == 3, [f.path for f in files]

    def test_auto_ccs_outputs_barcoded(self):
        input_file = op.join(
            TESTDATA, "auto_ccs_outputs_barcoded/file.datastore.json")
        args = [
            "python", "-m",
            "pbcoretools.tasks.auto_ccs_outputs_barcoded",
            input_file,
            "output.datastore.json"
        ]
        self._check_call(args)
        ds = DataStore.load_from_json("output.datastore.json")
        # 1 FASTQ, 1 FASTA, 2 ZIP
        assert len(ds.files) == 4
        self._check_datastore_files_exist("output.datastore.json")

    def test_auto_ccs_outputs(self):
        modes = ["consolidate", "fasta", "fastq"]
        for mode, output_file in zip(modes, self.OUTPUT_FILES):
            args = self._to_args(self.INPUT_FILE, mode)
            self._check_call(args)
            ds = DataStore.load_from_json("output.datastore.json")
            self._check_datastore_files(ds.files.values(), [output_file])
            self._check_datastore_files_exist("output.datastore.json")

    def test_auto_ccs_outputs_min_qv(self):
        args = self._to_args(self.INPUT_FILE, "fastq") + ["--min-qv", "10"]
        self._check_call(args)
        ds = DataStore.load_from_json("output.datastore.json")
        self._check_datastore_files(ds.files.values(),
                                    ["m54006_180707_211919.hifi.fastq.zip"])

    def test_auto_ccs_outputs_min_qv_2(self):
        ds_file = pbtestdata.get_file("rsii-ccs-multi-cell")
        args = [
            "python", "-m",
            "pbcoretools.tasks.auto_ccs_outputs",
            "fasta",
            ds_file,
            "output.datastore.json"
        ]
        self._check_call(args)
        ds = DataStore.load_from_json("output.datastore.json")
        assert sorted([op.basename(f.path) for f in ds.files.values()]) == [
            "multiple_movies.hifi.fasta.zip"]
        args.extend(["--min-qv", "0"])
        self._check_call(args)
        ds = DataStore.load_from_json("output.datastore.json")
        assert sorted([op.basename(f.path) for f in ds.files.values()]) == [
            "multiple_movies.hifi.fasta.zip"]
