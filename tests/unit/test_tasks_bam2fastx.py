"""
Tool contract wrapper tests of pbcoretools.tasks.bam2fast*.
"""

from zipfile import ZipFile
import subprocess
import tempfile
import logging
import shutil
import uuid
import os.path as op
import os
import pytest
import sys

from pbcore.io import (FastaReader, FastqReader, openDataSet,
                       SubreadSet, ConsensusReadSet, FastqWriter, FastqRecord,
                       TranscriptSet)
from pbcommand.testkit import PbIntegrationBase
from pbcommand.utils import which
from pbcommand.models.common import DataStore, DataStoreFile, FileTypes

import pbtestdata

from pbcoretools import pbvalidate

from base import get_temp_file
from test_file_utils import (validate_barcoded_datastore_files,
                             split_barcoded_dataset,
                             make_mock_laa_inputs,
                             make_fastq_inputs)

log = logging.getLogger(__name__)


def _get_zipped_fastx_file(zip_file):
    with ZipFile(zip_file, "r") as zip_in:
        file_name = zip_in.namelist()[0]
        dir_name = op.dirname(zip_file)
        zip_in.extractall(dir_name)
        return op.join(dir_name, file_name)


@pytest.mark.bam2fastx
class TestBam2Fasta(PbIntegrationBase):
    READER_CLASS = FastaReader
    EXTENSION = "fasta"

    def _get_counts(self, input_file, output_file):
        with SubreadSet(input_file) as ds:
            n_expected = len([rec for rec in ds])
        with self.READER_CLASS(output_file) as f:
            n_actual = len([rec for rec in f])
        return n_expected, n_actual

    def run_after(self, input_file, output_file, nrecords_expected=None):
        n_expected, n_actual = self._get_counts(input_file, output_file)
        assert n_actual == n_expected
        if nrecords_expected is not None:
            assert n_actual == nrecords_expected

    def run_and_check(self, args, input_file, output_file, nrecords_expected=None):
        self._check_call(args)
        self.run_after(input_file, output_file, nrecords_expected)

    def run_and_check_fastx(self, input_file, nrecords_expected=None):
        output_file = "exported.{e}".format(e=self.EXTENSION)
        args = [
            "python", "-m",
            "pbcoretools.tasks.bam2{e}".format(e=self.EXTENSION),
            input_file, output_file
        ]
        self.run_and_check(args, input_file, output_file, nrecords_expected)
        return output_file

    def test_bam2fastx(self):
        input_file = pbtestdata.get_file("subreads-xml")
        nrecords_expected = 117
        output_file = self.run_and_check_fastx(input_file, nrecords_expected)
        with SubreadSet(input_file) as ds_in:
            with self.READER_CLASS(output_file) as fa_out:
                for bam_rec, fa_rec in zip(ds_in, fa_out):
                    assert fa_rec.id == bam_rec.qName

    def test_bam2fastx_filtered(self):
        input_file = pbtestdata.get_file("subreads-xml")
        ds = SubreadSet(input_file, strict=True)
        ds.filters.addRequirement(length=[('>=', 1000)])
        input_tmp = get_temp_file(suffix=".subreadset.xml")
        ds.write(input_tmp)
        nrecords_expected = 13
        self.run_and_check_fastx(input_tmp, nrecords_expected)

    def test_bam2fastx_ignore_barcodes(self):
        input_file = pbtestdata.get_file("barcoded-subreadset")
        nrecords_expected = 3
        self.run_and_check_fastx(input_file, nrecords_expected)

    def test_bam2fastx_archive(self):
        input_file = pbtestdata.get_file("subreads-xml")
        nrecords_expected = 117
        output_file = "exported.{e}.zip".format(e=self.EXTENSION)
        args = args = [
            "python", "-m",
            "pbcoretools.tasks.bam2{e}_archive".format(e=self.EXTENSION),
            input_file, output_file
        ]
        self._check_call(args)
        output_file = _get_zipped_fastx_file(output_file)
        self.run_after(input_file, output_file, nrecords_expected)


class TestBam2Fastq(TestBam2Fasta):
    READER_CLASS = FastqReader
    EXTENSION = "fastq"


@pytest.mark.bam2fastx
class TestBam2FastxBarcoded(PbIntegrationBase):
    INPUT_FILE = pbtestdata.get_file("barcoded-subreadset")

    def _get_expected_file_names(self, extension):
        return [
            "subreads.lbc1--lbc1.{e}".format(e=extension),
            "subreads.lbc3--lbc3.{e}".format(e=extension),
            "subreads.unbarcoded.{e}".format(e=extension)
        ]

    def run_after(self, output_file, reader_class, extension):
        output_file = op.abspath(output_file)
        tmp_dir = tempfile.mkdtemp()
        _cwd = os.getcwd()
        try:
            os.chdir(tmp_dir)
            ZipFile(output_file, "r").extractall()
            file_names = sorted(os.listdir(tmp_dir))
            assert file_names == self._get_expected_file_names(extension)
            fastx_ids = ["m54008_160219_003234/74056024/3985_5421",  # bc 0
                         "m54008_160219_003234/28901719/5058_5262",  # bc 2
                         "m54008_160219_003234/4194401/236_10027"]  # bc -1
            for file_name, fastx_id in zip(file_names, fastx_ids):
                with reader_class(file_name) as f:
                    records = [rec.id for rec in f]
                    assert len(records) == 1
                    assert records[0] == fastx_id
        finally:
            os.chdir(_cwd)

    def _run_test_bam2fastx_barcoded(self, input_file, ext, reader_class):
        output_file = "subreads.{e}.zip".format(e=ext)
        args = args = [
            "python", "-m", "pbcoretools.tasks.bam2{e}_archive".format(e=ext),
            input_file, output_file
        ]
        self._check_call(args)
        self.run_after(output_file, reader_class, ext)

    def test_bam2fasta_barcoded(self):
        self._run_test_bam2fastx_barcoded(
            self.INPUT_FILE, "fasta", FastaReader)

    def test_bam2fastq_barcoded(self):
        self._run_test_bam2fastx_barcoded(
            self.INPUT_FILE, "fastq", FastqReader)


@pytest.mark.bam2fastx
class TestBam2FastxBarcodedNoLabels(TestBam2FastxBarcoded):
    INPUT_FILE = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name

    @classmethod
    def setup_class(cls):
        bam_files = []
        with SubreadSet(pbtestdata.get_file("barcoded-subreadset")) as ds_in:
            for er in ds_in.externalResources:
                bam_files.append(er.bam)
        with SubreadSet(*bam_files, strict=True) as ds_out:
            ds_out.write(cls.INPUT_FILE)

    def _get_expected_file_names(self, extension):
        return [
            "subreads.0--0.{e}".format(e=extension),
            "subreads.2--2.{e}".format(e=extension),
            "subreads.unbarcoded.{e}".format(e=extension)
        ]


@pytest.mark.bam2fastx
class TestBam2FastxTranscripts(PbIntegrationBase):
    INPUT_FILES = [
        "/pbi/dept/secondary/siv/testdata/isoseqs/TranscriptSet/polished.hq_tiny.transcriptset.xml",
        "/pbi/dept/secondary/siv/testdata/isoseqs/TranscriptSet/polished.lq_tiny.transcriptset.xml",
        pbtestdata.get_file("subreads-biosample-1")
    ]

    def run_after(self, hq_out, lq_out, ReaderClass):
        def _compare_records(bam_file, fx_records, transcript_type):
            with TranscriptSet(bam_file) as bam_reader:
                for bam_rec, fx_rec in zip(bam_reader.__iter__(), fx_records):
                    seqid = "UnnamedSample_{t}_{i}".format(t=transcript_type,
                                                   i=bam_rec.qName)
                    assert fx_rec.id == seqid
        with ReaderClass(hq_out) as hq_fastx:
            records = [rec for rec in hq_fastx]
            assert len(records) == 1
            _compare_records(self.INPUT_FILES[0], records, "HQ")
        with ReaderClass(lq_out) as lq_fastx:
            records = [rec for rec in lq_fastx]
            assert len(records) == 1
            _compare_records(self.INPUT_FILES[1], records, "LQ")

    def _get_args(self, fastx):
        return [
            "python", "-m",
            "pbcoretools.tasks.bam2{f}_transcripts".format(f=fastx),
            self.INPUT_FILES[0],
            self.INPUT_FILES[1],
            self.INPUT_FILES[2],
            "hq_transcripts.{f}".format(f=fastx),
            "lq_transcripts.{f}".format(f=fastx)
        ]

    def test_bam2fasta_transcripts(self):
        args = self._get_args("fasta")
        self._check_call(args)
        self.run_after(args[-2], args[-1], FastaReader)

    def test_bam2fastq_transcripts(self):
        args = self._get_args("fastq")
        self._check_call(args)
        self.run_after(args[-2], args[-1], FastqReader)
