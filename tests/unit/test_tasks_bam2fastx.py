
"""
Tool contract wrapper tests of pbcoretools.tasks.bam2fast*.
"""

from zipfile import ZipFile
import subprocess
import tempfile
import unittest
import logging
import shutil
import uuid
import os.path as op
import os
import sys

from pbcore.io import (FastaReader, FastqReader, openDataSet,
                       SubreadSet, ConsensusReadSet, FastqWriter, FastqRecord,
                       TranscriptSet)
import pbcommand.testkit
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


class Constants(object):
    BAM2FASTA = "bam2fasta"
    XMLLINT = "xmllint"


SIV_DATA_DIR = "/pbi/dept/secondary/siv/testdata"


def _to_skip_msg(exe):
    return "Missing {e} or {d}".format(d=SIV_DATA_DIR, e=exe)

# XXX hacks to make sure tools are actually available
HAVE_BAM2FASTX = which(Constants.BAM2FASTA) is not None
HAVE_DATA_DIR = op.isdir(SIV_DATA_DIR)
HAVE_XMLLINT = which(Constants.XMLLINT)

SKIP_MSG_BAM2FX = _to_skip_msg(Constants.BAM2FASTA)

skip_unless_bam2fastx = unittest.skipUnless(HAVE_BAM2FASTX, SKIP_MSG_BAM2FX)


class _BaseTestBam2Fasta(pbcommand.testkit.PbTestApp):
    INPUT_FILES = [get_temp_file(suffix=".subreadset.xml")]
    SRC_FILE = None # used to generate INPUT_FILES[0]
    MAX_NPROC = 24
    RESOLVED_NPROC = 1
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True
    READER_CLASS = FastaReader
    NRECORDS_EXPECTED = None

    def _get_output_file(self, rtc):
        return rtc.task.output_files[0]

    def _get_counts(self, rtc):
        with openDataSet(self.INPUT_FILES[0]) as ds:
            n_expected = len([rec for rec in ds])
        with self.READER_CLASS(self._get_output_file(rtc)) as f:
            n_actual = len([rec for rec in f])
        return n_expected, n_actual

    def run_after(self, rtc, output_dir):
        n_expected, n_actual = self._get_counts(rtc)
        self.assertEqual(n_actual, n_expected)
        if self.NRECORDS_EXPECTED is not None:
            self.assertEqual(n_actual, self.NRECORDS_EXPECTED)


@skip_unless_bam2fastx
class TestBam2Fasta(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta"
    NRECORDS_EXPECTED = 117
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta "
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2Fasta, cls).setUpClass()

    def run_after(self, rtc, output_dir):
        super(TestBam2Fasta, self).run_after(rtc, output_dir)
        with SubreadSet(rtc.task.input_files[0]) as ds_in:
            with FastaReader(rtc.task.output_files[0]) as fa_out:
                for bam_rec, fa_rec in zip(ds_in, fa_out):
                    self.assertEqual(fa_rec.id, bam_rec.qName)


@skip_unless_bam2fastx
class TestBam2FastaFiltered(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta "
    NRECORDS_EXPECTED = 13
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.filters.addRequirement(length=[('>=', 1000)])
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2FastaFiltered, cls).setUpClass()


@skip_unless_bam2fastx
class TestBam2FastaIgnoreBarcodes(_BaseTestBam2Fasta):
    """
    Make sure the base bam2fasta task always outputs a single FASTA file
    even when barcoding is present.
    """
    TASK_ID = "pbcoretools.tasks.bam2fasta"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta "
    SRC_FILE = pbtestdata.get_file("barcoded-subreadset")
    NRECORDS_EXPECTED = 3

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2FastaIgnoreBarcodes, cls).setUpClass()


@skip_unless_bam2fastx
class TestBam2Fastq(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fastq"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq "
    READER_CLASS = FastqReader
    NRECORDS_EXPECTED = 117
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2Fastq, cls).setUpClass()


@skip_unless_bam2fastx
class TestBam2FastqFiltered(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fastq"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq "
    READER_CLASS = FastqReader
    NRECORDS_EXPECTED = 13
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.filters.addRequirement(length=[('>=', 1000)])
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2FastqFiltered, cls).setUpClass()


def _get_zipped_fastx_file(zip_file):
    with ZipFile(zip_file, "r") as zip_in:
        file_name = zip_in.namelist()[0]
        dir_name = op.dirname(zip_file)
        zip_in.extractall(dir_name)
        return op.join(dir_name, file_name)


@skip_unless_bam2fastx
class TestBam2FastaArchive(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_archive"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta_archive"
    NRECORDS_EXPECTED = 117
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2FastaArchive, cls).setUpClass()

    def _get_output_file(self, rtc):
        return _get_zipped_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastqArchive(TestBam2Fastq):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fastq_archive --emit-tool-contract "
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fastq_archive --resolved-tool-contract "

    def _get_output_file(self, rtc):
        return _get_zipped_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastaCCS(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_ccs"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta_ccs"
    INPUT_FILES = [
        pbtestdata.get_file("rsii-ccs"),
        pbtestdata.get_file("subreads-sequel") # XXX NOT BARCODED!
    ]
    READER_CLASS = FastaReader
    NRECORDS_EXPECTED = None

    def _get_output_file(self, rtc):
        return _get_zipped_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastqCCS(TestBam2FastaCCS):
    TASK_ID = "pbcoretools.tasks.bam2fastq_ccs"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq_ccs"
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fastq_ccs --emit-tool-contract "
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fastq_ccs --resolved-tool-contract "
    READER_CLASS = FastqReader
    NRECORDS_EXPECTED = None

    def _get_output_file(self, rtc):
        return _get_zipped_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastaBarcoded(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.bam2fasta_archive"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta_archive"
    INPUT_FILES = [pbtestdata.get_file("barcoded-subreadset")]
    MAX_NPROC = 24
    RESOLVED_NPROC = 1
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True
    READER_CLASS = FastaReader
    EXT = "fasta"

    def _get_expected_file_names(self):
        return [
            "subreads.lbc1--lbc1.{e}".format(e=self.EXT),
            "subreads.lbc3--lbc3.{e}".format(e=self.EXT),
            "subreads.unbarcoded.{e}".format(e=self.EXT)
        ]

    def run_after(self, rtc, output_dir):
        tmp_dir = tempfile.mkdtemp()
        _cwd = os.getcwd()
        try:
            os.chdir(tmp_dir)
            ZipFile(rtc.task.output_files[0], "r").extractall()
            file_names = sorted(os.listdir(tmp_dir))
            self.assertEqual(file_names, self._get_expected_file_names())
            fastx_ids = ["m54008_160219_003234/74056024/3985_5421", # bc 0
                         "m54008_160219_003234/28901719/5058_5262", # bc 2
                         "m54008_160219_003234/4194401/236_10027" ] # bc -1
            for file_name, fastx_id in zip(file_names, fastx_ids):
                with self.READER_CLASS(file_name) as f:
                    records = [rec.id for rec in f]
                    self.assertEqual(len(records), 1)
                    self.assertEqual(records[0], fastx_id)
        finally:
            os.chdir(_cwd)


@skip_unless_bam2fastx
class TestBam2FastaCCSBarcoded(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_ccs"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta_ccs"
    INPUT_FILES = [
        pbtestdata.get_file("ccs-barcoded"),
        pbtestdata.get_file("barcoded-subreadset")
    ]
    READER_CLASS = FastaReader
    NRECORDS_EXPECTED = 2
    EXT = "fasta"

    def _get_output_file(self, rtc):
        return _get_zipped_fastx_file(rtc.task.output_files[0])

    def _get_expected_file_names(self):
        return [
            "ccs.Alice.lbc1--lbc1.{e}".format(e=self.EXT),
            "ccs.Charles.lbc3--lbc3.{e}".format(e=self.EXT)
        ]

    def run_after(self, rtc, output_dir):
        tmp_dir = tempfile.mkdtemp()
        _cwd = os.getcwd()
        try:
            os.chdir(tmp_dir)
            ZipFile(rtc.task.output_files[0], "r").extractall()
            file_names = sorted(os.listdir(tmp_dir))
            self.assertEqual(file_names, self._get_expected_file_names())
            fastx_ids = [
                "m54008_160219_003234/46727655/ccs", # bc 0
                "m54008_160219_003234/28901719/ccs", # bc 2
            ]
            for file_name, fastx_id in zip(file_names, fastx_ids):
                with self.READER_CLASS(file_name) as f:
                    records = [rec.id for rec in f]
                    self.assertEqual(len(records), 1)
                    self.assertEqual(records[0], fastx_id)
        finally:
            os.chdir(_cwd)


@skip_unless_bam2fastx
class TestBam2FastqCCSBarcoded(TestBam2FastaCCSBarcoded):
    TASK_ID = "pbcoretools.tasks.bam2fastq_ccs"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq_ccs"
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fastq_ccs --emit-tool-contract "
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fastq_ccs --resolved-tool-contract "
    READER_CLASS = FastqReader
    NRECORDS_EXPECTED = 2
    EXT = "fastq"

    def _get_output_file(self, rtc):
        return _get_zipped_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastaBarcodedNoLabels(TestBam2FastaBarcoded):
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name]

    @classmethod
    def setUpClass(cls):
        bam_files = []
        with SubreadSet(pbtestdata.get_file("barcoded-subreadset")) as ds_in:
            for er in ds_in.externalResources:
                bam_files.append(er.bam)
        with SubreadSet(*bam_files, strict=True) as ds_out:
            ds_out.write(cls.INPUT_FILES[0])
        super(TestBam2FastaBarcodedNoLabels, cls).setUpClass()

    def _get_expected_file_names(self):
        return [
            "subreads.0--0.{e}".format(e=self.EXT),
            "subreads.2--2.{e}".format(e=self.EXT),
            "subreads.unbarcoded.{e}".format(e=self.EXT)
        ]


@skip_unless_bam2fastx
class TestBam2FastqBarcoded(TestBam2FastaBarcoded):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fastq_archive --emit-tool-contract "
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fastq_archive --resolved-tool-contract "
    READER_CLASS = FastqReader
    EXT = "fastq"


@skip_unless_bam2fastx
class TestBam2FastqBarcodedNoLabels(TestBam2FastaBarcodedNoLabels):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fastq_archive --emit-tool-contract "
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fastq_archive --resolved-tool-contract "
    READER_CLASS = FastqReader
    EXT = "fastq"


def _setup_transcripts(hq_file, lq_file):
    from pbcoretools.tasks.filters import split_transcripts
    DS = "/pbi/dept/secondary/siv/testdata/isoseqs/TranscriptSet/polished.transcriptset.xml"
    split_transcripts(DS, hq_file, lq_file, 0.98)


@skip_unless_bam2fastx
class TestBam2FastaTranscripts(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.bam2fasta_transcripts"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fasta_transcripts"
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".transcriptset.xml").name,
        tempfile.NamedTemporaryFile(suffix=".transcriptset.xml").name,
        pbtestdata.get_file("subreads-biosample-1")
    ]
    READER_CLASS = FastaReader

    @classmethod
    def setUpClass(cls):
        _setup_transcripts(cls.INPUT_FILES[0], cls.INPUT_FILES[1])
        super(TestBam2FastaTranscripts, cls).setUpClass()

    def run_after(self, rtc, output_dir):
        def _compare_records(bam_file, fx_records, transcript_type):
            with TranscriptSet(bam_file) as bam_reader:
                for bam_rec, fx_rec in zip(bam_reader.__iter__(), fx_records):
                    seqid = "Alice_{t}_{i}".format(t=transcript_type,
                                                   i=bam_rec.qName)
                    self.assertEqual(fx_rec.id, seqid)
        with self.READER_CLASS(rtc.task.output_files[0]) as hq_fastx:
            records = [rec for rec in hq_fastx]
            self.assertEqual(len(records), 11701)
            _compare_records(rtc.task.input_files[0], records, "HQ")
        with self.READER_CLASS(rtc.task.output_files[1]) as lq_fastx:
            records = [rec for rec in lq_fastx]
            self.assertEqual(len(records), 44)
            _compare_records(rtc.task.input_files[1], records, "LQ")


@skip_unless_bam2fastx
class TestBam2FastqTranscripts(TestBam2FastaTranscripts):
    TASK_ID = "pbcoretools.tasks.bam2fastq_transcripts"
    DRIVER_BASE = "python -m pbcoretools.tasks.bam2fastq_transcripts"
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fastq_transcripts --emit-tool-contract"
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fastq_transcripts --resolved-tool-contract"
    READER_CLASS = FastqReader


# Ensure that pytest ignores the base-class.
del _BaseTestBam2Fasta
