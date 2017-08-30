
import subprocess
import tempfile
import unittest
import tarfile
import logging
import gzip
import os.path as op
import os
import sys

from pbcore.io import (FastaReader, FastqReader, openDataSet, HdfSubreadSet,
                       SubreadSet, ConsensusReadSet, FastqWriter, FastqRecord)
from pbcommand.testkit import PbTestApp
from pbcommand.utils import which

import pbtestdata

from pbcoretools.tasks.converters import split_laa_fastq, split_laa_fastq_archived, run_bam_to_bam
from pbcoretools import pbvalidate

from base import get_temp_file

log = logging.getLogger(__name__)


class Constants(object):
    BAX2BAM = "bax2bam"
    BAM2FASTA = "bam2fasta"
    BAM2BAM = "bam2bam"
    FASTA2REF = "fasta-to-reference"
    FASTA2GMAP = "fasta-to-gmap-reference"
    SLIMBAM = "slimbam"


SIV_DATA_DIR = "/pbi/dept/secondary/siv/testdata"


def _to_skip_msg(exe):
    return "Missing {e} or {d}".format(d=SIV_DATA_DIR, e=exe)

# XXX hacks to make sure tools are actually available
HAVE_BAX2BAM = which(Constants.BAX2BAM) is not None
HAVE_BAM2BAM = which(Constants.BAM2BAM) is not None
HAVE_BAM2FASTX = which(Constants.BAM2FASTA) is not None
HAVE_FASTA2REF = which(Constants.FASTA2REF) is not None
HAVE_FASTA2GMAP = which(Constants.FASTA2GMAP) is not None
HAVE_SLIMBAM = which(Constants.SLIMBAM) is not None
HAVE_DATA_DIR = op.isdir(SIV_DATA_DIR)
HAVE_DATA_AND_BAX2BAM = HAVE_BAX2BAM and HAVE_DATA_DIR
HAVE_DATA_AND_BAM2BAM = HAVE_BAM2BAM and HAVE_DATA_DIR

SKIP_MSG_BAX2BAM = _to_skip_msg(Constants.BAX2BAM)
SKIP_MSG_BAM2FX = _to_skip_msg(Constants.BAM2FASTA)
SKIP_MSG_BAM2BAM = _to_skip_msg(Constants.BAM2BAM)
SKIP_MSG_FASTA2REF = _to_skip_msg(Constants.FASTA2REF)
SKIP_MSG_FASTA2GMAP = _to_skip_msg(Constants.FASTA2GMAP)
SKIP_MSG_SLIMBAM = _to_skip_msg(Constants.SLIMBAM)

skip_unless_bax2bam = unittest.skipUnless(HAVE_DATA_AND_BAX2BAM, SKIP_MSG_BAX2BAM)
skip_unless_bam2fastx = unittest.skipUnless(HAVE_BAM2FASTX, SKIP_MSG_BAM2FX)
skip_unless_bam2bam = unittest.skipUnless(HAVE_DATA_AND_BAM2BAM, SKIP_MSG_BAM2BAM)
skip_unless_fasta2ref = unittest.skipUnless(HAVE_FASTA2REF, SKIP_MSG_FASTA2REF)
skip_unless_fasta2gmap = unittest.skipUnless(HAVE_FASTA2GMAP, SKIP_MSG_FASTA2GMAP)
skip_unless_slimbam = unittest.skipUnless(HAVE_SLIMBAM, SKIP_MSG_SLIMBAM)


def _get_bax2bam_inputs():
    """Little hackery to get the setup class Inputs and to avoid calls to
    setupclass if skiptest is used

    Nat: we want to test that this behaves properly when multiple movies are
    supplied as input, so we make an HdfSubreadSet on the fly from various
    bax files in testdata
    """
    if HAVE_DATA_AND_BAX2BAM:
        hdf_subread_xml = tempfile.NamedTemporaryFile(suffix=".hdfsubreadset.xml").name

        bax_files = (SIV_DATA_DIR + "/SA3-RS/lambda/2372215/0007_tiny/Analysis_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.bax.h5",
                     pbtestdata.get_file("rsii-bax-h5"))
        ds = HdfSubreadSet(*bax_files)
        ds.name = "lambda_rsii"
        assert len(set([f.movieName for f in ds.resourceReaders()])) == 2
        ds.write(hdf_subread_xml)
        return [hdf_subread_xml]
    else:
        # Assume the test data isn't found and the test won't be run
        return ["/path/to/this-test-should-be-skipped.txt"]


@skip_unless_bax2bam
class TestBax2Bam(PbTestApp):
    TASK_ID = "pbcoretools.tasks.h5_subreads_to_subread"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '

    # See comments above
    INPUT_FILES = _get_bax2bam_inputs()
    MAX_NPROC = 24

    RESOLVED_NPROC = 1
    RESOLVED_TASK_OPTIONS = {}
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True

    def run_after(self, rtc, output_dir):
        with SubreadSet(rtc.task.output_files[0]) as ds_out:
            self.assertEqual(len(ds_out.toExternalFiles()), 2)
            self.assertEqual(ds_out.name, "lambda_rsii")


class Bam2BamCore(object):
    SUBREADS = "/pbi/dept/secondary/siv/testdata/SA3-Sequel/phi29/315/3150101/r54008_20160219_002905/1_A01_micro/m54008_160219_003234_micro.subreadset.xml"
    BARCODES = "/pbi/dept/secondary/siv/barcodes/pacbio_barcodes_384/pacbio_barcodes_384.barcodeset.xml"

    def _validate_subreads(self, output_file):
        err, metrics = pbvalidate.validate_dataset(
            file_name=output_file,
            dataset_type="SubreadSet",
            quick=False,
            validate_index=True,
            strict=True)
        self.assertEqual(len(err), 0)
        with SubreadSet(output_file) as ds:
            self.assertEqual(len(ds.externalResources), 1)
            # make sure metadata are propagated
            md = ds.metadata
            self.assertEqual(
                md.collections.submetadata[0].attrib['InstrumentName'],
                "Inst54008")
            self.assertTrue(ds.externalResources[0].scraps is not None)
            self.assertEqual(ds.externalResources[0].barcodes, self.BARCODES)
            rr = ds.resourceReaders()[0]
            self.assertTrue(rr.pbi.hasBarcodeInfo)
            self.assertTrue("barcoded" in ds.name)
            self.assertTrue("barcoded" in ds.tags)
            #self.assertEqual(len(rr.pbi.bcReverse), 13194)


@skip_unless_bam2bam
class TestBam2Bam(unittest.TestCase, Bam2BamCore):

    def setUp(self):
        # FIXME workaround for 'nose' stupidity
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    def test_run_bam_to_bam(self):
        tmp_out = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        run_bam_to_bam(self.SUBREADS, self.BARCODES, tmp_out)
        self._validate_subreads(tmp_out)


@skip_unless_bam2bam
class TestBam2BamTask(PbTestApp, Bam2BamCore):
    TASK_ID = "pbcoretools.tasks.bam2bam_barcode"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [Bam2BamCore.SUBREADS, Bam2BamCore.BARCODES]
    MAX_NPROC = 8
    RESOLVED_NPROC = 8
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True

    def run_after(self, rtc, output_dir):
        self._validate_subreads(rtc.task.output_files[0])


class _BaseTestBam2Fasta(PbTestApp):
    TASK_ID = None # XXX override in subclasses
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
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
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fasta --emit-tool-contract"
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fasta --resolved-tool-contract"
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2Fasta, cls).setUpClass()


@skip_unless_bam2fastx
class TestBam2FastaFiltered(_BaseTestBam2Fasta):
    NRECORDS_EXPECTED = 13
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fasta --emit-tool-contract"
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fasta --resolved-tool-contract"
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
    DRIVER_EMIT = "python -m pbcoretools.tasks.bam2fasta --emit-tool-contract"
    DRIVER_RESOLVE = "python -m pbcoretools.tasks.bam2fasta --resolved-tool-contract"
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
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
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
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    READER_CLASS = FastqReader
    NRECORDS_EXPECTED = 13
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.filters.addRequirement(length=[('>=', 1000)])
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2FastqFiltered, cls).setUpClass()


def _get_tarred_fastx_file(tar_gz_file):
    file_name = tarfile.TarFile(fileobj=gzip.open(tar_gz_file)).getnames()[0]
    dir_name = op.dirname(tar_gz_file)
    return op.join(dir_name, file_name)


@skip_unless_bam2fastx
class TestBam2FastaArchive(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    NRECORDS_EXPECTED = 117
    SRC_FILE = pbtestdata.get_file("subreads-xml")

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(cls.SRC_FILE, strict=True)
        ds.write(cls.INPUT_FILES[0])
        super(TestBam2FastaArchive, cls).setUpClass()

    def _get_output_file(self, rtc):
        return _get_tarred_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastqArchive(TestBam2Fastq):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)

    def _get_output_file(self, rtc):
        return gzip.open(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastaCCS(_BaseTestBam2Fasta):
    TASK_ID = "pbcoretools.tasks.bam2fasta_ccs"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    INPUT_FILES = [pbtestdata.get_file("rsii-ccs")]
    READER_CLASS = FastaReader
    NRECORDS_EXPECTED = None

    def _get_output_file(self, rtc):
        return _get_tarred_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastqCCS(TestBam2FastaCCS):
    TASK_ID = "pbcoretools.tasks.bam2fastq_ccs"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    READER_CLASS = FastqReader
    NRECORDS_EXPECTED = None

    def _get_output_file(self, rtc):
        return _get_tarred_fastx_file(rtc.task.output_files[0])


@skip_unless_bam2fastx
class TestBam2FastaBarcoded(PbTestApp):
    TASK_ID = "pbcoretools.tasks.bam2fasta_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("barcoded-subreadset")]
    MAX_NPROC = 24
    RESOLVED_NPROC = 1
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True
    READER_CLASS = FastaReader
    EXT = "fasta"

    def run_after(self, rtc, output_dir):
        tmp_dir = tempfile.mkdtemp()
        _cwd = os.getcwd()
        try:
            os.chdir(tmp_dir)
            args = ["tar", "xzf", rtc.task.output_files[0]]
            self.assertEqual(subprocess.call(args), 0)
            file_names = sorted(os.listdir(tmp_dir))
            self.assertEqual(file_names, [
                             "reads.{e}.lbc1__lbc1.{e}".format(e=self.EXT),
                             "reads.{e}.lbc3__lbc3.{e}".format(e=self.EXT),
                             "reads.{e}.unbarcoded.{e}".format(e=self.EXT)
            ])
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
class TestBam2FastqBarcoded(TestBam2FastaBarcoded):
    TASK_ID = "pbcoretools.tasks.bam2fastq_archive"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    READER_CLASS = FastqReader
    EXT = "fastq"


@skip_unless_fasta2ref
class TestFastaToReference(PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta_to_reference"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".fasta").name]
    DATASET_TYPE = "ReferenceSet"

    @classmethod
    def setUpClass(cls):
        with open(cls.INPUT_FILES[0], "w") as fasta:
            fasta.write(">chr1\n")
            fasta.write("\n".join(["".join(["acgta"]*12)]*4))

    def run_after(self, rtc, output_dir):
        from pbcoretools.pbvalidate import validate_dataset
        e, m = validate_dataset(
            file_name=rtc.task.output_files[0],
            dataset_type=self.DATASET_TYPE,
            validate_index=True,
            strict=True)
        self.assertEqual(len(e), 0, str(e))


@skip_unless_fasta2gmap
class TestFastaToGmapReference(TestFastaToReference):
    TASK_ID = "pbcoretools.tasks.fasta_to_gmap_reference"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DATASET_TYPE = "GmapReferenceSet"


class TestFasta2Fofn(PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta2fofn"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("lambda-fasta")]
    IS_DISTRIBUTED = False
    RESOLVED_IS_DISTRIBUTED = False


class TestFasta2ReferenceSet(PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta2referenceset"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("lambda-fasta")]



def _get_fastq_records():
    return [
        FastqRecord("Barcode1--1_Cluster0_Phase0_NumReads91",
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        FastqRecord("Barcode1--1_Cluster0_Phase1_NumReads90",
                    "AAAAAAAGAAAAAAAAAAAAAAATAAAAAAAAAAAAAA",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        FastqRecord("Barcode2--2_Cluster0_Phase0_NumReads91",
                    "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        FastqRecord("Barcode2--2_Cluster0_Phase1_NumReads90",
                    "CAAAAAAGAAAAAAAAAAAAAAATAAAAAAAAAAAAAC",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    ]

def _get_chimera_records():
    return [
        FastqRecord("Barcode3--3_Cluster0_Phase0_NumReads10",
                    "AAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCC",
                    qualityString="&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"),
        FastqRecord("Barcode4--4_Cluster0_Phase0_NumReads11",
                    "AAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTT",
                    qualityString="&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    ]


def _make_fastq_inputs(records, ofn=None):
    if ofn is None:
        ofn = tempfile.NamedTemporaryFile(suffix=".fastq").name
    with FastqWriter(ofn) as fastq_out:
        for rec in records:
            fastq_out.writeRecord(rec)
    return ofn


class TestSplitLAA(unittest.TestCase):
    """
    Unit tests for LAA FASTQ splitter.
    """

    def setUp(self):
        # FIXME workaround for 'nose' stupidity
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self._records = _get_fastq_records()
        self.input_file_name = _make_fastq_inputs(self._records)

    def test_split_laa_fastq(self):
        ifn = self.input_file_name
        ofb = tempfile.NamedTemporaryFile().name
        ofs = split_laa_fastq(ifn, ofb)
        self.assertEqual(len(ofs), 2)
        for i, ofn in enumerate(ofs):
            with FastqReader(ofn) as fastq_in:
                recs = [rec for rec in fastq_in]
                for j in range(2):
                    self.assertEqual(str(recs[j]), str(self._records[(i*2)+j]))

    def test_split_laa_fastq_archived(self):
        ifn = self.input_file_name
        ofn = tempfile.NamedTemporaryFile(suffix=".gz").name
        rc = split_laa_fastq_archived(ifn, ofn)
        self.assertEqual(rc, 0)
        # now with a different extension
        ofn = tempfile.NamedTemporaryFile(suffix=".tar.gz").name
        rc = split_laa_fastq_archived(ifn, ofn)
        self.assertEqual(rc, 0)


class TestSplitLAATask(PbTestApp):
    TASK_ID = "pbcoretools.tasks.split_laa_fastq"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [
        tempfile.NamedTemporaryFile(suffix=".fastq").name,
        tempfile.NamedTemporaryFile(suffix=".fastq").name
    ]

    @classmethod
    def setUpClass(cls):
        _make_fastq_inputs(_get_fastq_records(), cls.INPUT_FILES[0])
        _make_fastq_inputs(_get_chimera_records(), cls.INPUT_FILES[1])


class TestContigSet2Fasta(PbTestApp):
    TASK_ID = "pbcoretools.tasks.contigset2fasta"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("contigset")]


@skip_unless_slimbam
class TestSlimbam(PbTestApp):
    TASK_ID = "pbcoretools.tasks.slimbam"
    DRIVER_EMIT = "python -m pbcoretools.tasks.converters emit-tool-contract {i} ".format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name]

    @classmethod
    def setUpClass(cls):
        ds_start = pbtestdata.get_file("internal-subreads")
        with SubreadSet(ds_start, strict=True) as ds_in:
            ds_in.makePathsAbsolute()
            ds_in.updateCounts()
            ds_in.write(cls.INPUT_FILES[0])

    def run_after(self, rtc, output_dir):
        errors, metrics = pbvalidate.validate_dataset(rtc.task.output_files[0],
            dataset_type="SubreadSet", validate_index=True, strict=True)
        self.assertEqual(len(errors), 0)
        with SubreadSet(rtc.task.input_files[0], strict=True) as ds_in:
            with SubreadSet(rtc.task.output_files[0], strict=True) as ds_out:
                self.assertEqual(ds_in.numRecords, ds_out.numRecords)
                self.assertEqual(ds_in.totalLength, ds_out.totalLength)
                bam_in = ds_in.externalResources[0].resourceId
                bam_out = ds_out.externalResources[0].resourceId
                factor = op.getsize(bam_in) / op.getsize(bam_out)
                self.assertTrue(factor >= 3, "File size larger than expected")
