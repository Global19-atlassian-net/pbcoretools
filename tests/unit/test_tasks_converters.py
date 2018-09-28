
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

from pbcore.io import (FastaReader, FastqReader, openDataSet, HdfSubreadSet,
                       SubreadSet, ConsensusReadSet, FastqWriter, FastqRecord,
                       TranscriptSet)
import pbcommand.testkit
from pbcommand.utils import which
from pbcommand.models.common import DataStore, DataStoreFile, FileTypes

import pbtestdata

from pbcoretools import pbvalidate
from pbcoretools.tasks.barcoding import _ds_to_datastore

from base import get_temp_file
from test_file_utils import (validate_barcoded_datastore_files,
                             split_barcoded_dataset,
                             make_mock_laa_inputs,
                             make_fastq_inputs)

log = logging.getLogger(__name__)


class Constants(object):
    BAX2BAM = "bax2bam"
    FASTA2REF = "fasta-to-reference"
    FASTA2GMAP = "fasta-to-gmap-reference"
    SLIMBAM = "slimbam"


SIV_DATA_DIR = "/pbi/dept/secondary/siv/testdata"


def _to_skip_msg(exe):
    return "Missing {e} or {d}".format(d=SIV_DATA_DIR, e=exe)

# XXX hacks to make sure tools are actually available
HAVE_BAX2BAM = which(Constants.BAX2BAM) is not None
HAVE_FASTA2REF = which(Constants.FASTA2REF) is not None
HAVE_FASTA2GMAP = which(Constants.FASTA2GMAP) is not None
HAVE_SLIMBAM = which(Constants.SLIMBAM) is not None
HAVE_DATA_DIR = op.isdir(SIV_DATA_DIR)
HAVE_DATA_AND_BAX2BAM = HAVE_BAX2BAM and HAVE_DATA_DIR

SKIP_MSG_BAX2BAM = _to_skip_msg(Constants.BAX2BAM)
SKIP_MSG_FASTA2REF = _to_skip_msg(Constants.FASTA2REF)
SKIP_MSG_FASTA2GMAP = _to_skip_msg(Constants.FASTA2GMAP)
SKIP_MSG_SLIMBAM = _to_skip_msg(Constants.SLIMBAM)

skip_unless_bax2bam = unittest.skipUnless(HAVE_DATA_AND_BAX2BAM, SKIP_MSG_BAX2BAM)
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
class TestBax2Bam(pbcommand.testkit.PbTestApp):
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


class _BaseTestBam2Fasta(pbcommand.testkit.PbTestApp):
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


@skip_unless_fasta2ref
class TestFastaToReference(pbcommand.testkit.PbTestApp):
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
        e, m = pbvalidate.validate_dataset(
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


class TestFasta2Fofn(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta2fofn"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("lambda-fasta")]
    IS_DISTRIBUTED = False
    RESOLVED_IS_DISTRIBUTED = False


class TestFasta2ReferenceSet(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.fasta2referenceset"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("lambda-fasta")]


class TestContigSet2Fasta(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.contigset2fasta"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [pbtestdata.get_file("contigset")]


@skip_unless_slimbam
class TestSlimbam(pbcommand.testkit.PbTestApp):
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


def get_f(filename):
    return pbtestdata.get_file(filename)

def generate_datastore(filename):
    # Generate datastore file on the fly
    input_dataset = get_f(filename)
    out_datastore_json = tempfile.NamedTemporaryFile(suffix='datastore.json').name
    _ds_to_datastore(input_dataset, out_datastore_json, 'generate_datastore')
    return out_datastore_json


class TestSubreadSetToDatastore(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.subreads_to_datastore"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.barcoding emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [get_f('subreads-xml')]


class TestCCSToDatastore(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.ccs_to_datastore"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.barcoding emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [get_f('ccs-barcoded')]


class TestTranscriptSetToDatastore(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.transcripts_to_datastore"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [get_f('transcripts-xml')]


class TestDatastoreToSubreadSet(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_subreads"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.barcoding emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [generate_datastore('subreads-xml')]


class TestDatastoreToCCS(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_ccs"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.barcoding emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.barcoding run-rtc '
    INPUT_FILES = [generate_datastore('ccs-barcoded')]


class TestDatastoreToTranscriptSet(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_transcripts"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [generate_datastore('transcripts-xml')]


class TestDatastoreToAlignmentSet(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_alignments"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [generate_datastore('aligned-ds-2')]


class TestDatastoreToCCSAlignmentSet(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_ccs_alignments"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [generate_datastore('ccs-xml-aligned')]


class TestDatastoreToTranscriptAlignmentSet(pbcommand.testkit.PbTestApp):
    TASK_ID = "pbcoretools.tasks.datastore_to_transcript_alignments"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.converters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.converters run-rtc '
    INPUT_FILES = [generate_datastore('transcript-alignments-xml')]
