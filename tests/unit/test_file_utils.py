
"""
Unit tests for pbcoretools.file_utils
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

from pbcore.io import (FastaReader, FastqReader, openDataSet, HdfSubreadSet,
                       SubreadSet, ConsensusReadSet, FastqWriter, FastqRecord,
                       TranscriptSet)
from pbcommand.models.common import DataStore, DataStoreFile, FileTypes
from pbcommand.utils import which

import pbtestdata

from pbcoretools.file_utils import (
    split_laa_fastq,
    split_laa_fastq_archived,
    get_ds_name,
    get_barcode_sample_mappings,
    make_barcode_sample_csv,
    make_combined_laa_zip,
    update_barcoded_sample_metadata,
    discard_bio_samples)


HAVE_XMLLINT = which("xmllint")

def _validate_dataset_xml(file_name):
    if HAVE_XMLLINT and "PB_DATASET_XSD" in os.environ:
        args = ["xmllint", "--schema", os.environ["PB_DATASET_XSD"], file_name]
        subprocess.check_output(args, stderr=subprocess.STDOUT)


def _get_fastq_records():
    # these correspond to barcodes in the barcoded-subreadset dataset
    # provided by pbtestdata
    return [
        FastqRecord("Barcodelbc1--lbc1_Cluster0_Phase0_NumReads91",
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        FastqRecord("Barcodelbc1--lbc1_Cluster0_Phase1_NumReads90",
                    "AAAAAAAGAAAAAAAAAAAAAAATAAAAAAAAAAAAAA",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        FastqRecord("Barcodelbc3--lbc3_Cluster0_Phase0_NumReads91",
                    "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        FastqRecord("Barcodelbc3--lbc3_Cluster0_Phase1_NumReads90",
                    "CAAAAAAGAAAAAAAAAAAAAAATAAAAAAAAAAAAAC",
                    qualityString="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    ]


def make_fastq_inputs(records=None, ofn=None):
    if records is None:
        records = _get_fastq_records()
    if ofn is None:
        ofn = tempfile.NamedTemporaryFile(suffix=".fastq").name
    with FastqWriter(ofn) as fastq_out:
        for rec in records:
            fastq_out.writeRecord(rec)
    return ofn


def make_mock_laa_inputs(fastq_file, csv_file):
    make_fastq_inputs(_get_fastq_records(), fastq_file)
    csv_tmp = open(csv_file, "w")
    csv_tmp.write("a,b\nc,d\ne,f")
    csv_tmp.close()


def split_barcoded_dataset(file_name, ext=".subreadset.xml"):
    from pbcoretools.bamSieve import filter_reads
    ds_in = openDataSet(file_name)
    ds_dir = tempfile.mkdtemp()
    ds_files = []
    for bc, label in zip([0, 2], ["lbc1--lbc1", "lbc3--lbc3"]):
        ds_tmp = op.join(ds_dir, "lima_output.{l}{e}".format(l=label, e=ext))
        filter_reads(
            input_bam=file_name,
            output_bam=ds_tmp,
            whitelist=[bc],
            use_barcodes=True)
        ds_files.append(DataStoreFile(uuid.uuid4(),
                                      "barcoding.tasks.lima-out-0",
                                      ds_in.datasetType,
                                      ds_tmp))
    return DataStore(ds_files)


def validate_barcoded_datastore_files(self, subreads, datastore):
    bio_sample_names = {
        "lbc1--lbc1": "Alice",
        "lbc3--lbc3": "Charles"
    }
    self.assertEqual(len(datastore.files), 2)
    ds_in = SubreadSet(subreads)
    for f in datastore.files.values():
        with SubreadSet(f.path) as ds:
            # FIXME need better testing here
            self.assertEqual(len(ds.filters), 1)
            bc_label = op.basename(f.path).split(".")[1]
            bio_name = bio_sample_names[bc_label]
            coll = ds.metadata.collections[0]
            self.assertEqual(len(coll.wellSample.bioSamples), 1)
            self.assertEqual(coll.wellSample.bioSamples[0].name, bio_name)
            self.assertEqual(ds.metadata.provenance.parentDataSet.uniqueId,
                             ds_in.uuid)
            self.assertEqual(ds.name, "{n} ({s})".format(n=ds_in.name,
                                                         s=bio_name))
            self.assertEqual(ds.uuid, f.uuid)
            md_tags = [r['tag'] for r in ds.metadata.record['children']]
            self.assertEqual(md_tags[0:4],
                             ["TotalLength", "NumRecords", "Provenance", "Collections"])
            _validate_dataset_xml(f.path)


class TestSplitLAA(unittest.TestCase):
    """
    Unit tests for LAA FASTQ splitter.
    """

    SUBREADS = pbtestdata.get_file("barcoded-subreadset")

    def setUp(self):
        # FIXME workaround for 'nose' stupidity
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self._records = _get_fastq_records()
        self.input_file_name = make_fastq_inputs(self._records)
        self._subreads = pbtestdata.get_file("barcoded-subreadset")

    def test_split_laa_fastq(self):
        ifn = self.input_file_name
        ofb = tempfile.NamedTemporaryFile().name
        ofs = split_laa_fastq(ifn, ofb, self._subreads)
        self.assertEqual(len(ofs), 2)
        suffixes = sorted([".".join(of.split('.')[1:]) for of in ofs])
        self.assertEqual(suffixes, ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq'])
        for i, ofn in enumerate(ofs):
            with FastqReader(ofn) as fastq_in:
                recs = [rec for rec in fastq_in]
                for j in range(2):
                    self.assertEqual(str(recs[j]), str(
                        self._records[(i * 2) + j]))

    def test_split_laa_fastq_archived(self):
        ifn = self.input_file_name
        ofn = tempfile.NamedTemporaryFile(suffix=".zip").name
        rc = split_laa_fastq_archived(ifn, ofn, self._subreads)
        self.assertEqual(rc, 0)
        with ZipFile(ofn, "r") as zip_in:
            files = zip_in.namelist()
            suffixes = sorted([".".join(of.split('.')[1:]) for of in files])
            self.assertEqual(suffixes, ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq'])

    def test_get_barcode_sample_mappings(self):
        with SubreadSet(self._subreads) as ds:
            # just double-checking that the XML defines more samples than are
            # actually present in the BAM
            assert len(ds.metadata.collections[0].wellSample.bioSamples) == 3
        samples = get_barcode_sample_mappings(SubreadSet(self._subreads))
        self.assertEqual(samples, {'lbc3--lbc3': 'Charles',
                                   'lbc1--lbc1': 'Alice'})

    def test_make_barcode_sample_csv(self):
        csv_file = tempfile.NamedTemporaryFile(suffix=".csv").name
        sample_mappings = make_barcode_sample_csv(self._subreads, csv_file)
        self.assertEqual(sample_mappings, {'lbc3--lbc3': 'Charles', 'lbc1--lbc1': 'Alice'})
        with open(csv_file) as f:
            self.assertEqual(
                f.read(), "Barcode Name,Bio Sample Name\nlbc1--lbc1,Alice\nlbc3--lbc3,Charles\n")

    def test_make_combined_laa_zip(self):
        INPUT_FILES = [
            tempfile.NamedTemporaryFile(suffix=".fastq").name,
            tempfile.NamedTemporaryFile(suffix=".csv").name,
            pbtestdata.get_file("barcoded-subreadset")
        ]
        make_mock_laa_inputs(INPUT_FILES[0], INPUT_FILES[1])
        zip_out = tempfile.NamedTemporaryFile(suffix=".zip").name
        rc = make_combined_laa_zip(INPUT_FILES[0], INPUT_FILES[1],
                                   INPUT_FILES[2], zip_out)
        self.assertEqual(rc, 0)
        self.assertTrue(op.getsize(zip_out) != 0)
        with ZipFile(zip_out, "r") as zip_file:
            file_names = set(zip_file.namelist())
            self.assertTrue("Barcoded_Sample_Names.csv" in file_names)
            self.assertTrue("consensus_sequence_statistics.csv" in file_names)


class TestBarcodeUtils(unittest.TestCase):

    def test_discard_bio_samples(self):
        ds = SubreadSet(pbtestdata.get_file("barcoded-subreadset"))
        discard_bio_samples(ds, "lbc1--lbc1")
        coll = ds.metadata.collections[0]
        self.assertEqual(len(coll.wellSample.bioSamples), 1)
        self.assertEqual(coll.wellSample.bioSamples[0].name, "Alice")
        # No matching BioSample records
        ds = SubreadSet(pbtestdata.get_file("barcoded-subreadset"))
        coll = ds.metadata.collections[0]
        coll.wellSample.bioSamples.pop(1)
        coll.wellSample.bioSamples.pop(1)
        bioSample = coll.wellSample.bioSamples[0]
        while len(bioSample.DNABarcodes) > 0:
            bioSample.DNABarcodes.pop(0)
        self.assertEqual(len(coll.wellSample.bioSamples), 1)
        discard_bio_samples(ds, "lbc1--lbc1")
        self.assertEqual(len(coll.wellSample.bioSamples), 1)
        self.assertEqual(coll.wellSample.bioSamples[0].name, "lbc1--lbc1")
        self.assertEqual(coll.wellSample.bioSamples[
                         0].DNABarcodes[0].name, "lbc1--lbc1")
        # no BioSample records
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))
        coll = ds.metadata.collections[0]
        self.assertEqual(len(coll.wellSample.bioSamples), 0)
        discard_bio_samples(ds, "lbc1--lbc1")
        self.assertEqual(len(coll.wellSample.bioSamples), 1)
        self.assertEqual(coll.wellSample.bioSamples[0].name, "lbc1--lbc1")
        self.assertEqual(coll.wellSample.bioSamples[
                         0].DNABarcodes[0].name, "lbc1--lbc1")

    def test_get_ds_name(self):
        ds = SubreadSet(pbtestdata.get_file("barcoded-subreadset"))
        name = get_ds_name(ds, "My Data", "My Barcode")
        self.assertEqual(name, "My Data (multiple samples)")
        for coll in ds.metadata.collections:
            while len(coll.wellSample.bioSamples) > 0:
                coll.wellSample.bioSamples.pop(0)
        name = get_ds_name(ds, "My Data", "My Barcode")
        self.assertEqual(name, "My Data (My Barcode)")
        ds = SubreadSet(pbtestdata.get_file("barcoded-subreadset"))
        for coll in ds.metadata.collections:
            while len(coll.wellSample.bioSamples) > 1:
                coll.wellSample.bioSamples.pop(1)
        name = get_ds_name(ds, "My Data", "My Barcode")
        expected = "My Data ({s})".format(
            s=ds.metadata.collections[0].wellSample.bioSamples[0].name)
        self.assertEqual(name, expected)
        ds = SubreadSet(ds.externalResources[0].bam)
        name = get_ds_name(ds, "My Data", "My Barcode")
        self.assertEqual(name, "My Data (My Barcode)")
        name = get_ds_name(ds, "My Data", None)
        self.assertEqual(name, "My Data (unknown sample)")

    def test_update_barcoded_sample_metadata(self):
        INPUT_FILES = [
            tempfile.NamedTemporaryFile(suffix=".datastore.json").name,
            pbtestdata.get_file("barcoded-subreadset"),
            pbtestdata.get_file("barcodeset")
        ]
        ds = split_barcoded_dataset(INPUT_FILES[1])
        ds.write_json(INPUT_FILES[0])
        base_dir = tempfile.mkdtemp()
        datastore = update_barcoded_sample_metadata(base_dir,
                                                    INPUT_FILES[0],
                                                    INPUT_FILES[1],
                                                    INPUT_FILES[2])
        validate_barcoded_datastore_files(self, INPUT_FILES[1], datastore)
