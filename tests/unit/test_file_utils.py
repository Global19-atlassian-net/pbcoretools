# -*- coding: utf-8 -*-

"""
Unit tests for pbcoretools.file_utils
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
from pbcore.io.dataset.DataSetMembers import DNABarcode
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
    mock_update_barcoded_sample_metadata,
    discard_bio_samples,
    set_bio_samples,
    parse_biosamples_csv,
    add_mock_collection_metadata,
    force_set_all_well_sample_names,
    force_set_all_bio_sample_names,
    sanitize_dataset_tags)


def _validate_dataset_xml(file_name):
    if "PB_DATASET_XSD" in os.environ:
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
    from pbcoretools.bamsieve import filter_reads
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
                                      "barcoding.tasks.lima-0",
                                      ds_in.datasetType,
                                      ds_tmp))
    return DataStore(ds_files)


def validate_barcoded_datastore_files(self, subreads, datastore,
                                      have_collection_metadata=True,
                                      use_barcode_uuids=True,
                                      number_of_expected_filters=1):
    """
    This is linked to the PacBioTestData file 'barcoded-subreadset', which
    has been manually edited in support of this test.
    """
    bio_sample_names = {
        "lbc1--lbc1": "Alice",
        "lbc3--lbc3": "Charles"
    }
    dna_bc_uuids = {
        "lbc1--lbc1": "dffb30e8-9243-4743-9980-468a20952167",
        "lbc3--lbc3": "eef1a8ea-c6a7-4233-982a-d426e1e7d8c9"
    }
    assert len(datastore.files) == 2
    ds_in = SubreadSet(subreads)
    for f in datastore.files.values():
        _validate_dataset_xml(f.path)
        with SubreadSet(f.path) as ds:
            # FIXME need better testing here
            assert len(ds.filters) == number_of_expected_filters
            assert ds.uuid == f.uuid
            bc_label = op.basename(f.path).split(".")[1]
            bio_name = bio_sample_names[bc_label]
            expected_tags = set(
                ["TotalLength", "NumRecords", "Provenance", "BioSamples"])
            assert len(ds.metadata.bioSamples) == 1
            assert ds.metadata.bioSamples[0].name == bio_name
            assert ds.metadata.provenance.parentDataSet.uniqueId == ds_in.uuid
            assert ds.name == "{n} ({s})".format(n=ds_in.name, s=bio_name)
            if use_barcode_uuids:
                assert ds.uuid == dna_bc_uuids[bc_label]
            else:
                assert ds.uuid != dna_bc_uuids[bc_label]
            md_tags = [r['tag'] for r in ds.metadata.record['children']]
            assert set(md_tags).issuperset(expected_tags)


class TestSplitLAA:
    """
    Unit tests for LAA FASTQ splitter.
    """

    SUBREADS = pbtestdata.get_file("barcoded-subreadset")

    def setup_method(self, method):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self._records = _get_fastq_records()
        self.input_file_name = make_fastq_inputs(self._records)
        self._subreads = pbtestdata.get_file("barcoded-subreadset")

    def test_split_laa_fastq(self):
        ifn = self.input_file_name
        ofb = tempfile.NamedTemporaryFile().name
        ofs = split_laa_fastq(ifn, ofb, self._subreads)
        assert len(ofs) == 2
        suffixes = sorted([".".join(of.split('.')[1:]) for of in ofs])
        assert suffixes == ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq']
        for i, ofn in enumerate(ofs):
            with FastqReader(ofn) as fastq_in:
                recs = [rec for rec in fastq_in]
                for j in range(2):
                    assert str(recs[j]) == str(
                        self._records[(i * 2) + j])

    def test_split_laa_fastq_archived(self):
        ifn = self.input_file_name
        ofn = tempfile.NamedTemporaryFile(suffix=".zip").name
        rc = split_laa_fastq_archived(ifn, ofn, self._subreads)
        assert rc == 0
        with ZipFile(ofn, "r") as zip_in:
            files = zip_in.namelist()
            suffixes = sorted([".".join(of.split('.')[1:]) for of in files])
            assert suffixes == ['Alice.lbc1--lbc1.fastq', 'Charles.lbc3--lbc3.fastq']

    def test_get_barcode_sample_mappings(self):
        with SubreadSet(self._subreads) as ds:
            # just double-checking that the XML defines more samples than are
            # actually present in the BAM
            assert len(ds.metadata.bioSamples) == 3
        samples = get_barcode_sample_mappings(SubreadSet(self._subreads))
        assert samples == {
            'lbc3--lbc3': 'Charles',
            'lbc1--lbc1': 'Alice'}

    def test_make_barcode_sample_csv(self):
        csv_file = tempfile.NamedTemporaryFile(suffix=".csv").name
        sample_mappings = make_barcode_sample_csv(self._subreads, csv_file)
        assert sample_mappings == {
            'lbc3--lbc3': 'Charles',
            'lbc1--lbc1': 'Alice'}
        with open(csv_file) as f:
            assert f.read() == "Barcode Name,Bio Sample Name\nlbc1--lbc1,Alice\nlbc3--lbc3,Charles\n"

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
        assert rc == 0
        assert op.getsize(zip_out) != 0
        with ZipFile(zip_out, "r") as zip_file:
            file_names = set(zip_file.namelist())
            assert "Barcoded_Sample_Names.csv" in file_names
            assert "consensus_sequence_statistics.csv" in file_names


class TestBarcodeUtils:
    SUBREADS = pbtestdata.get_file("barcoded-subreadset")

    def test_discard_bio_samples(self):
        ds = SubreadSet(self.SUBREADS)
        discard_bio_samples(ds, "lbc1--lbc1")
        coll = ds.metadata.collections[0]
        assert len(ds.metadata.bioSamples) == 1
        assert ds.metadata.bioSamples[0].name == "Alice"
        # No matching BioSample records
        ds = SubreadSet(self.SUBREADS)
        ds.metadata.bioSamples.pop(1)
        ds.metadata.bioSamples.pop(1)
        bioSample = ds.metadata.bioSamples[0]
        while len(bioSample.DNABarcodes) > 0:
            bioSample.DNABarcodes.pop(0)
        assert len(ds.metadata.bioSamples) == 1
        discard_bio_samples(ds, "lbc1--lbc1")
        assert len(ds.metadata.bioSamples) == 1
        assert ds.metadata.bioSamples[0].name == "lbc1--lbc1"
        assert ds.metadata.bioSamples[0].DNABarcodes[0].name == "lbc1--lbc1"
        # no BioSample records
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))
        assert len(ds.metadata.bioSamples) == 0
        discard_bio_samples(ds, "lbc1--lbc1")
        assert len(ds.metadata.bioSamples) == 1
        assert ds.metadata.bioSamples[0].name == "lbc1--lbc1"
        assert ds.metadata.bioSamples[0].DNABarcodes[0].name == "lbc1--lbc1"

    def test_get_ds_name(self):
        ds = SubreadSet(self.SUBREADS)
        name = get_ds_name(ds, "My Data", "My Barcode")
        assert name == "My Data (multiple samples)"
        while len(ds.metadata.bioSamples) > 0:
            ds.metadata.bioSamples.pop(0)
        name = get_ds_name(ds, "My Data", "My Barcode")
        assert name == "My Data (My Barcode)"
        ds = SubreadSet(self.SUBREADS)
        while len(ds.metadata.bioSamples) > 1:
            ds.metadata.bioSamples.pop(1)
        name = get_ds_name(ds, "My Data", "My Barcode")
        expected = "My Data ({s})".format(s=ds.metadata.bioSamples[0].name)
        assert name == expected
        ds = SubreadSet(ds.externalResources[0].bam)
        name = get_ds_name(ds, "My Data", "My Barcode")
        assert name == "My Data (My Barcode)"
        name = get_ds_name(ds, "My Data", None)
        assert name == "My Data (unknown sample)"

    @pytest.mark.xmllint
    @pytest.mark.constools
    def test_update_barcoded_sample_metadata(self):
        datastore_tmp = tempfile.NamedTemporaryFile(
            suffix=".datastore.json").name
        barcodes = pbtestdata.get_file("barcodeset")
        ds = split_barcoded_dataset(self.SUBREADS)
        ds.write_json(datastore_tmp)
        base_dir = tempfile.mkdtemp()
        datastore = update_barcoded_sample_metadata(base_dir,
                                                    datastore_tmp,
                                                    self.SUBREADS,
                                                    barcodes)
        validate_barcoded_datastore_files(self, self.SUBREADS, datastore)
        # now with use_barcode_uuids=False
        datastore = update_barcoded_sample_metadata(base_dir,
                                                    datastore_tmp,
                                                    self.SUBREADS,
                                                    barcodes,
                                                    use_barcode_uuids=False)
        validate_barcoded_datastore_files(self, self.SUBREADS, datastore,
                                          use_barcode_uuids=False)
        # test that it works with no collection metadata
        ss = SubreadSet(self.SUBREADS)
        ss.metadata.collections = None
        ss_tmp = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        ss.write(ss_tmp)
        ds = split_barcoded_dataset(ss_tmp)
        ds.write_json(datastore_tmp)
        base_dir = tempfile.mkdtemp()
        datastore = update_barcoded_sample_metadata(base_dir,
                                                    datastore_tmp,
                                                    self.SUBREADS,
                                                    barcodes)
        validate_barcoded_datastore_files(self, self.SUBREADS, datastore,
                                          have_collection_metadata=False)

    @pytest.mark.xmllint
    def test_mock_update_barcoded_sample_metadata(self):
        tmp_dir = tempfile.mkdtemp()
        datastore_tmp = op.join(tmp_dir, "lima.datastore.json")
        barcodeset = pbtestdata.get_file("barcodeset")
        barcodes = ["lbc1--lbc1", "lbc3--lbc3"]
        files = [
            op.join(tmp_dir, "lima.lbc1--lbc1.subreadset.xml"),
            op.join(tmp_dir, "lima.lbc3--lbc3.subreadset.xml")
        ]
        uuids = [uuid.uuid4() for fn in files]
        # XXX these are hardcoded to match the actual barcoded test input
        bc_uuids = [
            "dffb30e8-9243-4743-9980-468a20952167",
            "eef1a8ea-c6a7-4233-982a-d426e1e7d8c9"
        ]
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))

        def _add_barcoded_sample(sn, bn, id_):
            ds.metadata.bioSamples.addSample(sn)
            ds.metadata.bioSamples[-1].DNABarcodes.addBarcode(bn)
            ds.metadata.bioSamples[-1].DNABarcodes[-1].uniqueId = id_
        _add_barcoded_sample("Alice", "lbc1--lbc1", bc_uuids[0])
        _add_barcoded_sample("Charles", "lbc3--lbc3", bc_uuids[1])
        tmp_ds = op.join(tmp_dir, "input.subreadset.xml")
        ds.write(tmp_ds)
        for fn, bc, dsid in zip(files, barcodes, uuids):
            ds = SubreadSet(tmp_ds)
            ds.uuid = str(dsid)
            ds.name = ds.name + " ({b})".format(b=bc)
            ds.write(fn)
        ds_files = [DataStoreFile(dsid, "barcoding.tasks.lima-0",
                                  FileTypes.DS_SUBREADS.file_type_id, fn) for (dsid, fn) in zip(uuids, files)]
        ds = DataStore(ds_files)
        ds.write_json(datastore_tmp)
        base_dir = tempfile.mkdtemp()
        datastore = mock_update_barcoded_sample_metadata(base_dir,
                                                         datastore_tmp,
                                                         tmp_ds,
                                                         barcodeset)
        validate_barcoded_datastore_files(self, tmp_ds, datastore,
                                          number_of_expected_filters=0)

    def test_add_mock_collection_metadata(self):
        bam = pbtestdata.get_file("subreads-bam")
        ds = SubreadSet(bam)
        assert len(ds.metadata.collections) == 0
        add_mock_collection_metadata(ds)
        assert len(ds.metadata.collections) == 1
        assert ds.metadata.collections[0].context == "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0"

    def test_force_set_all_well_sample_names(self):
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))
        force_set_all_well_sample_names(ds, "My Test Sample")
        assert ds.metadata.collections[0].wellSample.name == "My Test Sample"

    def test_force_set_all_bio_sample_names(self):
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))
        force_set_all_bio_sample_names(ds, "My Test BioSample")
        assert ds.metadata.bioSamples[0].name == "My Test BioSample"

    def test_sanitize_dataset_tags(self):
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))
        base_name = ds.name
        ds.name = ds.name + " (filtered) (CCS)"
        ds.tags = "subreads,hidden,testdata,filtered "
        sanitize_dataset_tags(ds)
        assert ds.name == base_name + " (CCS)"
        assert ds.tags == "hidden,subreads,testdata"
        sanitize_dataset_tags(ds, remove_hidden=True)
        assert ds.name == base_name + " (CCS)"
        assert ds.tags == "subreads,testdata"
        ds.tags = ", hidden, ccs"
        sanitize_dataset_tags(ds, remove_hidden=True)
        assert ds.tags == "ccs"

    def test_parse_biosamples_csv(self):
        csv = "Barcode,BioSample Name\nlbc1--lbc1,Alice\nlbc2--lbc2,Bob"
        csv_tmp = tempfile.NamedTemporaryFile(suffix=".csv").name
        with open(csv_tmp, "w") as csv_out:
            csv_out.write(csv)
        records = parse_biosamples_csv(csv_tmp)
        assert records == [
            ("lbc1--lbc1", "Alice"),
            ("lbc2--lbc2", "Bob")]

    def test_parse_biosamples_csv_non_ascii(self):
        csv = "Barcode,BioSample Name\nlbc1--lbc1,Alice™\nlbc2--lbc2,Bob¢"
        csv_tmp = tempfile.NamedTemporaryFile(suffix=".csv").name
        with open(csv_tmp, "w") as csv_out:
            csv_out.write(csv)
        with pytest.raises(ValueError):
            parse_biosamples_csv(csv_tmp)

    def test_set_bio_samples(self):
        ds = SubreadSet(pbtestdata.get_file("subreads-sequel"))
        samples = [("lbc1--lbc1", "Alice"), ("lbc2--lbc2", "Bob")]
        set_bio_samples(ds, samples)
        samples_out = {
            s.DNABarcodes[0].name: s.name for s in ds.metadata.bioSamples}
        assert samples_out == dict(samples)
