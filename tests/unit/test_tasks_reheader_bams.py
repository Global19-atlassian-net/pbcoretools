"""
Unit and integration tests for pbcoretools.tasks.reheader_bams
"""

import random
import os.path as op
import os

from pbcommand.testkit import PbIntegrationBase
from pbcore.io import IndexedBamReader, openDataSet

from pbcoretools.tasks.reheader_bams import reheader_bam, reheader_dataset_bams

import pbtestdata

class TestReheaderBams(PbIntegrationBase):
    BIOSAMPLE_NAME = "mock-sm-{u:05d}".format(u=random.randint(1, 99999))
    LIBRARY_NAME = "mock-lb-{u:05d}".format(u=random.randint(1, 99999))

    def _validate_bam(self, bam_reader):
        for rg in bam_reader.readGroupTable:
            assert rg.SampleName == self.BIOSAMPLE_NAME
            assert rg.LibraryName == self.LIBRARY_NAME

    # this guards against overwriting the input dataset BAMs by mistake
    def _validate_input_bam(self, bam_reader):
        for rg in bam_reader.readGroupTable:
            assert rg.SampleName != self.BIOSAMPLE_NAME
            assert rg.LibraryName != self.LIBRARY_NAME

    def _validate_records(self, reader_in, reader_out):
        assert len(reader_in) == len(reader_out)
        for rec_in, rec_out in zip(reader_in, reader_out):
            assert rec_in.qName == rec_out.qName

    def test_reheader_bam(self):
        ofn = "subreads_out.bam"
        bam_file = pbtestdata.get_file("subreads-bam")
        reheader_bam(bam_file, ofn, self.BIOSAMPLE_NAME, self.LIBRARY_NAME)
        assert op.isfile(ofn) and op.isfile(ofn + ".pbi")
        with IndexedBamReader(ofn) as bam_out:
            self._validate_bam(bam_out)
            with IndexedBamReader(bam_file) as bam_in:
                self._validate_input_bam(bam_in)
                self._validate_records(bam_in, bam_out)

    def _validate_dataset(self, ds_out):
        for bam_res in ds_out.resourceReaders():
            self._validate_bam(bam_res)

    def _validate_input_dataset(self, ds_in):
        for bam_res in ds_in.resourceReaders():
            self._validate_input_bam(bam_res)

    def _run_reheader_dataset_bams(self, ds_file):
        with openDataSet(ds_file) as ds:
            ds_out = reheader_dataset_bams(ds,
                                           os.getcwd(),
                                           self.BIOSAMPLE_NAME,
                                           self.LIBRARY_NAME)
            self._validate_dataset(ds_out)
            self._validate_records(ds, ds_out)

    def test_reheader_dataset_bams_subreads(self):
        self._run_reheader_dataset_bams(pbtestdata.get_file("subreads-xml"))

    def test_reheader_dataset_bams_ccs(self):
        self._run_reheader_dataset_bams(pbtestdata.get_file("ccs-sequel"))

    def test_reheader_dataset_bams_ccs_barcoded(self):
        self._run_reheader_dataset_bams(pbtestdata.get_file("ccs-barcoded"))

    def _run_cli(self, ds_file, ds_out_file):
        args = [
            "python", "-m", "pbcoretools.tasks.reheader_bams",
            ds_file, ds_out_file,
            "--biosample-name", self.BIOSAMPLE_NAME,
            "--library-name", self.LIBRARY_NAME
        ]
        self._check_call(args)
        assert op.isfile(ds_out_file)
        with openDataSet(ds_out_file) as ds_out:
            self._validate_dataset(ds_out)
            with openDataSet(ds_file) as ds_in:
                self._validate_input_dataset(ds_in)
                self._validate_records(ds_in, ds_out)

    def test_integration_subreads(self):
        ds_out = "tmp_out.subreadset.xml"
        self._run_cli(pbtestdata.get_file("subreads-sequel"), ds_out)

    def test_integration_ccs(self):
        ds_out = "tmp_out.consensusreadset.xml"
        self._run_cli(pbtestdata.get_file("ccs-sequel"), ds_out)
