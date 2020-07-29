import tempfile

from pbcommand.testkit import PbIntegrationBase
from pbcore.io import openDataSet

from pbcoretools.tasks.memory.estimate_lima_memory import estimate_lima_memory
from pbcoretools.tasks.memory.get_dataset_size import get_dataset_size

import pbtestdata


class TestEstimateLimaMemory(PbIntegrationBase):
    TINY_DATA = pbtestdata.get_file("subreads-sequel")
    TINY_BARCODES = pbtestdata.get_file("barcodeset")
    BIG_BARCODES = "/pbi/dept/secondary/siv/barcodes/Sequel_RSII_384_barcodes_v1/Sequel_RSII_384_barcodes_v1.barcodeset.xml"
    BIG_DATA = "/pbi/dept/secondary/siv/testdata/Spider/all4mers/rSPOC1_20180629_223342/1_A01/mSPOC1_180629_223410.subreadset.xml"
    CCS_DATA = "/pbi/dept/secondary/siv/testdata/SA3-Sequel/bcol/m54119_161211_175055.consensusreadset.xml"

    def test_estimate_lima_memory(self):
        mem_gb = estimate_lima_memory(self.TINY_BARCODES, self.TINY_DATA, True)
        assert mem_gb == 2
        # this is silly of course.  but it's technically possible with the
        # Sequel II system, so we might as well just deal with it
        mem_gb = estimate_lima_memory(self.BIG_BARCODES, self.BIG_DATA, False)
        assert mem_gb == 551
        # this is a more realistic case - 147K barcode pairs but the BAM file
        # is small enough to fit in the default footprint
        mem_gb = estimate_lima_memory(self.BIG_BARCODES, self.CCS_DATA, False)
        assert mem_gb == 2

    def test_integration_tiny(self):
        args = [
            "python3", "-m", "pbcoretools.tasks.memory.estimate_lima_memory",
            self.TINY_BARCODES, self.TINY_DATA, "--symmetric"
        ]
        self._check_call(args)
        with open("lima_mem_gb.txt") as txt_out:
            assert txt_out.read() == "2"

    def test_integration_big(self):
        args = [
            "python3", "-m", "pbcoretools.tasks.memory.estimate_lima_memory",
            self.BIG_BARCODES, self.BIG_DATA, "--asymmetric"
        ]
        self._check_call(args)
        with open("lima_mem_gb.txt") as txt_out:
            assert txt_out.read() == "551"

    def test_defined_biosamples(self):
        # XXX awful dependency but it makes testing easier
        from pbcoretools.file_utils import set_bio_samples
        ds_tmp = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        bc = openDataSet(self.BIG_BARCODES)
        with openDataSet(self.BIG_DATA, trustCounts=True) as ds:
            bcs = [("bc1001--bc1{:03d}".format(x), "Sample {}".format(x))
                   for x in range(384)]
            set_bio_samples(ds, bcs)
            ds.write(ds_tmp)
        mem_gb = estimate_lima_memory(self.BIG_BARCODES, ds_tmp, False)
        assert mem_gb == 2


class TestGetDatasetSize(PbIntegrationBase):
    TINY_DATA = pbtestdata.get_file("subreads-sequel")
    BIG_DATA = TestEstimateLimaMemory.BIG_DATA
    TINY_REF = pbtestdata.get_file("lambdaNEB")
    BIG_REF = "/pbi/dept/secondary/siv/references/human_hs37d5/referenceset.xml"

    def test_get_dataset_size(self):
        tiny_xml = pbtestdata.get_file("subreads-sequel")
        m = get_dataset_size(tiny_xml, True, True)
        assert m.numRecords == 20
        assert m.totalLengthMb == 1
        assert m.indexSizeGb == 2
        m = get_dataset_size(tiny_xml, False, False)
        assert m.numRecords == 20
        assert m.totalLengthMb == 1
        assert m.indexSizeGb == 1
        m = get_dataset_size(self.BIG_DATA, True, True)
        assert m.numRecords == 805580876
        assert m.totalLengthMb == 271330
        assert m.indexSizeGb == 45
        m = get_dataset_size(self.TINY_REF, False, False)
        assert m.numRecords == 1
        assert m.totalLengthMb == 1
        m = get_dataset_size(self.BIG_REF, False, False)
        assert m.numRecords == 86
        assert m.totalLengthMb == 2993

    def _verify_outputs(self, numRecords, totalLengthMb, indexSizeGb):
        vals = [numRecords, totalLengthMb, indexSizeGb]
        files = ["numrecords.txt", "totallength.txt", "indexsize.txt"]
        for fn, val in zip(files, vals):
            with open(fn, "rt") as f:
                assert f.read() == str(val)

    def test_integration_tiny(self):
        args = [
            "python3", "-m", "pbcoretools.tasks.memory.get_dataset_size",
            pbtestdata.get_file("subreads-sequel"),
            "--skip-counts", "--get-index-size"
        ]
        self._check_call(args)
        self._verify_outputs(20, 1, 2)

    def test_integration_big(self):
        args = [
            "python3", "-m", "pbcoretools.tasks.memory.get_dataset_size",
            self.BIG_DATA,
            "--skip-counts", "--get-index-size"
        ]
        self._check_call(args)
        self._verify_outputs(805580876, 271330, 45)

    def test_integration_ref(self):
        args = [
            "python3", "-m", "pbcoretools.tasks.memory.get_dataset_size",
            self.TINY_REF
        ]
        self._check_call(args)
        self._verify_outputs(1, 1, 1)
