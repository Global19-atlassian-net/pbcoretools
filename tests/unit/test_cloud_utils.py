
import os.path as op
import os

from pbcoretools.cloud_utils import get_zmw_bgzf_borders, get_bam_offsets, split_bam, extract_bam_chunk, combine_with_header

from pbcommand.testkit import PbIntegrationBase
from pbcore.io import openDataSet, BamReader, IndexedBamReader, PacBioBamIndex
import pbtestdata


class TestCloudUtils(PbIntegrationBase):
    DS1 = pbtestdata.get_file("subreads-xml")
    DS2 = pbtestdata.get_file("subreads-sequel")

    def _get_bam_path(self, ds_path):
        with openDataSet(ds_path) as ds:
            return ds.resourceReaders()[0].filename

    def _remove_all(self):
        for file_name in os.listdir(os.getcwd()):
            if file_name.startswith("reads.chunk") and file_name.endswith(".bam"):
                os.remove(op.join(os.getcwd(), file_name))

    def test_split_bam(self):
        bam_file1 = self._get_bam_path(self.DS1)
        CHUNKS_IN = [1, 2, 3, 4]
        CHUNKS_OUT = [1, 2, 3, 3]
        for n_in, n_expected in zip(CHUNKS_IN, CHUNKS_OUT):
            nchunks = split_bam(bam_file1, n_in)
            assert nchunks == n_expected
            bam_in = IndexedBamReader(bam_file1)
            records_in = [rec.qName for rec in bam_in]
            records_out = []
            for i in range(n_expected):
                bam_out = BamReader("reads.chunk%d.bam" % i)
                records_out.extend([rec.qName for rec in bam_out])
            assert records_in == records_out
            self._remove_all()

    def test_get_zmw_bgzf_borders(self):
        bam_file = self._get_bam_path(self.DS1)
        pbi_file = bam_file + ".pbi"
        pbi = PacBioBamIndex(pbi_file)
        offsets = get_zmw_bgzf_borders(pbi)
        assert offsets == [(0, 1650, 396), (16, 7247, 26575), (48, 30983, 77209)]
        bam_file = self._get_bam_path(self.DS2)
        pbi_file = bam_file + ".pbi"
        pbi = PacBioBamIndex(pbi_file)
        offsets = get_zmw_bgzf_borders(pbi)
        assert offsets == [(0, 5177614, 447)]

    def test_get_bam_offsets(self):
        bam_file = self._get_bam_path(self.DS1)
        offsets = get_bam_offsets(bam_file, 4)
        assert offsets == [396, 26575, 77209]
        offsets = get_bam_offsets(bam_file, 3)
        assert offsets == [396, 26575, 77209]
        offsets = get_bam_offsets(bam_file, 2)
        assert offsets == [396, 77209]
        offsets = get_bam_offsets(bam_file, 1)
        assert offsets == [396]

    def test_combine_with_header(self):
        bam_file = self._get_bam_path(self.DS1)
        bam_size = op.getsize(bam_file)
        # see above - these are known boundaries for this particular input
        byte_ranges = [(396, 26575), (26575, 77209), (77209, bam_size)]
        with open(bam_file, "rb") as bam_in:
            with open("header.bam", "wb") as header_out:
                header_out.write(bam_in.read(396))
            for i, (start, end) in enumerate(byte_ranges):
                with open("tmp.chunk%d.bam" % i, "wb") as chunk_out:
                    bam_in.seek(start)
                    nbytes = end - start
                    chunk_out.write(bam_in.read(nbytes))
        for i in range(3):
            combine_with_header("header.bam", "tmp.chunk%d.bam" % i, "combined.chunk%d.bam" % i)
        bam_in = IndexedBamReader(bam_file)
        records_in = [rec.qName for rec in bam_in]
        records_out = []
        for i in range(3):
            bam_out = BamReader("combined.chunk%d.bam" % i)
            records_out.extend([rec.qName for rec in bam_out])
        assert records_in == records_out
