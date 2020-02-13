import tempfile
import sys

from pbcore.io.FastaIO import FastaRecord, FastaWriter
from pbcore.io import FastqWriter, ContigSet

from pbcoretools.chunking.chunk_utils import (write_pbcore_records,
                                              write_contigset_records,
                                              to_chunked_fasta_files,
                                              to_chunked_fastq_files,
                                              guess_optimal_max_nchunks_for_consensus)

import mock


class TestChunkUtils:

    def test_write_pbcore_records(self):
        records = [FastaRecord("chr1", "acgt"), FastaRecord("chr2", "tgca")]
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
        write_pbcore_records(FastaWriter, records, tmp_fasta)
        with open(tmp_fasta) as fasta_in:
            lines = fasta_in.read().splitlines()
            assert lines == [">chr1", "acgt", ">chr2", "tgca"]

    def test_write_contigset_records(self):
        records = [FastaRecord("chr1", "acgt"), FastaRecord("chr2", "tgca")]
        tmp_contigs = tempfile.NamedTemporaryFile(suffix=".contigset.xml").name
        write_contigset_records(FastaWriter, records, tmp_contigs)
        with ContigSet(tmp_contigs) as ds_in:
            rec2 = [(rec.id, rec.sequence) for rec in ds_in]
            assert rec2 == [("chr1", "acgt"), ("chr2", "tgca")]

    def test_to_chunked_fasta_files(self):
        records = mock.to_fasta_records(10)
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
        write_pbcore_records(FastaWriter, records, tmp_fasta)
        tmp_dir = tempfile.mkdtemp()
        chunks = list(to_chunked_fasta_files(
            tmp_fasta, 5, tmp_dir, "fasta_chunk", ".fasta"))
        assert len(chunks) == 5

    def test_to_chunked_fastq_files(self):
        records = mock.to_fastq_records(10)
        tmp_fastq = tempfile.NamedTemporaryFile(suffix=".fastq").name
        write_pbcore_records(FastqWriter, records, tmp_fastq)
        tmp_dir = tempfile.mkdtemp()
        chunks = list(to_chunked_fastq_files(
            tmp_fastq, 5, tmp_dir, "fastq_chunk", ".fastq"))
        assert len(chunks) == 5

    def test_guess_optimal_max_nchunks_for_consensus(self):
        assert guess_optimal_max_nchunks_for_consensus(400000) == 12
        assert guess_optimal_max_nchunks_for_consensus(4000000) == 19
        assert guess_optimal_max_nchunks_for_consensus(40000000) == 51
        assert guess_optimal_max_nchunks_for_consensus(400000000) == 83
        assert guess_optimal_max_nchunks_for_consensus(sys.maxsize) == 96
        assert guess_optimal_max_nchunks_for_consensus(400000000, 24) == 24
