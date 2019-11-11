
import tempfile
import unittest
import sys

from pbcore.io.FastaIO import FastaRecord, FastaWriter
from pbcore.io import FastqWriter, ContigSet

from pbcoretools.chunking.chunk_utils import (write_pbcore_records,
                                              write_contigset_records,
                                              to_chunked_fasta_files,
                                              to_chunked_fastq_files,
                                              guess_optimal_max_nchunks_for_consensus)

import mock

class TestChunkUtils(unittest.TestCase):

    def test_write_pbcore_records(self):
        records = [FastaRecord("chr1", "acgt"), FastaRecord("chr2", "tgca")]
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
        write_pbcore_records(FastaWriter, records, tmp_fasta)
        with open(tmp_fasta) as fasta_in:
            lines = fasta_in.read().splitlines()
            self.assertEqual(lines, [">chr1", "acgt", ">chr2", "tgca"])

    def test_write_contigset_records(self):
        records = [FastaRecord("chr1", "acgt"), FastaRecord("chr2", "tgca")]
        tmp_contigs = tempfile.NamedTemporaryFile(suffix=".contigset.xml").name
        write_contigset_records(FastaWriter, records, tmp_contigs)
        with ContigSet(tmp_contigs) as ds_in:
            rec2 = [ (rec.id, rec.sequence) for rec in ds_in ]
            self.assertEqual(rec2, [("chr1", "acgt"), ("chr2", "tgca")])

    def test_to_chunked_fasta_files(self):
        records = mock.to_fasta_records(10)
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
        write_pbcore_records(FastaWriter, records, tmp_fasta)
        tmp_dir = tempfile.mkdtemp()
        chunks = list(to_chunked_fasta_files(tmp_fasta, 5, tmp_dir, "fasta_chunk", ".fasta"))
        self.assertEqual(len(chunks), 5)

    def test_to_chunked_fastq_files(self):
        records = mock.to_fastq_records(10)
        tmp_fastq = tempfile.NamedTemporaryFile(suffix=".fastq").name
        write_pbcore_records(FastqWriter, records, tmp_fastq)
        tmp_dir = tempfile.mkdtemp()
        chunks = list(to_chunked_fastq_files(tmp_fastq, 5, tmp_dir, "fastq_chunk", ".fastq"))
        self.assertEqual(len(chunks), 5)

    def test_guess_optimal_max_nchunks_for_consensus(self):
        self.assertEqual(guess_optimal_max_nchunks_for_consensus(400000), 12)
        self.assertEqual(guess_optimal_max_nchunks_for_consensus(4000000), 19)
        self.assertEqual(guess_optimal_max_nchunks_for_consensus(40000000), 51)
        self.assertEqual(guess_optimal_max_nchunks_for_consensus(400000000), 83)
        self.assertEqual(guess_optimal_max_nchunks_for_consensus(sys.maxsize), 96)
        self.assertEqual(guess_optimal_max_nchunks_for_consensus(400000000, 24), 24)
