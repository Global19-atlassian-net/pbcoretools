
# FIXME too much copy-paste from test_tasks_scatter_gather.py

from zipfile import ZipFile
import subprocess
import tempfile
import unittest
import json
import os.path as op

from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter

from test_tasks_scatter_gather import MOCK_GFF_RECORDS, MOCK_VCF_RECORDS, MOCK_VCF_HEADER
from test_chunking_gather import create_zip


class GatherTextRecordsBase(object):
    RECORDS = []
    RECORD_HEADER = None
    EXTENSION = None

    @classmethod
    def setUpClass(cls):
        cls.INPUT_FILES = []
        base = tempfile.mkdtemp()
        for i in range(2):
            file_name = "%s.%d%s" % (base, i + 1, cls.EXTENSION)
            with open(file_name, 'w') as f:
                if cls.RECORD_HEADER is not None:
                    f.write(cls.RECORD_HEADER)
                f.write("\n".join(cls._get_chunk_records(i)))
                f.write("\n")  # XXX we need this for CSV gather
                cls.INPUT_FILES.append(file_name)

    @classmethod
    def _get_chunk_records(cls, i_chunk):
        return cls.RECORDS[i_chunk * 2:(i_chunk + 1) * 2]

    def _validate_result(self, gathered_file):
        base, ext = op.splitext(gathered_file)
        self.assertEqual(ext, "%s" % self.EXTENSION)
        with open(gathered_file) as f:
            lines_ = f.readlines()
            lines = self._get_lines(lines_)
            self.assertEqual(lines, self.RECORDS)
            self.validate_content(lines_)

    def validate_content(self, lines):
        pass

    def test_run_tool_and_validate(self):
        tmp_out = tempfile.NamedTemporaryFile().name + self.EXTENSION
        args = ["pbtools-gather", tmp_out] + self.INPUT_FILES
        subprocess.check_call(args)
        self._validate_result(tmp_out)


class TestGatherToolGff(GatherTextRecordsBase, unittest.TestCase):
    RECORDS = MOCK_GFF_RECORDS
    RECORD_HEADER = "##gff-version 3\n##source-id ipdSummary\n"
    EXTENSION = ".gff"

    @classmethod
    def _get_chunk_records(cls, i_chunk):
        if i_chunk == 0: return cls.RECORDS[2:]
        else: return cls.RECORDS[0:2]

    def _get_lines(self, lines):
        return [l.strip() for l in lines if l[0] != '#']

    def validate_content(self, lines):
        self.assertEqual(len(lines), 6)
        self.assertEqual(lines[1].strip(), "##source-id ipdSummary")


class TestGatherToolVcf(GatherTextRecordsBase, unittest.TestCase):
    RECORDS = MOCK_VCF_RECORDS
    RECORD_HEADER = MOCK_VCF_HEADER
    EXTENSION = ".vcf"

    @classmethod
    def _get_chunk_records(cls, i_chunk):
        if i_chunk == 0: return cls.RECORDS[2:]
        else: return cls.RECORDS[0:2]

    def _get_lines(self, lines):
        return [l.strip() for l in lines if l[0] != '#']

    def validate_content(self, lines):
        self.assertEqual(len(lines), 9)
        self.assertEqual(lines[3].strip(), "##reference=ecoliK12_pbi_March2013.fasta")


class TestGatherToolFasta(unittest.TestCase):
    CHUNK_CONTIGS = [
        ("lambda_NEB3011_0_30", "GGGCGGCGACCTCGCGGGTTTTCGCTATTT"),
        ("lambda_NEB3011_60_90", "CACTGAATCATGGCTTTATGACGTAACATC"),
        ("lambda_NEB3011_30_60", "GTGGACTCGGAGCAGTTCGGCAGCCAGCAG")
    ]
    EXTENSION = ".fasta"
    READER = FastaReader
    EXTRA_ARGS = []

    @classmethod
    def setUpClass(cls):
        cls.INPUT_FILES = []
        for header, seq in cls.CHUNK_CONTIGS:
            cls.INPUT_FILES.append(cls._write_fastx_file(header, seq))

    @classmethod
    def _write_fastx_file(cls, header, seq):
        fn = tempfile.NamedTemporaryFile(suffix=".fasta").name
        suffix = "|arrow"
        with FastaWriter(fn) as f:
            f.writeRecord("{h}{s}".format(h=header, s=suffix), seq)
        return fn

    def test_run_tool_and_validate(self):
        tmp_out = tempfile.NamedTemporaryFile().name + self.EXTENSION
        args = ["pbtools-gather", tmp_out] + self.INPUT_FILES + self.EXTRA_ARGS
        subprocess.check_call(args)
        self._validate_result(tmp_out)

    def _validate_result(self, gathered_file):
        with self.READER(gathered_file) as fastx_out:
            records = [rec for rec in fastx_out]
            self.assertEqual(len(records), 3)
            for (header, seq), rec2 in zip(self.CHUNK_CONTIGS, records):
                self.assertEqual(header+"|arrow", rec2.header)
                self.assertEqual(seq, rec2.sequence)


class TestGatherToolFastq(TestGatherToolFasta):
    EXTENSION = ".fastq"
    READER = FastqReader

    @classmethod
    def _write_fastx_file(cls, header, seq):
        fn = tempfile.NamedTemporaryFile(suffix=".fastq").name
        suffix = "|arrow"
        with FastqWriter(fn) as f:
            f.writeRecord("{h}{s}".format(h=header, s=suffix), seq, [35]*len(seq))
        return fn


class TestGatherToolFastaJoinContigs(TestGatherToolFasta):
    EXTRA_ARGS = ["--join-contigs"]

    def _validate_result(self, gathered_file):
        with self.READER(gathered_file) as fastx_out:
            records = [rec for rec in fastx_out]
            self.assertEqual(len(records), 1)
            combined_seq = "".join([x[1] for x in self.CHUNK_CONTIGS])
            self.assertEqual(len(combined_seq), len(records[0].sequence))


class TestGatherToolFastqJoinContigs(TestGatherToolFastaJoinContigs):
    EXTENSION = ".fastq"
    READER = FastqReader

    @classmethod
    def _write_fastx_file(cls, header, seq):
        fn = tempfile.NamedTemporaryFile(suffix=".fastq").name
        suffix = "|arrow"
        with FastqWriter(fn) as f:
            f.writeRecord("{h}{s}".format(h=header, s=suffix), seq, [35]*len(seq))
        return fn


class TestGatherToolZip(unittest.TestCase):

    def test_run_tool_and_validate(self):
        inputs = []
        for i in range(2):
            fn = tempfile.NamedTemporaryFile(suffix=".zip").name
            create_zip(fn)
            inputs.append(fn)
        tmp_out = tempfile.NamedTemporaryFile(suffix=".zip").name
        args = ["pbtools-gather", tmp_out] + inputs
        subprocess.check_call(args)
        uuids = set()
        with ZipFile(tmp_out, "r") as gathered:
            for member in gathered.namelist():
                d = json.loads(gathered.open(member).read())
                uuids.add(d["uuid"])
        self.assertEqual(len(uuids), 4)
