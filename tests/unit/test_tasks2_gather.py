
from zipfile import ZipFile
import tempfile
import textwrap
import json
import os.path as op

from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter
from pbcommand.testkit import PbIntegrationBase

from test_chunking_gather import create_zip
from utils import skip_if_no_internal_data


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
        self._check_call(args)
        self._validate_result(tmp_out)


class TestGatherToolCsv(GatherTextRecordsBase, PbIntegrationBase):
    RECORDS = [
        "contig1,3000000,170",
        "contig2,90000,180",
        "contig3,58000,159",
        "contig4,20000,160",
    ]
    RECORD_HEADER = "contig_id,length,coverage\n"
    EXTENSION = ".csv"

    def _get_lines(self, lines):
        return [l.strip() for l in lines[1:]]


MOCK_GFF_RECORDS = [
    "contig1\tkinModCall\tmodified_base\t1\t1\t31\t+\t.\tcoverage=169",
    "contig1\tkinModCall\tmodified_base\t2\t2\t41\t-\t.\tcoverage=170",
    "contig1\tkinModCall\tmodified_base\t3\t3\t51\t+\t.\tcoverage=168",
    "contig1\tkinModCall\tmodified_base\t4\t4\t60\t-\t.\tcoverage=173",
]


class TestGatherToolGff(GatherTextRecordsBase, PbIntegrationBase):
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


MOCK_VCF_RECORDS = textwrap.dedent('''\
    ecoliK12_pbi_March2013 84 . TG T 48 PASS DP=53
    ecoliK12_pbi_March2013 218 . GA G 47 PASS DP=58
    ecoliK12_pbi_March2013 1536 . G GC 47 PASS DP=91
    ''').rstrip().replace(' ', '\t').split('\n')
MOCK_VCF_HEADER = textwrap.dedent('''\
    ##fileformat=VCFv4.3
    ##fileDate=20170328
    ##source=GenomicConsensusV2.2.0
    ##reference=ecoliK12_pbi_March2013.fasta
    ##contig=<ID=ecoliK12_pbi_March2013,length=4642522>
    #CHROM POS ID REF ALT QUAL FILTER INFO
    ''')

class TestGatherToolVcf(GatherTextRecordsBase, PbIntegrationBase):
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


class TestGatherToolFasta(PbIntegrationBase):
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
        self._check_call(args)
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


class TestGatherToolZip(PbIntegrationBase):

    def test_run_tool_and_validate(self):
        inputs = []
        for i in range(2):
            fn = tempfile.NamedTemporaryFile(suffix=".zip").name
            create_zip(fn)
            inputs.append(fn)
        tmp_out = tempfile.NamedTemporaryFile(suffix=".zip").name
        args = ["pbtools-gather", tmp_out] + inputs
        self._check_call(args)
        uuids = set()
        with ZipFile(tmp_out, "r") as gathered:
            for member in gathered.namelist():
                d = json.loads(gathered.open(member).read())
                uuids.add(d["uuid"])
        self.assertEqual(len(uuids), 4)


@skip_if_no_internal_data
def test_gather_datastore_json():
    import subprocess
    from pbcommand.models import DataStore
    d = '/pbi/dept/secondary/siv/testdata/pbsvtools-unittest/data/test_scatter_align_datastore/'
    if1 = op.join(d, '1.aln.datastore.json')
    if2 = op.join(d, '2.aln.datastore.json')
    of = tempfile.NamedTemporaryFile(suffix=".datastore.json").name
    args = ['python', '-m', 'pbcoretools.tasks2.gather', of, if1, if2]
    subprocess.check_call(args)
    out_fns = DataStore.load_from_json(of).to_dict()['files']
    expected_bam_1 = op.join(d, '1.bam')
    expected_bam_2 = op.join(d, '2.bam')
    assert out_fns[0]['path'] == expected_bam_1
    assert out_fns[1]['path'] == expected_bam_2


class TestGatherToolBed(PbIntegrationBase):
    def test_gather_bed(self):
        if1 = "test_gather_bed_1.bed"
        with open(if1, 'w') as writer:
            writer.write("#chr\tstart\tend\n")
            writer.write("1\t2\t3\n")
            writer.write("2\t3\t4\n")
            writer.write("")
        if2 = "test_gather_bed_2.bed"
        with open(if2, 'w') as writer:
            writer.write("#chr\tstart\tend\n")
            writer.write("1\t2\t3\n")
            writer.write("")
        of = "test_gather_bed_out.bed"
        args = [
            "python", "-m", "pbcoretools.tasks2.gather",
            of, if1, if2
        ]
        self._check_call(args)
        out = open(of, 'r').readlines()
        expected = ['#chr\tstart\tend\n', '1\t2\t3\n', '2\t3\t4\n', '1\t2\t3\n']
        self.assertEqual(out, expected)
