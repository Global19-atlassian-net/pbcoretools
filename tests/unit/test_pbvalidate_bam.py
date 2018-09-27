
"""
Unit tests for BAM validation, using embedded SAM-format string that is
manipulated to trigger various errors.
"""

from cStringIO import StringIO
import subprocess
import unittest
import logging
import os.path
import time
op = os.path
import sys

import pysam

from pbcoretools.pbvalidate.core import ValidateFile, ValidateRecord, ValidateFileObject
from pbcoretools.pbvalidate.bam import ValidateReadGroup
from pbcoretools.pbvalidate import bam
import pbcore.io

TESTDATA = "/pbi/dept/secondary/siv/testdata"

# of course PacBio's real reads are MUCH longer
rec1 = "movie1/54130/0_10\t2\tecoliK12_pbi_March2013_2955000_to_2980000\t2\t10\t%(cigar)s\t*\t0\t0\tAATGAGGAGA\t*\tRG:Z:%(rg_id)s\tdq:Z:2222'$22'2\tdt:Z:NNNNAGNNGN\tip:B:C,255,2,0,10,22,34,0,2,3,0,16\tiq:Z:(+#1'$#*1&\tmq:Z:&1~51*5&~2\tnp:i:1\tqe:i:10\tqs:i:%(qs)d\trq:f:%(rq)s\tsn:B:f,%(sn)s\tsq:Z:<32<4<<<<3\tzm:i:54130\tAS:i:-3020\tNM:i:134\tcx:i:2"
sam_str_ = """\
@HD\tVN:1.5\tSO:%(so)s\tpb:3.0.1
@SQ\tSN:ecoliK12_pbi_March2013_2955000_to_2980000\tLN:25000\tM5:734d5f3b2859595f4bd87a2fe6b7389b
@RG\tID:%(rg_id)s%(pl)s\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;%(ipd)s=ip;FRAMERATEHZ=75.0;BASECALLERVERSION=%(bcv)s;BINDINGKIT=%(bk)s;SEQUENCINGKIT=%(sk)s\tPU:movie1
@PG\tID:bax2bam-0.0.2\tPN:bax2bam\tVN:0.0.2\tDS:bax2bam converts the legacy PacBio basecall format (bax.h5) into the BAM basecall format.\tCL:bax2bam in.bax.h5 out.bam
movie1/54130/0_10\t2\tecoliK12_pbi_March2013_2955000_to_2980000\t2\t10\t%(cigar)s\t*\t0\t0\tAATGAGGAGA\t*\tRG:Z:%(rg_id)s\tdq:Z:2222'$22'2\tdt:Z:NNNNAGNNGN\tip:B:%(ip)s\tiq:Z:(+#1'$#*1&\tmq:Z:&1~51*5&~2\tnp:i:1\tqe:i:10\tqs:i:%(qs)d\trq:f:%(rq)s\tsn:B:f,%(sn)s\tsq:Z:<32<4<<<<3\tzm:i:54130\tAS:i:-3020\tNM:i:134\tcx:i:2
movie1/54130/10_20\t2\tecoliK12_pbi_March2013_2955000_to_2980000\t12\t10\t%(cigar2)s\t*\t0\t0\tAATGAGGAGA\t*\tRG:Z:%(rg_id)s\tdq:Z:2222'$22'2\tdt:Z:NNNNAGNNGN\tip:B:C,255,2,0,10,22,34,0,2,3,0,16\tiq:Z:(+#1'$#*1&\tmq:Z:&1~51*5&~2\tnp:i:1\tqe:i:20\tqs:i:10\trq:f:0.854\tsn:B:f,2.0,2.0,2.0,2.0\tsq:Z:<32<4<<<<3\tzm:i:54130\tAS:i:-3020\tNM:i:134\tcx:i:2
%(qname)s\t%(flag)s\tecoliK12_pbi_March2013_2955000_to_2980000\t22\t10\t%(cigar)s\t*\t0\t0\tAATGAGGAGA\t%(qual)s\tRG:Z:%(rg_id)s\tdq:Z:2222'$22'2\tdt:Z:%(dt)s\tip:B:C,255,2,0,10,22,34,0,2,3,0,16\tiq:Z:(+#1'$#*1&\tmq:Z:&1~51*5&~2\tnp:i:1\tqe:i:%(qe)d\tqs:i:20\trq:f:0.854\tsn:B:f,2.0,2.0,2.0,2.0\tsq:Z:%(sq)s\tzm:i:54130\tAS:i:-3020\tNM:i:134\tcx:i:2"""


basic_tags = {
    "bcv": "2.1",
    "bk": "100356300",
    "cigar": "10=",
    "cigar2": "10=",
    "dt": "NNNNAGNNGN",
    "flag": "2",
    "ipd": "Ipd:CodecV1",
    "ip": "C,255,2,0,10,22,34,0,2,3,0,16",
    "pl": "\tPL:PACBIO",
    "qe": 30,
    "qname": "movie1/54130/20_30",
    "qs": 0,
    "qual": "*",
    "rg_id": "3f58e5b8",
    "rq": "0.854",
    "sk": "100356200",
    "sn": "2.0,2.0,2.0,2.0",
    "so": "coordinate",
    "sq": "<32<4<<<<3",
}
sam_str_good = sam_str_ % basic_tags

bad_tags = {
    "bcv": "9.8.7.6.5",
    "bk": "null",
    "cigar": "4M2I2D4M",
    "cigar2": "4=2I2D4=",
    "dt": "012345689",
    "flag": "4",
    "ipd": "Ipd",
    "ip": "C,255,2,0,10,22,34,0,2,3,0,16",
    "pl": "",
    "qe": 29,
    "qname": "movie1_54130_10-20",
    "qs": 10,
    "qual": "9999999999",
    "rg_id": "2f48aec3",
    "rq": "2001",
    "sk": "null",
    "sn": "-1,-1,-1,-1",
    "so": "unsorted",
    "sq": "!!!!!!!!!!",
}

unmapped_sam_str_ = """\
@HD\tVN:1.5\tSO:%(so)s\tpb:3.0.1
@RG\tID:b5482b33\tPL:%(pl)s\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;%(ipd)s=ip;BINDINGKIT=100356300;SEQUENCINGKIT=100356200;FRAMERATEHZ=75.0;BASECALLERVERSION=%(bcv)s;FRAMERATEHZ=75.000000\tPU:m140906_231018_42161_c100676332550000001823129611271486_s1_p0
@PG\tID:bax2bam-0.0.2\tPN:bax2bam\tVN:0.0.2\tDS:bax2bam converts the legacy PacBio basecall format (bax.h5) into the BAM basecall format.\tCL:bax2bam in.bax.h5 out.bam
%(qname)s\t%(flag)s\t%(rname)s\t%(pos)s\t255\t*\t*\t0\t0\tAAAGAGAGAG\t*\tRG:Z:b5482b33\tdq:Z:2222222222\tdt:Z:NNNNNNNNNN\tip:B:C,255,9,20,43,38,12,9,30,39,22\tiq:Z:,*11111001\tmq:Z:&47088')34\tnp:i:1\tqe:i:10\tqs:i:0\trq:f:0.811\tsn:B:f,%(sn)s%(sq)s\tzm:i:%(zm)s\tcx:i:2
m140906_231018_42161_c100676332550000001823129611271486_s1_p0/%(zm2)s/0_10\t4\t*\t0\t255\t*\t*\t0\t0\tAAAGAGAGAG\t*\tRG:Z:b5482b33\tdq:Z:2222222222\tdt:Z:NNNNNNNNNN\tip:B:C,255,9,20,43,38,12,9,30,39,22\tiq:Z:,*11111001\tmq:Z:&47088')34\tnp:i:1\tqe:i:10\tqs:i:0\trq:f:0.811\tsn:B:f,2.0,2.0,2.0,2.0\tsq:Z:8<4<:<6<0<\tzm:i:%(zm2)s\tcx:i:2"""

basic_tags2 = {
    "bcv": "2.1",
    "flag": "4",
    "ipd": "Ipd:CodecV1",
    "sq": "\tsq:Z:8<4<:<6<0<",
    "so": "queryname",
    "sn": "2.0,2.0,2.0,2.0",
    "pos": "0",
    "rname": "*",
    "zm": "8",
    "pl": "PACBIO",
    "qname": "m140906_231018_42161_c100676332550000001823129611271486_s1_p0/8/0_10",
    "zm2": "9",
}

bad_tags2 = {
    "bcv": "9.8.7.6.5",
    "flag": "1",
    "ipd": "Ipd",
    "sq": "",
    "so": "unsorted",
    "sn": "-1,-1,-1,-1",
    "pos": "1",
    "rname": "*",
    "zm": "9",
    "pl": "OTHER",
    "qname": "movie1/8/0_10",
    "zm2": "7",
}
# XXX attempting to screw with sort order, but this isn't currently checked...
bad_tags3 = dict(basic_tags2)
bad_tags3["zm2"] = "7"
bad_tags4 = dict(basic_tags)
bad_tags4["ip"] = "S,255,2,0,512,22,34,0,2,3,0,16"
sam_strings = [
    sam_str_ % basic_tags,
    (sam_str_ + "\n" + rec1) % bad_tags,
    unmapped_sam_str_ % basic_tags2,
    unmapped_sam_str_ % bad_tags2,
    #unmapped_sam_str_ % bad_tags3,
    sam_str_ % bad_tags4
]

DATA_DIR = op.join(op.dirname(op.dirname(__file__)), "data")


def generate_data_files(dir_name=None):
    import logging
    logging.basicConfig(level=logging.INFO)
    if dir_name is not None:
        os.chdir(dir_name)
    with open("tst1.fasta", "w") as f:
        f.write(">ecoliK12_pbi_March2013_2955000_to_2980000\n")
        f.write("AAAGAGAGAG" * 2500)
    pysam.samtools.faidx("tst1.fasta", catch_stdout=False)
    for i in range(len(sam_strings)):
        sam_file = "tst_%d_subreads.sam" % (i + 1)
        bam_file = "tst_%d_subreads.bam" % (i + 1)
        with open(sam_file, "w") as sam_out:
            sam_out.write(sam_strings[i])
        logging.info("Converting {s} to BAM".format(s=sam_file))
        # FIXME pysam is way broken - can't handle unmapped input?
        # convert to bam using pysam
        # with pysam.AlignmentFile(sam_file, "r", check_sq=False) as sam_in:
        #    with pysam.AlignmentFile(bam_file, "wb",
        #                             template=sam_in) as bam_out:
        #        for s in sam_in:
        #            bam_out.write(s)
        args = ["samtools", "view", "-b", "-o", bam_file, sam_file]
        assert subprocess.call(args) == 0, args
        os.remove(sam_file)
        # XXX don't create .pbi for this file, we want it to be absent
        if bam_file != "tst_2_subreads.bam":
            logging.info("Indexing {b}".format(b=bam_file))
            subprocess.call(["pbindex", bam_file])


def remove_data_files():
    os.remove("tst1.fasta")
    for i in range(len(sam_strings)):
        os.remove("tst_%d_subreads.sam" % (i + 1))
        os.remove("tst_%d_subreads.bam" % (i + 1))


class TestCase (unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        return

    def test_api_1(self):
        file_name = op.join(DATA_DIR, "tst_1_subreads.bam")
        _validators = bam.get_validators(validate_index=True)

        def _run_validators(f, expected_failures, validators=_validators):
            errors = []
            for v in validators:
                if isinstance(v, ValidateFileObject):
                    if not v.validate(f):
                        errors.extend(v.to_errors(f))
            for rg in f.peer.header["RG"]:
                for v in validators:
                    if isinstance(v, ValidateReadGroup):
                        if not v.validate(rg):
                            errors.extend(v.to_errors(rg))
            for aln in f:
                for v in validators:
                    if isinstance(v, ValidateRecord):
                        if not v.validate(aln):
                            errors.extend(v.to_errors(aln))
            found = sorted(list(set([type(e).__name__ for e in errors])))
            expected = sorted(list(set(expected_failures)))
            self.assertEqual(found, expected)
        bam_file = pbcore.io.BamReader(file_name)
        _run_validators(f=bam_file, expected_failures=[])
        # now a bad one
        file_name = op.join(DATA_DIR, "tst_2_subreads.bam")
        bam_file = pbcore.io.BamReader(file_name)
        _run_validators(f=bam_file, expected_failures=[
            'AlignmentCigarMatchError', 'AlignmentNotUniqueError',
            'AlignmentUnmappedError', 'BasecallerVersionError',
            'MissingCodecError', 'MissingIndexError', 'MissingPlatformError',
            'QnameFormatError', 'QnameRangeError',
            'ReadGroupChemistryError',
            'ReadGroupIdMismatchError', "ReadLengthError", 'TagValueError',
            'UninitializedSNRError', 'UnsortedError'])
        # a good unaligned file
        file_name = op.join(DATA_DIR, "tst_3_subreads.bam")
        bam_file = pbcore.io.BamReader(file_name)
        _run_validators(f=bam_file, expected_failures=[])
        # a bad unaligned file
        file_name = op.join(DATA_DIR, "tst_4_subreads.bam")
        bam_file = pbcore.io.BamReader(file_name)
        _run_validators(f=bam_file, expected_failures=[
            'BasecallerVersionError', 'MissingCodecError',
            'QnameHoleNumberError', 'QnameMovieError',
            'ReadGroupChemistryError', 'UninitializedSNRError',
            'UnmappedPropertiesError', 'UnsortedError',
            'WrongPlatformError'])
        _validators = bam.get_validators(aligned=True)
        _run_validators(f=bam_file, expected_failures=[
            'BasecallerVersionError',
            'FileNotAlignedError', 'MissingCodecError',
            'QnameHoleNumberError', 'QnameMovieError',
            'ReadGroupChemistryError',
            'UninitializedSNRError', 'UnsortedError',
            'WrongPlatformError'],
            validators=_validators)

    def test_1(self):
        file_name = op.join(DATA_DIR, "tst_1_subreads.bam")
        e, c = bam.validate_bam(file_name)
        self.assertEqual(len(e), 0)

    def test_1b(self):
        file_name = op.join(DATA_DIR, "tst_1_subreads.bam")
        e, c = bam.validate_bam(file_name, aligned=False, contents="CCS")
        errors = sorted([type(err).__name__ for err in e])
        self.assertEqual(errors,
                         ['FileAlignedError', 'FileContentMismatchError'])

    def test_1c_reference_fasta(self):
        file_name = op.join(DATA_DIR, "tst_1_subreads.bam")
        fasta_file = op.join(DATA_DIR, "tst1.fasta")
        e, c = bam.validate_bam(file_name, reference=fasta_file)
        self.assertEqual(len(e), 0)
        e, c = bam.validate_bam(file_name, aligned=False)

    def test_2(self):
        # now some errors
        file_name = op.join(DATA_DIR, "tst_2_subreads.bam")
        e, c = bam.validate_bam(file_name, validate_index=True)
        errors = sorted(list(set([type(err).__name__ for err in e])))
        self.assertEqual(errors,
                         ['AlignmentCigarMatchError',
                          'AlignmentUnmappedError',
                          'BasecallerVersionError',
                          'MissingCodecError', 'MissingIndexError',
                          'MissingPlatformError', 'QnameFormatError',
                          'QnameRangeError', 'ReadGroupChemistryError',
                          'ReadGroupIdMismatchError', "ReadLengthError",
                          'TagValueError',
                          'UninitializedSNRError', 'UnsortedError'])

    def test_3_unmapped(self):
        file_name = op.join(DATA_DIR, "tst_3_subreads.bam")
        e, c = bam.validate_bam(file_name)
        self.assertEqual(len(e), 0)

    def test_4_unmapped(self):
        file_name = op.join(DATA_DIR, "tst_4_subreads.bam")
        e, c = bam.validate_bam(file_name)
        errors1 = sorted([type(err).__name__ for err in e])
        self.assertEqual(errors1, ['BasecallerVersionError',
                                   'MissingCodecError',
                                   'QnameHoleNumberError', 'QnameMovieError',
                                   'ReadGroupChemistryError',
                                   'UninitializedSNRError',
                                   'UnmappedPropertiesError', 'UnsortedError',
                                   'WrongPlatformError'])
        e, c = bam.validate_bam(file_name, aligned=True)
        errors2 = sorted([type(err).__name__ for err in e])
        self.assertEqual(errors2,
                         ['BasecallerVersionError',
                          'FileNotAlignedError', 'MissingCodecError',
                          'QnameHoleNumberError', 'QnameMovieError',
                          'ReadGroupChemistryError',
                          'UninitializedSNRError', 'UnsortedError',
                          'WrongPlatformError'])
        # this should yield the same result as the first run
        e, c = bam.validate_bam(file_name, aligned=False)
        errors3 = sorted([type(err).__name__ for err in e])
        self.assertEqual(errors3, errors1)

    def test_bad_encoding(self):
        file_name = op.join(DATA_DIR, "tst_5_subreads.bam")
        e, c = bam.validate_bam(file_name)
        errors1 = sorted([type(err).__name__ for err in e])
        self.assertEqual(errors1, ["BadEncodingError"])

    def test_exit_code_0(self):
        file_name = op.join(DATA_DIR, "tst_1_subreads.bam")
        rc = subprocess.call(["pbvalidate", file_name])
        self.assertEqual(rc, 0)

    def test_exit_code_1(self):
        file_name = op.join(DATA_DIR, "tst_s_subreads.bam")
        rc = subprocess.call(["pbvalidate", file_name])
        self.assertEqual(rc, 1)

    @unittest.skipUnless(op.isdir(TESTDATA), "Testdata not found")
    def test_transcript_bam(self):
        BAM = "/pbi/dept/secondary/siv/testdata/isoseqs/TranscriptSet/unpolished.bam"
        e, c = bam.validate_bam(BAM)
        self.assertEqual(len(e), 0)

    @unittest.skipUnless(op.isdir(TESTDATA), "Testdata not found")
    def test_overlapping_alignments(self):
        BAM = "/pbi/dept/secondary/siv/testdata/pbreports-unittest/data/mapping_stats/pbmm2/aligned.bam"
        e, c = bam.validate_bam(BAM, aligned=True)
        errors2 = list(set(sorted([type(err).__name__ for err in e])))
        self.assertEqual(errors2, ["AlignmentNotUniqueError"])


if __name__ == "__main__":
    if "--make-files" in sys.argv:
        generate_data_files(DATA_DIR)
    else:
        unittest.main()
