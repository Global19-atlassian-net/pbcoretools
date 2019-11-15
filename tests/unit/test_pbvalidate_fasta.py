
import subprocess
import tempfile
import unittest
import shutil
import os.path as op
import os

import pbcore.io

from pbcoretools.pbvalidate.fasta import *
from pbcoretools.pbvalidate import fasta
from pbcoretools.pbvalidate.core import run_validators, ValidatorErrorContext

test_sequences = [
    # 1. this should work
    """\
>chr1 Jackalope chromosome 1;length=7
GATTACA
>chr2 Jackalope chromosome 2;length=7
TTACAGA""",
    # 2. empty line
    """\
>chr1
gattaca

>chr2
ttacaga""",
    # 3. whitespace
    """\
>chr1
 GATTACA
>chr2
TTACAGA """,
    # 4. no wrapping; this is actually okay now
    """\
>chr1
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>chr2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc""",
    # 5. inconsistent wrapping (x2)
    """\
>chr1
aaaaaaaaa
ccccccccccc
>chr2
ggggggggg
ttttttt
ccccccccc""",
    # 6. globally inconsistent wrapping
    """\
>chr1
aaaaaaaaaa
cccccccccc
ggggggg
>chr2
tttttttt
cccccccc
gggggggg
""",
    # 7. this should be okay (note use of 'u', which is allowed)
    """\
>chr1
aaaaaaaaaa
cccccccccc
ggggg
>chr2
tttttttttt
cccccccccc
uuuuuuuuuu
>chr3
atatatat""",
    # 8. bad sequence
    """\
>chr1
aaaa..aa
>chr2 this alone should be okay
GaTcUSN
>chr3
aaZgtc""",
    # 9. bad sequence in strict mode
    """\
>chr1
GATcUSN
>chr2
nnnnnAACTnn""",
    # 10. '>' in header
    """\
>chr1 >12345
ttttt
>chr2 ;!.+ this should be okay
aaaaaaaaa""",
    # 11. asterisk at start of identifier
    """\
>*chr1 this should fail
aaaaaaaaa
>chr2* this should work
tttttt""",
    # 12. other bad characters
    """\
>human:chr1
ttttttt
>human,chr2
aaaaaaaaaaaa
>human"chr3"
ccccccccccccccc
>human;chr3 this should be okay
gggggggggg""",
    # 13. blank identifier
    """\
> chr1
ttttttt
>chr2
aaaaa""",
    # 14. duplicate identifier
    """\
>Jackalope chromosome 1;length=7
GATTACA
>Jackalope chromosome 2;length=7
TTACAGA""",
    # 15. missing sequence
    """\
>chr1
cgta
>chr2""",
]


def validate_file(file_name, strict=False, validate_index=False):
    validators = fasta.get_validators(
        strict=strict,
        validate_index=validate_index)
    return fasta.run_validators(
        context_class=ValidatorErrorContext,
        path=file_name,
        reader_class=pbcore.io.FastaReader,
        validators=validators)


class TestCase (unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._cwd = os.getcwd()
        tmp_dir = tempfile.mkdtemp()
        os.chdir(tmp_dir)
        for i, seq in enumerate(test_sequences, start=1):
            open("test_%d.fa" % i, "w").write(seq)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._cwd)

    def test_allowed(self):
        e, c = validate_file("test_1.fa")
        self.assertEqual(len(e), 0)

    def test_index(self):
        e, c = validate_file("test_1.fa", validate_index=True)
        self.assertEqual(len(e), 1)

    def test_empty_line(self):
        e, c = validate_file("test_2.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], EmptyLineError)

    # pure formatting errors (independent of parser)
    def test_whitespace(self):
        e, c = validate_file("test_3.fa")
        # FIXME there are actually two errors corresponding to different lines
        # in the file, but the current structure of the code only allows us
        # to catch one of these
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], WhitespaceError)

    def test_wrapping(self):
        e, c = validate_file("test_4.fa")
        self.assertEqual(len(e), 0)
        #self.assertIsInstance(e[0], NoWrappingError)
        e, c = validate_file("test_5.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], SeqWrappingError)
        #e, c = validate_file("test_6.fa")
        #self.assertEqual(len(e), 1)
        #self.assertIsInstance(e[0], GlobalWrappingError)
        e, c = validate_file("test_7.fa")
        self.assertEqual(len(e), 0)

    # sequence errors
    def test_bad_nucleotide(self):
        e, c = validate_file("test_8.fa")
        self.assertEqual(len(e), 2)
        self.assertIsInstance(e[0], BadNucleotideError)
        e, c = validate_file("test_9.fa", strict=True)
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], BadNucleotideError)
        e, c = validate_file("test_9.fa", strict=False)
        self.assertEqual(len(e), 0)

    def test_extra_gt(self):
        e, c = validate_file("test_10.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], ExtraGTError)

    def test_identifier(self):
        e, c = validate_file("test_11.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], IdentifierAsteriskError)
        e, c = validate_file("test_12.fa")
        self.assertEqual(len(e), 3)
        self.assertIsInstance(e[0], BadIdentifierError)
        e, c = validate_file("test_13.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], BlankIdentifierError)
        e, c = validate_file("test_14.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], DuplicateIdError)

    def test_misc(self):
        e, c = validate_file("test_15.fa")
        self.assertEqual(len(e), 1)
        self.assertIsInstance(e[0], MissingSequenceError)

    def test_gzip(self):
        """Test that gzipped files are handled correctly"""
        file_name = op.join(op.dirname(op.dirname(__file__)), "data",
                            "tst2.fasta.gz")
        e, m = validate_file(file_name)
        self.assertEqual(len(e), 1)
        self.assertEqual(str(
            e[0]), "Inconsistent line wrapping for sequence 'ecoliK12_pbi_March2013_2955000_to_2980000'")

    def test_exit_code_0(self):
        rc = subprocess.call(["pbvalidate", "test_1.fa"])
        self.assertEqual(rc, 0)

    def test_exit_code_1(self):
        rc = subprocess.call(["pbvalidate", "test_2.fa"])
        self.assertEqual(rc, 1)

    @unittest.skip("No longer applicable in Python3???")
    def test_carriage_returns(self):
        file_name = op.join(op.dirname(op.dirname(__file__)), "data",
                            "bc_bad_returns.fasta")
        e, m = validate_file(file_name)
        self.assertEqual(len(e), 1)

    def test_dos(self):
        file_name = op.join(op.dirname(op.dirname(__file__)), "data",
                            "tst_dos.fasta")
        e, m = validate_file(file_name)
        self.assertEqual(len(e), 0)

    def test_fsa_extension(self):
        shutil.copyfile("test_1.fa", "test_1.fsa")
        rc = subprocess.call(["pbvalidate", "test_1.fsa"])
        self.assertEqual(rc, 0)


if __name__ == "__main__":
    unittest.main()
