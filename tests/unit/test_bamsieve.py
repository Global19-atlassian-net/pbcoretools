
import subprocess
import tempfile
import unittest
import os.path as op

from pbcore.io import openDataFile, openDataSet
import pbcore.data

from pbcoretools import bamSieve

DATA_DIR = op.join(op.dirname(op.dirname(__file__)), "data")
SUBREADS1 = op.join(DATA_DIR, "tst_1_subreads.bam")
DS1 = op.join(DATA_DIR, "tst_1.subreadset.xml")
SUBREADS2 = op.join(DATA_DIR, "tst_3_subreads.bam")
DS2 = op.join(DATA_DIR, "tst_3.subreadset.xml")
SUBREADS3 = pbcore.data.getUnalignedBam()
SUBREADS4 = pbcore.data.getBamAndCmpH5()[0]
CCS = pbcore.data.getCCSBAM()


class TestBamSieve(unittest.TestCase):

    def test_whitelist(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        WHITELIST = set([24962, 32901, 30983])

        def _run_with_whitelist(wl):
            rc = bamSieve.filter_reads(
                input_bam=SUBREADS3,
                output_bam=ofn,
                whitelist=wl)
            self.assertEqual(rc, 0)
            with openDataFile(ofn, strict=False) as bam_out:
                have_zmws = set([rec.HoleNumber for rec in bam_out])
                self.assertEqual(have_zmws, WHITELIST)
        _run_with_whitelist(WHITELIST)
        _run_with_whitelist(",".join([str(x) for x in list(WHITELIST)]))
        tmp_wl = tempfile.NamedTemporaryFile(suffix=".txt").name
        with open(tmp_wl, "w") as wl_out:
            wl_out.write("\n".join([str(x) for x in list(WHITELIST)]))
        _run_with_whitelist(tmp_wl)
        # now with a BAM file as whitelist
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            whitelist=SUBREADS4)
        with openDataFile(ofn, strict=False) as bam_out:
            self.assertEqual(117, len([rec for rec in bam_out]))

    def test_dataset_io(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        rc = bamSieve.filter_reads(
            input_bam=DS2,
            output_bam=ofn,
            whitelist="8")
        self.assertEqual(rc, 0)
        with openDataSet(ofn, strict=False) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(have_zmws, set([8]))

    def test_blacklist(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name

        def _run_with_blacklist(bl):
            rc = bamSieve.filter_reads(
                input_bam=SUBREADS2,
                output_bam=ofn,
                blacklist=bl)
            self.assertEqual(rc, 0)
            with openDataFile(ofn, strict=False) as bam_out:
                have_zmws = set([rec.HoleNumber for rec in bam_out])
                self.assertEqual(have_zmws, set([9]))
        _run_with_blacklist(set([8]))
        _run_with_blacklist("8,233")
        tmp_bl = tempfile.NamedTemporaryFile(suffix=".txt").name
        with open(tmp_bl, "w") as bl_out:
            bl_out.write("8\n233")
        _run_with_blacklist(tmp_bl)

    def test_percentage(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            percentage=50,
            seed=12345)
        self.assertEqual(rc, 0)
        with openDataFile(ofn, strict=False) as bam_out:
            n_records = len([rec for rec in bam_out])
            print n_records
            #self.assertEqual(n_records, 1)

    def test_error(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamSieve.filter_reads(
            input_bam=DS1,
            output_bam=ofn,
            whitelist=set([5, 6, 7, 8]),
            blacklist=set([1, 2, 3, 4]))
        self.assertEqual(rc, 1)
        rc = bamSieve.filter_reads(
            input_bam=DS1,
            output_bam=ofn,
            whitelist=set([5, 6, 7, 8]),
            percentage=50)
        self.assertEqual(rc, 1)
        rc = bamSieve.filter_reads(
            input_bam=DS1,
            output_bam=ofn,
            percentage=500)
        self.assertEqual(rc, 1)
        # dataset output, but BAM input
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS2,
            output_bam=ofn,
            percentage=50)
        self.assertEqual(rc, 1)

    def test_integration(self):
        args = ["bamSieve", "--help"]
        with tempfile.TemporaryFile() as stdout:
            with tempfile.TemporaryFile() as stderr:
                rc = subprocess.call(args, stdout=stdout, stderr=stderr)
                self.assertEqual(rc, 0)
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        args = [
            "bamSieve",
            "--log-level", "ERROR",
            "--whitelist", "8,233",
            SUBREADS2,
            ofn
        ]
        rc = subprocess.call(args)
        self.assertEqual(rc, 0)
        with openDataFile(ofn, strict=False) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(have_zmws, set([8]))


if __name__ == "__main__":
    unittest.main()
