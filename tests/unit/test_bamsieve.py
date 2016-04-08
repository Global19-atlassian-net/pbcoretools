
import subprocess
import tempfile
import unittest
import shutil
import os.path as op
import os

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
BARCODED = op.join(DATA_DIR, "barcoded.subreads.bam")
BARCODED_DS = op.join(DATA_DIR, "barcoded.subreadset.xml")

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
        # make sure paths are absolute
        tmpdir = tempfile.mkdtemp()
        ofn2 = op.join(tmpdir, op.basename(ofn))
        shutil.copyfile(ofn, ofn2)
        with openDataSet(ofn2, strict=False) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(have_zmws, set([8]))

    def test_dataset_relative_paths(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        basename = op.basename(ofn).split(".")[0]
        rc = bamSieve.filter_reads(
            input_bam=DS2,
            output_bam=ofn,
            whitelist="8",
            relative=True)
        self.assertEqual(rc, 0)
        # move everything to another directory and make sure paths still work
        tmpdir = tempfile.mkdtemp()
        for file_name in os.listdir(op.dirname(ofn)):
            if file_name.startswith(basename):
                shutil.move(op.join(op.dirname(ofn), file_name),
                            op.join(tmpdir, file_name))
        ofn2 = op.join(tmpdir, op.basename(ofn))
        with openDataSet(ofn2, strict=False) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(have_zmws, set([8]))

    def test_anonymize(self):
        ofn1 = tempfile.NamedTemporaryFile(suffix=".bam").name
        ofn2 = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn1,
            whitelist=set([24962]))
        self.assertEqual(rc, 0)
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn2,
            whitelist=set([24962]),
            anonymize=True)
        self.assertEqual(rc, 0)
        with openDataFile(ofn1, strict=False) as bam1:
            with openDataFile(ofn2, strict=False) as bam2:
                for rec1, rec2 in zip(bam1, bam2):
                    self.assertEqual(rec1.qName, rec2.qName)
                    self.assertNotEqual(rec1.peer.seq, rec2.peer.seq)

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

    def test_barcodes(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamSieve.filter_reads(
            input_bam=BARCODED,
            output_bam=ofn,
            whitelist=[0],
            use_barcodes=True)
        with openDataFile(ofn, strict=False) as bam_out:
            zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(len(zmws), 1)
            self.assertTrue(74056024 in zmws)

    def test_subreadset_scraps(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        rc = bamSieve.filter_reads(
            input_bam=BARCODED_DS,
            output_bam=ofn,
            whitelist=[74056024])
        self.assertEqual(rc, 0)
        def _verify():
            with openDataFile(ofn, strict=False) as ds_out:
                ext_res = ds_out.externalResources[0]
                for bam_file in [ext_res.bam, ext_res.scraps]:
                    with openDataFile(bam_file, strict=False) as bam:
                        zmws = set([rec.HoleNumber for rec in bam])
                        self.assertEqual(len(zmws), 1)
                        self.assertTrue(74056024 in zmws)
        _verify()
        rc = bamSieve.filter_reads(
            input_bam=BARCODED_DS,
            output_bam=ofn,
            blacklist=[28901719])
        self.assertEqual(rc, 0)

    def test_percentage(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            percentage=50,
            seed=12345)
        self.assertEqual(rc, 0)
        with openDataFile(ofn, strict=False) as bam_out:
            zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(len(zmws), 24)

    def test_count(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamSieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            count=1,
            seed=12345)
        self.assertEqual(rc, 0)
        with openDataFile(ofn, strict=False) as bam_out:
            zmws = set([rec.HoleNumber for rec in bam_out])
            self.assertEqual(len(zmws), 1)

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
        rc = bamSieve.filter_reads(
            input_bam=DS1,
            output_bam=ofn,
            percentage=50,
            count=1)
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
