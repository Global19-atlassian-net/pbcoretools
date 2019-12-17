import subprocess
import warnings
import tempfile
import shutil
import os.path as op
import os
import pytest

from pbcommand.models import FileTypes
from pbcore.io import openDataFile, BamReader, SubreadSet

import pbtestdata

from pbcoretools import bamsieve
from pbcoretools.bamsieve import UserError

DATA_DIR = op.join(op.dirname(op.dirname(__file__)), "data")
SUBREADS1 = op.join(DATA_DIR, "tst_1_subreads.bam")
DS1 = op.join(DATA_DIR, "tst_1.subreadset.xml")
SUBREADS2 = op.join(DATA_DIR, "tst_3_subreads.bam")
DS2 = op.join(DATA_DIR, "tst_3.subreadset.xml")
SUBREADS3 = pbtestdata.get_file("subreads-bam")
SUBREADS4 = pbtestdata.get_file("aligned-bam")
CCS = pbtestdata.get_file("ccs-bam")
BARCODED = pbtestdata.get_file("barcoded-subreads-bam")
BARCODED_DS = pbtestdata.get_file("barcoded-subreadset")
SUBREADS_STS = pbtestdata.get_file("subreads-sequel")


class TestBamSieve:

    def test_whitelist(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        WHITELIST = set([24962, 32901, 30983])

        def _run_with_whitelist(wl):
            rc = bamsieve.filter_reads(
                input_bam=SUBREADS3,
                output_bam=ofn,
                whitelist=wl)
            assert rc == 0
            with BamReader(ofn) as bam_out:
                have_zmws = set([rec.HoleNumber for rec in bam_out])
                assert have_zmws == WHITELIST
        _run_with_whitelist(WHITELIST)
        _run_with_whitelist(",".join([str(x) for x in list(WHITELIST)]))
        tmp_wl = tempfile.NamedTemporaryFile(suffix=".txt").name
        with open(tmp_wl, "w") as wl_out:
            wl_out.write("\n".join([str(x) for x in list(WHITELIST)]))
        _run_with_whitelist(tmp_wl)
        # now with a BAM file as whitelist
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            whitelist=SUBREADS4)
        with BamReader(ofn) as bam_out:
            assert 117 == len([rec for rec in bam_out])

    def test_subreads_whitelist(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        ofn2 = tempfile.NamedTemporaryFile(suffix=".bam").name
        WHITELIST = set(['m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/1650/1920_2155',
                         'm140905_042212_sidney_c100564852550000001823085912221377_s1_X0/7957/9554_9634',
                         'm140905_042212_sidney_c100564852550000001823085912221377_s1_X0/1650/2200_3298'])
        ZMWS = set([1650, 7957])

        def _run_with_whitelist(wl):
            rc = bamsieve.filter_reads(
                input_bam=SUBREADS3,
                output_bam=ofn,
                whitelist=wl,
                use_subreads=True)
            assert rc == 0
            with BamReader(ofn) as bam_out:
                have_zmws = set([rec.HoleNumber for rec in bam_out])
                assert have_zmws == ZMWS
                qnames = set([rec.qName for rec in bam_out])
                assert qnames == WHITELIST

        _run_with_whitelist(WHITELIST)
        _run_with_whitelist(",".join([str(x) for x in list(WHITELIST)]))
        tmp_wl = tempfile.NamedTemporaryFile(suffix=".txt").name
        with open(tmp_wl, "w") as wl_out:
            wl_out.write("\n".join([str(x) for x in list(WHITELIST)]))
        _run_with_whitelist(tmp_wl)
        # now with a BAM file as whitelist
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn2,
            use_subreads=True,
            whitelist=ofn)

        with BamReader(ofn) as bam_out:
            subreads = set([x.qName for x in bam_out])
        with BamReader(ofn2) as bam_out:
            subreads2 = set([x.qName for x in bam_out])
        assert subreads == subreads2

    def test_subreads_blacklist(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        ofn2 = tempfile.NamedTemporaryFile(suffix=".bam").name
        BLACKLIST = set(['m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/1650/1920_2155',
                         'm140905_042212_sidney_c100564852550000001823085912221377_s1_X0/7957/9554_9634',
                         'm140905_042212_sidney_c100564852550000001823085912221377_s1_X0/1650/2200_3298'])

        def _run_with_blacklist(bl):
            rc = bamsieve.filter_reads(
                input_bam=SUBREADS3,
                output_bam=ofn,
                blacklist=bl,
                use_subreads=True)
            assert rc == 0
            with BamReader(ofn) as bam_out:
                qnames = set([rec.qName for rec in bam_out])
                assert qnames & BLACKLIST == set()
                assert len([x for x in bam_out]) == 114

        _run_with_blacklist(BLACKLIST)
        _run_with_blacklist(",".join([str(x) for x in list(BLACKLIST)]))
        tmp_wl = tempfile.NamedTemporaryFile(suffix=".txt").name
        with open(tmp_wl, "w") as wl_out:
            wl_out.write("\n".join([str(x) for x in list(BLACKLIST)]))
        _run_with_blacklist(tmp_wl)

        # now with the BAM file we just made as blacklist
        EXPECTED_OUT = BLACKLIST
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn2,
            use_subreads=True,
            blacklist=ofn)

        with BamReader(ofn) as bam_out:
            subreads = set([x.qName for x in bam_out])
        with BamReader(ofn2) as bam_out:
            subreads2 = set([x.qName for x in bam_out])
        assert subreads & subreads2 == set()
        assert subreads2 == EXPECTED_OUT

        # now an integration test, because this is used in Cromwell workflow
        ofn3 = tempfile.NamedTemporaryFile(suffix=".subreads.bam").name
        args = ["bamsieve", "--subreads", "--blacklist", ofn, SUBREADS3, ofn3]
        rc = subprocess.check_call(args)
        with BamReader(ofn3) as bam_out:
            subreads3 = set([x.qName for x in bam_out])
            assert subreads & subreads3 == set()
            assert subreads3 == EXPECTED_OUT
        # and again, with a dataset as input
        ds_tmp = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        with SubreadSet(ofn) as ds:
            ds.write(ds_tmp)
        ofn4 = tempfile.NamedTemporaryFile(suffix=".subreads.bam").name
        args = ["bamsieve", "--subreads",
                "--blacklist", ds_tmp, SUBREADS3, ofn4]
        rc = subprocess.check_call(args)
        with BamReader(ofn4) as bam_out:
            subreads4 = set([x.qName for x in bam_out])
            assert subreads & subreads4 == set()
            assert subreads4 == EXPECTED_OUT

    def test_dataset_io(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        rc = bamsieve.filter_reads(
            input_bam=DS2,
            output_bam=ofn,
            whitelist="8")
        assert rc == 0
        with SubreadSet(ofn, strict=False) as bam_out:
            with SubreadSet(DS2) as ds_in:
                assert ds_in.uuid != bam_out.uuid
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            assert have_zmws == set([8])
        # make sure paths are absolute
        tmpdir = tempfile.mkdtemp()
        ofn2 = op.join(tmpdir, op.basename(ofn))
        shutil.copyfile(ofn, ofn2)
        with SubreadSet(ofn2, strict=False) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            assert have_zmws == set([8])

    def test_dataset_relative_paths(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        basename = op.basename(ofn).split(".")[0]
        rc = bamsieve.filter_reads(
            input_bam=DS2,
            output_bam=ofn,
            whitelist="8",
            relative=True)
        assert rc == 0
        # move everything to another directory and make sure paths still work
        tmpdir = tempfile.mkdtemp()
        for file_name in os.listdir(op.dirname(ofn)):
            if file_name.startswith(basename):
                shutil.move(op.join(op.dirname(ofn), file_name),
                            op.join(tmpdir, file_name))
        ofn2 = op.join(tmpdir, op.basename(ofn))
        with SubreadSet(ofn2, strict=False) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            assert have_zmws == set([8])

    def test_anonymize(self):
        ofn1 = tempfile.NamedTemporaryFile(suffix=".bam").name
        ofn2 = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn1,
            whitelist=set([24962]))
        assert rc == 0
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn2,
            whitelist=set([24962]),
            anonymize=True)
        assert rc == 0
        with openDataFile(ofn1) as bam1:
            with openDataFile(ofn2) as bam2:
                for rec1, rec2 in zip(bam1, bam2):
                    assert rec1.qName == rec2.qName
                    assert rec1.peer.seq != rec2.peer.seq

    def test_blacklist(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name

        def _run_with_blacklist(bl):
            rc = bamsieve.filter_reads(
                input_bam=SUBREADS2,
                output_bam=ofn,
                blacklist=bl)
            assert rc == 0
            with BamReader(ofn) as bam_out:
                have_zmws = set([rec.HoleNumber for rec in bam_out])
                assert have_zmws == set([9])
        _run_with_blacklist(set([8]))
        _run_with_blacklist("8,233")
        tmp_bl = tempfile.NamedTemporaryFile(suffix=".txt").name
        with open(tmp_bl, "w") as bl_out:
            bl_out.write("8\n233")
        _run_with_blacklist(tmp_bl)

    def test_barcodes(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamsieve.filter_reads(
            input_bam=BARCODED,
            output_bam=ofn,
            whitelist=[0],
            use_barcodes=True)
        with BamReader(ofn) as bam_out:
            zmws = set([rec.HoleNumber for rec in bam_out])
            assert len(zmws) == 1
            assert 74056024 in zmws

    def test_subreadset_scraps(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        rc = bamsieve.filter_reads(
            input_bam=BARCODED_DS,
            output_bam=ofn,
            whitelist=[74056024])
        assert rc == 0

        def _verify():
            with SubreadSet(ofn, strict=False) as ds_out:
                ext_res = ds_out.externalResources[0]
                assert ext_res.bam.endswith(".subreads.bam")
                assert ext_res.scraps.endswith(".scraps.bam")
                for bam_file in [ext_res.bam, ext_res.scraps]:
                    with BamReader(bam_file) as bam:
                        zmws = set([rec.HoleNumber for rec in bam])
                        assert len(zmws) == 1
                        assert 74056024 in zmws
        _verify()
        rc = bamsieve.filter_reads(
            input_bam=BARCODED_DS,
            output_bam=ofn,
            count=1,
            seed=1)
        _verify()
        rc = bamsieve.filter_reads(
            input_bam=BARCODED_DS,
            output_bam=ofn,
            blacklist=[28901719])
        assert rc == 0

    def test_percentage(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            percentage=50,
            seed=12345)
        assert rc == 0
        with BamReader(ofn) as bam_out:
            zmws = set([rec.HoleNumber for rec in bam_out])
            assert len(zmws) == 24

    def test_count(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS3,
            output_bam=ofn,
            count=1,
            seed=12345)
        assert rc == 0
        with BamReader(ofn) as bam_out:
            zmws = set([rec.HoleNumber for rec in bam_out])
            assert len(zmws) == 1

    def test_count_overflow(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        with warnings.catch_warnings(record=True) as w:
            rc = bamsieve.filter_reads(
                input_bam=SUBREADS3,
                output_bam=ofn,
                count=100000,
                seed=12345)
            assert rc == 0
            assert len(w) == 1
            with BamReader(ofn) as bam_out:
                zmws = set([rec.HoleNumber for rec in bam_out])
                assert len(zmws) == 48

    @pytest.mark.constools
    def test_sts_xml(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        rc = bamsieve.filter_reads(
            input_bam=SUBREADS_STS,
            output_bam=ofn,
            count=1,
            seed=12345,
            keep_original_uuid=True)
        assert rc == 0
        with SubreadSet(ofn, strict=True) as ds:
            with SubreadSet(SUBREADS_STS) as ds_in:
                assert ds_in.uuid == ds.uuid
            for er in ds.externalResources:
                if er.metaType == FileTypes.BAM_SUB.file_type_id:
                    assert er.sts is not None
                    assert os.path.exists(er.sts)
                    assert er.sts != ds_in.externalResources[0].sts
                    break
            else:
                self.fail("Can't find subreads BAM")

    def test_error(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        with pytest.raises(UserError) as exc:
            rc = bamsieve.filter_reads(
                input_bam=DS1,
                output_bam=ofn,
                whitelist=set([5, 6, 7, 8]),
                blacklist=set([1, 2, 3, 4]))
        with pytest.raises(UserError) as exc:
            rc = bamsieve.filter_reads(
                input_bam=DS1,
                output_bam=ofn,
                whitelist=set([5, 6, 7, 8]),
                percentage=50)
        with pytest.raises(UserError) as exc:
            rc = bamsieve.filter_reads(
                input_bam=DS1,
                output_bam=ofn,
                percentage=500)
        with pytest.raises(UserError) as exc:
            rc = bamsieve.filter_reads(
                input_bam=DS1,
                output_bam=ofn,
                percentage=50,
                count=1)
        # dataset output, but BAM input
        ofn = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        with pytest.raises(UserError) as exc:
            rc = bamsieve.filter_reads(
                input_bam=SUBREADS2,
                output_bam=ofn,
                percentage=50)

    def test_integration(self):
        args = ["bamsieve", "--help"]
        with tempfile.TemporaryFile() as stdout:
            with tempfile.TemporaryFile() as stderr:
                rc = subprocess.call(args, stdout=stdout, stderr=stderr)
                assert rc == 0
        ofn = tempfile.NamedTemporaryFile(suffix=".bam").name
        args = [
            "bamsieve",
            "--log-level", "ERROR",
            "--whitelist", "8,233",
            SUBREADS2,
            ofn
        ]
        rc = subprocess.call(args)
        assert rc == 0
        with BamReader(ofn) as bam_out:
            have_zmws = set([rec.HoleNumber for rec in bam_out])
            assert have_zmws == set([8])
