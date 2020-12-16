
import tempfile
import shutil
import json
import os.path as op
import os

import pytest

from pbcommand.testkit import PbIntegrationBase
from pbcommand.models.common import DataStore
from pbcore.io import BarcodeSet, SubreadSet

import pbtestdata

class TestRequireSinglePrimer(PbIntegrationBase):

    def _make_barcodeset(self, fasta_str):
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fasta").name
        with open(tmp_fasta, "w") as fasta_out:
            fasta_out.write(fasta_str)
        tmp_bc = tempfile.NamedTemporaryFile(suffix=".barcodeset.xml").name
        with BarcodeSet(tmp_fasta, generateIndices=True) as ds_bc:
            ds_bc.write(tmp_bc)
        return tmp_bc

    def test_integration_success(self):
        fasta_str = ">bc1\nacgtacgt"
        bc_single = self._make_barcodeset(fasta_str)
        args = [
            "python3", "-m", "pbcoretools.tasks.require_single_primer",
            bc_single
        ]
        self._check_call(args)
        assert not op.isfile("alarms.json")

    def test_integration_fail(self):
        assert not op.isfile("alarms.json")
        fasta_str = ">bc1\nacgtacgt\n>bc2\ntgcatgca"
        bc_multi = self._make_barcodeset(fasta_str)
        args = [
            "python3", "-m", "pbcoretools.tasks.require_single_primer",
            bc_multi
        ]
        with pytest.raises(Exception):
            self._check_call(args)
        assert op.isfile("alarms.json")


class TestDeleteBamResources(PbIntegrationBase):

    def test_integration(self):
        ds_in = pbtestdata.get_file("aligned-xml")
        prefix = op.basename(ds_in).split(".")[0]
        ds_dir = op.dirname(ds_in)
        for file_name in os.listdir(ds_dir):
            if file_name.startswith(prefix):
                shutil.copyfile(op.join(ds_dir, file_name), file_name)
        files = list(os.listdir(os.getcwd()))
        assert len(files) == 4
        args = [
            "python3", "-m", "pbcoretools.tasks.delete_bam_resources",
            op.basename(ds_in)
        ]
        self._check_call(args)
        files = list(os.listdir(os.getcwd()))
        n_removed = 0
        for fn in files:
            if fn.startswith(prefix) and not fn.endswith(".xml"):
                assert op.islink(fn) and not op.exists(op.realpath(fn))
                n_removed += 1
        assert n_removed == 3
        with open("removed_files.fofn") as fofn:
            assert len(fofn.readlines()) == 3


class TestConsolidateReadsBam(PbIntegrationBase):

    def test_integration_simple(self):
        ds_in = pbtestdata.get_file("ccs-sequel")
        args = [
            "python3", "-m", "pbcoretools.tasks.consolidate_reads_bam",
            ds_in
        ]
        self._check_call(args)
        assert op.isfile("reads.bam")


def _get_tmp_file_of_type(ft):
    ext = "." + ft.ext
    fn = tempfile.NamedTemporaryFile(suffix=ext).name
    with open(fn, "wt") as f_out:
        f_out.write("asdf")
    return fn


class TestCollectIsoseqFiles(PbIntegrationBase):

    def test_integration_simple_1(self):
        from pbcoretools.tasks.isoseq.collect_files import FILE_IDS_AND_NAMES
        args = ["python3", "-m" "pbcoretools.tasks.isoseq.collect_files"]
        for file_id, file_type, label in FILE_IDS_AND_NAMES:
            args.extend([
                "--{}".format(file_id.replace("_", "-")),
                _get_tmp_file_of_type(file_type)
            ])
        args.extend(["--single-sample", "--datastore", "datastore.json"])
        args.append(tempfile.NamedTemporaryFile(suffix=".bam").name)
        self._check_call(args)
        ds = DataStore.load_from_json("datastore.json")
        assert len(ds.files) == len(FILE_IDS_AND_NAMES)
