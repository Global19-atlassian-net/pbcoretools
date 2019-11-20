import xml.etree.cElementTree as ET
import subprocess
import itertools
import tempfile
import logging
import shutil
import pytest
import time
import re
import os

import numpy as np

from pbcore.io import PacBioBamIndex, IndexedBamReader
from pbcore.io import openIndexedAlignmentFile
from pbcore.io import (DataSet, SubreadSet, ReferenceSet, AlignmentSet,
                       openDataSet, DataSetMetaTypes,
                       ConsensusReadSet, ConsensusAlignmentSet)
from pbcore.io.dataset.DataSetIO import dsIdToSuffix
from pbcore.io.dataset.DataSetMembers import ExternalResource, Filters
from pbcore.io.dataset.DataSetWriter import toXml
from pbcore.io.dataset.DataSetValidator import validateFile
import pbcore.data.datasets as data
import pbcore.data as otherdata

import pbtestdata

log = logging.getLogger(__name__)


def _is_relative(xmlfile):
    for event, element in ET.iterparse(xmlfile, events=("start",)):
        if element.get("ResourceId", "noop").startswith(os.sep):
            return False
    return True


class TestDataSet:
    """Unit and integrationt tests for the DataSet class and \
    associated module functions"""

    def _check_cmd(self, cmd):
        log.debug(cmd)
        subprocess.check_call(cmd.split())

    def _run_cmd_with_output(self, cmd, outfile):
        self._check_cmd(cmd)
        assert os.path.exists(outfile)

    @pytest.mark.constools
    def test_split_cli(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset split --outdir {o} --contigs --chunks 2 {d}".format(
            o=outdir,
            d=data.getXml(7))
        self._check_cmd(cmd)
        assert os.path.exists(
            os.path.join(outdir, 'pbalchemysim0.chunk0.alignmentset.xml'))
        assert os.path.exists(
            os.path.join(outdir, 'pbalchemysim0.chunk1.alignmentset.xml'))

    def test_split_cli_targetsize(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset split --outdir {o} --zmws --targetSize 2 {d}".format(
            o=outdir,
            d=data.getXml(7))
        self._check_cmd(cmd)
        for i in range(4):
            assert os.path.exists(
                os.path.join(
                    outdir,
                    'pbalchemysim0.chunk{}.alignmentset.xml'.format(i)))

        # pbalchemysim0.chunk4.alignmentset.xml shouldn't exist
        assert not os.path.exists(
            os.path.join(
                outdir,
                'pbalchemysim0.chunk4.alignmentset.xml'))

    def test_copyTo_cli(self):
        # To a fname:
        # absolute:
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        cmd = "dataset copyto {i} {o}".format(i=data.getXml(7), o=fn)
        self._run_cmd_with_output(cmd, fn)
        sset = AlignmentSet(fn, strict=True)
        assert not _is_relative(fn)

    def test_copyTo_cli_relative(self):
        # relative:
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        cmd = "dataset copyto --relative {i} {o}".format(
            i=data.getXml(7), o=fn)
        self._run_cmd_with_output(cmd, fn)
        sset = AlignmentSet(fn, strict=True)
        assert _is_relative(fn)

    def test_copyTo_cli_absolute_dir(self):
        # to a directory:
        # absolute:
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, os.path.split(data.getXml(7))[1])
        cmd = "dataset copyto {i} {o}".format(i=data.getXml(7), o=outdir)
        self._run_cmd_with_output(cmd, fn)
        sset = AlignmentSet(fn, strict=True)
        assert not _is_relative(fn)

    def test_copyTo_cli_relative_dir(self):
        # relative:
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, os.path.split(data.getXml(7))[1])
        cmd = "dataset copyto --relative {i} {o}".format(i=data.getXml(7),
                                                         o=outdir)
        self._run_cmd_with_output(cmd, fn)
        sset = AlignmentSet(fn, strict=True)
        assert _is_relative(fn)

    def test_newUuid_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(7))
        aln.copyTo(fn)
        pre_uuid = AlignmentSet(fn).uuid
        cmd = "dataset newuuid {d}".format(d=fn)
        self._run_cmd_with_output(cmd, fn)
        post_uuid = AlignmentSet(fn).uuid
        assert pre_uuid != post_uuid
        cmd = "dataset newuuid --updateCounts {d}".format(d=fn)
        self._run_cmd_with_output(cmd, fn)
        post_uuid = AlignmentSet(fn).uuid
        assert pre_uuid != post_uuid

    def test_newUuid_clone_cli(self):
        fn_orig = data.getXml(7)
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, 'fn.alignmentset.xml')
        fn2 = os.path.join(outdir, 'fn2.alignmentset.xml')
        with AlignmentSet(fn_orig) as aln:
            aln.copyTo(fn)
            shutil.copy(fn, fn2)

        pre_uuid = AlignmentSet(fn).uuid
        pre_uuid2 = AlignmentSet(fn2).uuid
        assert pre_uuid == pre_uuid2

        cmd = "dataset newuuid {d}".format(d=fn)
        self._run_cmd_with_output(cmd, fn)

        cmd = "dataset newuuid {d}".format(d=fn2)
        self._run_cmd_with_output(cmd, fn2)

        post_uuid = AlignmentSet(fn).uuid
        post_uuid2 = AlignmentSet(fn2).uuid
        assert pre_uuid != post_uuid
        assert pre_uuid2 != post_uuid2
        # HASH, THEREFORE THESE ARE EQUAL:
        assert post_uuid == post_uuid2

    def test_newUuid_random_cli(self):
        fn_orig = data.getXml(7)
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, 'fn.alignmentset.xml')
        fn2 = os.path.join(outdir, 'fn2.alignmentset.xml')
        with AlignmentSet(fn_orig) as aln:
            aln.copyTo(fn)
            shutil.copy(fn, fn2)

        pre_uuid = AlignmentSet(fn).uuid
        pre_uuid2 = AlignmentSet(fn2).uuid
        assert pre_uuid == pre_uuid2

        cmd = "dataset newuuid --random {d}".format(d=fn)
        self._run_cmd_with_output(cmd, fn)

        cmd = "dataset newuuid --random {d}".format(d=fn2)
        self._run_cmd_with_output(cmd, fn2)

        post_uuid = AlignmentSet(fn).uuid
        post_uuid2 = AlignmentSet(fn2).uuid
        assert pre_uuid != post_uuid
        assert pre_uuid2 != post_uuid2
        # RANDOM, THEREFORE THESE ARE NOT EQUAL:
        assert post_uuid != post_uuid2

    def test_relativize_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(7))
        aln.copyTo(fn)
        assert not _is_relative(fn)
        cmd = "dataset relativize {d}".format(d=fn)
        self._run_cmd_with_output(cmd, fn)
        assert _is_relative(fn)

    def test_absolutize_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(7))
        aln.copyTo(fn, relative=True)
        assert _is_relative(fn)
        cmd = "dataset absolutize {d}".format(d=fn)
        self._run_cmd_with_output(cmd, fn)
        assert not _is_relative(fn)

    def test_absolutize_cli_2(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        outfn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(7))
        aln.copyTo(fn, relative=True)
        assert _is_relative(fn)
        cmd = "dataset absolutize {d} --outdir {o}".format(d=fn, o=outfn)
        self._run_cmd_with_output(cmd, fn)
        assert os.path.exists(outfn)
        assert _is_relative(fn)
        assert not _is_relative(outfn)

    def test_absolutize_cli_3(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, os.path.split(fn)[1])
        aln = AlignmentSet(data.getXml(7))
        aln.copyTo(fn, relative=True)
        assert _is_relative(fn)
        cmd = "dataset absolutize {d} --outdir {o}".format(d=fn, o=outdir)
        self._run_cmd_with_output(cmd, fn)
        assert os.path.exists(outfn)
        assert _is_relative(fn)
        assert not _is_relative(outfn)

    @pytest.mark.internal_data
    def test_loadmetadata_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        log.debug(fn)

        aln = AlignmentSet(data.getXml(7))
        aln.metadata.collections = None
        aln.copyTo(fn)
        aln.close()
        del aln
        assert os.path.exists(fn)

        aln = AlignmentSet(fn)
        assert not aln.metadata.collections

        cmd = "dataset loadmetadata {i} {m}".format(
            i=fn,
            m=("/pbi/dept/secondary/siv/testdata/"
               "SA3-Sequel/lambda/roche_SAT/"
               "m54013_151205_032353.run.metadata.xml"))
        self._check_cmd(cmd)
        aln = AlignmentSet(fn)
        assert aln.metadata.collections

    @pytest.mark.internal_data
    def test_loadmetadata_from_dataset_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        log.debug(fn)

        aln = AlignmentSet(data.getXml(7))
        aln.metadata.collections = None
        aln.copyTo(fn)
        aln.close()
        del aln
        assert os.path.exists(fn)

        aln = AlignmentSet(fn)
        assert not aln.metadata.collections

        cmd = "dataset loadmetadata {i} {m}".format(
            i=fn,
            m=("/pbi/dept/secondary/siv/testdata/"
               "SA3-Sequel/lambda/roche_SAT/"
               "m54013_151205_032353.subreadset.xml"))
        self._check_cmd(cmd)
        aln = AlignmentSet(fn)
        assert aln.metadata.collections

    @pytest.mark.internal_data
    def test_loadmetadata_from_dataset_create_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        fn2 = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        log.debug(fn)

        aln = AlignmentSet(data.getXml(7))
        aln.metadata.collections = None
        aln.copyTo(fn)
        aln.close()
        del aln
        assert os.path.exists(fn)

        aln = AlignmentSet(fn)
        assert not aln.metadata.collections

        cmd = "dataset create --metadata {m} {o} {i}".format(
            o=fn2,
            i=fn,
            m=("/pbi/dept/secondary/siv/testdata/"
               "SA3-Sequel/lambda/roche_SAT/"
               "m54013_151205_032353.subreadset.xml"))
        self._check_cmd(cmd)
        aln = AlignmentSet(fn2)
        assert aln.metadata.collections

    def test_contigset_split_cli(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset split --outdir {o} --chunks 2 {d}".format(
            o=outdir,
            d=data.getXml(8))
        self._check_cmd(cmd)
        assert os.path.exists(
            os.path.join(outdir,
                         "pbalchemysim0.chunk0.referenceset.xml"))
        assert os.path.exists(
            os.path.join(outdir,
                         "pbalchemysim0.chunk1.referenceset.xml"))

    def test_filter_cli(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, "filtered8.xml")
        log.debug(outfn)
        cmd = "dataset filter {i} {o} {f}".format(
            i=data.getXml(7),
            o=outfn,
            f="rname=E.faecalis.1")
        self._run_cmd_with_output(cmd, outfn)
        aln = AlignmentSet(data.getXml(7))
        aln.filters.addRequirement(rname=[('=', 'E.faecalis.1')])
        aln.updateCounts()
        dset = AlignmentSet(outfn)
        assert str(aln.filters) == str(dset.filters)
        assert aln.totalLength == dset.totalLength
        assert aln.numRecords == dset.numRecords

    def _get_mock_alignment_set_out(self, outdir):
        return os.path.join(outdir, os.path.basename(data.getXml(11)))

    def test_create_cli(self):
        log.debug("Absolute")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset create --type AlignmentSet {o} {i1} {i2}".format(
            o=os.path.join(outdir, 'pbalchemysim.alignmentset.xml'),
            i1=data.getXml(7), i2=data.getXml(10))
        self._check_cmd(cmd)
        assert os.path.exists(
            os.path.join(outdir, os.path.basename(data.getXml(11))))

    def test_create_cli_relative(self):
        log.debug("Relative")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = ("dataset create --relative --type AlignmentSet "
               "{o} {i1} {i2}".format(
                   o=ofn,
                   i1=data.getXml(7),
                   i2=data.getXml(10)))
        self._check_cmd(cmd)
        assert os.path.exists(ofn)

    def test_create_cli_noclobber(self):
        log.debug("Don't clobber existing")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        with open(ofn, "w") as xml_out:
            xml_out.write("<AlignmentSet/>")
        mtime = os.path.getmtime(ofn)
        # in case we do create a new file (bad), we need to detect that with a
        # different mtime (seconds from epoch)
        time.sleep(1)
        cmd = ("dataset create --relative --type AlignmentSet "
               "{o} {i1} {i2}".format(
                   o=os.path.join(outdir, 'pbalchemysim.alignmentset.xml'),
                   i1=data.getXml(7),
                   i2=data.getXml(10)))
        log.debug(cmd)
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.check_call(cmd.split())
        assert os.path.exists(ofn)
        assert mtime == os.path.getmtime(ofn)

    def test_create_cli_clobber(self):
        log.debug("Force clobber existing")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        with open(ofn, "w") as xml_out:
            xml_out.write("<AlignmentSet/>")
        mtime = os.path.getmtime(ofn)
        # We want to create a new file, we need to detect that with a
        # different mtime (seconds from epoch)
        time.sleep(1)
        cmd = ("dataset create --force --relative --type AlignmentSet "
               "{o} {i1} {i2}".format(
                   o=os.path.join(outdir, 'pbalchemysim.alignmentset.xml'),
                   i1=data.getXml(7),
                   i2=data.getXml(10)))
        self._run_cmd_with_output(cmd, ofn)
        assert mtime != os.path.getmtime(ofn)

    def test_create_cli_automatic_type(self):
        log.debug("No type specified")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = "dataset create {o} {i1} {i2}".format(
            o=ofn, i1=data.getXml(7), i2=data.getXml(10))
        self._run_cmd_with_output(cmd, ofn)
        aset = AlignmentSet(ofn)
        shutil.rmtree(outdir)

    @pytest.mark.constools
    def test_create_cli_generate_indices(self):
        log.debug("Generate existing indices")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = ("dataset create --type AlignmentSet "
               "--generateIndices {o} {i1} {i2}").format(
            o=ofn, i1=data.getXml(7), i2=data.getXml(10))
        self._run_cmd_with_output(cmd, ofn)
        aset = AlignmentSet(ofn, strict=True)
        shutil.rmtree(outdir)

    @pytest.mark.constools
    def test_create_cli_generate_indices_2(self):
        log.debug("Generate existing indices no type specified")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = ("dataset create "
               "--generateIndices {o} {i1} {i2}").format(
            o=ofn, i1=data.getXml(7), i2=data.getXml(10))
        self._run_cmd_with_output(cmd, ofn)
        aset = AlignmentSet(ofn, strict=True)
        shutil.rmtree(outdir)

    @pytest.mark.constools
    def test_create_cli_generate_indices_3(self):
        log.debug("Generate indices")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ifn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml",
                                          dir=outdir).name
        log.info(ifn)
        # copy the resources and xml
        aset = AlignmentSet(data.getXml(11))
        aset.copyTo(ifn)
        aset = AlignmentSet(ifn)
        # get the key resource filename
        ifn = aset.toExternalFiles()[0]
        log.info(ifn)
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = ("dataset create --type AlignmentSet "
               "--generateIndices {} {}").format(ofn, ifn)
        self._run_cmd_with_output(cmd, ofn)
        aset = AlignmentSet(ofn, strict=True)
        shutil.rmtree(outdir)

    @pytest.mark.constools
    def test_create_cli_generate_indices_4(self):
        log.debug("Generate indices, no type specified")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ifn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml",
                                          dir=outdir).name
        log.info(ifn)
        # copy the resources and xml
        aset = AlignmentSet(data.getXml(11))
        aset.copyTo(ifn)
        aset = AlignmentSet(ifn)
        # get the key resource filename
        ifn = aset.toExternalFiles()[0]
        log.info(ifn)
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = ("dataset create --generateIndices {} {}").format(ofn, ifn)
        self._run_cmd_with_output(cmd, ofn)
        aset = AlignmentSet(ofn, strict=True)
        shutil.rmtree(outdir)

    def test_create_cli_reference_fasta_fname(self):
        log.debug("Include reference fasta fname")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        ofn = self._get_mock_alignment_set_out(outdir)
        cmd = ("dataset create "
               "--reference-fasta-fname {r} {o} {i1} {i2}").format(
            o=ofn, i1=data.getXml(7), i2=data.getXml(10),
            r=otherdata.getFasta())

        def _run_and_validate(args, file_name):
            self._run_cmd_with_output(cmd, ofn)
            aset = AlignmentSet(ofn, strict=True)
            for res in aset.externalResources:
                assert res.reference == otherdata.getFasta()
        _run_and_validate(cmd, ofn)
        cmd2 = cmd + " --type AlignmentSet"
        os.remove(ofn)
        _run_and_validate(cmd2, ofn)
        shutil.rmtree(outdir)

    @pytest.mark.constools
    def test_create_cli_reference_fasta(self):
        tmp_dir = tempfile.mkdtemp(suffix="dataset-unittest")
        fasta = os.path.join(tmp_dir, "reference.fasta")
        with open(fasta, "w") as fasta_out:
            fasta_out.write(">chr1\nacgtacgtacgt")
        ref_xml = os.path.join(tmp_dir, "test.referenceset.xml")
        cmd = "dataset create {d} {f} --generateIndices --type ReferenceSet --name test_reference_name --organism test_reference_organism --ploidy octaploid".format(
            d=ref_xml, f=fasta)
        self._run_cmd_with_output(cmd, ref_xml)
        ref = ReferenceSet(ref_xml)
        assert ref.metadata.organism == "test_reference_organism"

    @pytest.mark.constools
    def test_loadstats_cli(self):
        outfile = os.path.join(
            tempfile.mkdtemp(suffix="dataset-unittest"),
            'withStats.alignmentset.xml')
        cmd = "dataset loadstats {d} {s} --outfile {o}".format(
            o=outfile, d=data.getXml(7), s=data.getStats())
        self._run_cmd_with_output(cmd, outfile)
        aln = AlignmentSet(outfile)
        assert aln.metadata.summaryStats

    def test_dataset_create_set_sample_names(self):
        sample_args = "--well-sample-name WELLSAMPLE --bio-sample-name BIOSAMPLE".split()
        outfile = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        cmd = " ".join(["dataset", "create", "--force", outfile,
                        pbtestdata.get_file("subreads-bam")] + sample_args)
        self._run_cmd_with_output(cmd, outfile)
        with SubreadSet(outfile) as ds:
            assert len(ds.metadata.collections) == 1
            assert ds.metadata.collections[0].wellSample.name == "WELLSAMPLE"
            assert ds.metadata.bioSamples[0].name == "BIOSAMPLE"
            assert len(ds.metadata.bioSamples) == 1
        # now with existing samples
        outfile = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        cmd = " ".join(["dataset", "create", "--force", outfile,
                        pbtestdata.get_file("barcoded-subreadset")] + sample_args)
        self._run_cmd_with_output(cmd, outfile)
        with SubreadSet(outfile) as ds:
            assert len(ds.metadata.collections) == 1
            assert ds.metadata.collections[0].wellSample.name == "WELLSAMPLE"
            biosamples = {s.name for s in ds.metadata.bioSamples}
            assert biosamples == {"BIOSAMPLE"}
