
import os
import re
import logging
import itertools
import tempfile

import numpy as np
import unittest
import shutil
from unittest.case import SkipTest

from pbcore.io import PacBioBamIndex, IndexedBamReader
from pbcore.io import openIndexedAlignmentFile
from pbcore.io.dataset.utils import BamtoolsVersion
from pbcore.io import (DataSet, SubreadSet, ReferenceSet, AlignmentSet,
                       openDataSet, DataSetMetaTypes, HdfSubreadSet,
                       ConsensusReadSet, ConsensusAlignmentSet)
from pbcore.io.dataset.DataSetIO import _dsIdToSuffix
from pbcore.io.dataset.DataSetMembers import ExternalResource, Filters
from pbcore.io.dataset.DataSetWriter import toXml
from pbcore.io.dataset.DataSetValidator import validateFile
from pbcore.util.Process import backticks
import pbcore.data.datasets as data
import pbcore.data as upstreamdata
import xml.etree.cElementTree as ET

log = logging.getLogger(__name__)

def _is_relative(xmlfile):
    for event, element in ET.iterparse(xmlfile, events=("start",)):
        if element.get("ResourceId", "noop").startswith(os.sep):
            return False
    return True

def _check_constools():
    cmd = "dataset"
    o, r, m = backticks(cmd)
    if r != 2:
        return False

    if not BamtoolsVersion().good:
        log.warn("Bamtools not found or out of date")
        return False

    cmd = "pbindex"
    o, r, m = backticks(cmd)
    if r != 1:
        return False

    cmd = "samtools"
    o, r, m = backticks(cmd)
    if r != 1:
        return False
    return True

def _internal_data():
    if os.path.exists("/pbi/dept/secondary/siv/testdata"):
        return True
    return False

class TestDataSet(unittest.TestCase):
    """Unit and integrationt tests for the DataSet class and \
    associated module functions"""

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_split_cli(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset split --outdir {o} --contigs --chunks 2 {d}".format(
            o=outdir,
            d=data.getXml(8))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(
            os.path.join(outdir, 'pbalchemysim0.chunk0.alignmentset.xml')))
        self.assertTrue(os.path.exists(
            os.path.join(outdir, 'pbalchemysim0.chunk1.alignmentset.xml')))


    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_copyTo_cli(self):
        # To a fname:
        # absolute:
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        cmd = "dataset copyto {i} {o}".format(i=data.getXml(8), o=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        sset = AlignmentSet(fn, strict=True)
        self.assertFalse(_is_relative(fn))

        # relative:
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        cmd = "dataset copyto --relative {i} {o}".format(i=data.getXml(8), o=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        sset = AlignmentSet(fn, strict=True)
        self.assertTrue(_is_relative(fn))

        # to a directory:
        # absolute:
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, os.path.split(data.getXml(8))[1])
        cmd = "dataset copyto {i} {o}".format(i=data.getXml(8), o=outdir)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        sset = AlignmentSet(fn, strict=True)
        self.assertFalse(_is_relative(fn))

        # relative:
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, os.path.split(data.getXml(8))[1])
        cmd = "dataset copyto --relative {i} {o}".format(i=data.getXml(8),
                                                         o=outdir)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        sset = AlignmentSet(fn, strict=True)
        self.assertTrue(_is_relative(fn))

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_newUuid_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(8))
        aln.copyTo(fn)
        pre_uuid = AlignmentSet(fn).uuid
        cmd = "dataset newuuid {d}".format(d=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        post_uuid = AlignmentSet(fn).uuid
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        self.assertNotEqual(pre_uuid, post_uuid)

    def test_newUuid_clone_cli(self):
        fn_orig = data.getXml(8)
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, 'fn.alignmentset.xml')
        fn2 = os.path.join(outdir, 'fn2.alignmentset.xml')
        with AlignmentSet(fn_orig) as aln:
            aln.copyTo(fn)
            shutil.copy(fn, fn2)

        pre_uuid = AlignmentSet(fn).uuid
        pre_uuid2 = AlignmentSet(fn2).uuid
        self.assertEqual(pre_uuid, pre_uuid2)

        cmd = "dataset newuuid {d}".format(d=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))

        cmd = "dataset newuuid {d}".format(d=fn2)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn2))

        post_uuid = AlignmentSet(fn).uuid
        post_uuid2 = AlignmentSet(fn2).uuid
        self.assertNotEqual(pre_uuid, post_uuid)
        self.assertNotEqual(pre_uuid2, post_uuid2)
        # HASH, THEREFORE THESE ARE EQUAL:
        self.assertEqual(post_uuid, post_uuid2)

    def test_newUuid_random_cli(self):
        fn_orig = data.getXml(8)
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        fn = os.path.join(outdir, 'fn.alignmentset.xml')
        fn2 = os.path.join(outdir, 'fn2.alignmentset.xml')
        with AlignmentSet(fn_orig) as aln:
            aln.copyTo(fn)
            shutil.copy(fn, fn2)

        pre_uuid = AlignmentSet(fn).uuid
        pre_uuid2 = AlignmentSet(fn2).uuid
        self.assertEqual(pre_uuid, pre_uuid2)

        cmd = "dataset newuuid --random {d}".format(d=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))

        cmd = "dataset newuuid --random {d}".format(d=fn2)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn2))

        post_uuid = AlignmentSet(fn).uuid
        post_uuid2 = AlignmentSet(fn2).uuid
        self.assertNotEqual(pre_uuid, post_uuid)
        self.assertNotEqual(pre_uuid2, post_uuid2)
        # RANDOM, THEREFORE THESE ARE NOT EQUAL:
        self.assertNotEqual(post_uuid, post_uuid2)

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_relativize_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(8))
        aln.copyTo(fn)
        self.assertFalse(_is_relative(fn))
        cmd = "dataset relativize {d}".format(d=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        self.assertTrue(_is_relative(fn))

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_absolutize_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(8))
        aln.copyTo(fn, relative=True)
        self.assertTrue(_is_relative(fn))
        cmd = "dataset absolutize {d}".format(d=fn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        self.assertFalse(_is_relative(fn))

        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        outfn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        aln = AlignmentSet(data.getXml(8))
        aln.copyTo(fn, relative=True)
        self.assertTrue(_is_relative(fn))
        cmd = "dataset absolutize {d} --outdir {o}".format(d=fn, o=outfn)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        self.assertTrue(os.path.exists(outfn))
        self.assertTrue(_is_relative(fn))
        self.assertFalse(_is_relative(outfn))

        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, os.path.split(fn)[1])
        aln = AlignmentSet(data.getXml(8))
        aln.copyTo(fn, relative=True)
        self.assertTrue(_is_relative(fn))
        cmd = "dataset absolutize {d} --outdir {o}".format(d=fn, o=outdir)
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(fn))
        self.assertTrue(os.path.exists(outfn))
        self.assertTrue(_is_relative(fn))
        self.assertFalse(_is_relative(outfn))

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_loadmetadata_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        log.debug(fn)

        aln = AlignmentSet(data.getXml(8))
        aln.metadata.collections = None
        aln.copyTo(fn)
        aln.close()
        del aln
        self.assertTrue(os.path.exists(fn))

        aln = AlignmentSet(fn)
        self.assertFalse(aln.metadata.collections)

        cmd = "dataset loadmetadata {i} {m}".format(
            i=fn,
            m=("/pbi/dept/secondary/siv/testdata/"
               "SA3-Sequel/lambda/roche_SAT/"
               "m54013_151205_032353.run.metadata.xml"))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0, m)
        aln = AlignmentSet(fn)
        self.assertTrue(aln.metadata.collections)

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_loadmetadata_from_dataset_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        log.debug(fn)

        aln = AlignmentSet(data.getXml(8))
        aln.metadata.collections = None
        aln.copyTo(fn)
        aln.close()
        del aln
        self.assertTrue(os.path.exists(fn))

        aln = AlignmentSet(fn)
        self.assertFalse(aln.metadata.collections)

        cmd = "dataset loadmetadata {i} {m}".format(
            i=fn,
            m=("/pbi/dept/secondary/siv/testdata/"
               "SA3-Sequel/lambda/roche_SAT/"
               "m54013_151205_032353.subreadset.xml"))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0, m)
        aln = AlignmentSet(fn)
        self.assertTrue(aln.metadata.collections)

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_loadmetadata_from_dataset_create_cli(self):
        fn = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        fn2 = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
        log.debug(fn)

        aln = AlignmentSet(data.getXml(8))
        aln.metadata.collections = None
        aln.copyTo(fn)
        aln.close()
        del aln
        self.assertTrue(os.path.exists(fn))

        aln = AlignmentSet(fn)
        self.assertFalse(aln.metadata.collections)

        cmd = "dataset create --metadata {m} {o} {i}".format(
            o=fn2,
            i=fn,
            m=("/pbi/dept/secondary/siv/testdata/"
               "SA3-Sequel/lambda/roche_SAT/"
               "m54013_151205_032353.subreadset.xml"))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0, m)
        aln = AlignmentSet(fn2)
        self.assertTrue(aln.metadata.collections)

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_contigset_split_cli(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset split --outdir {o} --chunks 2 {d}".format(
            o=outdir,
            d=data.getXml(9))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(
            os.path.join(outdir,
                         "pbalchemysim0.chunk0.referenceset.xml")))
        self.assertTrue(os.path.exists(
            os.path.join(outdir,
                         "pbalchemysim0.chunk1.referenceset.xml")))

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_filter_cli(self):
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, "filtered8.xml")
        log.debug(outfn)
        cmd = "dataset filter {i} {o} {f}".format(
            i=data.getXml(8),
            o=outfn,
            f="rname=E.faecalis.1")
        log.debug(cmd)
        o, r, m = backticks(cmd)
        if r != 0:
            log.debug(m)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(outfn))
        aln = AlignmentSet(data.getXml(8))
        aln.filters.addRequirement(rname=[('=', 'E.faecalis.1')])
        aln.updateCounts()
        dset = AlignmentSet(outfn)
        self.assertEqual(str(aln.filters), str(dset.filters))
        self.assertEqual(aln.totalLength, dset.totalLength)
        self.assertEqual(aln.numRecords, dset.numRecords)

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_create_cli(self):
        log.debug("Absolute")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = "dataset create --type AlignmentSet {o} {i1} {i2}".format(
            o=os.path.join(outdir, 'pbalchemysim.alignmentset.xml'),
            i1=data.getXml(8), i2=data.getXml(11))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(
            os.path.join(outdir, os.path.basename(data.getXml(12)))))

        log.debug("Relative")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        cmd = ("dataset create --relative --type AlignmentSet "
               "{o} {i1} {i2}".format(
                   o=os.path.join(outdir, 'pbalchemysim.alignmentset.xml'),
                   i1=data.getXml(8),
                   i2=data.getXml(11)))
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(
            os.path.join(outdir, os.path.basename(data.getXml(12)))))

    @unittest.skipIf(not _check_constools(),
                     "bamtools or pbindex not found, skipping")
    def test_loadstats_cli(self):
        outfile = os.path.join(
            tempfile.mkdtemp(suffix="dataset-unittest"),
            'withStats.alignmentset.xml')
        cmd = "dataset loadstats {d} {s} --outfile {o}".format(
            o=outfile, d=data.getXml(8), s=data.getStats())
        log.debug(cmd)
        o, r, m = backticks(cmd)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(outfile))
        aln = AlignmentSet(outfile)
        self.assertTrue(aln.metadata.summaryStats)
