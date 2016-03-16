
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

log = logging.getLogger(__name__)

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
