import xml.etree.ElementTree as ET
import logging
import tempfile
import os
import itertools
import pytest

from urllib.parse import urlparse

from pbcore.io.dataset.utils import _infixFname
from pbcore.io import (DataSet, SubreadSet, ConsensusReadSet,
                       ReferenceSet, ContigSet, AlignmentSet,
                       FastaReader, FastaWriter, IndexedFastaReader,
                       ConsensusAlignmentSet,
                       openDataFile, FastaWriter)
import pbcore.data.datasets as data
from pbcore.io.dataset.DataSetValidator import validateXml

import pbtestdata

log = logging.getLogger(__name__)


class TestDataSet:
    """Unit and integrationt tests for the DataSet class and \
    associated module functions"""

    @pytest.mark.constools
    def test_alignmentset_consolidate(self):

        log.debug("Test through API")
        aln = AlignmentSet(pbtestdata.get_file("aligned-ds-2"))
        assert len(aln.toExternalFiles()) == 2
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, 'merged.bam')
        aln.consolidate(outfn)
        assert os.path.exists(outfn)
        assert len(aln.toExternalFiles()) == 1
        nonCons = AlignmentSet(pbtestdata.get_file("aligned-ds-2"))
        assert len(nonCons.toExternalFiles()) == 2
        for read1, read2 in zip(sorted(list(aln)), sorted(list(nonCons))):
            assert read1 == read2
        assert len(aln) == len(nonCons)

        # Test that it is a valid xml:
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        datafile = os.path.join(outdir, "apimerged.bam")
        xmlfile = os.path.join(outdir, "apimerged.xml")
        log.debug(xmlfile)
        aln.write(xmlfile)

        log.debug("Test with cheap filter")
        aln = AlignmentSet(pbtestdata.get_file("aligned-ds-2"))
        assert len(list(aln)) == 21
        aln.filters.addRequirement(length=[(">=", 10000)])
        assert len(list(aln)) == 10
        assert len(aln.toExternalFiles()) == 2
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, 'merged.bam')
        aln.consolidate(outfn)
        assert os.path.exists(outfn)
        assert len(aln.toExternalFiles()) == 1
        nonCons = AlignmentSet(pbtestdata.get_file("aligned-ds-2"))
        nonCons.filters.addRequirement(length=[(">=", 10000)])
        assert len(nonCons.toExternalFiles()) == 2
        for read1, read2 in zip(sorted(list(aln)), sorted(list(nonCons))):
            assert read1 == read2
        assert len(list(aln)) == len(list(nonCons))

        log.debug("Test with not refname filter")
        # This isn't trivial with bamtools
        """
        aln = AlignmentSet(data.getXml(11))
        assert len(list(aln)) == 177
        aln.filters.addRequirement(rname=[('!=', 'B.vulgatus.5')])
        assert len(list(aln)) == 7
        assert len(aln.toExternalFiles()) == 2
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, 'merged.bam')
        aln.consolidate(outfn)
        assert os.path.exists(outfn)
        assert len(aln.toExternalFiles()) == 1
        nonCons = AlignmentSet(data.getXml(11))
        nonCons.filters.addRequirement(rname=[('!=', 'B.vulgatus.5')])
        assert len(nonCons.toExternalFiles()) == 2
        for read1, read2 in zip(sorted(list(aln)), sorted(list(nonCons))):
            assert read1 == read2
        assert len(list(aln)) == len(list(nonCons))
        """

        log.debug("Test with expensive filter")
        aln = AlignmentSet(data.getXml(11))
        assert len(list(aln)) == 177
        aln.filters.addRequirement(accuracy=[('>', '.85')])
        assert len(list(aln)) == 174
        assert len(aln.toExternalFiles()) == 2
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, 'merged.bam')
        aln.consolidate(outfn)
        assert os.path.exists(outfn)
        assert len(aln.toExternalFiles()) == 1
        nonCons = AlignmentSet(data.getXml(11))
        nonCons.filters.addRequirement(accuracy=[('>', '.85')])
        assert len(nonCons.toExternalFiles()) == 2
        for read1, read2 in zip(sorted(list(aln)), sorted(list(nonCons))):
            assert read1 == read2
        assert len(list(aln)) == len(list(nonCons))

        log.debug("Test cli")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        datafile = os.path.join(outdir, "merged.bam")
        xmlfile = os.path.join(outdir, "merged.xml")
        cmd = "dataset consolidate {i} {d} {x}".format(i=data.getXml(11),
                                                       d=datafile,
                                                       x=xmlfile)
        log.debug(cmd)
        subprocess.check_call(cmd.split())

    @pytest.mark.skip(reason="DISABLED FOR AUTOMATED TESTING")
    @pytest.mark.constools
    def test_alignmentset_partial_consolidate(self):
        testFile = ("/pbi/dept/secondary/siv/testdata/SA3-DS/"
                    "lambda/2372215/0007_tiny/Alignment_"
                    "Results/m150404_101626_42267_c10080"
                    "7920800000001823174110291514_s1_p0."
                    "all.alignmentset.xml")
        aln = AlignmentSet(testFile)
        nonCons = AlignmentSet(testFile)
        assert len(aln.toExternalFiles()) == 3
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        outfn = os.path.join(outdir, 'merged.bam')
        aln.consolidate(outfn, numFiles=2)
        assert not os.path.exists(outfn)
        assert os.path.exists(_infixFname(outfn, "0"))
        assert os.path.exists(_infixFname(outfn, "1"))
        assert len(aln.toExternalFiles()) == 2
        assert len(nonCons.toExternalFiles()) == 3
        for read1, read2 in zip(sorted(list(aln)), sorted(list(nonCons))):
            assert read1 == read2
        assert len(aln) == len(nonCons)

        log.debug("Test cli")
        outdir = tempfile.mkdtemp(suffix="dataset-unittest")
        datafile = os.path.join(outdir, "merged.bam")
        xmlfile = os.path.join(outdir, "merged.xml")
        cmd = "dataset consolidate --numFiles 2 {i} {d} {x}".format(
            i=testFile, d=datafile, x=xmlfile)
        log.debug(cmd)
        subprocess.check_call(cmd.split())
