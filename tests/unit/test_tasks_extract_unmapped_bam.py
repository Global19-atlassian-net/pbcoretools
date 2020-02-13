import pytest
import tempfile
import re

from pbcommand.models.common import DataStore
from pbcommand.testkit import PbIntegrationBase
from pbcore.io import AlignmentSet, SubreadSet

from pbcoretools import bamsieve

import pbtestdata


def assert_no_reads_in_common(self, alignment_file, output_file):
    with AlignmentSet(alignment_file) as mapped:
        mapped_subreads = set(
            zip(mapped.index.holeNumber, mapped.index.qStart))
        with SubreadSet(output_file) as unmapped:
            unmapped_subreads = set(
                zip(unmapped.index.holeNumber, unmapped.index.qStart))
            assert len(unmapped_subreads) > 0
            assert len(mapped_subreads & unmapped_subreads) == 0


def _make_filtered(ds_file):
    tmp_file = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
    bamsieve.filter_reads(
        input_bam=ds_file,
        output_bam=tmp_file,
        blacklist={49050})
    return tmp_file


@pytest.mark.constools
class TestExtractUnmappedBam(PbIntegrationBase):

    def test_run_bamsieve_extract_unmapped(self):
        mapped = _make_filtered(pbtestdata.get_file("aligned-xml"))
        subreads = pbtestdata.get_file("subreads-xml")
        args = [
            "bamsieve",
            "--subreads",
            "--blacklist",
            mapped,
            subreads,
            "unmapped.subreads.bam"
        ]
        self._check_call(args)
        assert_no_reads_in_common(self, mapped, "unmapped.subreads.bam")
