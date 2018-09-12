
import tempfile
import unittest
import re

from pbcommand.testkit.core import PbTestApp
from pbcommand.models.common import DataStore
from pbcore.io import AlignmentSet, SubreadSet

from pbcoretools.tasks.extract_unmapped_bam import (run_extract_unmapped,
                                                    make_unmapped_bam)
from pbcoretools import bamSieve

import pbtestdata


def assert_empty_datastore(self, file_name):
    ds = DataStore.load_from_json(file_name)
    self.assertEqual(len(ds.files), 0)


def validate_outputs(self, output_file, alignment_file):
    ds = DataStore.load_from_json(output_file)
    self.assertEqual(len(ds.files), 1)
    unmapped_file = list(ds.files.values())[0].path
    return assert_no_reads_in_common(self, alignment_file, unmapped_file)


def assert_no_reads_in_common(self, alignment_file, output_file):
    with AlignmentSet(alignment_file) as mapped:
        mapped_zmws = set(mapped.index.holeNumber)
        with SubreadSet(output_file) as unmapped:
            unmapped_zmws = set(unmapped.index.holeNumber)
            self.assertTrue(len(unmapped_zmws) > 0)
            self.assertEqual(len(mapped_zmws & unmapped_zmws), 0)


def _make_filtered(ds_file):
    tmp_file = tempfile.NamedTemporaryFile(suffix=".alignmentset.xml").name
    bamSieve.filter_reads(
        input_bam=ds_file,
        output_bam=tmp_file,
        blacklist={49050})
    return tmp_file


class TestExtractUnmappedBam(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mapped = _make_filtered(pbtestdata.get_file("aligned-xml"))

    def test_get_blacklist(self):
        from pbcoretools.tasks.extract_unmapped_bam import _get_blacklist, _get_blacklist_pbi
        subreads = pbtestdata.get_file("subreads-sequel")
        pbi = re.sub(".subreadset.xml", ".subreads.bam.pbi", subreads)
        ds = SubreadSet(subreads, skipCounts=True)
        blacklist = sorted(_get_blacklist(ds))
        self.assertEqual(blacklist, [(-2081539485, 5177614), (-2081539485, 6160775)])
        blacklist = sorted(_get_blacklist_pbi(pbi))
        self.assertEqual(blacklist, [(-2081539485, 5177614), (-2081539485, 6160775)])

    def test_make_unmapped_bam(self):
        subreads = pbtestdata.get_file("subreads-xml")
        output_bam = tempfile.NamedTemporaryFile(suffix=".subreads.bam").name
        make_unmapped_bam(self.mapped, subreads, output_bam)
        assert_no_reads_in_common(self, self.mapped, output_bam)

    def test_run_extract_unmapped(self):
        subreads = pbtestdata.get_file("subreads-xml")
        output_ds = tempfile.NamedTemporaryFile(suffix=".datastore.json").name
        run_extract_unmapped(self.mapped, subreads, output_ds)
        validate_outputs(self, output_ds, self.mapped)

    def test_run_extract_unmapped_no_output(self):
        subreads = pbtestdata.get_file("subreads-xml")
        output_ds = tempfile.NamedTemporaryFile(suffix=".datastore.json").name
        run_extract_unmapped(self.mapped, subreads, output_ds, False)
        assert_empty_datastore(self, output_ds)


class TestExtractUnalignedTCDefaults(PbTestApp):
    DRIVER_BASE = "python -m pbcoretools.tasks.extract_unmapped_bam"
    INPUT_FILES = [
        _make_filtered(pbtestdata.get_file("aligned-xml")),
        pbtestdata.get_file("subreads-xml")
    ]

    def run_after(self, rtc, output_dir):
        assert_empty_datastore(self, rtc.task.output_files[0])


class TestExtractUnalignedTC(TestExtractUnalignedTCDefaults):
    TASK_OPTIONS = {
        "pbcoretools.task_options.output_unaligned_bam": True
    }

    def run_after(self, rtc, output_dir):
        validate_outputs(self, rtc.task.output_files[0], self.INPUT_FILES[0])
