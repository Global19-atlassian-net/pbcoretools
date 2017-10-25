
import logging
import os.path as op

from pbcore.io import (FastaReader, FastqReader, openDataSet, HdfSubreadSet,
                       SubreadSet, ConsensusReadSet)
import pbcore.data.datasets as data
from pbcommand.testkit import PbTestApp

from base import get_temp_file

import pbtestdata

log = logging.getLogger(__name__)

DATA = op.join(op.dirname(__file__), "data")

SIV_DATA_DIR = "/pbi/dept/secondary/siv/testdata"

def _to_skip_msg(exe):
    return "Missing {e} or {d}".format(d=SIV_DATA_DIR, e=exe)

class TestFilterDataSet(PbTestApp):
    TASK_ID = "pbcoretools.tasks.filterdataset"
    DRIVER_EMIT = 'python -m pbcoretools.tasks.filters emit-tool-contract {i} '.format(i=TASK_ID)
    DRIVER_RESOLVE = 'python -m pbcoretools.tasks.filters run-rtc '
    TASK_OPTIONS = {"pbcoretools.task_options.other_filters": "length <= 1400"}
    RESOLVED_TASK_OPTIONS = {
        "pbcoretools.task_options.other_filters": "length <= 1400"}
    INPUT_FILES = [get_temp_file(suffix=".subreadset.xml")]
    MAX_NPROC = 24
    RESOLVED_NPROC = 1
    IS_DISTRIBUTED = True
    RESOLVED_IS_DISTRIBUTED = True
    READER_CLASS = SubreadSet
    N_EXPECTED = 18
    EXPECTED_FILTER_STR = "( length <= 1400 )"

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(data.getXml(10), strict=True)
        ds.write(cls.INPUT_FILES[0])

    def _get_counts(self, rtc):
        n_expected = self.N_EXPECTED
        with self.READER_CLASS(rtc.task.output_files[0]) as f:
            n_actual = len(f)
        return n_expected, n_actual

    def _get_filters(self, rtc):
        with self.READER_CLASS(rtc.task.output_files[0]) as f:
            return str(f.filters)

    def run_after(self, rtc, output_dir):
        n_expected, n_actual = self._get_counts(rtc)
        self.assertEqual(self._get_filters(rtc), self.EXPECTED_FILTER_STR)
        self.assertEqual(n_actual, n_expected)
        ds = openDataSet(rtc.task.output_files[0])
        self.assertTrue(ds.name.endswith("(filtered)"))
        self.assertTrue("filtered" in ds.tags)


class TestFilterDataSetBq(TestFilterDataSet):

    TASK_OPTIONS = {"pbcoretools.task_options.other_filters": "length >= 10 AND bq >= 10"}
    RESOLVED_TASK_OPTIONS = TASK_OPTIONS
    N_EXPECTED = 2
    EXPECTED_FILTER_STR = "( length >= 10 AND bq >= 10 )"

    @classmethod
    def setUpClass(cls):
        ds = SubreadSet(pbtestdata.get_file("barcoded-subreadset"),
                        strict=True)
        ds.filters.addRequirement(bq=[('>', 30)])
        assert len(ds) == 1
        ds.write(cls.INPUT_FILES[0])
