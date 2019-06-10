
import unittest
import tempfile
import logging
import uuid
import os.path as op

from pbcore.io import (FastaReader, FastqReader, openDataSet,
                       SubreadSet, ConsensusReadSet, TranscriptSet)
import pbcore.data.datasets as data

from pbcoretools.filters import combine_filters, run_filter_dataset, sanitize_read_length

from base import get_temp_file, IntegrationBase

import pbtestdata

log = logging.getLogger(__name__)

class TestFilterDataSet(IntegrationBase):
    BASE_ARGS = [
        "python", "-m", "pbcoretools.tasks2.dataset_filter"
    ]
    READER_CLASS = SubreadSet

    def _set_up_basic(self):
        input_file = get_temp_file(suffix=".subreadset.xml")
        ds = SubreadSet(data.getXml(10), strict=True)
        ds.metadata.addParentDataSet(uuid.uuid4(),
                                     ds.datasetType,
                                     createdBy="AnalysisJob",
                                     timeStampedName="")
        ds.write(input_file)
        return input_file, len(ds)

    def _set_up_combine_filters(self):
        ds_in = get_temp_file(suffix=".subreadset.xml")
        with SubreadSet(pbtestdata.get_file("subreads-xml"), strict=True) as ds:
            assert len(ds) == 117 and len(ds.filters) == 0
            ds.filters.addRequirement(length=[('>=', 1000)])
            assert len(ds) == 13
            ds.write(ds_in)
        return ds_in, 13

    def _get_counts(self, output_file):
        with self.READER_CLASS(output_file) as f:
            return len(f)

    def _get_filters(self, output_file):
        with self.READER_CLASS(output_file) as f:
            return str(f.filters)

    def run_after(self, output_file, n_expected, expected_filter_str):
        n_actual = self._get_counts(output_file)
        self.assertEqual(self._get_filters(output_file), expected_filter_str)
        self.assertEqual(n_actual, n_expected)
        ds = openDataSet(output_file)
        self.assertEqual(len(ds.metadata.provenance), 0)
        self.assertTrue(ds.name.endswith("(filtered)"))
        self.assertTrue("filtered" in ds.tags)
        return ds

    def test_filter_dataset(self):
        ds_in, n_input = self._set_up_basic()
        ds_out = get_temp_file(suffix=".subreadset.xml")
        args = self.BASE_ARGS + [ds_in, ds_out, "length <= 1400"]
        self._check_call(args)
        n_expected = 18
        expected_filter_str = "( length <= 1400 )"
        self.run_after(ds_out, n_expected, expected_filter_str)

    def test_filter_dataset_nofilter(self):
        ds_in, n_input = self._set_up_basic()
        ds_out = get_temp_file(suffix=".subreadset.xml")
        args = self.BASE_ARGS + [ds_in, ds_out, ""]
        self._check_call(args)
        n_expected = n_input
        expected_filter_str = ""
        self.run_after(ds_out, n_expected, expected_filter_str)

    def test_filter_dataset_bq(self):
        ds_in = get_temp_file(suffix=".subreadset.xml")
        ds = SubreadSet(pbtestdata.get_file("barcoded-subreadset"),
                        strict=True)
        ds.filters.addRequirement(bq=[('>=', 31)])
        assert len(ds) == 1
        ds.write(ds_in)
        ds_out = get_temp_file(suffix=".subreadset.xml")
        args = self.BASE_ARGS + [ds_in, ds_out, "length >= 10 AND bq >= 10"]
        self._check_call(args)
        n_expected = 2
        expected_filter_str = "( bq >= 10 AND length >= 10 )"
        self.run_after(ds_out, n_expected, expected_filter_str)

    def test_filter_dataset_downsample(self):
        ds_in = get_temp_file(suffix=".subreadset.xml")
        with SubreadSet(pbtestdata.get_file("subreads-xml"), strict=True) as ds:
            assert len(ds) == 117 and len(ds.filters) == 0
            ds.write(ds_in)
        ds_out = get_temp_file(suffix=".subreadset.xml")
        args = self.BASE_ARGS + [ds_in, ds_out, "", "--downsample", "2"]
        self._check_call(args)
        n_expected = 54
        expected_filter_str = "( Uint32Cast(zm) % 2 == 0 )"
        ds = self.run_after(ds_out, n_expected, expected_filter_str)
        self.assertTrue("downsampled" in ds.tags)

    def test_filter_dataset_combine_filters(self):
        ds_in, n_input = self._set_up_combine_filters()
        ds_out = get_temp_file(suffix=".subreadset.xml")
        args = self.BASE_ARGS + [ds_in, ds_out, "rq >= 0.901"]
        self._check_call(args)
        n_expected = 18
        expected_filter_str = "( length >= 1000 AND rq >= 0.901 )"
        self.run_after(ds_out, 12, expected_filter_str)

    def test_combine_filters(self):
        ds_in, n_input = self._set_up_combine_filters()
        with openDataSet(ds_in, strict=True) as ds:
            filters = {"rq": [(">=", 0.901)]}
            combine_filters(ds, filters)
            ds.reFilter(light=False)
            self.assertEqual(len(ds), 12)
            filters = {"rq": [(">=", 0.8)], "length": [(">=", 500)]}
            combine_filters(ds, filters)
            ds.reFilter(light=False)
            self.assertEqual(len(ds), 48)
            ds.filters = None
            filters = {"zm": [("==", "0", 2)]}
            combine_filters(ds, filters)
            ds.reFilter(light=False)
            self.assertEqual(len(ds), 54)

    def test_combine_filters_run_filter_dataset(self):
        ds_in, n_input = self._set_up_combine_filters()
        ds_out = get_temp_file(suffix=".subreadset.xml")
        my_filters = "rq >= 0.901"
        run_filter_dataset(ds_in, ds_out, 0, my_filters)
        with openDataSet(ds_out, strict=True) as ds:
            self.assertEqual(len(ds), 12)
        my_filters = "rq >= 0.8"
        run_filter_dataset(ds_in, ds_out, 500, my_filters)
        with openDataSet(ds_out, strict=True) as ds:
            self.assertEqual(len(ds), 48)
