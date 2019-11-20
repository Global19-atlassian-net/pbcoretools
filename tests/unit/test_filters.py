import subprocess
import tempfile
import os.path as op
import pytest

from pbcoretools.filters import (run_filter_dataset,
                                 sanitize_read_length)
from pbcoretools.DataSetEntryPoints import parse_filter_list
from pbcore.io import openDataFile, openDataSet
import pbcore.data.datasets as data

from pbcoretools import bamsieve


class TestFilterDataSet:

    def test_dataset_io_sanitizing(self):
        ssfn = data.getXml(7)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name

        # some smoke tests:
        run_filter_dataset(ssfn, ofn, "", "")
        run_filter_dataset(ssfn, ofn, 0, "")
        run_filter_dataset(ssfn, ofn, 0, "None")
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, "None", "None")
        run_filter_dataset(ssfn, ofn, None, None)
        run_filter_dataset(ssfn, ofn, 0, 0)
        run_filter_dataset(ssfn, ofn, False, False)
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, True, False)
        run_filter_dataset(ssfn, ofn, False, True)
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, True, True)
        run_filter_dataset(ssfn, ofn, 1, True)
        run_filter_dataset(ssfn, ofn, "0", "None")
        run_filter_dataset(ssfn, ofn, " 0 ", "None")
        run_filter_dataset(ssfn, ofn, "-1", "none")
        run_filter_dataset(ssfn, ofn, " -1 ", "none")
        run_filter_dataset(ssfn, ofn, -1, "none")
        run_filter_dataset(ssfn, ofn, "-1", "blah")
        run_filter_dataset(ssfn, ofn, "-1", "None AND None")
        run_filter_dataset(ssfn, ofn, "1.0", "None AND None")
        run_filter_dataset(ssfn, ofn, " 1.0 ", "None AND None")
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, " 1. 0 ", "None AND None")
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, "1. 0", "None AND None")
        run_filter_dataset(ssfn, ofn, "10.0", "None AND None")
        run_filter_dataset(ssfn, ofn, "0.01", "None AND None")
        run_filter_dataset(ssfn, ofn, "100000.01", "None AND None")
        run_filter_dataset(ssfn, ofn, ".00", "None AND None")
        run_filter_dataset(ssfn, ofn, 1.0, "None AND None")
        run_filter_dataset(ssfn, ofn, 10.0, "None AND None")
        run_filter_dataset(ssfn, ofn, 0.01, "None AND None")
        run_filter_dataset(ssfn, ofn, 100000.01, "None AND None")
        run_filter_dataset(ssfn, ofn, 100000, "None AND None")
        run_filter_dataset(ssfn, ofn, .00, "None AND None")
        run_filter_dataset(ssfn, ofn, 1., "None AND None")
        run_filter_dataset(ssfn, ofn, "1.", "None AND None")
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, "None.1", "None AND None")
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, "None1", "None AND None")
        with pytest.raises(ValueError):
            run_filter_dataset(ssfn, ofn, "1.1None", "None AND None")

        assert sanitize_read_length("") is None
        assert sanitize_read_length(0) is None
        assert sanitize_read_length(0) is None
        with pytest.raises(ValueError):
            sanitize_read_length("None")
        assert sanitize_read_length(None) is None
        assert sanitize_read_length(0) is None
        assert sanitize_read_length(False) is None
        with pytest.raises(ValueError):
            sanitize_read_length(True)
        assert sanitize_read_length(False) is None
        with pytest.raises(ValueError):
            sanitize_read_length(True)
        assert 1 == sanitize_read_length(1)
        assert 0 == sanitize_read_length("0")
        assert 0 == sanitize_read_length(" 0 ")
        assert -1 == sanitize_read_length("-1")
        assert -1 == sanitize_read_length(" -1 ")
        assert -1 == sanitize_read_length(-1)
        assert 1 == sanitize_read_length("1.0")
        assert 1 == sanitize_read_length(" 1.0 ")
        assert 1 == sanitize_read_length(1.0)
        with pytest.raises(ValueError):
            sanitize_read_length(" 1. 0 ")
        with pytest.raises(ValueError):
            sanitize_read_length("1. 0")
        assert 10 == sanitize_read_length("10.0")
        assert 10 == sanitize_read_length(10.0)
        assert 0 == sanitize_read_length("0.01")
        assert 0 == sanitize_read_length(0.01)
        assert 100000 == sanitize_read_length("100000.01")
        assert 100000 == sanitize_read_length(100000.01)
        assert 100000 == sanitize_read_length(100000)
        # These two are somewhat annoying:
        assert 0 == sanitize_read_length(".00")
        assert sanitize_read_length(.00) is None
        assert 1 == sanitize_read_length(1.)
        assert 1 == sanitize_read_length("1.")
        with pytest.raises(ValueError):
            sanitize_read_length("None.1")
        with pytest.raises(ValueError):
            sanitize_read_length("None1")
        with pytest.raises(ValueError):
            sanitize_read_length("1.1None")

    def test_filter_application(self):
        ssfn = data.getXml(7)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name

        # some smoke tests:
        run_filter_dataset(ssfn, ofn, "0", "None")
        ds = openDataSet(ofn)
        assert len(ds) == 92
        run_filter_dataset(ssfn, ofn, "10000", "None")
        ds = openDataSet(ofn)
        assert len(ds) == 0

        run_filter_dataset(ssfn, ofn, "100", "rq > .7")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( rq > .7 AND length >= 100 )"

        # AND conjunction
        run_filter_dataset(ssfn, ofn, "100", "length < 5000 AND rq > .7")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( length < 5000 AND rq > .7 AND length >= 100 )"

        # semicolon conjunction

        run_filter_dataset(ssfn, ofn, "100", "length < 5000; rq > .7")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( length < 5000 AND rq > .7 AND length >= 100 )"

        run_filter_dataset(ssfn, ofn, 0,
                           "length >= 1000; length <= 5000; rq >= .7")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( length >= 1000 AND length <= 5000 AND rq >= .7 )"

        run_filter_dataset(ssfn, ofn, 0,
                           "length gte 1000; length lte 5000; rq >= .7")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( length gte 1000 AND length lte 5000 AND rq >= .7 )"

    def test_filter_comma_raises(self):
        with pytest.raises(ValueError):
            ssfn = data.getXml(7)
            ofn = tempfile.NamedTemporaryFile(suffix=".xml").name
            run_filter_dataset(ssfn, ofn, "100", "rq > .7, length < 5000")

    def test_filter_more(self):
        ssfn = data.getXml(7)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name

        # zm=[3,4,5] condition
        run_filter_dataset(ssfn, ofn, 0,
                           "zm=[3,4,5] AND length >= 1000")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( zm = [3,4,5] AND length >= 1000 )"

        # zm=[3,4,5] condition
        run_filter_dataset(ssfn, ofn, 0, "zm=[3,4,5]; length >= 1000")
        ds = openDataSet(ofn)
        assert str(ds.filters) == "( zm = [3,4,5] AND length >= 1000 )"

        # zm=[3,4,5] condition by itself
        run_filter_dataset(ssfn, ofn, 0,
                           "zm=[3,4,5]")
        ds = openDataSet(ofn)
        assert str(ds.filters) == '( zm = [3,4,5] )'

    def test_datset_name(self):
        ssfn = data.getXml(7)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name
        run_filter_dataset(ssfn, ofn, "0", "None")
        ds = openDataSet(ofn)
        assert ds.name.endswith("(filtered)")
        assert "filtered" in ds.tags

    def test_parse_filter_list(self):
        f1 = ["rq >= 0.7 AND length <= 10000"]
        f2 = parse_filter_list(f1)
        assert f2 == {"rq": [(">=", "0.7")], "length": [("<=", "10000")]}
