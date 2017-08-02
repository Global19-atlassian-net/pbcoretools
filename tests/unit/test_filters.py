
import subprocess
import tempfile
import unittest
import os.path as op

from pbcoretools.tasks.filters import (run_filter_dataset,
                                       sanitize_read_length)
from pbcore.io import openDataFile, openDataSet
import pbcore.data.datasets as data

from pbcoretools import bamSieve

class TestFilterDataSet(unittest.TestCase):

    def test_dataset_io_sanitizing(self):
        ssfn = data.getXml(8)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name

        # some smoke tests:
        run_filter_dataset(ssfn, ofn, "", "")
        run_filter_dataset(ssfn, ofn, 0, "")
        run_filter_dataset(ssfn, ofn, 0, "None")
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, "None", "None")
        run_filter_dataset(ssfn, ofn, None, None)
        run_filter_dataset(ssfn, ofn, 0, 0)
        run_filter_dataset(ssfn, ofn, False, False)
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, True, False)
        run_filter_dataset(ssfn, ofn, False, True)
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, True, True)
        run_filter_dataset(ssfn, ofn, 1, True)
        run_filter_dataset(ssfn, ofn, "0", "None")
        run_filter_dataset(ssfn, ofn, " 0 ", "None")
        run_filter_dataset(ssfn, ofn, "-1", "none")
        run_filter_dataset(ssfn, ofn, " -1 ", "none")
        run_filter_dataset(ssfn, ofn, -1, "none")
        run_filter_dataset(ssfn, ofn, "-1", "blah")
        run_filter_dataset(ssfn, ofn, "-1", "None, None")
        run_filter_dataset(ssfn, ofn, "1.0", "None, None")
        run_filter_dataset(ssfn, ofn, " 1.0 ", "None, None")
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, " 1. 0 ", "None, None")
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, "1. 0", "None, None")
        run_filter_dataset(ssfn, ofn, "10.0", "None, None")
        run_filter_dataset(ssfn, ofn, "0.01", "None, None")
        run_filter_dataset(ssfn, ofn, "100000.01", "None, None")
        run_filter_dataset(ssfn, ofn, ".00", "None, None")
        run_filter_dataset(ssfn, ofn, 1.0, "None, None")
        run_filter_dataset(ssfn, ofn, 10.0, "None, None")
        run_filter_dataset(ssfn, ofn, 0.01, "None, None")
        run_filter_dataset(ssfn, ofn, 100000.01, "None, None")
        run_filter_dataset(ssfn, ofn, 100000, "None, None")
        run_filter_dataset(ssfn, ofn, .00, "None, None")
        run_filter_dataset(ssfn, ofn, 1., "None, None")
        run_filter_dataset(ssfn, ofn, "1.", "None, None")
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, "None.1", "None, None")
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, "None1", "None, None")
        with self.assertRaises(ValueError):
            run_filter_dataset(ssfn, ofn, "1.1None", "None, None")

        self.assertEqual(None, sanitize_read_length(""))
        self.assertEqual(None, sanitize_read_length(0))
        self.assertEqual(None, sanitize_read_length(0))
        with self.assertRaises(ValueError):
            sanitize_read_length("None")
        self.assertEqual(None, sanitize_read_length(None))
        self.assertEqual(None, sanitize_read_length(0))
        self.assertEqual(None, sanitize_read_length(False))
        with self.assertRaises(ValueError):
            sanitize_read_length(True)
        self.assertEqual(None, sanitize_read_length(False))
        with self.assertRaises(ValueError):
            sanitize_read_length(True)
        self.assertEqual(1, sanitize_read_length(1))
        self.assertEqual(0, sanitize_read_length("0"))
        self.assertEqual(0, sanitize_read_length(" 0 "))
        self.assertEqual(-1, sanitize_read_length("-1"))
        self.assertEqual(-1, sanitize_read_length(" -1 "))
        self.assertEqual(-1, sanitize_read_length(-1))
        self.assertEqual(1, sanitize_read_length("1.0"))
        self.assertEqual(1, sanitize_read_length(" 1.0 "))
        self.assertEqual(1, sanitize_read_length(1.0))
        with self.assertRaises(ValueError):
            sanitize_read_length(" 1. 0 ")
        with self.assertRaises(ValueError):
            sanitize_read_length("1. 0")
        self.assertEqual(10, sanitize_read_length("10.0"))
        self.assertEqual(10, sanitize_read_length(10.0))
        self.assertEqual(0, sanitize_read_length("0.01"))
        self.assertEqual(0, sanitize_read_length(0.01))
        self.assertEqual(100000, sanitize_read_length("100000.01"))
        self.assertEqual(100000, sanitize_read_length(100000.01))
        self.assertEqual(100000, sanitize_read_length(100000))
        # These two are somewhat annoying:
        self.assertEqual(0, sanitize_read_length(".00"))
        self.assertEqual(None, sanitize_read_length(.00))
        self.assertEqual(1, sanitize_read_length(1.))
        self.assertEqual(1, sanitize_read_length("1."))
        with self.assertRaises(ValueError):
            sanitize_read_length("None.1")
        with self.assertRaises(ValueError):
            sanitize_read_length("None1")
        with self.assertRaises(ValueError):
            sanitize_read_length("1.1None")

    def test_filter_application(self):
        ssfn = data.getXml(8)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name

        # some smoke tests:
        run_filter_dataset(ssfn, ofn, "0", "None")
        ds = openDataSet(ofn)
        self.assertEqual(len(ds), 92)
        run_filter_dataset(ssfn, ofn, "10000", "None")
        ds = openDataSet(ofn)
        self.assertEqual(len(ds), 0)

        run_filter_dataset(ssfn, ofn, "100", "rq > .7")
        ds = openDataSet(ofn)
        self.assertEqual(str(ds.filters),
                         "( rq > .7 AND length >= 100 )")

        # semicolon:
        run_filter_dataset(ssfn, ofn, "100", "rq > .7; length < 5000")
        ds = openDataSet(ofn)
        self.assertEqual(
            str(ds.filters),
            '( length < 5000 AND rq > .7 AND length >= 100 )')

        # comma
        run_filter_dataset(ssfn, ofn, "100", "rq > .7, length < 5000")
        ds = openDataSet(ofn)
        self.assertEqual(
            str(ds.filters),
            '( length < 5000 AND rq > .7 AND length >= 100 )')

        run_filter_dataset(ssfn, ofn, 0,
                           "rq>=.7, length >= 1000, length <= 5000")
        ds = openDataSet(ofn)
        self.assertEqual(
            str(ds.filters),
            '( length >= 1000 AND length <= 5000 AND rq >= .7 )')

        run_filter_dataset(ssfn, ofn, 0,
                           "rq>=.7, length gte 1000, length lte 5000")
        ds = openDataSet(ofn)
        self.assertEqual(
            str(ds.filters),
            '( length gte 1000 AND length lte 5000 AND rq >= .7 )')


    def test_datset_name(self):
        ssfn = data.getXml(8)
        ofn = tempfile.NamedTemporaryFile(suffix=".xml").name
        run_filter_dataset(ssfn, ofn, "0", "None")
        ds = openDataSet(ofn)
        self.assertTrue(ds.name.endswith("(filtered)"))
        self.assertTrue("filtered" in ds.tags)


if __name__ == "__main__":
    unittest.main()
