
import subprocess
import unittest
import tempfile
import logging
import os

from pbcoretools.tasks2.chunk_reads_for_lima import split_reads

import pbtestdata

log = logging.getLogger(__name__)

class TestChunkReadsForLima(unittest.TestCase):

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp_dir = tempfile.mkdtemp()
        log.info("temp dir is %s", self._tmp_dir)
        os.chdir(self._tmp_dir)

    def tearDown(self):
        os.chdir(self._cwd)

    def test_split_reads(self):
        ds_file = pbtestdata.get_file("subreads-sequel")
        nchunks = split_reads(ds_file, 0, 1, False)
        self.assertEqual(nchunks, 2)
        nchunks = split_reads(ds_file, 1, 1, False)
        self.assertEqual(nchunks, 1)

    def test_split_reads_peek_guess(self):
        ds_file = pbtestdata.get_file("subreads-sequel")
        nchunks = split_reads(ds_file, 0, 1, True)
        self.assertEqual(nchunks, 1)
        nchunks = split_reads(ds_file, 2, 1, True)
        self.assertEqual(nchunks, 1)

    def _check_call(self, args):
        with open("stdout", "w") as stdout:
            with open("stderr", "w") as stderr:
                return subprocess.check_call(args, stdout=stdout, stderr=stderr)

    def test_integration(self):
        args = [
            "python", "-m", "pbcoretools.tasks2.chunk_reads_for_lima",
            pbtestdata.get_file("subreads-sequel"),
            "--targetSize", "1"
        ]
        self._check_call(args)
        files = [f for f in os.listdir(".") if f.endswith("set.xml")]
        self.assertEqual(len(files), 2)

    def test_integration_peek_guess(self):
        args = [
            "python", "-m", "pbcoretools.tasks2.chunk_reads_for_lima",
            pbtestdata.get_file("subreads-sequel"),
            "--targetSize", "1"
        ]
        args.append("--lima-peek-guess")
        self._check_call(args)
        files = [f for f in os.listdir(".") if f.endswith("set.xml")]
        self.assertEqual(len(files), 1)
