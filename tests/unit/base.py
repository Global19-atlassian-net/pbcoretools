
import subprocess
import unittest
import tempfile
import logging
import os.path as op
import os

TESTDATA = "/pbi/dept/secondary/siv/testdata/pbcoretools-unittest/data"

skip_if_no_testdata = unittest.skipUnless(op.isdir(TESTDATA), "Testdata not found")

log = logging.getLogger(__name__)

def _get_temp_file(suffix, dir_):
    t = tempfile.NamedTemporaryFile(suffix=suffix, delete=False, dir=dir_)
    t.close()
    return t.name


def get_temp_file(suffix="", dir_=None):
    return _get_temp_file(suffix, dir_=dir_)

def get_temp_dir(suffix=""):
    """This will make subdir in the root tmp dir"""
    return tempfile.mkdtemp(dir=None, suffix=suffix)


class IntegrationBase(unittest.TestCase):

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp_dir = tempfile.mkdtemp()
        os.chdir(self._tmp_dir)

    def tearDown(self):
        os.chdir(self._cwd)

    def _check_call(self, args):
        log.info("Writing logs to subprocess.std* in %s", self._tmp_dir)
        with open("subprocess.stdout", "w") as stdout:
            with open("subprocess.stderr", "w") as stderr:
                return subprocess.check_call(args, stdout=stdout, stderr=stderr)
