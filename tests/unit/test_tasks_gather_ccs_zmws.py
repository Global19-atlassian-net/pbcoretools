
import unittest
import tempfile
import gzip
import json

from pbcommand.testkit.core import PbIntegrationBase

from pbcoretools.tasks.gather_ccs_zmws import gather_chunks, FIELDS


class TestGatherCcsZmws(PbIntegrationBase, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        d1 = {"zmws": [{FIELDS[x]:x for x in range(5)}]}
        d2 = {"zmws": [{FIELDS[x]:x+5 for x in range(5)}]}
        j1 = tempfile.NamedTemporaryFile(suffix=".json.gz").name
        j2 = tempfile.NamedTemporaryFile(suffix=".json.gz").name
        for d, fn in zip([d1, d2], [j1, j2]):
            with gzip.open(fn, mode="wt") as gz_out:
                gz_out.write(json.dumps(d))
        cls.INFILES = [j1, j2]

    def _validate_output(self, file_name):
        with gzip.open(file_name, mode="rt") as gz_in:
            d = json.loads(gz_in.read())
            assert d == {"zmws": [
                {FIELDS[x]:x for x in range(5)},
                {FIELDS[x]:x+5 for x in range(5)}
            ]}

    def test_gather(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".json.gz").name
        gather_chunks(self.INFILES, ofn)
        self._validate_output(ofn)

    def test_integration(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".json.gz").name
        args = [
            "python3", "-m", "pbcoretools.tasks.gather_ccs_zmws",
             ofn
        ] + self.INFILES
        self._check_call(args)
        self._validate_output(ofn)
