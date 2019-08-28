
import subprocess
import unittest
import tempfile
import logging
import os

from pbcoretools.tasks2.check_reads_for_lima import is_ccs_demultiplexed

import pbtestdata

log = logging.getLogger(__name__)

class TestCheckReadsForLima(unittest.TestCase):

    def test_input_is_not_demultiplexed(self):
        ds_file = pbtestdata.get_file("ccs-sequel")
        is_demultiplexed = is_ccs_demultiplexed(ds_file)
        self.assertEqual(is_demultiplexed, False)

    def test_input_is_demultiplexed(self):
        ds_file = pbtestdata.get_file("ccs-barcoded")
        is_demultiplexed = is_ccs_demultiplexed(ds_file)
        self.assertEqual(is_demultiplexed, True)
