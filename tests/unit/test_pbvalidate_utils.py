
from cStringIO import StringIO
import xml.dom.minidom
import unittest
import os.path as op
import os

import pbcoretools.pbvalidate.utils
import pbcoretools.pbvalidate.main

DATA_DIR = op.join(op.dirname(op.dirname(__file__)), "data")


class MiscTests(unittest.TestCase):

    def test_junit_output(self):
        file_name_1 = os.path.join(DATA_DIR, "tst_1_subreads.bam")
        file_name_2 = os.path.join(DATA_DIR, "tst_2_subreads.bam")
        results = []
        for file_name in [file_name_1, file_name_2]:
            p = pbcoretools.pbvalidate.main.get_parser()
            args = p.parse_args([file_name, "--quiet"])
            results.append(pbcoretools.pbvalidate.main.run_validator(args))
        out = StringIO()
        pbcoretools.pbvalidate.utils.generate_multiple_file_junit_report(
            results=results,
            xml_out=out,
            skipped_files=[os.path.join(DATA_DIR, "tst_3_subreads.bam")])
        dom = xml.dom.minidom.parseString(out.getvalue())
        tests = dom.getElementsByTagName("testcase")
        self.assertEqual(len(tests), 3)
        failures = dom.getElementsByTagName("failure")
        self.assertEqual(len(failures), 1)
        skipped = dom.getElementsByTagName("skipped")
        self.assertEqual(len(skipped), 1)

    def test_single_file_xunit_output(self):
        file_name = os.path.join(DATA_DIR, "tst_2_subreads.bam")
        p = pbcoretools.pbvalidate.main.get_parser()
        args = p.parse_args([file_name, "--quiet", "--index"])
        result = pbcoretools.pbvalidate.main.run_validator(args)
        dom = result.to_xml()
        failures = dom.getElementsByTagName("failure")
        self.assertEqual(len(failures), 16)


if __name__ == "__main__":
    unittest.main()
