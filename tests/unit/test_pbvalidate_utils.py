import xml.dom.minidom
import tempfile
import os.path as op
import os

import pbcoretools.pbvalidate.utils
import pbcoretools.pbvalidate.main

DATA_DIR = op.join(op.dirname(op.dirname(__file__)), "data")


class MiscTests:

    def test_junit_output(self):
        file_name_1 = os.path.join(DATA_DIR, "tst_1_subreads.bam")
        file_name_2 = os.path.join(DATA_DIR, "tst_2_subreads.bam")
        results = []
        for file_name in [file_name_1, file_name_2]:
            p = pbcoretools.pbvalidate.main.get_parser()
            args = p.parse_args([file_name, "--quiet"])
            results.append(pbcoretools.pbvalidate.main.run_validator(args))
        tmp_out = tempfile.NamedTemporaryFile(suffix=".xml").name
        with open(tmp_out, "wt") as out:
            pbcoretools.pbvalidate.utils.generate_multiple_file_junit_report(
                results=results,
                xml_out=out,
                skipped_files=[os.path.join(DATA_DIR, "tst_3_subreads.bam")])
        with open(tmp_out, "rt") as out:
            dom = xml.dom.minidom.parseString(out.read())
            tests = dom.getElementsByTagName("testcase")
            assert len(tests) == 3
            failures = dom.getElementsByTagName("failure")
            assert len(failures) == 1
            skipped = dom.getElementsByTagName("skipped")
            assert len(skipped) == 1

    def test_single_file_xunit_output(self):
        file_name = os.path.join(DATA_DIR, "tst_2_subreads.bam")
        p = pbcoretools.pbvalidate.main.get_parser()
        args = p.parse_args([file_name, "--quiet", "--index"])
        result = pbcoretools.pbvalidate.main.run_validator(args)
        dom = result.to_xml()
        failures = dom.getElementsByTagName("failure")
        assert len(failures) == 14
