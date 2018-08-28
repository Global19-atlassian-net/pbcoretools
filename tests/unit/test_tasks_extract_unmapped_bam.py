
from pbcommand.testkit.core import PbTestApp
from pbcommand.models.common import DataStore
from pbcore.io import AlignmentSet, SubreadSet

import pbtestdata


class TestExtractUnalignedTCDefaults(PbTestApp):
    DRIVER_BASE = "python -m pbcoretools.tasks.extract_unmapped_bam"
    INPUT_FILES = [
        pbtestdata.get_file("aligned-xml"),
        pbtestdata.get_file("subreads-xml")
    ]

    def run_after(self, rtc, output_dir):
        ds = DataStore.load_from_json(rtc.task.output_files[0])
        self.assertEqual(len(ds.files), 0)


class TestExtractUnalignedTC(TestExtractUnalignedTCDefaults):
    TASK_OPTIONS = {
        "pbcoretools.task_options.output_unaligned_bam": True
    }

    def run_after(self, rtc, output_dir):
        ds = DataStore.load_from_json(rtc.task.output_files[0])
        self.assertEqual(len(ds.files), 1)
        with AlignmentSet(self.INPUT_FILES[0]) as mapped:
            mapped_zmws = set(mapped.index.holeNumber)
            with SubreadSet(list(ds.files.values())[0].path) as unmapped:
                unmapped_zmws = set(unmapped.index.holeNumber)
                self.assertEqual(len(mapped_zmws & unmapped_zmws), 0)
