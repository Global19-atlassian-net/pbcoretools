#!/usr/bin/env python
"""
Consolidate TranscriptSet into bam files, allowing adding a prefix
string (e.g., 'mysample_HQ_') to every transcript names.
"""

import sys
import logging
import pysam

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import FileTypes, get_pbparser, ResourceTypes
from pysam import AlignmentFile  # pylint: disable=no-member, no-name-in-module

from pbcore.io import TranscriptSet
from pbcoretools.tasks.consolidate_alignments import Constants as BaseConstants, run_consolidate
from pbcoretools.file_utils import get_prefixes

logger = logging.getLogger()


def get_consolidate_parser(tool_id, file_type, driver_exe, version, description):
    """
    Input:
        idx - 0 SubreadSet
        idx - 1 HQ TranscriptSet
        idx - 2 LQ TranscriptSet
    Output:
        idx - 0 HQ TranscriptSet, of which read names have biosample_HQ prefix
        idx - 1 LQ TranscriptSet, of which read names have biosample_LQ prefix
        idx - 2 HQ DataStore of output TranscriptSet
        idx - 3 LQ DataStore of output TranscriptSet
    """
    ds_type = file_type.file_type_id.split(".")[-1]
    p = get_pbparser(tool_id,
                     version,
                     "{t} consolidate".format(t=ds_type),
                     description,
                     driver_exe,
                     is_distributed=True,
                     resource_types=(ResourceTypes.TMP_DIR,))
    p.add_input_file_type(FileTypes.DS_SUBREADS,
                          "subreads",
                          "Input SubreadSet",
                          "SubreadSet with biosample metadata.")
    p.add_input_file_type(file_type,
                          "hq_ds_in",
                          "Input High Quality {t}".format(t=ds_type),
                          "Gathered {t} to consolidate".format(t=ds_type))
    p.add_input_file_type(file_type,
                          "lq_ds_in",
                          "Input Low Quality {t}".format(t=ds_type),
                          "Gathered {t} to consolidate".format(t=ds_type))
    p.add_output_file_type(file_type,
                           "hq_ds_out",
                           "Output High Quality ",
                           description="Output {t} of consolidated bam files".format(t=ds_type),
                           default_name="combined.hq")
    p.add_output_file_type(file_type,
                           "lq_ds_out",
                           "Output Low Quality ",
                           description="Output {t} of consolidated bam files".format(t=ds_type),
                           default_name="combined.lq")
    p.add_output_file_type(FileTypes.JSON,
                           "hq_datastore",
                           "JSON Datastore",
                           description="Datastore containing High Quality {t}".format(t=ds_type),
                           default_name="resources.hq")
    p.add_output_file_type(FileTypes.JSON,
                           "lq_datastore",
                           "JSON Datastore",
                           description="Datastore containing Low Quality {t}".format(t=ds_type),
                           default_name="resources.lq")
    return p


class Constants(BaseConstants):
    TOOL_ID = "pbcoretools.tasks.consolidate_transcripts"
    INPUT_FILE_TYPE = FileTypes.DS_TRANSCRIPT
    TOOL_DESC = __doc__
    DRIVER = "python -m {} --resolved-tool-contract ".format(TOOL_ID)


def consolidate_transcripts(ds_in, prefix):
    """Return a function which
    - must take (new_resource_file, numFiles, useTmp) as input,
    - should consolidate ds_in (input transcripset)
    - should add biosample prefix to transcript read names
    """
    def _consolidate_transcripts_f(new_resource_file, numFiles, useTmp,
                                   perfix=prefix, ds_in=ds_in):
        external_files = ds_in.toExternalFiles()
        assert len(external_files) >= 1, "{!r} must contain one or more bam files".format(ds_in)
        header = AlignmentFile(external_files[0], 'rb', check_sq=False).header
        with AlignmentFile(new_resource_file, 'wb', header=header) as writer:
            for external_file in external_files:
                with AlignmentFile(external_file, 'rb', check_sq=False) as reader:
                    for record in reader:
                        record.query_name = prefix + record.query_name
                        writer.write(record)
        ds_in = TranscriptSet(new_resource_file)  # override ds_in
    return _consolidate_transcripts_f


def __runner(ds_items):
    for ds_in, ds_out, datastore, prefix in ds_items:
        def func(ds_in):
            return consolidate_transcripts(ds_in, prefix=prefix)
        run_consolidate(dataset_file=ds_in,
                        output_file=ds_out,
                        datastore_file=datastore,
                        consolidate=True,
                        n_files=1,
                        task_id=Constants.TOOL_ID,
                        consolidate_f=func)
        # At this piont datastore contains paths to bam/bai/pbi files, now override
        # datastore with TranscriptSet
        from pbcoretools.datastore_utils import dataset_to_datastore
        dataset_to_datastore(ds_out, datastore, source_id=Constants.TOOL_ID)
    return 0


def args_runner(args):
    hq_prefix, lq_prefix = get_prefixes(args.subreads)
    ds_items = [
        (args.hq_ds_in, args.hq_ds_out, args.hq_datastore, hq_prefix),
        (args.lq_ds_in, args.lq_ds_out, args.lq_datastore, lq_prefix)
    ]
    return __runner(ds_items)


def rtc_runner(rtc):
    hq_prefix, lq_prefix = get_prefixes(rtc.task.input_files[0])
    ds_items = [
        (rtc.task.input_files[1], rtc.task.output_files[0], rtc.task.output_files[2], hq_prefix),
        (rtc.task.input_files[2], rtc.task.output_files[1], rtc.task.output_files[3], lq_prefix)
    ]
    return __runner(ds_items)


def main(argv=sys.argv):
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger()
    parser = get_consolidate_parser(Constants.TOOL_ID, Constants.INPUT_FILE_TYPE,
                                    Constants.DRIVER, Constants.VERSION, Constants.TOOL_DESC)
    return pbparser_runner(argv[1:],
                           parser,
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
