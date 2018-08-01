"""
Tool contract wrappers for miscellaneous quick functions involving filtration
"""

import logging
import sys
import re

from pbcoretools.DataSetEntryPoints import parse_filter_list
from pbcore.io import openDataSet, TranscriptSet
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, OutputFileType

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_NAMESPACE = 'pbcoretools'
    DRIVER_BASE = "python -m pbcoretools.tasks.filters "

    # Iso-Seq
    TRANSCRIPT_QV_CUTOFF = 0.99


registry = registry_builder(Constants.TOOL_NAMESPACE, Constants.DRIVER_BASE)


rl_opt = QuickOpt(0, "Minimum subread length",
                  "Minimum length of subreads")

filters_opt = QuickOpt(
    "",
    "Filters to add to the DataSet",
    ("A semicolon or comma separated list of other filters "
     "to add to the DataSet"))

downsample_opt = QuickOpt(0, "Downsampling Factor", "If other than 0 or 1, this will add a filter to the input dataset to sample a random selection of ZMWs instead of running over the full dataset.  For example, a downsample factor of 10 means that 1/10th (10%) of ZMWs will be used.")

subreads_file_type = OutputFileType(FileTypes.DS_SUBREADS.file_type_id,
                                    "SubreadSet", "Filtered SubreadSet XML",
                                    "Filtered SubreadSet XML", "filtered")

def sanitize_read_length(read_length):
    if read_length:
        if not re.search('^-?\d*(\.\d*)?$', str(read_length).strip()):
            raise ValueError('read_length filter value "{v}" is not a '
                             'number'.format(v=read_length))
        try:
            return int(read_length)
        except ValueError:
            return int(float(read_length))


def combine_filters(ds, filters):
    if len(ds.filters) > 0:
        for old_filter in ds.filters:
            log.info("Combining user-supplied filters with existing filter '%s", old_filter)
            for name, options in filters.items():
                for option_params in options:
                    oper, value = option_params[0:2]
                    modulo = None
                    if len(option_params) > 2:
                        modulo = option_params[2]
                    have_req = False
                    for old_prop in old_filter.plist:
                        if old_prop.name == name and old_prop.operator == oper:
                            old_prop.value = value
                            have_req = True
                    if not have_req:
                        old_filter.addRequirement(name, oper, value,
                                                  modulo=modulo)
        ds.filters._runCallbacks()
    else:
        ds.filters.addFilter(**filters)


def run_filter_dataset(in_file, out_file, read_length, other_filters,
                       downsample_factor=0):
    dataSet = openDataSet(in_file)
    dataSet.updateCounts() # just in case
    rlen = sanitize_read_length(read_length)
    filters = {}
    if other_filters and other_filters != "None":
        if ' AND ' in str(other_filters):
            filters = parse_filter_list(str(other_filters).split(' AND '))
        else:
            filters = parse_filter_list(str(other_filters).split(','))
        log.info("{i} other filters will be added".format(i=len(filters)))
    combine_filters(dataSet, filters)
    tags = {t.strip() for t in dataSet.tags.strip().split(",")}
    if rlen:
        combine_filters(dataSet, {'length': [('>=', rlen)]})
    if not downsample_factor in [0, 1]:
        combine_filters(dataSet, {'zm': [("==", "0", downsample_factor)]})
        tags.add("downsampled")
    dataSet.updateCounts()
    tags.add("filtered")
    dataSet.tags = ",".join(list(tags))
    if not "(filtered)" in dataSet.name:
        dataSet.name = dataSet.name + " (filtered)"
    dataSet.newUuid()
    dataSet.write(out_file)
    return 0


@registry("filterdataset", "0.3.0",
          FileTypes.DS_SUBREADS,
          subreads_file_type, is_distributed=True, nproc=1,
          options={"read_length":rl_opt,
                   "downsample_factor":downsample_opt,
                   "other_filters":filters_opt})
def run_filterDataSet(rtc):
    return run_filter_dataset(
        rtc.task.input_files[0], rtc.task.output_files[0],
        rtc.task.options["pbcoretools.task_options.read_length"],
        rtc.task.options["pbcoretools.task_options.other_filters"],
        rtc.task.options["pbcoretools.task_options.downsample_factor"])


def _split_transcripts(transcripts, hq_file, lq_file, cutoff):
    with TranscriptSet(transcripts, strict=True) as ds_in:
        ds_hq = ds_in.copy()
        ds_lq = ds_in.copy()
        combine_filters(ds_hq, {'rq': [(">=", cutoff)]})
        combine_filters(ds_lq, {'rq': [("<", cutoff)]})
        ds_hq.updateCounts()
        ds_lq.updateCounts()
        ds_hq.write(hq_file)
        ds_lq.write(lq_file)
    return 0


hq_file_type = OutputFileType(FileTypes.DS_TRANSCRIPT.file_type_id,
                              "hq_transcripts",
                              "HQ TranscriptSet",
                              "Hiqh-Quality TranscriptSet XML",
                              "hq_transcripts")
lq_file_type = OutputFileType(FileTypes.DS_TRANSCRIPT.file_type_id,
                              "lq_transcripts",
                              "LQ TranscriptSet",
                              "Low-Quality TranscriptSet XML",
                              "lq_transcripts")
hq_qv_cutoff = QuickOpt(Constants.TRANSCRIPT_QV_CUTOFF,
                        "QV cutoff for HQ transcripts",
                        "Minimum read quality required for a transcript to be considered 'high-quality'")

@registry("split_transcripts", "0.1.1",
          FileTypes.DS_TRANSCRIPT,
          (hq_file_type, lq_file_type),
          is_distributed=True,
          nproc=1,
          options={"hq_qv_cutoff": hq_qv_cutoff})
def _run_split_transcripts(rtc):
    cutoff = rtc.task.options["pbcoretools.task_options.hq_qv_cutoff"]
    return _split_transcripts(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        rtc.task.output_files[1],
        cutoff)


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
