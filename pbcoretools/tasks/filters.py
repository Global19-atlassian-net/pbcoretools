
"""
Tool contract wrappers for miscellaneous quick functions involving filtration
"""

import functools
import tempfile
import logging
import shutil
import gzip
import re
import os
import sys

from pbcoretools.DataSetEntryPoints import parse_filter_list
from pbcore.io import (SubreadSet, HdfSubreadSet, FastaReader, FastaWriter,
                       FastqReader, FastqWriter, openDataSet)
from pbcommand.engine import run_cmd
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes

log = logging.getLogger(__name__)

TOOL_NAMESPACE = 'pbcoretools'
DRIVER_BASE = "python -m pbcoretools.tasks.filters "

registry = registry_builder(TOOL_NAMESPACE, DRIVER_BASE)


rl_opt = QuickOpt(0, "Minimum subread length",
    "Minimum length of subreads")

filters_opt = QuickOpt("", "Filters to add to the DataSet",
    "A comma separated list of other filters to add to the DataSet")


def run_filter_dataset(in_file, out_file, read_length, other_filters):
    dataSet = openDataSet(in_file)
    if read_length:
        dataSet.filters.addRequirement(length=[('>=', read_length)])
        log.info("Readlength filter added")
    if other_filters:
        filters = parse_filter_list(str(other_filters).split(','))
        dataSet.filters.addRequirement(**filters)
        log.info("{i} other filters added".format(i=len(filters)))
    if read_length or other_filters:
        dataSet.updateCounts()
    dataSet.write(out_file)
    return 0

@registry("filterdataset", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.DS_SUBREADS, is_distributed=True, nproc=1,
          options={"read_length":rl_opt,
                   "other_filters":filters_opt})
def run_filterDataSet(rtc):
    return run_filter_dataset(
        rtc.task.input_files[0], rtc.task.output_files[0],
        rtc.task.options["pbcoretools.task_options.read_length"],
        rtc.task.options["pbcoretools.task_options.other_filters"])

if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
