"""
Dataset filtering tool that incorporates downsampling; replaces old pbsmrtpipe
task pbcoretools.tasks.filter_dataset.
"""

import logging
import re
import os.path as op
import sys


from pbcoretools.DataSetEntryPoints import parse_filter_list
from pbcore.io import openDataSet

log = logging.getLogger(__name__)
__version__ = "0.1"


def sanitize_read_length(read_length):
    if read_length:
        if not re.search(r'^-?\d*(\.\d*)?$', str(read_length).strip()):
            raise ValueError('read_length filter value "{v}" is not a '
                             'number'.format(v=read_length))
        try:
            return int(read_length)
        except ValueError:
            return int(float(read_length))


def combine_filters(ds, filters):
    if len(ds.filters) > 0:
        for old_filter in ds.filters:
            log.info(
                "Combining user-supplied filters with existing filter '%s", old_filter)
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
                       downsample_factor=0,
                       min_rq=None):
    rlen = sanitize_read_length(read_length)
    filters = {}
    if other_filters and other_filters != "None":
        filters = parse_filter_list([str(other_filters)])
        log.info("{i} other filters will be added".format(i=len(filters)))
    tags = set()
    if rlen or min_rq is not None or len(filters) > 0 or not downsample_factor in [0, 1]:
        dataSet = openDataSet(in_file)
        orig_uuid = dataSet.uuid
        dataSet.updateCounts()  # just in case
        combine_filters(dataSet, filters)
        tags.update({t.strip() for t in dataSet.tags.strip().split(",")})
        if rlen:
            combine_filters(dataSet, {'length': [('>=', rlen)]})
        if min_rq is not None and min_rq > 0:
            combine_filters(dataSet, {'rq': [('>=', min_rq)]})
        if not downsample_factor in [0, 1]:
            combine_filters(dataSet, {'zm': [("==", "0", downsample_factor)]})
            tags.add("downsampled")
        dataSet.updateCounts()
        # XXX note we do *not* set a new UUID in case we want to keep a parent-
        # child relationship to the input dataset.  since the filtered dataset
        # will not be imported back into SMRT Link it is ok to keep the
        # original UUID
        dataSet.uuid = orig_uuid
    else:
        # if we're not actually changing anything, don't load indices
        dataSet = openDataSet(in_file, skipCounts=True)
    tags.add("filtered")
    dataSet.tags = ",".join(list(tags))
    if not "(filtered)" in dataSet.name:
        dataSet.name = dataSet.name + " (filtered)"
    if len(dataSet.metadata.provenance) > 0:
        log.warning("Removing existing provenance record: %s",
                    dataSet.metadata.provenance)
        dataSet.metadata.provenance = None
    dataSet.write(out_file)
    return 0
