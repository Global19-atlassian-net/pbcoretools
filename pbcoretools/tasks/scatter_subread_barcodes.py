
"""
Scatter subreads by barcodes, used for input to Long Amplicon Analysis
processing.
"""

import functools
import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

from pbcoretools.tasks.scatter_subread_zmws import args_runner_impl, rtc_runner_impl, add_base_subread_scatter_options, ScatterSubreadsConstants
import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.subreadset_barcode_scatter"


class Constants(ScatterSubreadsConstants):
    TOOL_ID = "pbcoretools.tasks.subreadset_barcode_scatter"
    DRIVER_EXE = "python -m pbcoretools.tasks.scatter_subread_barcodes --resolved-tool-contract "


def get_contract_parser_impl(C):
    p = get_scatter_pbparser(C.TOOL_ID, "0.1.3",
        "%sSet barcode scatter" % C.READ_TYPE,
        "Scatter %s DataSet by barcodes" % C.READ_TYPE, C.DRIVER_EXE,
        C.CHUNK_KEYS, is_distributed=True)
    return add_base_subread_scatter_options(C, p)


get_contract_parser = functools.partial(get_contract_parser_impl, Constants)

def run_main(chunk_output_json, dataset_xml, max_nchunks, output_dir):
    return CU.write_subreadset_barcode_chunks_to_file(
        chunk_file=chunk_output_json,
        dataset_path=dataset_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk_dataset",
        chunk_ext=FileTypes.DS_SUBREADS.ext)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           functools.partial(args_runner_impl, run_main),
                           functools.partial(rtc_runner_impl, run_main),
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
