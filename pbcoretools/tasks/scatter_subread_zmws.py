
"""
Scatter subreads by ZMW range, used for input to ConsensusRead processing.
"""

import functools
import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.subreadset_zmw_scatter"


class ScatterSubreadsConstants(object):
    DEFAULT_NCHUNKS = 5
    DATASET_TYPE = FileTypes.DS_SUBREADS
    CHUNK_KEYS = ("$chunk.subreadset_id", )
    READ_TYPE = "Subread"
    READ_TYPE_ABBREV = "subread"


class Constants(ScatterSubreadsConstants):
    TOOL_ID = "pbcoretools.tasks.subreadset_zmw_scatter"
    DRIVER_EXE = "python -m pbcoretools.tasks.scatter_subread_zmws --resolved-tool-contract "


def get_contract_parser_impl(C):
    p = get_scatter_pbparser(C.TOOL_ID, "0.1.3",
        "%sSet ZMW scatter" % C.READ_TYPE,
        "Scatter %s DataSet by ZMWs" % C.READ_TYPE, C.DRIVER_EXE,
        C.CHUNK_KEYS, is_distributed=True)
    return add_base_subread_scatter_options(C, p)


def add_base_subread_scatter_options(C, p):
    p.add_input_file_type(C.DATASET_TYPE,
                          "dataset",
                          "%sSet" % C.READ_TYPE,
                          "Pac Bio Fasta format")

    p.add_output_file_type(FileTypes.CHUNK,
                           "chunk_report_json",
                           "Chunk %sSet" % C.READ_TYPE,
                           "PacBio Chunked JSON %sSet" % C.READ_TYPE,
                           "%sset_chunked" % C.READ_TYPE_ABBREV)

    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.scatter_subread_max_nchunks",
              "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")

    p.add_str("pbcoretools.task_options.scatter_subreadset_chunk_key",
              "chunk_key", "$chunk:subreadset_id", "Chunk key",
              "Chunk key to use (format $chunk:{chunk-key}")
    return p

get_contract_parser = functools.partial(get_contract_parser_impl, Constants)

def run_main(chunk_output_json, dataset_xml, max_nchunks, output_dir):
    return CU.write_subreadset_zmw_chunks_to_file(
        chunk_file=chunk_output_json,
        dataset_path=dataset_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk_dataset",
        chunk_ext=FileTypes.DS_SUBREADS.ext)


def args_runner_impl(args):
    return run_main(args.chunk_report_json,
                    args.subreadset,
                    args.max_nchunks,
                    output_dir=os.path.dirname(args.chunk_report_json))


def rtc_runner_impl(run_func, rtc):
    return run_func(rtc.task.output_files[0],
                    rtc.task.input_files[0],
                    rtc.task.max_nchunks,
                    output_dir=os.path.dirname(rtc.task.output_files[0]))


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
