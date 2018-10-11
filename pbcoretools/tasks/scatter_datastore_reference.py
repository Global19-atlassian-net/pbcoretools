
"""
Scatter datastore json to multiple datastore json.
Input datastore json contains an input dataset which must be one of (SubreadSet, ConsensusReadSet,
TranscriptSet). `scatter_datastore` scatters the input dataset by ZMW range,
and writes each scattered dataset into an output datastore json file.

Used for input to pbmm2 align (Subreads, CCS, TranscriptSet) to reference.
"""

import functools
import logging
import os
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.scatter_datastore_reference"
    DRIVER_EXE = "python -m {} --resolved-tool-contract ".format(TOOL_ID)
    DEFAULT_NCHUNKS = 5
    CHUNK_KEYS = (CU.Constants.CHUNK_KEY_DATASTORE_JSON,)

    # Only supports datastore json of following types:
    ALLOWED_TYPES = (
        FileTypes.DS_SUBREADS,
        FileTypes.DS_CCS,
        FileTypes.DS_TRANSCRIPT
    )
    READ_TYPE = ','.join([str(t.file_type_id) for t in ALLOWED_TYPES])


def get_contract_parser_impl(C):
    p = get_scatter_pbparser(C.TOOL_ID, "0.0.1", "DataStore scatter",
                             "Scatter DataStore of {} reads by ZMWs".format(C.READ_TYPE),
                             C.DRIVER_EXE, C.CHUNK_KEYS, is_distributed=True)
    p.add_input_file_type(FileTypes.JSON, "in_datastore_json",
                          "DataStore", "PacBio DataStore Json")
    p.add_input_file_type(FileTypes.DS_REF, "in_reference_xml",
                          "ReferenceSet", "Pac Bio ReferenceSet XML")
    p.add_output_file_type(FileTypes.CHUNK, "out_chunk_json",
                           "Chunk datastore of %s" % C.READ_TYPE,
                           "PacBio Chunked JSON", "datastore_chunked")
    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.scatter_datastore_max_nchunks",
              "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p

get_contract_parser = functools.partial(get_contract_parser_impl, Constants)


def run_main(in_datastore_json, in_ref_xml, out_chunk_json, max_nchunks, out_dir):
    """
    Convert an input datastore json file to a dataset.
    """
    CU.write_datastore_chunks_to_file(out_chunk_json,
                                      in_datastore_json,
                                      in_ref_xml,
                                      max_nchunks,
                                      out_dir,
                                      "chunk_datastore",
                                      'datastore.json')
    return 0


def args_runner_impl(args):
    return run_main(args.in_datastore_json,
                    args.in_reference_xml,
                    args.out_chunk_json,
                    args.max_nchunks,
                    out_dir=os.path.dirname(args.chunk_datastore_json))


def rtc_runner_impl(rtc):
    return run_main(rtc.task.input_files[0], # datastore json
                    rtc.task.input_files[1], # reference xml
                    rtc.task.output_files[0], # chunk json
                    rtc.task.max_nchunks,
                    out_dir=os.path.dirname(rtc.task.output_files[0]))


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner_impl,
                           rtc_runner_impl,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
