
"""
Scatter processed RNA transcripts by pseudo-ZMW range, also passing along the
original input SubreadSet.
"""

import logging
import os
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.scatter_transcripts"
    VERSION = "0.1"
    DEFAULT_NCHUNKS = 2
    DRIVER_EXE = "python -m pbcoretools.tasks.scatter_transcripts --resolved-tool-contract "
    DATASET_TYPE = FileTypes.DS_TRANSCRIPT
    CHUNK_KEYS = (CU.Constants.CHUNK_KEY_TRANSCRIPT,
                  CU.Constants.CHUNK_KEY_SUBSET)
    READ_TYPE = "Transcript"
    READ_TYPE_ABBREV = "transcript"


def _get_parser():
    p = get_scatter_pbparser(Constants.TOOL_ID, "0.1.3",
        "TranscriptSet ZMW scatter",
        "Scatter Transcript DataSet by ZMWs",
        Constants.DRIVER_EXE,
        Constants.CHUNK_KEYS,
        is_distributed=True)
    p.add_input_file_type(FileTypes.DS_TRANSCRIPT,
                          "transcripts",
                          "TranscriptSet",
                          "Pac Bio Transcript DataSet XML format")
    p.add_input_file_type(FileTypes.DS_SUBREADS,
                          "subreads",
                          "SubreadSet",
                          "Pac Bio Subread DataSet XML format")
    p.add_output_file_type(FileTypes.CHUNK,
                           "chunk_report_json",
                           "Chunk TranscriptSet",
                           "PacBio Chunked JSON TranscriptSet",
                           "transcriptset_chunked")
    return p


def run_main(chunk_output_json, dataset_xml, subreads_in, max_nchunks, output_dir):
    return CU.write_transcriptset_zmw_chunks_to_file(
        chunk_file=chunk_output_json,
        dataset_path=dataset_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk_transcriptset",
        chunk_ext=FileTypes.DS_TRANSCRIPT.ext,
        extra_chunk_keys={CU.Constants.CHUNK_KEY_SUBSET: subreads_in})


def _args_runner(args):
    return run_main(args.chunk_report_json, args.transcripts, args.subreads,
                    args.max_nchunks, os.path.dirname(args.chunk_report_json))


def _rtc_runner(rtc):
    return run_main(rtc.task.output_files[0],
                    rtc.task.input_files[0],
                    rtc.task.input_files[1],
                    rtc.task.max_nchunks,
                    os.path.dirname(rtc.task.output_files[0]))


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           _get_parser(),
                           _args_runner,
                           _rtc_runner,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
