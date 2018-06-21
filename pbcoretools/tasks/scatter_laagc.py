"""
Scatter SubreadSet by barcode with companion ReferenceSet for LAAgc
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
    TOOL_ID = "pbcoretools.tasks.scatter_laagc"
    VERSION = "0.1"
    CHUNK_KEYS = ('$chunk.subreadset_id', "$chunk.reference_id")
    DRIVER_BASE = "python -m pbcoretools.tasks.scatter_laagc --resolved-tool-contract "
    OPT_CHUNK_KEY = "pbcoretools.task_options.dev_scatter_chunk_key"
    OPT_MAX_NCHUNKS = 'pbcoretools.task_options.scatter_laagc_max_chunks'
    DEFAULT_NCHUNKS = 12


def get_contract_parser():
    p = get_scatter_pbparser(Constants.TOOL_ID,
                             Constants.VERSION,
                             "Scatter SubreadSet",
                             "Pacbio DataSet SubreadSet",
                             Constants.DRIVER_BASE,
                             Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_REF,
                          "ds_reference",
                          "ReferenceSet",
                          "Pac Bio ReferenceSet XML")

    p.add_input_file_type(FileTypes.DS_SUBREADS,
                          "subreads",
                          "SubreadSet",
                          "Pacbio DataSet SubreadSet")

    p.add_output_file_type(FileTypes.CHUNK,
                           "cjson_out",
                           "Chunk JSON",
                           "Chunked JSON",
                           "subreads_reference.chunked")

    # max nchunks for this specific task
    p.add_int(Constants.OPT_MAX_NCHUNKS,
              "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p


def run_main(ds_xml, reference_set_xml, output_json, max_nchunks, output_dir):
    extra_chunk_keys = {"$chunk.reference_id": reference_set_xml}
    CU.write_subreadset_barcode_chunks_to_file(
        chunk_file=output_json,
        dataset_path=ds_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk.subreadset",
        chunk_ext=FileTypes.DS_SUBREADS.ext,
        extra_chunk_keys=extra_chunk_keys)
    return 0


def args_runner(args):
    # FIXME. The chunk needs to be passed directly the func
    chunk_key = args.chunk_key
    output_dir = os.path.dirname(args.cjson_out)
    return run_main(args.subreads, args.ds_reference, args.cjson_out, args.max_nchunks, output_dir)


def rtc_runner(rtc):
    return run_main(rtc.task.input_files[1],
                    rtc.task.input_files[0],
                    rtc.task.output_files[0],
                    rtc.task.max_nchunks,
                    os.path.dirname(rtc.task.output_files[0]))


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
