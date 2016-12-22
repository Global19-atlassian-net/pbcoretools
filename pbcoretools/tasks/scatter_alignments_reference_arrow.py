"""
Scatter AlignmentSet, ReferenceSet -> Chunk.JSON

This is similar to scatter_alignments_refernce.py but splits by ZMW instead of by contig.

"""
import logging
import os
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes
from pbcoretools.tasks import scatter_alignments_reference


import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.alignment_zmw_scatter"
MODULE_NAME = "pbcoretools.tasks.scatter_alignments_reference_arrow"


class Constants(object):
    DEFAULT_NCHUNKS = 16
    CHUNK_KEYS = ('$chunk.alignmentset_id', "$chunk.reference_id")
    DRIVER_BASE = "python -m {module} --resolved-tool-contract "
    OPT_CHUNK_KEY = "pbcoretools.task_options.dev_scatter_chunk_key"
    OPT_MAX_NCHUNKS = 'pbcoretools.task_options.scatter_alignments_reference_max_nchunks'


def run_main(ds_xml, reference_set_xml, output_json, max_nchunks, output_dir):
    CU.write_alignmentset_chunks_to_file(output_json,
                                         ds_xml,
                                         reference_set_xml,
                                         max_nchunks,
                                         output_dir,
                                         "chunk_alignmentset",
                                         FileTypes.DS_ALIGN.ext,
                                         by_zmw = True)
    return 0


def args_runner(args):
    # FIXME. The chunk needs to be passed directly the func
    chunk_key = args.chunk_key
    #chunk_key = "alignmentset_id"
    output_dir = os.path.dirname(args.cjson_out)
    return run_main(args.alignment_ds, args.ds_reference, args.cjson_out, args.max_nchunks, output_dir)


def rtc_runner(rtc):
    return run_main(rtc.task.input_files[0],
                    rtc.task.input_files[1],
                    rtc.task.output_files[0],
                    rtc.task.max_nchunks,
                    os.path.dirname(rtc.task.output_files[0]))


def main(argv=sys.argv):
    mp = scatter_alignments_reference.get_contract_parser(
        tool_id=TOOL_ID,
        module_name=MODULE_NAME)
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
