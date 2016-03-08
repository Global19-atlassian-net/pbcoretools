
import logging
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import run_main_gather_alignmentset

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_alignmentset"
    CHUNK_KEY = "$chunk:alignmentset_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_alignments --resolved-tool-contract "
    CONSOLIDATE_ID = "pbcoretools.task_options.consolidate_aligned_bam"
    N_FILES_ID = "pbcoretools.task_options.consolidate_n_files"


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev Alignments Gather",
                            "General Chunk Alignments Gather",
                            Constants.DRIVER,
                            is_distributed=False)

    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with AlignmentSet chunk key")

    p.add_output_file_type(FileTypes.DS_ALIGN,
                           "ds_out",
                           "AlignmentSet",
                           "Gathered AlignmentSet",
                           default_name="gathered")
    # XXX moved this functionality to pbalign
#    p.add_boolean(Constants.CONSOLIDATE_ID, "consolidate",
#        default=False,
#        name="Consolidate .bam",
#        description="Merge chunked/gathered .bam files")
#    p.add_int(Constants.N_FILES_ID, "consolidate_n_files",
#        default=1,
#        name="Number of .bam files",
#        description="Number of .bam files to create in consolidate mode")
    return p


def args_runner(args):
    return run_main_gather_alignmentset(
        chunk_input_json=args.cjson_in,
        output_file=args.ds_out,
        chunk_key=args.chunk_key,
        consolidate=False)#args.consolidate,
        #consolidate_n_files=args.consolidate_n_files)


def rtc_runner(rtc):
    return run_main_gather_alignmentset(
        chunk_input_json=rtc.task.input_files[0],
        output_file=rtc.task.output_files[0],
        chunk_key=rtc.task.chunk_key,
        consolidate=False)#rtc.task.options[Constants.CONSOLIDATE_ID],
        #consolidate_n_files=rtc.task.options[Constants.N_FILES_ID])


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
