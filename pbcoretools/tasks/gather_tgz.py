import logging
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import run_main_gather_tgz

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_tgz"
    CHUNK_KEY = "$chunk.tgz_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_tgz --resolved-tool-contract "
    OPT_CHUNK_KEY = 'pbcoretools.task_options.gather_tgz_chunk_key'


def _get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev TGZ Gather",
                            "General Chunk TGZ Gather",
                            Constants.DRIVER,
                            is_distributed=True)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with TGZ chunk key")

    p.add_output_file_type(FileTypes.TGZ, "tgz_out",
                           "TGZ",
                           "Gathered TGZ", "gathered")

    # Only need to add to argparse layer for the commandline
    p.arg_parser.add_str(Constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         "$chunk.tgz_id",
                         "Chunk key",
                         "Chunk key to use (format $chunk.{chunk-key}")

    return p


def _args_runner(args):
    return run_main_gather_tgz(args.cjson_in, args.tgz_out, args.chunk_key)


def _rtc_runner(rtc):
    return run_main_gather_tgz(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.chunk_key)


def _main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           _get_parser(),
                           _args_runner,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(_main())
