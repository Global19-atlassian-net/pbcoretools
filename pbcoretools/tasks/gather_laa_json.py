
import logging
import json
import sys

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)

class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_laa_json"
    CHUNK_KEY = "$chunk.json_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_laa_json --resolved-tool-contract "
    OPT_CHUNK_KEY = 'pbcoretools.task_options.gather_json_chunk_key'


def gather_laa_json(input_files, output_file):
    d = {}
    for file_name in input_files:
        with open(file_name) as chunk_in:
            d.update(json.load(chunk_in))
    with open(output_file, "w") as json_out:
        json.dump(d, json_out)
    return 0


def run(chunk_input_json, output_file, chunk_key):
    chunks = load_pipeline_chunks_from_json(chunk_input_json)
    chunked_files = []
    for chunk in chunks:
        if chunk_key in chunk.chunk_keys:
            chunked_files.append(chunk.chunk_d[chunk_key])
        else:
            raise KeyError("Unable to find chunk key '{i}' in {p}".format(i=chunk_key, p=chunk))
    return gather_laa_json(chunked_files, output_file)


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "LAA JSON gather",
                            "Gather task for LAA subreads JSON output",
                            Constants.DRIVER,
                            is_distributed=False)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with JSON chunk key")
    p.add_output_file_type(FileTypes.JSON, "json_out", "JSON",
                           "Gathered JSON",
                           default_name="amplicon_analysis_subreads")
    return p


def args_runner(args):
    return run(args.cjson_in, args.json_out, Constants.CHUNK_KEY)


def rtc_runner(rtc):
    return run(rtc.task.input_files[0], rtc.task.output_files[0],
               rtc.task.chunk_key)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
