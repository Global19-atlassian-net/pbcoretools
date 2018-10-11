import logging
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.models.common import DataStore
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import get_datum_from_chunks_by_chunk_key
from pbcoretools.datastore_utils import datastore_to_datastorefile_objs


log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_datastore"
    CHUNK_KEY = "$chunk.datastore_id"
    VERSION = "0.1.0"
    DRIVER = "python -m {} --resolved-tool-contract ".format(TOOL_ID)
    ALLOWED_TYPES = (
        FileTypes.DS_SUBREADS,
        FileTypes.DS_CCS,
        FileTypes.DS_TRANSCRIPT,
        FileTypes.DS_ALIGN,
        FileTypes.DS_ALIGN_CCS,
        FileTypes.DS_ALIGN_TRANSCRIPT
    )


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev DataStore Gather",
                            "General Chunk DataStore Gather",
                            Constants.DRIVER,
                            is_distributed=True)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "CHUNK Json",
                          "Gathered CHUNK Json with DataStore chunk key")
    p.add_output_file_type(FileTypes.JSON, "out", "DataStore",
                           "Gathered DataStore", default_name="gathered")
    return p


def run_main(chunk_json, datastore_json, chunk_key):
    """run main"""
    chunks = load_pipeline_chunks_from_json(chunk_json)
    # Allow looseness
    if not chunk_key.startswith('$chunk.'):
        chunk_key = '$chunk.' + chunk_key
        log.warn("Prepending chunk key with '$chunk.' to '%s'", str(chunk_key))

    datastore_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    log.debug("Chunked DataStore files are %s.", (', '.join(datastore_files)))
    log.info("Gather chunked DataStore files to %s.", datastore_json)

    if len(datastore_files) == 0:  # Sanity check chunked datastore files are not empty
        raise ValueError("Expected one or more chunked datastore files!")

    # Sanity Check all DataStoreFile objs in datastore files must share the same file type and
    # the file type must be allowed
    objs, dataset_type_id, readcls, ext = datastore_to_datastorefile_objs(
        in_datastore_json=datastore_files[0], allowed_types=Constants.ALLOWED_TYPES)
    for chunked_file in datastore_files[1:]:
        _objs, _dataset_type_id, _readcls, _ext = datastore_to_datastorefile_objs(
            in_datastore_json=chunked_file, allowed_types=Constants.ALLOWED_TYPES)
        if _dataset_type_id != dataset_type_id:
            raise ValueError("Could not gather datastore files of mixed types {} and {}"
                             .format(_dataset_type_id, dataset_type_id))
        objs.extend(_objs)
    ds_out = DataStore(objs)
    ds_out.write_json(datastore_json)
    return 0


def args_runner(args):
    return run_main(args.cjson_in, args.out, Constants.CHUNK_KEY)

def rtc_runner(rtc):
    return run_main(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.chunk_key)

def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
