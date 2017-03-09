
"""
Scatter subreads by barcodes, used for input to Long Amplicon Analysis
processing.
"""

import functools
import logging
import os
import sys

from pbcore.io import ConsensusReadSet, ConsensusAlignmentSet, ReferenceSet
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes, PipelineChunk

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.scatter_ccs_aligned_barcodes"


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.scatter_ccs_aligned_barcodes"
    DEFAULT_NCHUNKS = 5
    DRIVER_EXE = "python -m pbcoretools.tasks.scatter_ccs_aligned_barcodes --resolved-tool-contract "
    DATASET_TYPE = FileTypes.DS_ALIGN_CCS
    CHUNK_KEYS = (CU.Constants.CHUNK_KEY_CCS_ALNSET,
                  CU.Constants.CHUNK_KEY_REF,
                  CU.Constants.CHUNK_KEY_CCSSET)
    READ_TYPE = "ConsensusAlignment"
    READ_TYPE_ABBREV = "ccsalignments"

def get_contract_parser_impl(C):
    p = get_scatter_pbparser(C.TOOL_ID, "0.1.3",
        "%sSet barcode scatter" % C.READ_TYPE,
        "Scatter %s DataSet by barcodess" % C.READ_TYPE, C.DRIVER_EXE,
        C.CHUNK_KEYS, is_distributed=True)

    p.add_input_file_type(FileTypes.DS_ALIGN_CCS,
                          "alignments",
                          "ConsensusAlignmentSet",
                          "ConsensusAlignmentSet XML")
    p.add_input_file_type(FileTypes.DS_REF,
                          "reference",
                          "ReferenceSet",
                          "ReferenceSet XML")
    p.add_input_file_type(FileTypes.DS_CCS,
                          "ccs",
                          "ConsensusReadSet",
                          "ConsensusReadSet XML")
    p.add_output_file_type(FileTypes.CHUNK,
                           "chunk_report_json",
                           "Chunk %sSet" % C.READ_TYPE,
                           "PacBio Chunked JSON %sSet" % C.READ_TYPE,
                           "%sset_chunked" % C.READ_TYPE_ABBREV)

    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.scatter_ccs_aln_max_nchunks",
              "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p

get_contract_parser = functools.partial(get_contract_parser_impl, Constants)


def to_chunked_ccsaln_files(alignments_path, reference_path, ccs_path,
                            max_total_nchunks, chunk_key, dir_name,
                            base_name, ext):
    dset = ConsensusAlignmentSet(alignments_path, strict=True)
    dset_chunks = dset.split(barcodes=True, chunks=max_total_nchunks)
    # sanity checking
    reference_set = ReferenceSet(reference_path, strict=True)
    ccs_set = ConsensusReadSet(ccs_path, strict=True)
    d = {}
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        d[Constants.CHUNK_KEYS[1]] = reference_path
        d[Constants.CHUNK_KEYS[2]] = ccs_path
        c = PipelineChunk(chunk_id, **d)
        yield c


def write_barcode_chunks_to_file(chunk_file, alignments_path, reference_path,
                                 ccs_path, max_total_chunks, dir_name,
                                 chunk_base_name, chunk_ext):
    chunks = list(to_chunked_ccsaln_files(alignments_path,
                                          reference_path,
                                          ccs_path,
                                          max_total_chunks,
                                          Constants.CHUNK_KEYS[0],
                                          dir_name, chunk_base_name,
                                          chunk_ext))
    CU.write_chunks_to_json(chunks, chunk_file)
    return 0


def run_main(chunk_output_json, alignments_xml, reference_xml, ccs_xml, max_nchunks, output_dir):
    return write_barcode_chunks_to_file(
        chunk_file=chunk_output_json,
        alignments_path=alignments_xml,
        reference_path=reference_xml,
        ccs_path=ccs_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk_dataset",
        chunk_ext=FileTypes.DS_ALIGN_CCS.ext)


def _args_runner(args):
    return run_main(args.chunk_report_json, args.alignments, args.reference,
                    args.ccs, args.max_nchunks, os.path.dirname(args.chunk_report_json))


def _rtc_runner(rtc):
    output_dir = os.path.dirname(rtc.task.output_files[0])
    max_nchunks = rtc.task.max_nchunks
    return run_main(rtc.task.output_files[0], rtc.task.input_files[0],
                    rtc.task.input_files[1], rtc.task.input_files[2],
                    max_nchunks, output_dir)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_runner,
                           _rtc_runner,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
