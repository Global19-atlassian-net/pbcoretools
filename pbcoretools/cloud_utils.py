"""
Utilities for running PacBio workflows on cloud services
"""

import logging
import math
import os.path as op

from pbcore.io import PacBioBamIndex

log = logging.getLogger(__name__)

# https://sourceforge.net/p/samtools/mailman/message/28413844/
BGZF_TERM = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'


def _write_bam_chunk(bam_in, bam_out, header_bytes, offset, record_n_bytes):
    bam_out.write(header_bytes)
    bam_in.seek(offset)
    bam_out.write(bam_in.read(record_n_bytes))
    bam_out.write(BGZF_TERM)


def extract_bam_chunk(bam_in,
                      output_file_name,
                      header_n_bytes,
                      offset,
                      record_n_bytes):
    with open(output_file_name, "wb") as bam_out:
        bam_in.seek(0)
        header_bytes = bam_in.read(header_n_bytes)
        _write_bam_chunk(bam_in, bam_out, header_bytes, offset, record_n_bytes)


def get_zmw_bgzf_borders(pbi):
    offsets = []
    byteOffsets = pbi.virtualFileOffset >> 16
    for i, (zmw, offset) in enumerate(zip(pbi.holeNumber, byteOffsets)):
        if i == 0 or (zmw != pbi.holeNumber[i-1] and offset != byteOffsets[i-1]):
            offsets.append((i, zmw, offset))
    return offsets


def get_bam_offsets(file_name, nchunks):
    assert nchunks >= 1
    pbi_file = file_name + ".pbi"
    if not op.exists(pbi_file):
        raise IOError("BAM input must be accompanied by PacBio index (.pbi)")
    pbi = PacBioBamIndex(pbi_file)
    header_n_bytes = pbi.virtualFileOffset[0] >> 16
    offsets = get_zmw_bgzf_borders(pbi)
    nchunks = min(nchunks, len(offsets))
    per_chunk = math.ceil(len(offsets) / nchunks)
    chunk_offsets = []
    i = 0
    while i < nchunks:
        k = i * per_chunk
        if k >= len(offsets):
            break
        (i_zmw, zmw, offset) = offsets[k]
        chunk_offsets.append(offset)
        i += 1
    return chunk_offsets


def split_bam(file_name, nchunks, prefix="reads"):
    """
    Given a BAM file name and target number of chunks, write out up to nchunks
    BAM files, each split at a ZMW boundary.
    """
    bam_size = op.getsize(file_name)
    offsets = get_bam_offsets(file_name, nchunks)
    header_n_bytes = offsets[0]
    log.info("Output will be split into {n} chunks".format(n=len(offsets)))
    with open(file_name, "rb") as bam_in:
        for i, offset in enumerate(offsets):
            if i < len(offsets) - 1:
                record_n_bytes = offsets[i+1] - offset
            else:
                record_n_bytes = bam_size - offset
            bam_out = "{p}.chunk{i}.bam".format(p=prefix, i=i)
            log.info("Writing chunk {i} to {f}".format(i=i, f=bam_out))
            extract_bam_chunk(bam_in, bam_out, header_n_bytes,
                              offset, record_n_bytes)
    return len(offsets)
