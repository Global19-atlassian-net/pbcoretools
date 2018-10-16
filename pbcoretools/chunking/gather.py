
from collections import defaultdict, namedtuple, OrderedDict
from functools import partial as P
from zipfile import ZipFile
import itertools
import argparse
import tarfile
import logging
import shutil
import json
import math
import os.path as op
import sys

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.validators import fofn_to_files, validate_file
from pbcommand.pb_io.report import load_report_from_json
from pbcommand.common_options import add_debug_option
from pbcommand.cli.utils import main_runner_default, subparser_builder
from pbcommand.cli import get_default_argparser
from pbcommand.models.report import Report

from pbcore.io import (SubreadSet, ContigSet, AlignmentSet, ConsensusReadSet,
                       ConsensusAlignmentSet, TranscriptSet, TranscriptAlignmentSet)
from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbcore.io.GffIO import merge_gffs_sorted
from pbcore.io.VcfIO import merge_vcfs_sorted

log = logging.getLogger(__name__)

__version__ = '1.1'


def _validate_chunk_json_file(path):
    chunks = load_pipeline_chunks_from_json(path)
    return path


def __gather_fastx(fastx_reader, fastx_writer, fastx_files, output_file):

    n = 0
    with fastx_writer(output_file) as writer:
        for fastx_file in fastx_files:
            with fastx_reader(fastx_file) as reader:
                for record in reader:
                    n += 1
                    writer.writeRecord(record)

    log.info("Completed gathering {n} files (with {x} records) to {f}".format(
        n=len(fastx_files), f=output_file, x=n))


gather_fasta = P(__gather_fastx, FastaReader, FastaWriter)
gather_fastq = P(__gather_fastx, FastqReader, FastqWriter)
gather_gff = merge_gffs_sorted
gather_vcf = merge_vcfs_sorted


def _read_header(csv_file):
    with open(csv_file, 'r') as f:
        header = f.readline()

    if ',' in header:
        return header.rstrip().split(',')
    else:
        return None


def __has_header(handle):
    next(handle)


def __has_header_and_one_record(handle):
    # has header
    __has_header(handle)
    # has at least one record
    next(handle)


def __csv_inspector(func, csv_file):

    is_empty = False
    with open(csv_file, 'r') as f:
        try:
            func(f)
        except StopIteration:
            is_empty = True

    return is_empty


_csv_is_empty = P(__csv_inspector, __has_header_and_one_record)
_csv_has_header = P(_csv_is_empty, __has_header)


def get_datum_from_chunks_by_chunk_key(chunks, chunk_key):
    log.info(
        "extracting datum from chunks using chunk-key '{c}'".format(c=chunk_key))
    datum = []
    for chunk in chunks:
        if chunk_key in chunk.chunk_keys:
            value = chunk.chunk_d[chunk_key]
            datum.append(value)
        else:
            raise KeyError("Unable to find chunk key '{i}' in {p}".format(
                i=chunk_key, p=chunk))

    return datum


def gather_csv(csv_files, output_file, skip_empty=True):
    """

    :param csv_files:
    :param output_file:
    :param skip_empty: Emtpy files with or without the header will be skipped

    :type skip_empty: bool
    :return:
    """

    header = _read_header(csv_files[0])
    #nfields = 0 if header is None else len(header)

    with open(output_file, 'w') as writer:
        if header is not None:
            writer.write(",".join(header) + "\n")
        for csv_file in csv_files:
            if not _csv_is_empty(csv_file):
                header = _read_header(csv_file)
                # should do a comparison of the headers to make sure they
                # have the same number of fields
                with open(csv_file, 'r') as f:
                    # skip header
                    _ = f.readline()
                    for record in f:
                        writer.write(record)

    log.info("successfully merged {n} files to {f}".format(
        n=len(csv_files), f=output_file))

    return output_file


def gather_report(json_files, output_file):
    """
    Combines statistics (usually raw counts) stored as JSON files.
    Data models: pbcommand.models.report
    """
    reports = [load_report_from_json(fn) for fn in json_files]
    merged = Report.merge(reports)
    with open(output_file, "w") as writer:
        writer.write(merged.to_json())
    return output_file


def gather_txt(input_files, output_file, skip_empty=False):
    """
    Very simple concatenation of text files labeled by file name.
    """
    lines = []
    for input_file in input_files:
        with open(input_file, "r") as txt:
            lines.append("### FILE {f}:\n{t}".format(f=input_file,
                                                     t=txt.read()))
    with open(output_file, "w") as out:
        out.write("\n\n\n".join(lines))
    return output_file


def gather_fofn(input_files, output_file, skip_empty=True):
    """
    This should be better spec'ed and impose a tighter constraint on the FOFN

    :param input_files: List of file paths
    :param output_file: File Path
    :param skip_empty: Ignore empty files

    :return: Output file

    :rtype: str
    """

    all_files = []
    for input_file in input_files:
        file_names = fofn_to_files(input_file)
        all_files.extend(file_names)

    with open(output_file, 'w') as f:
        f.write("\n".join(all_files))

    return output_file


def _sanitize_dataset_tags(dset):
    tags = set(dset.tags.split(","))
    if "chunked" in tags or "filtered" in tags:
        x.discard("chunked")
        x.discard("filtered")
        dset.tags = ",".join(sorted(list(dset.tags)))


def __gather_contigset(resource_file_extension, input_files, output_file,
                       new_resource_file=None,
                       skip_empty=True):
    """
    :param input_files: List of file paths
    :param output_file: File Path
    :param new_resource_file: the path of the file to which the other contig
                              files are consolidated
    :param skip_empty: Ignore empty files (doesn't do much yet)

    :return: Output file

    :rtype: str
    """
    if skip_empty:
        _input_files = []
        for file_name in input_files:
            cs = ContigSet(file_name)
            if len(cs.toExternalFiles()) > 0:
                _input_files.append(file_name)
        input_files = _input_files
    tbr = ContigSet(*input_files)
    if not new_resource_file:
        if output_file.endswith('xml'):
            new_resource_file = output_file[:-3] + resource_file_extension
    tbr.consolidate(new_resource_file)
    tbr.newUuid()
    _sanitize_dataset_tags(tbr)
    tbr.write(output_file)
    return output_file


gather_contigset = P(__gather_contigset, "fasta")


def gather_fastq_contigset(input_files, output_file):
    if len(input_files) == 1:
        shutil.copyfile(input_files[0], output_file)
    else:
        contigset_name = op.splitext(output_file)[0] + ".contigset.xml"
        __gather_contigset("fastq", input_files, contigset_name,
                           new_resource_file=output_file)
        assert op.isfile(output_file)
        return output_file


def _uniqueify_metadata(ds):
    if ds.metadata is not None and len(ds.metadata.collections) > 1:
        have_collections = set()
        k = 0
        while k < len(ds.metadata.collections):
            if ds.metadata.collections[k].uniqueId in have_collections:
                ds.metadata.collections.pop(k)
            else:
                have_collections.add(ds.metadata.collections[k].uniqueId)
                k += 1


def __gather_readset(dataset_type, input_files, output_file, skip_empty=True,
                     consolidate=False, consolidate_n_files=1):
    """
    :param input_files: List of file paths
    :param output_file: File Path
    :param skip_empty: Ignore empty files (doesn't do much yet)

    :return: Output file

    :rtype: str
    """
    tbr = dataset_type(*input_files)
    _uniqueify_metadata(tbr)
    if consolidate:
        new_resource_file = output_file[:-4] + ".bam"
        tbr.consolidate(new_resource_file, numFiles=consolidate_n_files)
        tbr.induceIndices()
    tbr.newUuid()
    _sanitize_dataset_tags(tbr)
    tbr.write(output_file)
    return output_file


gather_subreadset = P(__gather_readset, SubreadSet)
gather_alignmentset = P(__gather_readset, AlignmentSet)
gather_ccsset = P(__gather_readset, ConsensusReadSet)
gather_ccs_alignmentset = P(__gather_readset, ConsensusAlignmentSet)
gather_transcripts = P(__gather_readset, TranscriptSet)
gather_transcript_alignmentset = P(__gather_readset, TranscriptAlignmentSet)


def gather_bigwig(input_files, output_file):
    import pyBigWig
    chr_lengths = {}
    FileInfo = namedtuple("FileInfo", ("file_name", "file_id", "seqids"))
    files_info = []
    for i, file_name in enumerate(input_files):
        log.info("Reading header info from {f}...".format(f=file_name))
        if op.getsize(file_name) == 0:
            continue
        try:
            file_id = int(op.dirname(file_name).split("-")[-1])
        except ValueError:
            file_id = i
        bw_chunk = pyBigWig.open(file_name)
        seqids = []
        for (seqid, length) in bw_chunk.chroms().iteritems():
            chr_lengths.setdefault(seqid, 0)
            chr_lengths[seqid] = max(length, chr_lengths[seqid])
            seqids.append(seqid)
        files_info.append(FileInfo(file_name, file_id, seqids))
        bw_chunk.close()
    if len(files_info) == 0:
        with open(output_file, "wb") as f:
            return output_file
    bw = pyBigWig.open(output_file, "w")
    files_info.sort(lambda a, b: cmp(a.file_id, b.file_id))
    regions = OrderedDict()
    seqid_files = defaultdict(list)
    for f in files_info:
        for seqid in f.seqids:
            log.debug("{f} ({i}): {s} {l}".format(f=f.file_name,
                                                  i=f.file_id, s=seqid, l=chr_lengths[seqid]))
            regions[seqid] = chr_lengths[seqid]
            seqid_files[seqid].append(f)
    bw.addHeader([(k, v) for k, v in regions.iteritems()])
    seq_chunk = namedtuple("SeqChunk", ("starts", "ends", "values"))
    for (seqid, length) in regions.iteritems():
        log.info("Collecting values for {i}...".format(i=seqid))
        chunks = []
        k = 0
        for file_info in seqid_files[seqid]:
            log.info("Reading values from {f}".format(f=file_info.file_name))
            bw_chunk = pyBigWig.open(file_info.file_name)
            starts, ends, values = [], [], []
            chr_max = bw_chunk.chroms()[seqid]
            for i, val in enumerate(bw_chunk.values(seqid, 0, chr_max)):
                if not math.isnan(val):
                    starts.append(i)
                    ends.append(i + 1)
                    values.append(val)
                    k += 1
            chunks.append(seq_chunk(starts, ends, values))
            bw_chunk.close()
        chunks.sort(lambda a, b: cmp(a.starts[0], b.starts[0]))
        starts = list(itertools.chain(*[x.starts for x in chunks]))
        ends = list(itertools.chain(*[x.ends for x in chunks]))
        values = list(itertools.chain(*[x.values for x in chunks]))
        seqids = [seqid] * len(starts)
        log.info("Adding {i}:{s}-{e}".format(i=seqid, s=starts[0], e=ends[-1]))
        bw.addEntries(seqids, starts, ends=ends, values=values)
    bw.close()
    return output_file


def gather_tgz(input_files, output_file):
    with tarfile.open(output_file, mode="w:gz") as tgz_out:
        for tgz_file in input_files:
            with tarfile.open(tgz_file, mode="r:gz") as tgz_in:
                for member in tgz_in.getmembers():
                    tgz_out.addfile(member, tgz_in.extractfile(member.name))
    return output_file


def gather_zip(input_files, output_file):
    with ZipFile(output_file, "w", allowZip64=True) as zip_out:
        for zip_file in input_files:
            with ZipFile(zip_file, "r") as zip_in:
                for member in zip_in.namelist():
                    f = zip_in.open(member)
                    zip_out.writestr(member, f.read())
    return output_file
