"""
Demultiplex output files from joint clustering of multiple Iso-Seq samples,
where all transcripts or isoforms may belong to more than one sample.
"""
# FIXME all this should be rewritten in C++

from collections import defaultdict
import logging
import csv
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcommand.models.common import DataStore
from pbcore.io import ConsensusReadSet, TranscriptSet, FastaReader, FastaWriter

from pbcoretools.utils import get_base_parser
from pbcoretools.tasks.isoseq.collect_files import (FILE_IDS_AND_NAMES,
                                                    to_datastore_file)

log = logging.getLogger(__name__)

########################################################################
# FIXME this is copy-paste from pbreports.io.isoseq


def read_isoseq3_refine_flnc(flnc_ccs):
    """
    Read the header and index of the full-length non-chimeric (FLNC) CCS reads
    output by the 'isoseq3 refine' command, and return a dictionary mapping CCS
    read IDs to sample names.
    """
    flnc_read_samples = {}
    ds = ConsensusReadSet(flnc_ccs, strict=True)
    for rg in ds.readGroupTable:
        sel = ds.index.qId == rg.ID
        zmws = ds.index.holeNumber[sel]
        for zmw in zmws:
            qname = "{m}/{z}/ccs".format(m=rg.MovieName, z=zmw)
            assert not qname in flnc_read_samples, qname
            flnc_read_samples[qname] = rg.SampleName
    return flnc_read_samples


def read_isoseq3_collapsed_isoform_reads(tsv_file):
    """
    Parse the file collapse_isoforms.read_stat.txt from 'isoseq3 collapse' and
    return a dict of lists of FLNC read IDs keyed by isoform ID.
    """
    isoform_reads = defaultdict(list)
    with open(tsv_file, mode="rt") as tsv_in:
        reader = csv.reader(tsv_in, delimiter="\t")
        header = next(reader)
        assert header == "id length pbid".split()
        for rec in reader:
            isoform_reads[rec[2]].append(rec[0])
    return isoform_reads


def read_isoseq3_cluster_report_csv(csv_file):
    """
    Parse the CSV file (unpolished.cluster_report.csv) written by the
    'isoseq3 cluster' command, and return a dictionary mapping transcript IDs
    to CCS read qNames.
    """
    transcript_reads = defaultdict(list)
    with open(csv_file, mode="rt") as csv_in:
        reader = csv.reader(csv_in, delimiter=",")
        header = next(reader)
        assert header == ["cluster_id", "read_id", "read_type"], header
        for (cluster_id, read_id, read_type) in reader:
            transcript_reads[cluster_id].append(read_id)
    return transcript_reads
# end copy-paste
########################################################################


def _get_sample_to_read_mapping(read_stats, flnc_samples):
    isoform_reads = read_isoseq3_collapsed_isoform_reads(read_stats)
    sample_isoforms = defaultdict(list)
    for isoform_id, read_ids in isoform_reads.items():
        isoform_samples = set([])
        for read_id in read_ids:
            isoform_samples.add(flnc_samples[read_id])
        for sample_name in isoform_samples:
            sample_isoforms[sample_name].append(isoform_id)
    return sample_isoforms


def demultiplex_transcripts(transcripts_file,
                            transcript_reads,
                            flnc_samples,
                            sample_names,
                            transcript_type):
    """
    Split a TranscriptSet into per-sample FASTA files.
    """
    # this stores the record index of transcripts associated with each
    # sample
    sample_transcripts = defaultdict(list)
    sample_files = {}
    with TranscriptSet(transcripts_file, strict=True) as ds:
        log.info("collecting transcript record IDs per-sample")
        for i, zmw in enumerate(ds.index.holeNumber):
            transcript_id = "transcript/{z}".format(z=zmw)
            read_names = transcript_reads[transcript_id]
            samples = {flnc_samples[r] for r in read_names}
            for sample in samples:
                log.debug("Associated transcript %s with sample %s",
                          transcript_id, sample)
                sample_transcripts[sample].append(i)
        #header = dict(ds.resourceReaders()[0].peer.header)
        _cache = {}
        for j, sample in enumerate(sample_names, start=1):
            fasta_file = "unpolished-{j}.{t}.fasta".format(
                j=j, t=transcript_type)
            log.info("writing per-sample %s transcripts for '%s' to %s",
                     transcript_type, sample, fasta_file)
            with FastaWriter(fasta_file) as fasta_out:
                for i in sample_transcripts[sample]:
                    if i in _cache:
                        qname, sequence = _cache[i]
                    else:
                        rec = ds[i]
                        qname, sequence = rec.qName, rec.peer.query_sequence
                        _cache[i] = (qname, sequence)
                    header = "{q} sample:{s}".format(q=qname, s=sample)
                    fasta_out.writeRecord(header, sequence)
            sample_files[sample] = fasta_file
    return sample_files


def demultiplex_collapsed_isoforms(fasta_file,
                                   sample_isoforms,
                                   sample_names):
    """
    Split a collapsed isoforms FASTA into per-sample files.
    """
    collapsed_sequences = {}
    with FastaReader(fasta_file) as fasta_in:
        for rec in fasta_in:
            isoform_id = rec.id.split("|")[0]
            collapsed_sequences[isoform_id] = (rec.header, rec.sequence)
    sample_files = {}
    for j, sample in enumerate(sample_names, start=1):
        fasta_file = "collapsed-{j}.fasta".format(j=j)
        log.info("writing per-sample collapsed isoforms for '%s' to %s",
                 sample, fasta_file)
        with FastaWriter(fasta_file) as fasta_out:
            for isoform_id in sorted(sample_isoforms[sample]):
                header, sequence = collapsed_sequences[isoform_id]
                fasta_out.writeRecord(header, sequence)
        sample_files[sample] = fasta_file
    return sample_files


def demultiplex_transcripts_hqlq(transcripts_hq,
                                 transcripts_lq,
                                 cluster_report_csv,
                                 flnc_samples,
                                 sample_names):
    transcript_reads = read_isoseq3_cluster_report_csv(cluster_report_csv)
    fasta_hq = demultiplex_transcripts(transcripts_hq,
                                       transcript_reads,
                                       flnc_samples,
                                       sample_names,
                                       "hq")
    fasta_lq = demultiplex_transcripts(transcripts_lq,
                                       transcript_reads,
                                       flnc_samples,
                                       sample_names,
                                       "lq")
    return fasta_hq, fasta_lq


def run_args(args):
    files = {}
    log.info("Getting sample names from FLNC reads")
    flnc_samples = read_isoseq3_refine_flnc(args.flnc_ccs)
    sample_names = sorted(set(flnc_samples.values()))
    if args.transcripts_hq:
        assert args.cluster_csv is not None, "--cluster-csv argument required"
        files["hq_fasta"], files["lq_fasta"] = \
            demultiplex_transcripts_hqlq(args.transcripts_hq,
                                         args.transcripts_lq,
                                         args.cluster_csv,
                                         flnc_samples,
                                         sample_names)
    sample_isoforms = None
    if args.read_stats is not None:
        sample_isoforms = _get_sample_to_read_mapping(
            args.read_stats, flnc_samples)
    if args.collapse_fasta:
        assert sample_isoforms is not None, "--read-stats argument required"
        files["collapse_fasta"] = demultiplex_collapsed_isoforms(
            args.collapse_fasta,
            sample_isoforms,
            sample_names)
    for i, sample in enumerate(sample_names, start=1):
        datastore_files = []
        for file_id, file_type, label in FILE_IDS_AND_NAMES:
            if file_id in files:
                f = to_datastore_file(
                    file_name=files[file_id][sample],
                    file_id=file_id,
                    file_type=file_type,
                    label="{l} ({s})".format(l=label, s=sample))
                datastore_files.append(f)
        ds_file = "{p}-{i}.datastore.json".format(
            p=args.datastore_prefix, i=i)
        DataStore(datastore_files).write_json(ds_file)
        log.info("Wrote files for '%s' to %s", sample, ds_file)
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("flnc_ccs", help="Full-Length Non-Chimeric CCS Reads")
    p.add_argument("--cluster-csv", default=None, help="Cluster Report CSV")
    p.add_argument("--transcripts-hq", default=None, help="HQ TranscriptSet")
    p.add_argument("--transcripts-lq", default=None, help="LQ TranscriptSet")
    p.add_argument("--collapse-fasta", default=None, help="Collapsed Isoforms")
    #p.add_argument("--collapse-gff", default=None, help="Collapsed Isoform GFF")
    p.add_argument("--read-stats", default=None,
                   help="Collapsed Isoform Read Info")
    p.add_argument("--datastore-prefix",
                   default="sample",
                   help="Prefix for output per-sample datastore.json files")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
