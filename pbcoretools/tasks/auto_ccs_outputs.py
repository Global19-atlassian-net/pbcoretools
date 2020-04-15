"""
Generate single-file CCS BAM and FASTQ outputs from a ConsensusReadSet.
"""

import tempfile
import logging
import uuid
import sys
import os.path as op
import re

import numpy as np

from pbcommand.models import FileTypes, DataStoreFile, DataStore
from pbcommand.cli import pacbio_args_runner, get_default_argparser_with_base_opts
from pbcommand.utils import setup_log
from pbcore.io import ConsensusReadSet
from pbcore.util.statistics import accuracy_as_phred_qv, phred_qv_as_accuracy

from pbcoretools.bam2fastx import run_bam_to_fastq, run_bam_to_fasta
from pbcoretools.filters import combine_filters
from pbcoretools.utils import get_base_parser
from pbcoretools import __VERSION__

log = logging.getLogger(__name__)


class Constants:
    HIFI_RQ = 0.99
    BAM_ID = "ccs_bam_out"
    FASTA_ID = "ccs_fasta_out"
    FASTQ_ID = "ccs_fastq_out"
    FASTA2_ID = "ccs_fasta_lq_out"
    FASTQ2_ID = "ccs_fastq_lq_out"
    FASTA_FILE_IDS = [FASTA_ID, FASTA2_ID]
    FASTQ_FILE_IDS = [FASTQ_ID, FASTQ2_ID]


# SL-5606: if it's really HiFi we should label it accordingly, otherwise
# indicate QV in file name
def _get_suffix_and_label(min_rq):
    if min_rq >= Constants.HIFI_RQ:
        return ("hifi_reads", "HiFi Reads")
    else:
        qv = "Q{q}".format(q=int(np.round(accuracy_as_phred_qv(min_rq))))
        return ("{q}_reads".format(q=qv.lower()), "≥{q} Reads".format(q=qv))


def _to_datastore_file(file_name, file_id, file_type, description):
    return DataStoreFile(uuid.uuid4(),
                         file_id,
                         file_type.file_type_id,
                         op.abspath(file_name),
                         name=op.basename(file_name),
                         description=description)


def consolidate_bam(base_dir, file_prefix, dataset,
                    min_rq=Constants.HIFI_RQ):
    ds_out = dataset.copy()
    suffix, label = _get_suffix_and_label(min_rq)
    bam_file_name = op.join(base_dir, ".".join([file_prefix, suffix, "bam"]))
    combine_filters(ds_out, {"rq": [('>=', min_rq)]})
    ds_out.consolidate(bam_file_name)
    return _to_datastore_file(file_name=bam_file_name,
                              file_id=Constants.BAM_ID,
                              file_type=FileTypes.BAM_CCS,
                              description="{l} BAM file".format(l=label))


def _run_bam2fastx(file_type, dataset_file, fastx_file):
    if file_type == FileTypes.FASTQ:
        return run_bam_to_fastq(dataset_file, fastx_file)
    else:
        return run_bam_to_fasta(dataset_file, fastx_file)


def to_fastx_files(file_type,
                   ds,
                   ccs_dataset_file,
                   file_ids,
                   base_dir,
                   file_prefix,
                   min_rq=Constants.HIFI_RQ,
                   no_zip=False):

    def run_to_output(ccs_file, file_id, rq):
        suffix, label = _get_suffix_and_label(rq)
        base_name = ".".join([file_prefix, suffix, file_type.ext])
        fastx_file = op.join(base_dir, base_name)
        if not no_zip:
            fastx_file += ".zip"
        _run_bam2fastx(file_type, ccs_file, fastx_file)
        desc = "{l} ({t})".format(l=label, t=file_type.ext.upper())
        return _to_datastore_file(fastx_file, file_id, file_type, desc)

    ccs_hifi = ccs_dataset_file
    is_all_hifi = np.all(ds.index.readQual >= min_rq)
    if not is_all_hifi:
        ccs_hq, ccs_lq = ds.copy(), ds.copy()
        ccs_hifi = tempfile.NamedTemporaryFile(
            suffix=".consensusreadset.xml").name
        combine_filters(ccs_hq, {"rq": [('>=', min_rq)]})
        ccs_hq.write(ccs_hifi)
    return [run_to_output(ccs_hifi, file_ids[0], min_rq)]


def get_prefix_and_bam_file_name(ds, is_barcoded=False):
    bam_file_name = file_prefix = None
    if is_barcoded:
        assert len(ds.externalResources) == 1
        bam = ds.resourceReaders()[0]
        barcodes = list(set(zip(bam.pbi.bcForward, bam.pbi.bcReverse)))
        assert len(barcodes) == 1, "Multiple barcodes found in {f}: {b}".format(
            f=ds.fileNames[0], b=", ".join([str(b) for b in barcodes]))
        bam_file_name = ds.externalResources[0].bam
        log.info("Found a single barcoded BAM %s", bam_file_name)
        # we need to handle both '.bam' and '.ccs.bam'
        file_prefix = re.sub(".bam$", "",
                             re.sub(".ccs.bam$", "", op.basename(bam_file_name)))
    else:
        movies = sorted(list({rg.MovieName for rg in ds.readGroupTable}))
        file_prefix = "_".join(movies)
        if len(movies) > 1:
            log.warning("Multiple movies found: %s", movies)
            file_prefix = "multiple_movies"
    return bam_file_name, file_prefix


def run_ccs_bam_fastq_exports(ccs_dataset_file, base_dir, is_barcoded=False,
                              min_rq=Constants.HIFI_RQ, no_zip=False):
    """
    Take a ConsensusReadSet and write BAM/FASTQ files to the output
    directory.  If this is a demultiplexed dataset, it is assumed to have
    a single BAM file within a dataset that is already imported in SMRT Link.
    Note that this function runs the exports serially, and is therefore no
    longer used in this specific task, but rather in the barcoded version that
    runs in parallel.
    """
    datastore_files = []
    with ConsensusReadSet(ccs_dataset_file, strict=True) as ds:
        bam_file_name, file_prefix = get_prefix_and_bam_file_name(
            ds, is_barcoded)
        if bam_file_name is None:
            datastore_files.append(consolidate_bam(base_dir, file_prefix, ds, min_rq))
        fasta_file_ids = [Constants.FASTA_ID, Constants.FASTA2_ID]
        fastq_file_ids = [Constants.FASTQ_ID, Constants.FASTQ2_ID]
        datastore_files.extend(
            to_fastx_files(FileTypes.FASTA,
                           ds,
                           ccs_dataset_file,
                           fasta_file_ids,
                           base_dir,
                           file_prefix,
                           min_rq=min_rq,
                           no_zip=no_zip))
        datastore_files.extend(
            to_fastx_files(FileTypes.FASTQ,
                           ds,
                           ccs_dataset_file,
                           fastq_file_ids,
                           base_dir,
                           file_prefix,
                           min_rq=min_rq,
                           no_zip=no_zip))
    return datastore_files


def run_args(args):
    datastore_out = op.abspath(args.datastore_out)
    base_dir = op.dirname(datastore_out)
    datastore_files = []
    with ConsensusReadSet(args.dataset_file, strict=True) as ds:
        bam_file_name, file_prefix = get_prefix_and_bam_file_name(
            ds, is_barcoded=False)
        if args.mode == "fasta":
            datastore_files.extend(to_fastx_files(
                FileTypes.FASTA, ds, args.dataset_file, Constants.FASTA_FILE_IDS, base_dir, file_prefix, args.min_rq, no_zip=args.no_zip))
        elif args.mode == "fastq":
            datastore_files.extend(to_fastx_files(
                FileTypes.FASTQ, ds, args.dataset_file, Constants.FASTQ_FILE_IDS, base_dir, file_prefix, args.min_rq, no_zip=args.no_zip))
        elif args.mode == "consolidate":
            if bam_file_name is None:
                datastore_files.append(
                    consolidate_bam(base_dir, file_prefix, ds))
    DataStore(datastore_files).write_json(datastore_out)
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__VERSION__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("mode", choices=["consolidate", "fasta", "fastq"])
    p.add_argument("dataset_file")
    p.add_argument("datastore_out")
    p.add_argument("--min-rq",
                   dest="min_rq",
                   type=float,
                   default=Constants.HIFI_RQ,
                   help="Sets RQ cutoff for splitting output")
    p.add_argument("--min-qv",
                   dest="min_rq",
                   type=lambda arg: phred_qv_as_accuracy(int(arg)),
                   help="Alternative to --min-rq, on Phred scale (0-60)")
    p.add_argument("--no-zip",
                   action="store_true",
                   help="Disable ZIP output")
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
