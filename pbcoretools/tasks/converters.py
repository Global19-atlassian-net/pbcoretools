
"""
Tool contract wrappers for miscellaneous quick functions.
"""

from collections import defaultdict
from zipfile import ZipFile
import subprocess
import itertools
import functools
import tempfile
import logging
import shutil
import gzip
import copy
import re
import os.path as op
import os
import sys

from pbcore.io import (SubreadSet, HdfSubreadSet, FastaReader, FastaWriter,
                       FastqReader, FastqWriter, BarcodeSet, ExternalResource,
                       ExternalResources, openDataSet, ContigSet, ReferenceSet,
                       GmapReferenceSet)
from pbcommand.engine import run_cmd
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, SymbolTypes, OutputFileType, DataStore, DataStoreFile
from pbcommand.utils import walker

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_NAMESPACE = 'pbcoretools'
    DRIVER_BASE = "python -m pbcoretools.tasks.converters "

    # For large references setting this to give more overhead for
    # memory usage. The underlying tool should have a well defined
    # memory usage that is independent of reference size (if possible)
    DEFAULT_FASTA_CONVERT_MAX_NPROC = 4

    # default filter applied to output of 'lima'
    BARCODE_QUALITY_GREATER_THAN = 26
    ALLOWED_BC_TYPES = set([f.file_type_id for f in
                            [FileTypes.DS_SUBREADS, FileTypes.DS_CCS]])


registry = registry_builder(Constants.TOOL_NAMESPACE, Constants.DRIVER_BASE)


def _run_bax_to_bam(input_file_name, output_file_name):
    base_name = ".".join(output_file_name.split(".")[:-2])
    input_file_name_tmp = input_file_name
    # XXX bax2bam won't write an hdfsubreadset unless the input is XML too
    if input_file_name.endswith(".bax.h5"):
        input_file_name_tmp = tempfile.NamedTemporaryFile(
            suffix=".hdfsubreadset.xml").name
        ds_tmp = HdfSubreadSet(input_file_name)
        ds_tmp.write(input_file_name_tmp)
    args =[
        "bax2bam",
        "--subread",
        "-o", base_name,
        "--output-xml", output_file_name,
        "--xml", input_file_name_tmp
    ]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args),
                     stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code
    with SubreadSet(output_file_name) as ds:
        ds.assertIndexed()
    return 0


def run_bax_to_bam(input_file_name, output_file_name):
    with HdfSubreadSet(input_file_name) as ds_in:
        movies = set()
        for rr in ds_in.resourceReaders():
            movies.add(rr.movieName)
        if len(movies) > 1:
            out_dir = os.path.dirname(output_file_name)
            ds_out_files = []
            for bax_file in ds_in.toExternalFiles():
                output_file_name_tmp = os.path.join(out_dir, ".".join(
                    os.path.basename(bax_file).split(".")[:-2]) +
                    ".hdfsubreadset.xml")
                rc = _run_bax_to_bam(bax_file, output_file_name_tmp)
                if rc != 0:
                    log.error("bax2bam failed")
                    return rc
                ds_out_files.append(output_file_name_tmp)
            ds = SubreadSet(*ds_out_files)
            ds.name = ds_in.name
            if 'Description' in ds_in.objMetadata:
                ds.objMetadata['Description'] = ds_in.objMetadata['Description']
                ds.metadata.merge(ds_in.metadata)
            ds.write(output_file_name)
        else:
            return _run_bax_to_bam(input_file_name, output_file_name)
    return 0


# XXX no longer used
def _run_bam_to_bam(subread_set_file, barcode_set_file, output_file_name,
                    nproc=1, score_mode="symmetric"):
    """
    Run legacy barcoding with bam2bam (requires separate installation).
    """
    if not score_mode in ["symmetric", "asymmetric", "tailed"]:
        raise ValueError("Unrecognized score mode '{m}'".format(m=score_mode))
    bc = BarcodeSet(barcode_set_file)
    if len(bc.resourceReaders()) > 1:
        raise NotImplementedError("Multi-FASTA BarcodeSet input is not supported.")
    new_prefix = re.sub(".subreadset.xml$", "", output_file_name)
    args = [
        "bam2bam",
        "-j", str(nproc),
        "-b", str(nproc),
        "-o", new_prefix,
        "--barcodes", barcode_set_file,
        "--scoreMode", score_mode,
        subread_set_file
    ]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args),
                     stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code
    assert op.isfile(output_file_name)
    tmp_out = op.join(op.dirname(output_file_name),
                      "tmp_" + op.basename(output_file_name))
    shutil.move(output_file_name, tmp_out)
    with SubreadSet(tmp_out, strict=True) as ds:
        with SubreadSet(subread_set_file) as ds_in:
            ds.metadata = ds_in.metadata
            ds_name_out = ds_in.name
            if not "(barcoded)" in ds_name_out:
                ds_name_out = ds_name_out + " (barcoded)"
            ds.name = ds_name_out
            if ds.tags.strip() == "":
                ds.tags = "barcoded"
            elif not "barcoded" in ds.tags:
                ds.tags = ",".join(ds.tags.strip().split(",") + ["barcoded"])
        ds.updateCounts()
        ds.newUuid()
        ds.write(output_file_name)
    return 0


def _ungzip_fastx(gzip_file_name, fastx_file_name):
    """
    Decompress an output from bam2fastx.
    """
    with gzip.open(gzip_file_name, "rb") as gz_in:
        with open(fastx_file_name, "wb") as fastx_out:
            def _fread():
                return gz_in.read(1024)
            for chunk in iter(_fread, ''):
                fastx_out.write(chunk)


def archive_files(input_file_names, output_file_name, remove_path=True):
    """
    Create a zipfile from a list of input files.

    :param remove_path: if True, the directory will be removed from the input
                        file names before archiving.  All inputs and the output
                        file must be in the same directory for this to work.
    """
    archive_file_names = input_file_names
    if remove_path:
        archive_file_names = [op.basename(fn) for fn in archive_file_names]
    log.info("Creating zip file %s", output_file_name)
    with ZipFile(output_file_name, "w", allowZip64=True) as zip_out:
        for file_name, archive_file_name in zip(input_file_names,
                                                archive_file_names):
            zip_out.write(file_name, archive_file_name)
    return 0


def _run_bam_to_fastx(program_name, fastx_reader, fastx_writer,
                     input_file_name, output_file_name, tmp_dir=None):
    """
    Converts a dataset to a set of fastx file, possibly archived.
    Can take a subreadset or consensusreadset as input.
    Will convert to either fasta or fastq.
    If the dataset is barcoded, it will split the fastx files per-barcode.
    If the output file is .zip, the fastx file(s) will be archived accordingly.
    """
    assert isinstance(program_name, basestring)
    barcode_mode = False
    barcode_sets = set()
    if output_file_name.endswith(".zip"):
        with openDataSet(input_file_name) as ds_in:
            barcode_mode = ds_in.isBarcoded
            if barcode_mode:
                # attempt to collect the labels of barcodes used on this
                # dataset.  assumes that all BAM files used the same barcodes
                for bam in ds_in.externalResources:
                    if bam.barcodes is not None:
                        barcode_sets.add(bam.barcodes)
    barcode_labels = []
    if barcode_mode:
        if len(barcode_sets) == 1:
            bc_file = list(barcode_sets)[0]
            log.info("Reading barcode labels from %s", bc_file)
            try:
                with BarcodeSet(bc_file) as bc_in:
                    for bc in bc_in:
                        barcode_labels.append(bc.id)
            except IOError as e:
                log.error("Can't read %s", bc_file)
                log.error(e)
        elif len(barcode_sets) > 1:
            log.warn("Multiple barcode sets used for this SubreadSet:")
            for fn in sorted(list(barcode_sets)):
                log.warn("  %s", fn)
        else:
            log.info("No barcode labels available")
    base_ext = "." + re.sub("bam2", "", program_name)
    suffix = "{f}.gz".format(f=base_ext)
    tmp_out_prefix = tempfile.NamedTemporaryFile(dir=tmp_dir).name
    tmp_out_dir = op.dirname(tmp_out_prefix)
    args = [
        program_name,
        "-o", tmp_out_prefix,
        input_file_name,
    ]
    if barcode_mode:
        args.insert(1, "--split-barcodes")
    log.info(" ".join(args))
    result = run_cmd(" ".join(args),
                     stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    def _is_fastx_file(fn):
        return fn.startswith(tmp_out_prefix) and fn.endswith(suffix)
    try:
        if result.exit_code != 0:
            return result.exit_code
        else:
            if output_file_name.endswith(".zip"):
                tc_out_dir = op.dirname(output_file_name)
                fastx_file_names = []
                # find the barcoded FASTX files and un-gzip them to the same
                # output directory and file prefix as the ultimate output
                for fn in walker(tmp_out_dir, _is_fastx_file):
                    if barcode_mode:
                        # bam2fastx outputs files with the barcode indices
                        # encoded in the file names; here we attempt to
                        # translate these to barcode labels, falling back on
                        # the original indices if necessary
                        bc_fwd_rev = fn.split(".")[-3].split("_")
                        bc_label = "unbarcoded"
                        if (bc_fwd_rev != ["65535", "65535"] and
                            bc_fwd_rev != ["-1", "-1"]):
                            def _label_or_none(x):
                                try:
                                    bc = int(x)
                                    if bc < 0:
                                        return "none"
                                    elif bc < len(barcode_labels):
                                        return barcode_labels[bc]
                                except ValueError as e:
                                    pass
                                return x
                            bc_fwd_label = _label_or_none(bc_fwd_rev[0])
                            bc_rev_label = _label_or_none(bc_fwd_rev[1])
                            bc_label = "{f}__{r}".format(f=bc_fwd_label,
                                                         r=bc_rev_label)
                        suffix2 = ".{l}{t}".format(l=bc_label, t=base_ext)
                    else:
                        suffix2 = '.' + base_ext
                    fn_out = re.sub(".zip$", "", op.basename(output_file_name))
                    if fn_out.endswith(base_ext):
                        fn_out = re.sub(base_ext, suffix2, fn_out)
                    fastx_out = op.join(tc_out_dir, fn_out)
                    _ungzip_fastx(fn, fastx_out)
                    fastx_file_names.append(fastx_out)
                assert len(fastx_file_names) > 0
                return archive_files(fastx_file_names, output_file_name)
            else:
                tmp_out = "{p}{b}.gz".format(p=tmp_out_prefix, b=base_ext)
                _ungzip_fastx(tmp_out, output_file_name)
                os.remove(tmp_out)
    finally:
        for fn in walker(tmp_out_dir, _is_fastx_file):
            os.remove(fn)
    return 0


def split_laa_fastq(input_file_name, output_file_base):
    """
    Split an LAA FASTQ file into one file per barcode.
    """
    if op.getsize(input_file_name) == 0:
        return []
    records = defaultdict(list)
    with FastqReader(input_file_name) as fastq_in:
        for rec in fastq_in:
            bc_id = rec.id.split("_")[0]
            records[bc_id].append(rec)
    outputs = []
    for bc_id in sorted(records.keys()):
        ofn = "{b}.{i}.fastq".format(b=output_file_base, i=bc_id)
        with FastqWriter(ofn) as fastq_out:
            for rec in records[bc_id]:
                fastq_out.writeRecord(rec)
        outputs.append(ofn)
    return outputs


def split_laa_fastq_archived(input_file_name, output_file_name):
    """
    Split an LAA FASTQ file into one file per barcode and package as zip.
    """
    base, ext = op.splitext(output_file_name)
    assert (ext == ".zip")
    fastq_files = list(split_laa_fastq(input_file_name, base))
    if len(fastq_files) == 0: # workaround for empty input
        with ZipFile(output_file_name, "w", allowZip64=True) as zip_out:
            return 0
    return archive_files(fastq_files, output_file_name)


def run_fasta_to_fofn(input_file_name, output_file_name):
    args = ["echo", input_file_name, ">", output_file_name]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh = sys.stdout,
                     stderr_fh=sys.stderr)
    return result.exit_code


def run_fasta_to_referenceset(input_file_name, output_file_name):
    args = ["dataset create", "--type ReferenceSet", "--generateIndices",
            output_file_name, input_file_name]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh = sys.stdout,
                     stderr_fh=sys.stderr)
    return result.exit_code


def __run_fasta_to_reference(program_name, dataset_class,
                             input_file_name, output_file_name,
                             organism=None, reference_name=None,
                             ploidy="haploid"):
    if reference_name is None or reference_name == "":
        reference_name = op.splitext(op.basename(input_file_name))[0]
    ds_in = ContigSet(input_file_name)

    # For historical reasons, there's some munging between the ReferenceSet
    # name and the sub directories that are created within the job dir.
    rx = re.compile('[^A-Za-z0-9_]')
    sanitized_name = re.sub(rx, '_', reference_name)

    if len(ds_in.externalResources) > 1:
        raise TypeError("Only a single FASTA file is supported as input.")
    fasta_file_name = ds_in.externalResources[0].resourceId
    output_dir_name = op.dirname(output_file_name)
    args = [
        program_name,
        "--organism", str(organism) if organism != "" else "unknown",
        "--ploidy", str(ploidy) if ploidy != "" else "unknown",
        "--debug",
        fasta_file_name,
        output_dir_name,
        reference_name
    ]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh=sys.stdout, stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code

    # By convention, the dataset XML is written to this path
    dataset_xml = op.join(output_dir_name, sanitized_name,
                          "{t}.xml".format(t=dataset_class.__name__.lower()))

    log.info("Looking for DataSet XML `{}`".format(dataset_xml))
    assert op.isfile(dataset_xml), dataset_xml
    with dataset_class(dataset_xml, strict=True) as ds_ref:
        ds_ref.makePathsAbsolute()
        log.info("saving final {t} to {f}".format(
                 f=output_file_name, t=dataset_class.__name__))
        ds_ref.write(output_file_name)
    return 0


run_fasta_to_reference = functools.partial(__run_fasta_to_reference,
    "fasta-to-reference", ReferenceSet)
run_fasta_to_gmap_reference = functools.partial(__run_fasta_to_reference,
    "fasta-to-gmap-reference", GmapReferenceSet)


run_bam_to_fasta = functools.partial(_run_bam_to_fastx, "bam2fasta",
    FastaReader, FastaWriter)
run_bam_to_fastq = functools.partial(_run_bam_to_fastx, "bam2fastq",
    FastqReader, FastqWriter)


subreads_from_h5_file_type = OutputFileType(FileTypes.DS_SUBREADS.file_type_id,
                                            "Subreads", "Subread data in XML dataset",
                                            "Imported SubreadSet", "subreads")
subreads_barcoded_file_type = OutputFileType(FileTypes.DS_SUBREADS.file_type_id,
                                             "SubreadSet",
                                             "Barcoded Subreads",
                                             "Barcoded Subreads DataSet XML",
                                             "subreads_barcoded")

@registry("h5_subreads_to_subread", "0.1.0",
          FileTypes.DS_SUBREADS_H5,
          subreads_from_h5_file_type, is_distributed=True, nproc=1)
def run_bax2bam(rtc):
    return run_bax_to_bam(rtc.task.input_files[0], rtc.task.output_files[0])


fasta_file_type = OutputFileType(FileTypes.FASTA.file_type_id, "fasta", "FASTA file",
                                 "Reads in FASTA format", "reads")
fastq_file_type = OutputFileType(FileTypes.FASTQ.file_type_id, "fastq", "FASTQ file",
                                 "Reads in FASTQ format", "reads")


@registry("bam2fastq", "0.1.0",
          FileTypes.DS_SUBREADS,
          fastq_file_type, is_distributed=True, nproc=1)
def run_bam2fastq(rtc):
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


fofn_file_type = OutputFileType(FileTypes.FOFN.file_type_id, "FOFN file",
                                "FOFN file", "List of input files", "files")

@registry("fasta2fofn", "0.1.0",
          FileTypes.FASTA,
          fofn_file_type, is_distributed=False, nproc=1)
def run_fasta2fofn(rtc):
    return run_fasta_to_fofn(rtc.task.input_files[0], rtc.task.output_files[0])


ref_file_type = OutputFileType(FileTypes.DS_REF.file_type_id, "ReferenceSet",
                               "Reference Dataset",
                               "PacBio Reference DataSet XML", "reference")

@registry("fasta2referenceset", "0.1.0",
          FileTypes.FASTA,
          ref_file_type,
          is_distributed=True,
          nproc=Constants.DEFAULT_FASTA_CONVERT_MAX_NPROC)
def run_fasta2referenceset(rtc):
    return run_fasta_to_referenceset(rtc.task.input_files[0],
                                     rtc.task.output_files[0])


@registry("fasta_to_reference", "0.1.0",
          FileTypes.FASTA,
          ref_file_type,
          is_distributed=True,
          nproc=Constants.DEFAULT_FASTA_CONVERT_MAX_NPROC,
          options={
                "organism": "",
                "ploidy": "haploid",
                "reference_name":""
          })
def run_fasta_to_reference_pbscala(rtc):
    return run_fasta_to_reference(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        reference_name=rtc.task.options["pbcoretools.task_options.reference_name"],
        organism=rtc.task.options["pbcoretools.task_options.organism"],
        ploidy=rtc.task.options["pbcoretools.task_options.ploidy"])


gmap_ref_file_type = OutputFileType(FileTypes.DS_GMAP_REF.file_type_id, "GmapReferenceSet",
                               "GmapReferenceSet XML",
                               "PacBio GMAP Reference DataSet XML", "reference")

@registry("fasta_to_gmap_reference", "0.1.0",
          FileTypes.FASTA,
          gmap_ref_file_type, is_distributed=True,
          nproc=Constants.DEFAULT_FASTA_CONVERT_MAX_NPROC,
          options={
                "organism": "",
                "ploidy": "haploid",
                "reference_name":""
          })
def _run_fasta_to_gmap_reference(rtc):
    return run_fasta_to_gmap_reference(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        reference_name=rtc.task.options["pbcoretools.task_options.reference_name"],
        organism=rtc.task.options["pbcoretools.task_options.organism"],
        ploidy=rtc.task.options["pbcoretools.task_options.ploidy"])


consensus_zip_ftype = OutputFileType(FileTypes.ZIP.file_type_id,
                                     "fastq_split_zip",
                                     "Consensus Amplicons",
                                     "Consensus amplicons in FASTQ format, split by barcode",
                                     "consensus_fastq")
chimera_zip_ftype = OutputFileType(FileTypes.ZIP.file_type_id,
                                   "fastq_split_zip",
                                   "Chimeric/Noise Sequences by barcode",
                                   "Chimeric and noise sequences in FASTQ format, split by barcode",
                                   "chimera_fastq")

@registry("split_laa_fastq", "0.3.0",
          (FileTypes.FASTQ, FileTypes.FASTQ),
          (consensus_zip_ftype, chimera_zip_ftype),
          is_distributed=True, nproc=1)
def _run_split_laa_fastq(rtc):
    # XXX a bit of a hack to support unique file names for the FASTQ tarballs
    return max(split_laa_fastq_archived(rtc.task.input_files[0],
                                        rtc.task.output_files[0]),
               split_laa_fastq_archived(rtc.task.input_files[1],
                                        rtc.task.output_files[1]))


@registry("contigset2fasta", "0.1.0",
          FileTypes.DS_CONTIG,
          FileTypes.FASTA,
          is_distributed=True, nproc=1)
def contigset_to_fasta(rtc):
    with ContigSet(rtc.task.input_files[0]) as ds_in:
        if len(ds_in.externalResources) != 1:
            raise ValueError("This task assumes that the ContigSet contains "+
                             "only a single FASTA file.")
        file_name = ds_in.externalResources[0].resourceId
        os.symlink(file_name, rtc.task.output_files[0])
    return 0


def _run_slimbam(ext_res, nproc=1, exe="slimbam"):
    bam_in = ext_res.resourceId
    base, ext = op.splitext(op.splitext(op.basename(bam_in))[0])
    base_out = "{b}_slimmed".format(b=base)
    bam_out = base_out + ext + ".bam"
    args = [exe, "-o", base_out, "-j", str(nproc), bam_in]
    log.debug("COMMAND: " + " ".join(args))
    assert subprocess.call(args) == 0
    assert op.isfile(bam_out)
    ext_res.resourceId = bam_out
    log.info("Updated path: {o}".format(o=bam_out))
    for fi in ext_res.indices:
        if fi.metaType == FileTypes.I_PBI.file_type_id:
            args = ["pbindex", bam_out]
            log.debug("COMMAND: " + " ".join(args))
            assert subprocess.call(args) == 0
            pbi_out = bam_out + ".pbi"
            assert op.isfile(pbi_out)
            fi.resourceId = pbi_out
            log.info("Updated path: {o}".format(o=pbi_out))
        elif fi.metaType == FileTypes.I_BAI.file_type_id:
            args = ["samtools", "index", bam_out]
            log.debug("COMMAND: " + " ".join(args))
            assert subprocess.call(args) == 0
            bai_out = bam_out + ".bai"
            assert op.isfile(bai_out)
            fi.resourceId = bai_out
            log.info("Updated path: {o}".format(o=bai_out))
        else:
            local_file = op.basename(ext_res.resourceId)
            if not op.exists(local_file):
                shutil.copyfile(ext_res.resourceId, local_file)
                ext_res.resourceId = local_file
                log.info("Updated path: {o}".format(o=local_file))
    for er in ext_res.externalResources:
        if er.resourceId.endswith(".bam"):
            _run_slimbam(er, nproc=nproc, exe=exe)
    return ext_res


@registry("slimbam", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.DS_SUBREADS,
          is_distributed=True,
          nproc=SymbolTypes.MAX_NPROC,
          options={"slimbam_exe": "slimbam"})
def run_slimbam(rtc):
    log.warn("This task is for internal testing only; please do not use in customer-facing pipelines.")
    os.chdir(op.dirname(rtc.task.output_files[0]))
    with SubreadSet(rtc.task.input_files[0], strict=True) as ds:
        for er in ds.externalResources:
            if er.resourceId.endswith(".bam"):
                _run_slimbam(er, nproc=rtc.task.nproc, exe=rtc.task.options['pbcoretools.task_options.slimbam_exe'])
    ds.newUuid()
    ds.updateCounts()
    ds.write(rtc.task.output_files[0])
    return 0


def _iterate_datastore_read_set_files(datastore_file):
    """
    Iterate over SubreadSet or ConsensusReadSet files listed in a datastore JSON.
    """
    ds = DataStore.load_from_json(datastore_file)
    files = ds.files.values()
    for f in files:
        if f.file_type_id in Constants.ALLOWED_BC_TYPES:
            yield f


@registry("datastore_to_subreads", "0.2.0",
          FileTypes.DATASTORE,
          FileTypes.DS_SUBREADS,
          is_distributed=False,
          nproc=1)
def run_datastore_to_subreads(rtc):
    datasets = list(_iterate_datastore_read_set_files(rtc.task.input_files[0]))
    if len(datasets) > 0:
        with SubreadSet(*[f.path for f in datasets], strict=True) as ds:
            ds.newUuid()
            ds.write(rtc.task.output_files[0])
    else:
        raise ValueError("Expected one or more SubreadSets in datastore")
    return 0


def discard_bio_samples(subreads, barcode_label):
    """
    Remove any BioSample records from a SubreadSet that are not associated
    with the specified barcode.
    """
    for collection in subreads.metadata.collections:
        deletions = []
        for k, bio_sample in enumerate(collection.wellSample.bioSamples):
            barcodes = set([bc.name for bc in bio_sample.DNABarcodes])
            if barcode_label in barcodes:
                continue
            if len(barcodes) == 0:
                log.warn("No barcodes defined for sample %s", bio_sample.name)
            deletions.append(k)
        for k in reversed(deletions):
            collection.wellSample.bioSamples.pop(k)
        if len(collection.wellSample.bioSamples) == 0:
            log.warn("Collection %s has no BioSamples", collection.context)
            log.warn("Will create new BioSample and DNABarcode records")
            collection.wellSample.bioSamples.addSample(barcode_label)
            collection.wellSample.bioSamples[0].DNABarcodes.addBarcode(barcode_label)


def get_ds_name(ds, base_name, barcode_label):
    """
    Given the base (parent) dataset name, add a suffix indicating sample
    """
    suffix = "(unknown sample)"
    try:
        collection = ds.metadata.collections[0]
        n_samples = len(collection.wellSample.bioSamples)
        if n_samples == 1:
            suffix = "(%s)" % collection.wellSample.bioSamples[0].name
        elif n_samples > 1:
            suffix = "(multiple samples)"
        else:
            raise IndexError("No BioSample records present")
    except IndexError:
        if barcode_label is not None:
            suffix = "({l})".format(l=barcode_label)
    return "{n} {s}".format(n=base_name, s=suffix)


def update_barcoded_sample_metadata(base_dir,
                                    datastore_file,
                                    input_reads,
                                    barcode_set,
                                    isoseq_mode=False):
    """
    Given a datastore JSON of SubreadSets produced by barcoding, apply the
    following updates to each:
    1. Include only the BioSample(s) corresponding to its barcode
    2. Add the BioSample name to the dataset name
    3. Add a ParentDataSet record in the Provenance section.
    """
    datastore_files = []
    barcode_names = []
    with BarcodeSet(barcode_set) as bc_in:
        for rec in bc_in:
            barcode_names.append(rec.id)
    parent_ds = openDataSet(input_reads)
    for f in _iterate_datastore_read_set_files(datastore_file):
        ds_out = op.join(base_dir, op.basename(f.path))
        with openDataSet(f.path, strict=True) as ds:
            assert ds.datasetType in Constants.ALLOWED_BC_TYPES, ds.datasetType
            barcode_label = None
            ds_barcodes = sorted(list(set(zip(ds.index.bcForward, ds.index.bcReverse))))
            if isoseq_mode:
                ds_barcodes = sorted(list(set([tuple(sorted(bcs)) for bcs in ds_barcodes])))
            if len(ds_barcodes) == 1:
                bcf, bcr = ds_barcodes[0]
                barcode_label = "{f}--{r}".format(f=barcode_names[bcf],
                                                  r=barcode_names[bcr])
                try:
                    discard_bio_samples(ds, barcode_label)
                except Exception as e:
                    log.error(e)
                    log.warn("Continuing anyway, but results may not be "
                             "displayed correctly in SMRT Link")
            else:
                raise IOError(
                    "The file {f} contains multiple barcodes: {b}".format(
                    f=f.path, b="; ".join([str(bc) for bc in ds_barcodes])))
            ds.metadata.addParentDataSet(parent_ds.uuid,
                                         parent_ds.datasetType,
                                         createdBy="AnalysisJob",
                                         timeStampedName="")
            ds.name = get_ds_name(ds, parent_ds.name, barcode_label)
            ds.filters.addRequirement(
                bq=[('>', Constants.BARCODE_QUALITY_GREATER_THAN)])
            ds.newUuid()
            ds.write(ds_out)
            f_new = copy.deepcopy(f)
            f_new.path = ds_out
            f_new.uuid = ds.uuid
            datastore_files.append(f_new)
    return DataStore(datastore_files)


@registry("update_barcoded_sample_metadata", "0.3.0",
          (FileTypes.JSON, FileTypes.DS_SUBREADS, FileTypes.DS_BARCODE),
          FileTypes.DATASTORE,
          is_distributed=False,
          nproc=1)
def _run_update_barcoded_sample_metadata(rtc):
    base_dir = op.dirname(rtc.task.output_files[0])
    datastore = update_barcoded_sample_metadata(
        base_dir=op.dirname(rtc.task.output_files[0]),
        datastore_file=rtc.task.input_files[0],
        input_reads=rtc.task.input_files[1],
        barcode_set=rtc.task.input_files[2],
        isoseq_mode=False)
    datastore.write_json(rtc.task.output_files[0])
    return 0


@registry("update_barcoded_sample_metadata_ccs", "0.1.0",
          (FileTypes.JSON, FileTypes.DS_CCS, FileTypes.DS_BARCODE),
          FileTypes.DATASTORE,
          is_distributed=False,
          nproc=1)
def _run_update_barcoded_sample_metadata(rtc):
    base_dir = op.dirname(rtc.task.output_files[0])
    datastore = update_barcoded_sample_metadata(
        base_dir=op.dirname(rtc.task.output_files[0]),
        datastore_file=rtc.task.input_files[0],
        input_reads=rtc.task.input_files[1],
        barcode_set=rtc.task.input_files[2],
        isoseq_mode=True)
    datastore.write_json(rtc.task.output_files[0])
    return 0


ds_name_opt = QuickOpt("", "Name of Output Data Set",
                       "Name of new demultiplexed data set as it appears in "+
                       "SMRT Link")

@registry("reparent_subreads", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.DS_SUBREADS,
          is_distributed=False,
          nproc=1,
          options={"new_dataset_name":ds_name_opt})
def _run_reparent_subreads(rtc):
    NAME_OPT_ID = "pbcoretools.task_options.new_dataset_name"
    if rtc.task.options[NAME_OPT_ID].strip() == "":
        raise ValueError("New dataset name is required")
    with SubreadSet(rtc.task.input_files[0], strict=True) as ds_in:
        if len(ds_in.metadata.provenance) > 0:
            log.warn("Removing existing provenance record: %s",
                     ds_in.metadata.provenance)
            ds_in.metadata.provenance = None
        ds_in.name = rtc.task.options[NAME_OPT_ID]
        ds_in.newUuid(random=True)
        ds_in.write(rtc.task.output_files[0])
    return 0


def _ds_to_datastore(dataset_file, datastore_file,
                     source_id="pbcoretools.tasks.converters-out-0"):
    with openDataSet(dataset_file, strict=True) as ds:
        ds_file = DataStoreFile(ds.uniqueId, source_id, ds.datasetType, dataset_file)
        ds_out = DataStore([ds_file])
        ds_out.write_json(datastore_file)
    return 0


@registry("subreads_to_datastore", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.JSON,
          is_distributed=False,
          nproc=1)
def _run_subreads_to_datastore(rtc):
    return _ds_to_datastore(rtc.task.input_files[0],
                            rtc.task.output_files[0],
                            source_id=rtc.task.task_id + "-out-0")


@registry("ccs_to_datastore", "0.1.0",
          FileTypes.DS_CCS,
          FileTypes.JSON,
          is_distributed=False,
          nproc=1)
def _run_ccs_to_datastore(rtc):
    return _ds_to_datastore(rtc.task.input_files[0],
                            rtc.task.output_files[0],
                            source_id=rtc.task.task_id + "-out-0")


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
