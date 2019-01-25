
"""
Tool contract wrappers for miscellaneous file conversions.
"""

import subprocess
import tempfile
import functools
import logging
import shutil
import re
import os.path as op
import os
import sys

from pbcore.io import (SubreadSet, ContigSet, ConsensusReadSet, openDataSet,
                       ReferenceSet, GmapReferenceSet, HdfSubreadSet, TranscriptSet,
                       AlignmentSet, ConsensusAlignmentSet, TranscriptAlignmentSet)
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, SymbolTypes, OutputFileType, DataStore, DataStoreFile
from pbcommand.engine import run_cmd

from pbcoretools.tasks.barcoding import _ds_to_datastore
from pbcoretools.file_utils import iterate_datastore_read_set_files

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
    args = [
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


def _run_bax_to_bam_multi_dataset(input_file_name, output_file_name):
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
                ds.objMetadata['Description'] = ds_in.objMetadata[
                    'Description']
                ds.metadata.merge(ds_in.metadata)
            ds.write(output_file_name)
        else:
            return _run_bax_to_bam(input_file_name, output_file_name)
    return 0


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
def _run_bax2bam(rtc):
    return _run_bax_to_bam_multi_dataset(rtc.task.input_files[0], rtc.task.output_files[0])


fasta_file_type = OutputFileType(FileTypes.FASTA.file_type_id, "fasta", "FASTA file",
                                 "Reads in FASTA format", "reads")
fastq_file_type = OutputFileType(FileTypes.FASTQ.file_type_id, "fastq", "FASTQ file",
                                 "Reads in FASTQ format", "reads")


def _run_fasta_to_fofn(input_file_name, output_file_name):
    args = ["echo", input_file_name, ">", output_file_name]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    return result.exit_code


def _run_fasta_to_referenceset(input_file_name, output_file_name):
    args = ["dataset create", "--type ReferenceSet", "--generateIndices",
            output_file_name, input_file_name]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh=sys.stdout,
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
    rx = re.compile('[^A-Za-z0-9_\-\.]')
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
    result = run_cmd(" ".join(args), stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
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


fofn_file_type = OutputFileType(FileTypes.FOFN.file_type_id, "FOFN file",
                                "FOFN file", "List of input files", "files")


@registry("fasta2fofn", "0.1.0",
          FileTypes.FASTA,
          fofn_file_type, is_distributed=False, nproc=1)
def _run_fasta2fofn(rtc):
    return _run_fasta_to_fofn(rtc.task.input_files[0], rtc.task.output_files[0])


ref_file_type = OutputFileType(FileTypes.DS_REF.file_type_id, "ReferenceSet",
                               "Reference Dataset",
                               "PacBio Reference DataSet XML", "reference")


@registry("fasta2referenceset", "0.1.0",
          FileTypes.FASTA,
          ref_file_type,
          is_distributed=True,
          nproc=Constants.DEFAULT_FASTA_CONVERT_MAX_NPROC)
def _run_fasta2referenceset(rtc): # pragma: no cover
    return _run_fasta_to_referenceset(rtc.task.input_files[0],
                                      rtc.task.output_files[0])


@registry("fasta_to_reference", "0.1.0",
          FileTypes.FASTA,
          ref_file_type,
          is_distributed=True,
          nproc=Constants.DEFAULT_FASTA_CONVERT_MAX_NPROC,
          options={
              "organism": "",
              "ploidy": "haploid",
              "reference_name": ""
          })
def _run_fasta_to_reference_pbscala(rtc): # pragma: no cover
    return run_fasta_to_reference(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        reference_name=rtc.task.options[
            "pbcoretools.task_options.reference_name"],
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
              "reference_name": ""
          })
def _run_fasta_to_gmap_reference(rtc): # pragma: no cover
    return run_fasta_to_gmap_reference(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        reference_name=rtc.task.options[
            "pbcoretools.task_options.reference_name"],
        organism=rtc.task.options["pbcoretools.task_options.organism"],
        ploidy=rtc.task.options["pbcoretools.task_options.ploidy"])


@registry("contigset2fasta", "0.1.0",
          FileTypes.DS_CONTIG,
          FileTypes.FASTA,
          is_distributed=True, nproc=1)
def _contigset_to_fasta(rtc): # pragma: no cover
    with ContigSet(rtc.task.input_files[0]) as ds_in:
        if len(ds_in.externalResources) != 1:
            raise ValueError("This task assumes that the ContigSet contains " +
                             "only a single FASTA file.")
        file_name = ds_in.externalResources[0].resourceId
        os.symlink(file_name, rtc.task.output_files[0])
    return 0


def _run_slimbam(ext_res, nproc=1, exe="slimbam"): # pragma: no cover
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
def _run_slimbam_rtc(rtc): # pragma: no cover
    log.warn(
        "This task is for internal testing only; please do not use in customer-facing pipelines.")
    os.chdir(op.dirname(rtc.task.output_files[0]))
    with SubreadSet(rtc.task.input_files[0], strict=True) as ds:
        for er in ds.externalResources:
            if er.resourceId.endswith(".bam"):
                _run_slimbam(er, nproc=rtc.task.nproc, exe=rtc.task.options[
                             'pbcoretools.task_options.slimbam_exe'])
    ds.newUuid()
    ds.updateCounts()
    ds.write(rtc.task.output_files[0])
    return 0


def __run_datastore_to_dataset(rtc, filetype, readcls):
    datasets = list(iterate_datastore_read_set_files(rtc.task.input_files[0], [filetype]))
    if len(datasets) > 0:
        with readcls(*[f.path for f in datasets], strict=True, skipCounts=True) as ds:
            ds.newUuid()
            ds.write(rtc.task.output_files[0])
    else:
        raise ValueError("Expected one or more {} in datastore".format(readcls.__name__))
    return 0


# internal only
@registry("datastore_to_transcripts", "0.2.1",
          FileTypes.JSON,
          FileTypes.DS_TRANSCRIPT,
          is_distributed=False,
          nproc=1)
def run_datastore_to_transcripts(rtc):
    return __run_datastore_to_dataset(rtc, FileTypes.DS_TRANSCRIPT.file_type_id,
                                      TranscriptSet)


# internal only
@registry("datastore_to_alignments", "0.2.1",
          FileTypes.JSON,
          FileTypes.DS_ALIGN,
          is_distributed=False,
          nproc=1)
def run_datastore_to_alignments(rtc):
    return __run_datastore_to_dataset(rtc, FileTypes.DS_ALIGN.file_type_id,
                                      AlignmentSet)


# internal only
@registry("datastore_to_ccs_alignments", "0.2.1",
          FileTypes.JSON,
          FileTypes.DS_ALIGN_CCS,
          is_distributed=False,
          nproc=1)
def run_datastore_to_ccs_alignments(rtc):
    return __run_datastore_to_dataset(rtc, FileTypes.DS_ALIGN_CCS.file_type_id,
                                      ConsensusAlignmentSet)


# internal only
@registry("datastore_to_transcript_alignments", "0.2.1",
          FileTypes.JSON,
          FileTypes.DS_ALIGN_TRANSCRIPT,
          is_distributed=False,
          nproc=1)
def run_datastore_to_transcript_alignments(rtc):
    return __run_datastore_to_dataset(rtc, FileTypes.DS_ALIGN_TRANSCRIPT.file_type_id,
                                      TranscriptAlignmentSet)


@registry("transcripts_to_datastore", "0.1.2",
          FileTypes.DS_TRANSCRIPT,
          FileTypes.JSON,
          is_distributed=False,
          nproc=1)
def _run_transcripts_to_datastore(rtc):
    return _ds_to_datastore(rtc.task.input_files[0],
                            rtc.task.output_files[0],
                            source_id=rtc.task.task_id + "-out-0")


@registry("update_consensus_reads", "0.1.1",
          FileTypes.DS_CCS,
          FileTypes.DS_CCS,
          is_distributed=False, # requires skipCounts=True
          nproc=1,
          options=dict(use_run_design_uuid=False))
def _run_update_consensus_reads(rtc):
    with ConsensusReadSet(rtc.task.input_files[0], skipCounts=True) as ds:
        run_design_uuid = None
        if rtc.task.options["pbcoretools.task_options.use_run_design_uuid"]:
            uuids = set([])
            for collection in ds.metadata.collections:
                if collection.consensusReadSetRef is not None:
                    uuids.add(collection.consensusReadSetRef.uuid)
            if len(uuids) == 1:
                run_design_uuid = list(uuids)[0]
            elif len(uuids) == 0:
                log.warn("No pre-defined ConsensusReadSetRef UUID found")
            else:
                log.warn("Multiple ConsensusReadSetRef UUIDs found")
        if run_design_uuid is not None:
            ds.uuid = run_design_uuid
        else:
            ds.newUuid()
        ds.name = re.sub(" (filtered)", "", ds.name)
        tags = [t.strip() for t in ds.tags.split(",")]
        if "hidden" in tags:
            tags.remove("hidden")
        ds.tags = ",".join(tags)
        ds.write(rtc.task.output_files[0])
    return 0


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
