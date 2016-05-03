
"""
Tool contract wrappers for miscellaneous quick functions.
"""

import functools
import tempfile
import logging
import shutil
import gzip
import re
import os.path as op
import os
import sys

from pbcore.io import (SubreadSet, HdfSubreadSet, FastaReader, FastaWriter,
                       FastqReader, FastqWriter, BarcodeSet, ExternalResource,
                       ExternalResources, openDataSet, ContigSet, ReferenceSet)
from pbcommand.engine import run_cmd
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, SymbolTypes

log = logging.getLogger(__name__)

TOOL_NAMESPACE = 'pbcoretools'
DRIVER_BASE = "python -m pbcoretools.tasks.converters "

registry = registry_builder(TOOL_NAMESPACE, DRIVER_BASE)

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


def add_subread_resources(self, subreads, scraps=None, barcodes=None):
    assert subreads.endswith(".subreads.bam")
    ext_res_new = ExternalResource()
    ext_res_new.resourceId = subreads
    ext_res_new.metaType = 'PacBio.SubreadFile.SubreadBamFile'
    ext_res_new.addIndices([subreads + ".pbi"])
    if scraps is not None or barcodes is not None:
        ext_res_inner = ExternalResources()
        if scraps is not None:
            assert scraps.endswith(".scraps.bam")
            ext_res_scraps = ExternalResource()
            ext_res_scraps.resourceId = scraps
            ext_res_scraps.metaType = 'PacBio.SubreadFile.ScrapsBamFile'
            ext_res_scraps.addIndices([scraps + ".pbi"])
            ext_res_inner.append(ext_res_scraps)
        if barcodes is not None:
            assert barcodes.endswith(".barcodeset.xml"), barcodes
            ext_res_barcode = ExternalResource()
            ext_res_barcode.resourceId = barcodes
            ext_res_barcode.metaType = "PacBio.DataSet.BarcodeSet"
            ext_res_inner.append(ext_res_barcode)
            ext_res_new.append(ext_res_inner)
    self.externalResources.append(ext_res_new)


def run_bam_to_bam(subread_set_file, barcode_set_file, output_file_name,
                   nproc=1, score_mode="symmetric"):
    if not score_mode in ["asymmetric", "symmetric"]:
        raise ValueError("Unrecognized score mode '{m}'".format(m=score_mode))
    bc = BarcodeSet(barcode_set_file)
    if len(bc.resourceReaders()) > 1:
        raise NotImplementedError("Multi-FASTA BarcodeSet input is not supported.")
    barcode_fasta = bc.toExternalFiles()[0]
    with SubreadSet(subread_set_file) as ds:
        ds_new = SubreadSet(strict=True)
        for ext_res in ds.externalResources:
            subreads_bam = ext_res.bam
            scraps_bam = ext_res.scraps
            assert subreads_bam is not None
            if scraps_bam is None:
                raise TypeError("The input SubreadSet must include scraps.")
            new_prefix = op.join(op.dirname(output_file_name),
                re.sub(".subreads.bam", "_barcoded", op.basename(subreads_bam)))
            if not op.isabs(subreads_bam):
                subreads_bam = op.join(op.dirname(subread_set_file),
                    subreads_bam)
            if not op.isabs(scraps_bam):
                scraps_bam = op.join(op.dirname(subread_set_file), scraps_bam)
            args = [
                "bam2bam",
                "-j", str(nproc),
                "-b", str(nproc),
                "-o", new_prefix,
                "--barcodes", barcode_fasta,
                "--scoreMode", score_mode,
                subreads_bam, scraps_bam
            ]
            log.info(" ".join(args))
            result = run_cmd(" ".join(args),
                             stdout_fh=sys.stdout,
                             stderr_fh=sys.stderr)
            if result.exit_code != 0:
                return result.exit_code
            subreads_bam = new_prefix + ".subreads.bam"
            scraps_bam = new_prefix + ".scraps.bam"
            assert op.isfile(subreads_bam), "Missing {f}".format(f=subreads_bam)
            add_subread_resources(ds_new,
                subreads=subreads_bam,
                scraps=scraps_bam,
                barcodes=barcode_set_file)
        ds._filters.clearCallbacks()
        ds_new._filters = ds._filters
        ds_new._populateMetaTypes()
        ds_new.metadata = ds.metadata
        ds_new.name = ds.name + " (barcoded)"
        ds_new.updateCounts()
        ds_new.newUuid()
        ds_new.write(output_file_name)
    return 0


def _unzip_fastx(gzip_file_name, fastx_file_name):
    with gzip.open(gzip_file_name, "rb") as gz_in:
        with open(fastx_file_name, "wb") as fastx_out:
            fastx_out.write(gz_in.read())


def _run_bam_to_fastx(program_name, fastx_reader, fastx_writer,
                     input_file_name, output_file_name):
    assert isinstance(program_name, basestring)
    barcode_mode = False
    if output_file_name.endswith(".gz"):
        with openDataSet(input_file_name) as ds_in:
            barcode_mode = ds_in.isBarcoded
    tmp_out_prefix = tempfile.NamedTemporaryFile().name
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
    if result.exit_code != 0:
        return result.exit_code
    else:
        base_ext = re.sub("bam2", "", program_name) 
        if not barcode_mode:
            tmp_out = "{p}.{b}.gz".format(p=tmp_out_prefix, b=base_ext)
            assert os.path.isfile(tmp_out), tmp_out
            if output_file_name.endswith(".gz"):
                log.info("cp {t} {f}".format(t=tmp_out, f=output_file_name))
                shutil.copyfile(tmp_out, output_file_name)
            else:
                _unzip_fastx(tmp_out, output_file_name)
            os.remove(tmp_out)
        else:
            suffix = "{f}.gz".format(f=base_ext)
            tmp_out_dir = op.dirname(tmp_out_prefix)
            tc_out_dir = op.dirname(output_file_name)
            barcoded_file_names = []
            # find the barcoded FASTX files and unzip them to the same
            # output directory and file prefix as the ultimate output
            for fn in os.listdir(tmp_out_dir):
                fn = op.join(tmp_out_dir, fn)
                if fn.startswith(tmp_out_prefix) and fn.endswith(suffix):
                    bc_fwd_rev = fn.split(".")[-3].split("_")
                    suffix2 = ".{f}_{r}.{t}".format(
                        f=bc_fwd_rev[0], r=bc_fwd_rev[1], t=base_ext)
                    assert fn == tmp_out_prefix + suffix2 + ".gz"
                    fn_out = re.sub(".gz$", suffix2, output_file_name)
                    fastx_out = op.join(tc_out_dir, fn_out)
                    _unzip_fastx(fn, fastx_out)
                    barcoded_file_names.append(op.basename(fn_out))
                    os.remove(fn)
            assert len(barcoded_file_names) > 0
            # now make a gzipped tarball
            args = ["tar", "-czf", output_file_name] + barcoded_file_names
            log.info("Running '{a}'".format(a=" ".join(args)))
            _cwd = os.getcwd()
            try:
                # we want the files to have no leading path
                os.chdir(op.dirname(output_file_name))
                result = run_cmd(" ".join(args),
                                 stdout_fh=sys.stdout,
                                 stderr_fh=sys.stderr)
            except Exception:
                raise
            else:
                if result.exit_code != 0:
                    return result.exit_code
            finally:
                os.chdir(_cwd)
            assert op.isfile(output_file_name)
    return 0


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
    # the '.py' name difference will be resolved in pbdataset/pbcoretools, but
    # for now, work with either
    if result.exit_code == 127:
        args = ["dataset.py create", "--type ReferenceSet",
                "--generateIndices",
                output_file_name, input_file_name]
        log.info(" ".join(args))
        result = run_cmd(" ".join(args), stdout_fh = sys.stdout,
                         stderr_fh=sys.stderr)
    return result.exit_code


def run_fasta_to_reference(input_file_name, output_file_name,
                           organism=None, reference_name=None,
                           ploidy="haploid"):
    if reference_name is None or reference_name == "":
        reference_name = op.splitext(op.basename(input_file_name))[0]
    ds_in = ContigSet(input_file_name)
    if len(ds_in.externalResources) > 1:
        raise TypeError("Only a single FASTA file is supported as input.")
    fasta_file_name = ds_in.externalResources[0].resourceId
    output_dir_name = op.dirname(output_file_name)
    args = [
        "fasta-to-reference",
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
    ref_file = op.join(output_dir_name, reference_name, "referenceset.xml")
    assert op.isfile(ref_file)
    with ReferenceSet(ref_file, strict=True) as ds_ref:
        ds_ref.makePathsAbsolute()
        log.info("saving final ReferenceSet to {f}".format(f=output_file_name))
        ds_ref.write(output_file_name)
    return 0


run_bam_to_fasta = functools.partial(_run_bam_to_fastx, "bam2fasta",
    FastaReader, FastaWriter)
run_bam_to_fastq = functools.partial(_run_bam_to_fastx, "bam2fastq",
    FastqReader, FastqWriter)


@registry("h5_subreads_to_subread", "0.1.0",
          FileTypes.DS_SUBREADS_H5,
          FileTypes.DS_SUBREADS, is_distributed=True, nproc=1)
def run_bax2bam(rtc):
    return run_bax_to_bam(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2bam_barcode", "0.1.0",
          (FileTypes.DS_SUBREADS, FileTypes.DS_BARCODE),
          FileTypes.DS_SUBREADS,
          is_distributed=True,
          nproc=SymbolTypes.MAX_NPROC,
          options={"score_mode":"symmetric"})
def run_bam2bam(rtc):
    return run_bam_to_bam(
        subread_set_file=rtc.task.input_files[0],
        barcode_set_file=rtc.task.input_files[1],
        output_file_name=rtc.task.output_files[0],
        nproc=rtc.task.nproc,
        score_mode=rtc.task.options["pbcoretools.task_options.score_mode"])


@registry("bam2fastq", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.FASTQ, is_distributed=True, nproc=1)
def run_bam2fastq(rtc):
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fasta", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.FASTA, is_distributed=True, nproc=1)
def run_bam2fasta(rtc):
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fasta_archive", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.GZIP, is_distributed=True, nproc=1)
def run_bam2fasta_archive(rtc):
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fastq_archive", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.GZIP, is_distributed=True, nproc=1)
def run_bam2fasta_archive(rtc):
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("fasta2fofn", "0.1.0",
          FileTypes.FASTA,
          FileTypes.FOFN, is_distributed=False, nproc=1)
def run_fasta2fofn(rtc):
    return run_fasta_to_fofn(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("fasta2referenceset", "0.1.0",
          FileTypes.FASTA,
          FileTypes.DS_REF, is_distributed=True, nproc=1)
def run_fasta2referenceset(rtc):
    return run_fasta_to_referenceset(rtc.task.input_files[0],
                                     rtc.task.output_files[0])


@registry("fasta_to_reference", "0.1.0",
          FileTypes.DS_CONTIG,
          FileTypes.DS_REF, is_distributed=True, nproc=1,
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


@registry("bam2fastq_ccs", "0.1.0",
          FileTypes.DS_CCS,
          FileTypes.GZIP, is_distributed=True, nproc=1)
def run_bam2fastq_ccs(rtc):
    """
    Duplicate of run_bam2fastq, but with ConsensusReadSet as input.
    """
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fasta_ccs", "0.1.0",
          FileTypes.DS_CCS,
          FileTypes.GZIP, is_distributed=True, nproc=1)
def run_bam2fasta_ccs(rtc):
    """
    Duplicate of run_bam2fasta, but with ConsensusReadSet as input.
    """
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0])


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
