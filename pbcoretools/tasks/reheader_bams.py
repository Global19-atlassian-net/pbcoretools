"""
Utility task to set new SM and/or LB tags in the headers of every BAM resource
in a dataset.
"""

import subprocess
import logging
import os.path as op
import os
import sys

import pysam

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcore.io import openDataSet

from pbcoretools import file_utils
from pbcoretools import __VERSION__

log = logging.getLogger(__name__)


def reheader_bam(bam_file_in,
                 bam_file_out,
                 biosample_name=None,
                 library_name=None):
    """
    Write a new BAM file identical to the input except for substitution or
    addition of SM and/or LB tags in the @RG header.  If the tags are already
    present and current no file will be written.

    :return: True if header was changed, False if header is already current
    """
    # XXX https://github.com/pysam-developers/pysam/issues/939
    pysam.set_verbosity(0)  # pylint: disable=no-member
    was_changed = False
    with pysam.AlignmentFile(bam_file_in, "rb", check_sq=False) as bam_in:  # pylint: disable=no-member
        header = dict(bam_in.header)
        for rg in header["RG"]:
            if biosample_name:
                if rg.get("SM", None) != biosample_name:
                    was_changed = True
                rg["SM"] = biosample_name
            if library_name:
                if rg.get("LB", None) != library_name:
                    was_changed = True
                rg["LB"] = library_name
        if not was_changed:
            return False
        log.debug("Writing modified header and records to %s", bam_file_out)
        with pysam.AlignmentFile(bam_file_out,  # pylint: disable=no-member
                                 "wb",
                                 header=header) as bam_out:
            for rec in bam_in:
                bam_out.write(rec)
        log.debug("Running pbindex")
        subprocess.check_call(["samtools", "index", bam_file_out])
        subprocess.check_call(["pbindex", bam_file_out])
    return True


def reheader_dataset_bams(ds,
                          output_dir,
                          biosample_name=None,
                          library_name=None):
    ds = ds.copy()
    have_files = set(os.listdir(output_dir))

    def _get_ofn(fn):
        ofn = op.basename(fn)
        k = 1
        while ofn in have_files:
            ofn = op.splitext(ofn)[0] + "-{k}.bam".format(k=k)
            k += 1
        have_files.add(ofn)
        return op.join(output_dir, ofn)

    ds.close()
    for ext_res in ds.externalResources:
        if ext_res.bam:
            ofn = _get_ofn(ext_res.bam)
            log.info("Updating BAM file %s to %s", ext_res.bam, ofn)
            was_changed = reheader_bam(ext_res.bam, ofn,
                                       biosample_name=biosample_name,
                                       library_name=library_name)
            if was_changed:
                ext_res.resourceId = ofn
                ext_res.pbi = ext_res.bam + ".pbi"
                ext_res.bai = ext_res.bam + ".bai"
            else:
                log.warn("Skipped update because headers are already current")
    ds.updateCounts()
    ds.newUuid(random=True)
    return ds


def _run_args(args):
    if not args.biosample_name and not args.library_name:
        log.error("No biosample or library name specified")
        return 1
    ds_out_file = op.abspath(args.output_file)
    ds_out = reheader_dataset_bams(args.dataset,
                                   op.dirname(ds_out_file),
                                   args.biosample_name,
                                   args.library_name)
    if args.biosample_name:
        file_utils.force_set_all_bio_sample_names(ds_out, args.biosample_name)
    if args.library_name:
        file_utils.force_set_all_well_sample_names(ds_out, args.library_name)
    ds_out.write(ds_out_file)
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__VERSION__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("dataset", type=openDataSet, help="Path to input dataset")
    p.add_argument("output_file", help="Name of output dataset file")
    p.add_argument("--biosample-name", default=None, help="New BioSample Name")
    p.add_argument("--library-name", default=None, help="New Library Name")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=_run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
