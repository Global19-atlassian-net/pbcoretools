# FIXME consolidate with bamSieve

"""
Write a dataset containing only unmapped reads.
"""

import subprocess
import logging
import uuid
import os.path as op
import os
import sys

try:
    from pysam.calignmentfile import AlignmentFile # pylint: disable=no-name-in-module, import-error, fixme, line-too-long
except ImportError:
    from pysam.libcalignmentfile import AlignmentFile # pylint: disable=no-name-in-module, import-error, fixme, line-too-long

from pbcommand.models import get_pbparser, FileTypes, DataStore, DataStoreFile
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log, pool_map
from pbcore.io import openDataSet, SubreadSet, PacBioBamIndex

log = logging.getLogger(__name__)


class Constants(object):
    TASK_ID = "pbcoretools.tasks.extract_unmapped_bam"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.extract_unmapped_bam --resolved-tool-contract "
    OPT_OUTPUT_ID = "pbcoretools.task_options.output_unaligned_bam"


def _get_blacklist_pbi(pbi_file):
    pbi = PacBioBamIndex(pbi_file)
    return set(zip(pbi.qId, pbi.holeNumber))


# NOTE this ignores filters entirely - this works in practice because our
# AlignmentSet XMLs don't contain filters, but if that assumption changes we
# will need to use the combined dataset index instead of .pbi files
def _get_blacklist(ds, nproc=1):
    pbi_files = [ext_res.pbi for ext_res in ds.externalResources]
    results = pool_map(_get_blacklist_pbi, pbi_files, nproc)
    blacklist = set()
    for chunk in results:
        blacklist.update(chunk)
    return blacklist


def make_unmapped_bam(alignment_file, subread_file, output_bam, nproc=1):
    log.info("Reading blacklist from %s", alignment_file)
    ds_mapped = openDataSet(alignment_file,
                            strict=True,
                            skipCounts=True)
    blacklist = _get_blacklist(ds_mapped)
    n_rec = 0
    unmapped = set()
    log.info("Reading from %s", subread_file)
    with openDataSet(subread_file, strict=True) as ds_in:
        bam_template = ds_in.resourceReaders()[0].peer
        with AlignmentFile(output_bam, 'wb', template=bam_template) as out:
            for i_rec, (qId, zmw) in enumerate(zip(ds_in.index.qId,
                                                   ds_in.index.holeNumber)):
                if not (qId, zmw) in blacklist:
                    out.write(ds_in[i_rec].peer)
                    unmapped.add((qId, zmw))
                    n_rec += 1
    log.info("Wrote %d records from %d ZMWs", n_rec, len(unmapped))
    # the index is slightly superfluous but useful for testing
    return subprocess.check_call(["pbindex", output_bam])


def run_extract_unmapped(alignment_file,
                         subread_file,
                         datastore_file,
                         extract_unaligned=True,
                         nproc=1):
    datastore_files = []
    if extract_unaligned:
        base_dir = op.dirname(datastore_file)
        output_bam = op.join(base_dir, "unaligned.subreads.bam")
        make_unmapped_bam(alignment_file, subread_file, output_bam, nproc)
        datastore_files.append(
            DataStoreFile(uuid.uuid4(),
                          Constants.TASK_ID + "-out-1",
                          FileTypes.BAM_SUB.file_type_id,
                          output_bam))
    else:
        log.warn("extract_unaligned=False, will write empty datastore")
    datastore = DataStore(datastore_files)
    datastore.write_json(datastore_file)
    return 0


def _get_parser():  # pragma: no cover
    p = get_pbparser(Constants.TASK_ID,
                     Constants.VERSION,
                     "Extract unmapped BAM",
                     __doc__,
                     Constants.DRIVER,
                     is_distributed=True,
                     nproc=8)

    p.add_input_file_type(FileTypes.DS_ALIGN,
                          "align_in",
                          "Input AlignmentSet",
                          "Gathered AlignmentSet")
    p.add_input_file_type(FileTypes.DS_SUBREADS,
                          "subreads_in",
                          "Input SubreadSet",
                          "Input (unmapped) SubreadSet")
    p.add_output_file_type(FileTypes.DATASTORE,
                           "datastore_out",
                           "JSON datastore",
                           description="Unaligned BAM datastore",
                           default_name="unaligned")
    p.tool_contract_parser.add_boolean(Constants.OPT_OUTPUT_ID,
                                       "extract_unaligned",
                                       default=False,
                                       name="Output unaligned .bam",
                                       description="Output unaligned .bam")
    return p


def args_runner(args):  # pragma: no cover
    return run_extract_unmapped(
        args.align_in,
        args.subreads_in,
        args.datastore_out)


def rtc_runner(rtc):  # pragma: no cover
    return run_extract_unmapped(
        rtc.task.input_files[0],
        rtc.task.input_files[1],
        rtc.task.output_files[0],
        extract_unaligned=rtc.task.options[Constants.OPT_OUTPUT_ID],
        nproc=rtc.task.nproc)


def main(argv=sys.argv):  # pragma: no cover
    return pbparser_runner(argv[1:],
                           _get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
