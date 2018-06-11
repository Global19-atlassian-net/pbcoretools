
"""
Tool contract wrappers for utility tasks in LAA pipelines.
"""

import logging
import os.path as op
import os
import sys

from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, OutputFileType

from pbcoretools.file_utils import (split_laa_fastq_archived,
                                    make_combined_laa_zip)


log = logging.getLogger(__name__)


class Constants(object):
    TOOL_NAMESPACE = 'pbcoretools'
    DRIVER_BASE = "python -m pbcoretools.tasks.laa "


registry = registry_builder(Constants.TOOL_NAMESPACE, Constants.DRIVER_BASE)


consensus_zip_ftype = OutputFileType(
    FileTypes.ZIP.file_type_id,
    "fastq_split_zip",
    "Consensus Amplicons",
    "Consensus amplicons in FASTQ format, split by barcode",
    "consensus")
chimera_zip_ftype = OutputFileType(
    FileTypes.ZIP.file_type_id,
    "fastq_split_zip",
    "Chimeric/Noise Sequences by barcode",
    "Chimeric and noise sequences in FASTQ format, split by barcode",
    "chimera")


@registry("split_laa_fastq", "0.4",
          (FileTypes.FASTQ, FileTypes.FASTQ, FileTypes.DS_SUBREADS),
          (consensus_zip_ftype, chimera_zip_ftype),
          is_distributed=True, nproc=1)
def _run_split_laa_fastq(rtc):
    # XXX a bit of a hack to support unique file names for the FASTQ tarballs
    return max(split_laa_fastq_archived(rtc.task.input_files[0],
                                        rtc.task.output_files[0],
                                        rtc.task.input_files[2]),
               split_laa_fastq_archived(rtc.task.input_files[1],
                                        rtc.task.output_files[1],
                                        rtc.task.input_files[2]))


combined_zip_ftype = OutputFileType(FileTypes.ZIP.file_type_id,
                                    "consensus_combined_zip",
                                    "Consensus Sequences Summary",
                                    "Consensus Sequences Summary ZIP file",
                                    "consensus_sequences_summary")


@registry("make_combined_laa_zip", "0.1.3",
          (FileTypes.FASTQ, FileTypes.CSV, FileTypes.DS_SUBREADS),
          combined_zip_ftype,
          is_distributed=True,
          nproc=1)
def _run_make_combined_laa_zip(rtc):
    return make_combined_laa_zip(
        fastq_file=rtc.task.input_files[0],
        summary_csv=rtc.task.input_files[1],
        input_subreads=rtc.task.input_files[2],
        output_file_name=rtc.task.output_files[0])


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
