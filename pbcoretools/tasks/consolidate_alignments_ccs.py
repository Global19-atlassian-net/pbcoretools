
"""
Consolidate ConsensusAlignmentSet .bam files
"""

import logging
import sys

from pbcommand.models import FileTypes
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.tasks.consolidate_alignments import Constants as BaseConstants, get_parser, args_runner, rtc_runner


class Constants(BaseConstants):
    TOOL_ID = "pbcoretools.tasks.consolidate_alignments_ccs"
    VERSION = "0.2.0"
    DRIVER = "python -m pbcoretools.tasks.consolidate_alignments_ccs --resolved-tool-contract "


def main(argv=sys.argv):
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger()
    p = get_parser(
        tool_id=Constants.TOOL_ID,
        file_type=FileTypes.DS_ALIGN_CCS,
        driver_exe=Constants.DRIVER,
        version=Constants.VERSION,
        description=__doc__)
    return pbparser_runner(argv[1:],
                           p,
                           lambda rtc: args_runner(rtc, Constants.TOOL_ID),
                           lambda rtc: rtc_runner(rtc, Constants.TOOL_ID),
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
