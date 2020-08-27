"""
Restrict Trim Adapters workflow to a single primer (for now).
"""

import logging
import sys

from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log
from pbcommand.models.common import PacBioAlarm
from pbcore.io import BarcodeSet

from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


def _run_args(args):
    ds = BarcodeSet(args.barcodeset, strict=True)
    if len(ds) > 1:
        alarm = PacBioAlarm(
            exception=None,
            info=None,
            message="This application currently only supports a single PCR primer.  To trim adapters when multiple primers are used, please run lima on the command line.",
            name="Multiple Primers",
            severity=logging.ERROR,
            owner="python")
        alarm.to_json("alarms.json")
        return 1
    else:
        return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("barcodeset", help="BarcodeSet XML")
    return p


def _main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=_run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(_main(sys.argv))
