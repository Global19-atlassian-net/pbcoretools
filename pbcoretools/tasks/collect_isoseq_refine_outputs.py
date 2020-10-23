"""
Collect the outputs of the 'isoseq refine' command and write a Report JSON.
"""

import logging
import json
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.models.report import Report, Attribute
from pbcommand.utils import setup_log
from pbcore.io import ConsensusReadSet

log = logging.getLogger(__name__)
__version__ = "0.1"


def run_args(args):
    with ConsensusReadSet(args.flnc_ccs) as flnc:
        samples = {rg.SampleName for rg in flnc.readGroupTable}
        samples = sorted(list(samples))
        if len(samples) > 1:
            log.warning("Expected only one sample")
        sample = ",".join(samples)
        attributes = [
            Attribute("sample_name", value=sample)
        ]
        with open(args.summary_json, "rt") as json_in:
            rpt_d = json.loads(json_in.read())
            for k, v in rpt_d.items():
                attributes.append(Attribute(k, value=v))
        r = Report("isoseq_refine",
                   attributes=attributes,
                   dataset_uuids=[flnc.uuid])
        r.write_json(args.output)
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("flnc_ccs", help="FLNC CCS Reads")
    p.add_argument("summary_json", help="Summary JSON")
    p.add_argument("-o", "--output", default="isoseq_refine.report.json",
                   help="Name of output file")
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
