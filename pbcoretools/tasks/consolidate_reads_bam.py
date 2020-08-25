
import logging
import os.path
import sys

from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log
from pbcore.io import ConsensusReadSet

from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


def _run_args(args):
    ds = ConsensusReadSet(args.ccsxml, strict=True)
    orig_uuid = ds.uuid
    ds.consolidate("reads.bam", useTmp=False)
    bam_res = ds.externalResources[0]
    if args.zmws_json:
        bam_res._setSubResByMetaType("PacBio.FileTypes.json",
                                     args.zmws_json)
    if args.report_ccs_processing:
        bam_res._setSubResByMetaType("PacBio.FileTypes.JsonReport",
                                     args.report_ccs_processing)
    ds.uuid = orig_uuid
    ds.write("final.consensusreadset.xml")
    with open("reads.fofn", "wt") as fofn:
        fofn.write(os.path.abspath("reads.bam"))
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("ccsxml", help="Path to input dataset")
    p.add_argument("--zmws-json", action="store", default=None,
                   help="Path to zmws.json.gz")
    p.add_argument("--report-ccs-processing", action="store", default=None,
                   help="Path to ccs_processing.report.json")
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
