#!/usr/bin/env python

"""
Convert a mapped .pbi file to CSV for loading into R.  This is a workaround for
a bug in some unrolled alignment files that crashes pbbamr.
"""

import logging
import csv
import os
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcore.io import openDataFile

log = logging.getLogger(__name__)
__version__ = "0.1"


HEADERS = ["idx","file","hole","qstart","qend","qual","offset","flag","ref","tstart","tend","astart","aend","rc","matches","mismatches","mq","inserts","dels"]
HEADERS_BC = ["bcf", "bcr", "bcqual"]
HEADERS_SNR = ["snrA", "snrC", "snrG", "snrT"]
HEADERS_NPASSES = ["np"]

def _write_csv(rows, csv_file, headers=HEADERS):
    with open(csv_file, "wt") as csv_out:
        writer = csv.writer(csv_out,
                            delimiter=",",
                            lineterminator=os.linesep,
                            quotechar='"',
                            quoting=csv.QUOTE_MINIMAL)
        writer.writerow(headers)
        for i, row in enumerate(rows, start=1):
            row = [str(i)] + list(row)
            writer.writerow([str(item) for item in row])


def run_args(args):
    ds = openDataFile(args.dataset)
    get_full_bam = args.load_snr or args.load_numpasses
    is_barcoded = ds.isBarcoded
    headers = list(HEADERS)
    if is_barcoded:
        headers += HEADERS_BC
    if args.load_snr:
        headers += HEADERS_SNR
    if args.load_numpasses:
        headers += HEADERS_NPASSES
    rows = []
    for rr in ds.resourceReaders():
        identity = rr.pbi.identity
        for i, holeNumber in enumerate(rr.pbi.holeNumber):
            reference = rr.referenceInfo(rr.pbi.tId[i])[2]
            aLen = rr.pbi.aEnd[i] - rr.pbi.aStart[i]
            if aLen <= 0 or identity[i] < 0:
                log.warning("ZMW %s has negative-length alignment or negative computed identity, skipping", holeNumber)
                continue
            rc = "FALSE"
            if rr.pbi.isReverseStrand[i]:
                rc = "TRUE"
            row = [rr.filename,
                   holeNumber,
                   rr.pbi.qStart[i],
                   rr.pbi.qEnd[i],
                   rr.pbi.readQual[i],
                   rr.pbi.virtualFileOffset[i],
                   rr.pbi.contextFlag[i],
                   reference,
                   rr.pbi.tStart[i],
                   rr.pbi.tEnd[i],
                   rr.pbi.aStart[i],
                   rr.pbi.aEnd[i],
                   rc,
                   rr.pbi.nM[i],
                   rr.pbi.nMM[i],
                   rr.pbi.mapQV[i],
                   rr.pbi.nIns[i],
                   rr.pbi.nDel[i]]
            if is_barcoded:
                row.extend([
                    rr.pbi.bcForward[i],
                    rr.pbi.bcReverse[i],
                    rr.pbi.bcQual[i]
                ])
            if get_full_bam:
                rec = rr[i]
                if args.load_snr:
                    snr = rec.peer.get_tag("sn")
                    row.extend(snr)
                if args.load_numpasses:
                    row.append(rec.peer.get_tag("np"))
            rows.append(row)
    _write_csv(rows, args.csv_out, headers=headers)
    log.info("Wrote %s", args.csv_out)
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("dataset", help="PacBio dataset XML (AlignmentSet or ConsensusAlignmentSet)")
    p.add_argument("csv_out", help="CSV output file")
    p.add_argument("--load-snr", action="store_true", default=False,
                   help="Include per-read SNRs")
    p.add_argument("--load-numpasses", action="store_true", default=False,
                   help="Include numPasses (CCS only)")
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
