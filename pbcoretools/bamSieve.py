
"""
Tool for subsetting a BAM or PacBio DataSet file based on either a whitelist of
hole numbers or a percentage of reads to be randomly selected.
"""

from __future__ import division
from collections import defaultdict
import subprocess
import logging
import random
import os.path as op
import sys

from pysam import AlignmentFile

from pbcommand.common_options import add_log_quiet_option
from pbcommand.cli import pacbio_args_runner, \
    get_default_argparser_with_base_opts
from pbcommand.utils import setup_log
from pbcore.io import openDataFile, openDataSet

VERSION = "0.1"

log = logging.getLogger(__name__)


def _process_zmw_list(zmw_list):
    zmws = set()
    if zmw_list is None:
        return zmws
    elif isinstance(zmw_list, set):
        return zmw_list
    elif isinstance(zmw_list, (list, tuple)):
        return set(zmw_list)
    elif op.isfile(zmw_list):
        base, ext = op.splitext(zmw_list)
        if ext in [".bam", ".xml"]:
            with openDataFile(zmw_list) as ds_zmw:
                for f in ds_zmw.resourceReaders():
                    zmws.update(set(list(f.holeNumber)))
        else:
            with open(zmw_list) as f:
                lines = f.read().splitlines()
                zmws.update(set([int(x) for x in lines]))
    else:
        zmws.update(set([int(x) for x in zmw_list.split(",")]))
    return zmws


def filter_reads(input_bam,
                 output_bam,
                 whitelist=None,
                 blacklist=None,
                 percentage=None,
                 count=None,
                 seed=None):
    if output_bam is None:
        log.error("Must specify output file")
        return 1
    n_specified = 4 - [whitelist, blacklist, percentage, count].count(None)
    if n_specified != 1:
        log.error("You must choose one and only one of the following "+
                  "options: --whitelist, --blacklist, --count, --percentage")
        return 1
    if seed is not None:
        random.seed(seed)
    if whitelist is None and blacklist is None:
        if not 0 < percentage < 100 and not count > 0:
            log.error("No reads selected for output.")
            return 1
    output_ds = None
    if output_bam.endswith(".xml"):
        if not input_bam.endswith(".xml"):
            print "DataSet output only supported for DataSet inputs."
            return 1
        output_ds = output_bam
        output_bam = op.splitext(output_ds)[0] + ".bam"
    if output_bam == input_bam:
        log.error("Input and output files must not be the same path")
        return 1
    elif not output_bam.endswith(".bam"):
        log.error("Output file name must end in either '.bam' or '.xml'")
        return 1
    n_file_reads = 0
    have_zmws = set()
    with openDataFile(input_bam) as ds_in:
        if not ds_in.isIndexed:
            log.error("Input BAM must have accompanying .pbi index")
            return 1
        f1 = ds_in.resourceReaders()[0]
        if percentage is not None or count is not None:
            with AlignmentFile(output_bam, 'wb', template=f1.peer) as bam_out:
                zmw_dict = defaultdict(list)
                for i_file, bam in enumerate(ds_in.resourceReaders()):
                    for i_read, zmw in enumerate(bam.holeNumber):
                        movie = bam.readGroupInfo(bam.qId[i_read]).MovieName
                        zmw_dict[(movie, zmw, i_file)].append(i_read)
                zmws = zmw_dict.keys()
                n_zmws_start = len(zmws)
                n_zmws_out = count
                if percentage is not None:
                    n_zmws_out = int(n_zmws_start * percentage / 100.0)
                log.info("Will random select {n} / {d} ZMWs".format(
                    n=n_zmws_out, d=n_zmws_start))
                have_reads = set()
                while True:
                    i_zmw = random.randint(0, len(zmws) - 1)
                    if not zmws[i_zmw] in have_zmws:
                        movie, zmw, i_file = zmws[i_zmw]
                        bam = ds_in.resourceReaders()[i_file]
                        for i_read in zmw_dict[zmws[i_zmw]]:
                            assert not (i_file, i_read) in have_reads
                            bam_out.write(bam[i_read].peer)
                            have_reads.add((i_file, i_read))
                            n_file_reads += 1
                        have_zmws.add(zmws[i_zmw])
                    if len(have_zmws) == n_zmws_out:
                        break
        else:
            # convert these to Python sets
            _whitelist = _process_zmw_list(whitelist)
            _blacklist = _process_zmw_list(blacklist)
            have_zmws = set()
            n_reads = 0
            with AlignmentFile(output_bam, 'wb', template=f1.peer) as bam_out:
                for f in ds_in.resourceReaders():
                    for i_zmw, zmw in enumerate(f.holeNumber):
                        if ((len(_whitelist) > 0 and zmw in _whitelist) or
                                (len(_blacklist) > 0 and not zmw in _blacklist)):
                            rec = f[i_zmw]
                            bam_out.write(rec.peer)
                            have_zmws.add(zmw)
                            n_file_reads += 1
    if n_file_reads == 0:
        log.error("No reads written")
        return 1
    log.info("{n} reads from {z} ZMWs written".format(
        n=n_file_reads, z=len(have_zmws)))
    try:
        rc = subprocess.call(["pbindex", output_bam])
    except OSError as e:
        if e.errno == 2:
            log.warn("pbindex not present, will not create .pbi file")
        else:
            raise
    if output_ds is not None:
        with openDataSet(input_bam) as ds_in:
            ds_out = ds_in.__class__(output_bam)
            ds_out.write(output_ds)
            log.info("wrote {t} XML to {x}".format(
                t=ds_in.__class__.__name__, x=output_ds))
    return 0


def show_zmws(input_file):
    zmws = []
    with openDataFile(input_file) as ds_in:
        is_indexed = ds_in.isIndexed
        if not is_indexed:
            log.warning("Unindexed file(s), this may be very slow")
        for rr in ds_in.resourceReaders():
            if is_indexed:
                zmws.extend(list([int(x) for x in rr.holeNumber]))
            else:
                zmws.extend([int(rec.HoleNumber) for rec in rr])
    print "\n".join([str(x) for x in sorted(list(set(zmws)))])


def run(args):
    if args.show_zmws:
        if [args.whitelist, args.blacklist, args.percentage].count(None) != 3:
            log.warning("Ignoring unused filtering arguments")
        show_zmws(args.input_bam)
        return 0
    return filter_reads(
        input_bam=args.input_bam,
        output_bam=args.output_bam,
        whitelist=args.whitelist,
        blacklist=args.blacklist,
        percentage=args.percentage,
        count=args.count,
        seed=args.seed)


def get_parser():
    p = add_log_quiet_option(get_default_argparser_with_base_opts(
        version=VERSION,
        description=__doc__,
        default_level="WARN"))
    p.add_argument("input_bam",
                   help="Input BAM or DataSet from which reads will be read")
    p.add_argument("output_bam", nargs='?', default=None,
                   help="Output BAM or DataSet to which filtered reads will "
                        "be written")
    p.add_argument("--show-zmws", action="store_true", default=False,
                   help="Print a list of ZMWs and exit")
    p.add_argument("--whitelist", action="store", default=None,
                   help="Comma-separated list of ZMWs, or file containing " +
                        "whitelist of one hole number per line, or " +
                        "BAM/DataSet file from which to extract ZMWs")
    p.add_argument("--blacklist", action="store", default=None,
                   help="Opposite of --whitelist, specifies ZMWs to discard")
    p.add_argument("--percentage", action="store", type=float, default=None,
                   help="If you prefer to recover a percentage of a SMRTcell "
                        "rather than a specific list of reads specify that "
                        "percentage (range 0-100) here")
    p.add_argument("-n", "--count", action="store", type=int, default=None,
                   help="Recover a specific number of ZMWs picked at random")
    p.add_argument("-s", "--seed", action="store", type=int, default=None,
                   help="Random seed for selecting a percentage of reads")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=run,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
