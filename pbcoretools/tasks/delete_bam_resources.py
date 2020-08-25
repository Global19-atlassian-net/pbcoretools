"""
Delete all BAM and associated resources in chunked dataset files, leaving
behind broken symlinks to removed_file placeholders.
"""
from pathlib import Path
import logging
import os.path as op
import os
import sys

from pbcommand.cli import pacbio_args_runner
from pbcommand.utils import setup_log
from pbcore.io import openDataSet

from pbcoretools.utils import get_base_parser

log = logging.getLogger(__name__)


def _replace_with_link(fn):
    cwd = os.getcwd()
    base_dir = op.dirname(fn)
    try:
        placeholder_fn = op.join(base_dir, "removed_file")
        os.chdir(base_dir)
        os.remove(fn)
        Path(placeholder_fn).touch()
        os.symlink(op.basename(placeholder_fn), op.basename(fn))
        return placeholder_fn
    except (IOError, OSError) as e:
        log.error(e)
        return None
    finally:
        os.chdir(cwd)


def delete_bam_resources(datasets):
    placeholder_files = set()
    removed_files = []
    for dataset in datasets:
        dataset_files = set() #[dataset])
        ds = openDataSet(dataset, skipCounts=True)
        for resource in ds.externalResources:
            for fn in [resource.bam, resource.pbi, resource.bai]:
                if fn is not None and op.isfile(fn):
                    dataset_files.add(fn)
        with open("removed_files.fofn", mode="wt") as fofn:
            for fn in list(dataset_files):
                log.info("Deleting %s", fn)
                link_file = _replace_with_link(fn)
                if link_file is not None:
                    placeholder_files.add(link_file)
                    removed_files.append(fn)
                    fofn.write(fn + "\n")
    log.info("Cleaning up placeholder files to leave broken links")
    for fn in placeholder_files:
        os.remove(fn)
    return len(removed_files)


def _run_args(args):
    delete_bam_resources(args.datasets)
    return 0


def _get_parser():
    p = get_base_parser(__doc__)
    p.add_argument("datasets", nargs="+", help="Paths to chunked datasets")
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
