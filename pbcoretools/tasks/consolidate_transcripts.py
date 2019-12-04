#!/usr/bin/env python
"""
Consolidate TranscriptSet into bam files, allowing adding a prefix
string (e.g., 'mysample_HQ_') to every transcript names.
"""

import sys
import os.path as op
import logging
import subprocess

from pysam import AlignmentFile  # pylint: disable=no-member, no-name-in-module

from pbcommand.utils import setup_log
from pbcommand.cli import pacbio_args_runner
from pbcommand.models import FileTypes, DataStore, DataStoreFile
from pbcore.io import ConsensusAlignmentSet, TranscriptAlignmentSet, TranscriptSet, openDataSet

from pbcoretools.file_utils import get_prefixes
from pbcoretools.datastore_utils import dataset_to_datastore
from pbcoretools.utils import get_base_parser


def get_parser():
    """
    Input:
        idx - 0 SubreadSet
        idx - 1 HQ TranscriptSet
        idx - 2 LQ TranscriptSet
    Output:
        idx - 0 HQ TranscriptSet, of which read names have biosample_HQ prefix
        idx - 1 LQ TranscriptSet, of which read names have biosample_LQ prefix
        idx - 2 HQ DataStore of output TranscriptSet BAM file
        idx - 3 LQ DataStore of output TranscriptSet BAM file
    """
    p = get_base_parser(__doc__)
    p.add_argument("subreads", help="SubreadSet with biosample metadata.")
    p.add_argument("hq_ds_in", help="Gathered HQ transcripts")
    p.add_argument("lq_ds_in", help="Gathered LQ transcripts")
    p.add_argument("hq_ds_out", help="Output HQ transcripts")
    p.add_argument("lq_ds_out", help="Output LQ transcripts")
    p.add_argument(
        "hq_datastore", help="Datastore containing HQ transcripts BAM")
    p.add_argument(
        "lq_datastore", help="Datastore containing LQ transcripts BAM")
    return p


class Constants:
    TOOL_ID = "consolidate_transcripts"
    BAI_FILE_TYPES = {
        FileTypes.BAMBAI.file_type_id,
        FileTypes.I_BAI.file_type_id
    }


def consolidate_transcripts(ds_in, prefix):
    """Return a function which
    - must take (new_resource_file, numFiles, useTmp) as input,
    - should consolidate ds_in (input transcripset)
    - should add biosample prefix to transcript read names
    """
    def _consolidate_transcripts_f(new_resource_file, numFiles, useTmp,
                                   perfix=prefix, ds_in=ds_in):
        external_files = ds_in.toExternalFiles()
        assert len(
            external_files) >= 1, "{!r} must contain one or more bam files".format(ds_in)
        header = AlignmentFile(external_files[0], 'rb', check_sq=False).header
        with AlignmentFile(new_resource_file, 'wb', header=header) as writer:
            for external_file in external_files:
                with AlignmentFile(external_file, 'rb', check_sq=False) as reader:
                    for record in reader:
                        record.query_name = prefix + record.query_name
                        writer.write(record)
        # create pbi and bai index files for new_resource_file
        subprocess.check_call(["pbindex", new_resource_file])
        ds_in = TranscriptSet(new_resource_file)  # override ds_in
    return _consolidate_transcripts_f


def bam_of_dataset(dataset_fn):
    return op.splitext(dataset_fn)[0] + ".bam"


def get_reads_name(ds_in):
    if isinstance(ds_in, TranscriptAlignmentSet):
        return 'Aligned transcripts'
    if isinstance(ds_in, ConsensusAlignmentSet):
        return 'Aligned consensus reads'
    return 'Aligned reads'


def run_consolidate(dataset_file, output_file, datastore_file,
                    consolidate, n_files,
                    consolidate_f=lambda ds: ds.consolidate):
    datastore_files = []
    with openDataSet(dataset_file) as ds_in:
        if consolidate:
            if len(ds_in.toExternalFiles()) <= 0:
                raise ValueError(
                    "DataSet {} must contain one or more files!".format(dataset_file))
            new_resource_file = bam_of_dataset(output_file)
            consolidate_f(ds_in)(new_resource_file,
                                 numFiles=n_files, useTmp=False)
            # always display the BAM/BAI if consolidation is enabled
            # XXX there is no uniqueness constraint on the sourceId, but this
            # seems sloppy nonetheless - unfortunately I don't know how else to
            # make view rule whitelisting work
            reads_name = get_reads_name(ds_in)
            for ext_res in ds_in.externalResources:
                if ext_res.resourceId.endswith(".bam"):
                    ds_file = DataStoreFile(
                        ext_res.uniqueId,
                        Constants.TOOL_ID + "-out-2",
                        ext_res.metaType,
                        ext_res.bam,
                        name=reads_name,
                        description=reads_name)
                    datastore_files.append(ds_file)
                    # Prevent duplicated index files being added to datastore, since consolidated
                    # dataset may contain multiple indices pointing to the same physical file
                    added_resources = set()
                    for index in ext_res.indices:
                        if (index.metaType in Constants.BAI_FILE_TYPES and
                                index.resourceId not in added_resources):
                            added_resources.add(index.resourceId)
                            ds_file = DataStoreFile(
                                index.uniqueId,
                                Constants.TOOL_ID + "-out-3",
                                index.metaType,
                                index.resourceId,
                                name="Index of {}".format(reads_name.lower()),
                                description="Index of {}".format(reads_name.lower()))
                            datastore_files.append(ds_file)
        ds_in.newUuid()
        ds_in.write(output_file)
    datastore = DataStore(datastore_files)
    datastore.write_json(datastore_file)
    return 0


def __runner(ds_items):
    for ds_in, ds_out, datastore, prefix in ds_items:
        def func(ds_in):
            return consolidate_transcripts(ds_in, prefix=prefix)
        run_consolidate(dataset_file=ds_in,
                        output_file=ds_out,
                        datastore_file=datastore,
                        consolidate=True,
                        n_files=1,
                        consolidate_f=func)
        # At this piont, ds_out is the same as ds_in, override ds_out with
        # newly created, read name modified TranscriptSet
        new_resource_file = bam_of_dataset(ds_out)
        _ds_out = TranscriptSet(new_resource_file)
        _ds_out.newUuid()
        _ds_in = TranscriptSet(ds_in)
        _ds_out.tags = _ds_in.tags
        _ds_out.name = _ds_in.name
        _ds_out.write(ds_out)
        # At this piont datastore contains paths to bam/bai/pbi files, now override
        # datastore with TranscriptSet
        dataset_to_datastore(ds_out, datastore, source_id=Constants.TOOL_ID)
    return 0


def args_runner(args):
    hq_prefix, lq_prefix = get_prefixes(args.subreads)
    ds_items = [
        (args.hq_ds_in, args.hq_ds_out, args.hq_datastore, hq_prefix),
        (args.lq_ds_in, args.lq_ds_out, args.lq_datastore, lq_prefix)
    ]
    return __runner(ds_items)


def main(argv=sys.argv):
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger()
    return pacbio_args_runner(argv[1:],
                              get_parser(),
                              args_runner,
                              log,
                              setup_log)


if __name__ == '__main__':
    sys.exit(main())
