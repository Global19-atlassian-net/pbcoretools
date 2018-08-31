#/usr/bin/env python

from __future__ import print_function

from collections import defaultdict
from functools import reduce
import argparse
import logging
import string
import re
import os

from pbcore.io import DataSet, ContigSet, openDataSet, openDataFile
from pbcore.io.dataset import InvalidDataSetIOError
from pbcore.io.dataset.utils import _swapPath
from pbcore.io.dataset.DataSetMembers import Filters, OPMAP
from pbcore.io.dataset.DataSetValidator import validateFile
from pbcommand.validators import validate_output_dir

from pbcoretools.file_utils import (add_mock_collection_metadata,
                                    force_set_all_well_sample_names,
                                    force_set_all_bio_sample_names)

log = logging.getLogger(__name__)


def show_sample_names_if_defined(ds):
    if ds.metadata.collections:
        bio_samples = set()
        well_samples = set()
        for coll in ds.metadata.collections:
            well_samples.add(coll.wellSample.name)
            bio_samples.update({s.name for s in coll.wellSample.bioSamples})
        well_samples = sorted(list(well_samples))
        bio_samples = sorted(list(bio_samples))
        if not bio_samples:
            bio_samples = ["unknown"]
        print("Well sample(s)        : {s}".format(s=", ".join(well_samples)))
        print("Biological sample(s)  : {s}".format(s=", ".join(bio_samples)))


def summarizeXml(args):
    dset = openDataSet(args.infile, strict=args.strict)

    # check to see if there was an error updating the dataset length:
    numFlag = ""
    if dset.numRecords == 0:
        dset.updateCounts()
        if not dset._countsUpdated:
            numFlag = " Unable to update counts!"
    print("DataSet Type          : {f}".format(f=dset.datasetType))
    print("Name                  : {f}".format(f=dset.name))
    print("Id                    : {f}".format(f=dset.uuid))
    print("Number of records     : {r}{f}".format(r=dset.numRecords,
                                                  f=numFlag))
    print("Total number of bases : {r}{f}".format(r=dset.totalLength,
                                                  f=numFlag))
    print("# of Resources        : {r}".format(r=len(dset.toExternalFiles())))
    print("Filters               : {r}".format(r=str(dset.filters) if
                                               dset.filters else "None"))
    show_sample_names_if_defined(dset)
    if args.show_chemistry:
        print("Sequencing Chemistry  : {c}".format(
            c=", ".join(dset.sequencingChemistry)))
    for fname in dset.toExternalFiles():
        print(fname)
    return 0


def summarize_options(parser):
    parser.description = "Print basic information about a DataSet XML file"
    parser.add_argument("infile", type=str,
                        help="The xml file to summarize")
    parser.add_argument("--show-chemistry", action="store_true",
                        help="Show the sequencing chemistries deduced from BAM headers (e.g. 'P6-C4')")
    parser.set_defaults(func=summarizeXml)


def createXml(args):
    if os.path.exists(args.outfile) and not args.force:
        raise IOError("Output file {} already exists. Use --force to "
                      "clobber".format(args.outfile))
    if args.dsType is None:
        dset = openDataFile(*args.infile, strict=args.strict,
                            skipCounts=args.skipCounts,
                            generateIndices=args.generateIndices)
    else:
        dsTypes = DataSet.castableTypes()
        dset = dsTypes[args.dsType](*args.infile, strict=args.strict,
                                    skipCounts=args.skipCounts,
                                    generateIndices=args.generateIndices)
    if args.dsName != '':
        dset.name = args.dsName
    if args.metadata:
        dset.loadMetadata(args.metadata)
    if args.well_sample_name or args.bio_sample_name:
        if args.metadata:
            log.warn(
                "Setting the WellSample or BioSample name will overwrite fields pulled from %s", args.metadata)
        n_new_collections = add_mock_collection_metadata(dset)
        if n_new_collections > 0:
            log.warn(
                "Created new CollectionMetadata from blank template for %d movies", n_new_collections)
        if args.well_sample_name:
            force_set_all_well_sample_names(dset, args.well_sample_name)
        if args.bio_sample_name:
            force_set_all_bio_sample_names(dset, args.bio_sample_name)
    log.debug("Dataset created")
    dset.newUuid()
    dset.write(args.outfile, validate=args.novalidate, relPaths=args.relative)
    log.debug("Dataset written")
    return 0


def create_options(parser):
    parser.description = ('Create an XML file from fofn, BAM, or DataSet XML '
                          'files.')
    parser.add_argument("outfile", type=str, help="The XML file to create")
    parser.add_argument("infile", type=str, nargs='+',
                        help=("The fofn, BAM, or XML file(s) to make into "
                              "a DataSet XML"))
    parser.add_argument("--force", default=False, action='store_true',
                        help=("Clobber output file if it already exists"))
    parser.add_argument("--type", type=str,
                        dest='dsType',
                        choices=DataSet.castableTypes(),
                        help=("The type of XML to create (may be inferred "
                              "if not provided). "))
    parser.add_argument("--name", type=str, default='',
                        dest='dsName',
                        help="The name (in metadata) of the new DataSet")
    parser.add_argument("--generateIndices", action='store_true',
                        default=False,
                        help=("Generate index files (.pbi and .bai for BAM, "
                              ".fai for FASTA).  Requires pysam and "
                              "pbindex."))
    parser.add_argument("--metadata", type=str, default=None,
                        help=("A Sequel metadata.xml file (or DataSet) to "
                              "supply metadata"))
    parser.add_argument("--novalidate", action='store_false', default=True,
                        help=("Don't validate the resulting XML, don't modify "
                              "paths"))
    parser.add_argument("--relative", action='store_true', default=False,
                        help=("Make the included paths relative instead of "
                              "absolute (not compatible with --novalidate)"))
    parser.add_argument("--well-sample-name", action="store", default=None,
                        help=("Set the WellSample name for all movies (will generate new CollectionMetadata from blank template for any movies that are not already represented)."))
    parser.add_argument("--bio-sample-name", action="store", default=None,
                        help=("Set the BioSample name for all movies (will generate new CollectionMetadata from blank template for any movies that are not already represented)."))
    parser.set_defaults(func=createXml)

def pad_separators(base_set):
    """e.g. 'gt' will hit 'length', so we pad it to ' gt '"""
    padded=[]
    for sep in base_set:
        new_sep=[]
        if sep[0] in string.ascii_lowercase:
            new_sep.append(" ")
        new_sep.append(sep)
        if sep[-1] in string.ascii_lowercase:
            new_sep.append(" ")
        padded.append(''.join(new_sep))
    return padded

def parse_filter_list(filtStrs):
    filters=defaultdict(list)
    # pull them from the filter parser:
    separators=OPMAP.keys()
    # pad the ones that start and end with letters
    separators=pad_separators(separators)
    for filt in filtStrs:
        for sep in separators:
            if sep in filt:
                param, condition=filt.split(sep)
                condition=(sep.strip(), condition.strip())
                filters[param.strip()].append(condition)
                break
    return filters

def filterXml(args):
    if args.infile.endswith('xml'):
        dataSet=openDataSet(args.infile, strict=args.strict)
        filters=parse_filter_list(args.filters)
        dataSet.filters.addRequirement(**filters)
        dataSet.updateCounts()
        log.info("{i} filters added".format(i=len(filters)))
        dataSet.write(args.outfile)
    else:
        raise IOError("No files found/found to be compatible")
    return 0

def filter_options(parser):
    pbiFilterOptions=set(Filters()._pbiMappedVecAccMap().keys())
    bamFilterOptions=set(Filters()._bamAccMap.keys())
    parser.description=('Add filters to an XML file. Suggested fields: '
                          '{f}. More expensive fields: {b}.\nMultiple filters '
                          'of different names will be ANDed together, '
                          'multiple filters of the same name will be ORed '
                          'together, duplicating existing requirements'.format(
        f=sorted(list(pbiFilterOptions)),
        b=sorted(list(bamFilterOptions - pbiFilterOptions))))
    # parser.add_argument("infile", type=validate_file,
    parser.add_argument("infile", type=str,
                        help="The XML file to filter")
    parser.add_argument("outfile", type=str,
                        help="The resulting DataSet XML file")
    parser.add_argument("filters", type=str, nargs='+',
                        help=("The parameters, operators and values to filter "
                              "(e.g. 'rq>0.85')"))
    parser.set_defaults(func=filterXml)

def splitXml(args):
    log.debug("Starting split")
    dataSet=openDataSet(args.infile, strict=args.strict)
    chunks=len(args.outfiles)
    if args.chunks:
        chunks=args.chunks
    if isinstance(dataSet, ContigSet):
        dss=dataSet.split(chunks)
    else:
        dss=dataSet.split(chunks=chunks,
                            ignoreSubDatasets=(not args.subdatasets),
                            contigs=args.contigs,
                            maxChunks=args.maxChunks,
                            breakContigs=args.breakContigs,
                            targetSize=args.targetSize,
                            zmws=args.zmws,
                            barcodes=args.barcodes,
                            byRecords=(not args.byRefLength),
                            updateCounts=(not args.noCounts))
    log.debug("Splitting into {i} chunks".format(i=len(dss)))
    infix='chunk{i}'
    chNums=range(len(dss))
    if args.barcodes:
        infix='{i}'
        chNums=['_'.join(ds.barcodes).replace(
            '[', '').replace(']', '').replace(', ', '-') for ds in dss]
    nSuf=-2 if re.search(r".+\.\w+set\.xml", args.infile) else -1
    if not args.outfiles:
        if not args.outdir:
            args.outfiles=['.'.join(args.infile.split('.')[:nSuf] +
                                      [infix.format(i=chNum)] +
                                      args.infile.split('.')[nSuf:])
                             for chNum in chNums]
        else:
            args.outfiles=['.'.join(args.infile.split('.')[:nSuf] +
                                      [infix.format(i=chNum)] +
                                      args.infile.split('.')[nSuf:])
                             for chNum in chNums]
            args.outfiles=[os.path.join(args.outdir,
                                          os.path.basename(outfn))
                             for outfn in args.outfiles]
            num=len(dss)
            end=''
            if num > 5:
                num=5
                end='...'
            log.debug("Emitting {f} {e}".format(
                f=', '.join(args.outfiles[:num]),
                e=end))
    log.debug("Finished splitting, now writing")
    for out_fn, dset in zip(args.outfiles, dss):
        dset.write(out_fn)
    log.debug("Done writing files")
    return 0

def split_options(parser):
    parser.description="Split the DataSet"
    parser.add_argument("infile", type=str,
                        help="The DataSet XML file to split")
    parser.add_argument("--contigs", default=False, action='store_true',
                        help="Split on contigs")
    parser.add_argument("--barcodes", default=False, action='store_true',
                        help="Split on barcodes")
    parser.add_argument("--zmws", default=False, action='store_true',
                        help="Split on ZMWs")
    parser.add_argument("--byRefLength", default=True, action='store_true',
                        help="Split contigs by contig length")
    parser.add_argument("--noCounts", default=False, action='store_true',
                        help="Update DataSet counts after split")
    parser.add_argument("--chunks", default=0, type=int,
                        help="Split contigs into <chunks> total windows")
    parser.add_argument("--maxChunks", default=0, type=int,
                        help="Split contigs into at most <chunks> groups")
    parser.add_argument("--targetSize", default=5000, type=int,
                        help="Target number of records per chunk")
    parser.add_argument("--breakContigs", default=False, action='store_true',
                        help="Break contigs to get closer to maxCounts")
    parser.add_argument("--subdatasets", default=False, action='store_true',
                        help="Split on SubDataSets")
    parser.add_argument("--outdir", default=False, type=validate_output_dir,
                        help="Specify an output directory")
    parser.add_argument("outfiles", nargs=argparse.REMAINDER,
                        type=str, help="The resulting XML files (optional)")
    parser.set_defaults(func=splitXml)

def mergeXml(args):
    dss=[openDataSet(infn, strict=args.strict) for infn in args.infiles]
    allds=reduce(lambda ds1, ds2: ds1 + ds2, dss)
    if not allds is None:
        allds.updateCounts()
        allds.write(args.outfile)
    else:
        raise InvalidDataSetIOError("Merge failed, likely due to "
                                    "conflicting Filters")
    return 0

def merge_options(parser):
    parser.description=('Combine DataSet XML files without touching '
                          'resource files (e.g. BAM, Fasta, etc.)')
    parser.add_argument("outfile", type=str,
                        help="The resulting XML file")
    # parser.add_argument("infiles", type=validate_file, nargs='+',
    parser.add_argument("infiles", type=str, nargs='+',
                        help="The XML files to merge")
    parser.set_defaults(func=mergeXml)


def relativizeXml(args):
    dss=openDataSet(args.infile, strict=args.strict)
    dss.write(args.infile, relPaths=True)
    return 0

def relativize_options(parser):
    # doesn't make sense to have an outdir, that would screw with the relative
    # paths...
    parser.description='Make the paths in a DataSet XML file relative'
    parser.add_argument("infile", type=str,
                        help="The XML file to relativize")
    parser.set_defaults(func=relativizeXml)

def absolutizeXml(args):
    dss=openDataSet(args.infile, strict=args.strict)
    outfn=args.infile
    if args.outdir:
        if os.path.isdir(args.outdir):
            outfn=_swapPath(args.outdir, args.infile)
        else:
            outfn=args.outdir
    if args.updateCounts:
        dss.updateCounts()
    dss.write(outfn, relPaths=False)
    return 0

def absolutize_options(parser):
    parser.description=('Make the paths in an DataSet XML file absolute, '
                          'optionally writing the new DataSet to a new '
                          'location at the same time')
    parser.add_argument("infile", type=str,
                        help="The XML file to absolutize")
    parser.add_argument("--outdir", default=None, type=str,
                        help="Specify an optional output directory")
    parser.add_argument("--update", default=False, action="store_true",
                        dest="updateCounts",
                        help="Update dataset metadata")
    parser.set_defaults(func=absolutizeXml)

def copyToXml(args):
    dss=openDataSet(args.infile, strict=args.strict)
    outfn=args.outdir
    if os.path.isdir(args.outdir):
        outfn=_swapPath(args.outdir, args.infile)
    dss.copyTo(os.path.split(outfn)[0])
    dss.write(outfn, relPaths=args.relative)
    return 0

def copyTo_options(parser):
    parser.description=('Copy a DataSet XML and external resources to a new '
                          'location')
    parser.add_argument("infile", type=str,
                        help="The XML file to copy")
    parser.add_argument("outdir", type=str,
                        help="The copy to directory")
    parser.add_argument("--relative", action='store_true', default=False,
                        help=("Make the included paths relative instead of "
                              "absolute"))
    parser.set_defaults(func=copyToXml)

def newUuidXml(args):
    dss=openDataSet(args.infile, strict=args.strict)
    dss.newUuid(random=args.random)
    dss.write(args.infile, validate=False)
    return 0

def newUniqueId_options(parser):
    parser.description="Refresh a DataSet's UniqueId"
    parser.add_argument("infile", type=str,
                        help="The XML file to refresh")
    parser.add_argument("--random", action='store_true', default=False,
                        help=("Generate a random UUID, instead of a hash"))
    parser.set_defaults(func=newUuidXml)

def loadStatsXml(args):
    dset=openDataSet(args.infile, strict=args.strict)
    if len(dset.externalResources) > 1:
        log.info("More than one ExternalResource found, adding the "
                 "sts.xml nested external resource to the first one")
    dset.externalResources[0].sts=args.statsfile
    if args.outfile:
        dset.write(args.outfile, validate=False)
    else:
        dset.write(args.infile, validate=False)
    return 0

def loadStatsXml_options(parser):
    parser.description=("Add an sts.xml file external resource to a "
                          "DataSet XML file")
    parser.add_argument("infile", type=str,
                        help="The DataSet XML file to modify")
    parser.add_argument("statsfile", type=str,
                        help="The sts.xml file to load")
    parser.add_argument("--outfile", type=str, default=None,
                        help="The DataSet XML file to output")
    parser.set_defaults(func=loadStatsXml)

def loadMetadataXml(args):
    dset=openDataSet(args.infile, strict=args.strict)
    dset.loadMetadata(args.metadata)
    if args.outfile:
        dset.write(args.outfile, validate=False)
    else:
        dset.write(args.infile, validate=False)
    return 0

def loadMetadataXml_options(parser):
    parser.description=('Copy the contents of a Sequel metadata.xml file '
                          'into a DataSet XML file')
    parser.add_argument("infile", type=str,
                        help="The DataSet XML file to modify")
    parser.add_argument("metadata", type=str,
                        help=("The Sequel metadata.xml file to load (or "
                              "DataSet XML to borrow from)"))
    parser.add_argument("--outfile", type=str, default=None,
                        help="The DataSet XML file to output")
    parser.set_defaults(func=loadMetadataXml)

def validateXml(args):
    validateFile(args.infile, args.skipFiles)
    print("{f} is valid DataSet XML with valid ResourceId "
          "references".format(f=args.infile))
    return 0

def validate_options(parser):
    parser.description=('Validate DataSet XML and ResourceId files '
                          '(XML validation '
                          'only available when PyXB is installed)')
    parser.add_argument("infile", type=str,
                        help="The DataSet XML file to validate")
    parser.add_argument("--skipFiles",
                        default=False, action='store_true',
                        help="Skip validating ResourceIds")
    parser.set_defaults(func=validateXml)

def consolidateXml(args):
    """Combine BAMs and apply the filters described in the XML file, producing
    one consolidated XML"""
    dset=openDataSet(args.infile)
    dset.consolidate(args.datafile, numFiles=args.numFiles, useTmp=(not
                     args.noTmp))
    dset.write(args.xmlfile)
    return 0

def consolidate_options(parser):
    parser.description=('Combine the resource files (BAM, fasta, etc.) '
                          'and apply the filters described in a DataSet XML '
                          'file')
    # parser.add_argument("infile", type=validate_file,
    parser.add_argument("--numFiles", type=int, default=1,
                        help="The number of data files to produce")
    parser.add_argument("--noTmp", default=False, action='store_true',
                        help="Don't copy to a tmp location to ensure local"
                             " disk use")
    parser.add_argument("infile", type=str,
                        help="The DataSet XML file to consolidate")
    parser.add_argument("datafile", type=str,
                        help=("The resulting data file (name will be modified "
                              "if more than one output data file is requested "
                              "using numFiles)"))
    parser.add_argument("xmlfile", type=str,
                        help="The resulting DataSet XML file")
    parser.set_defaults(func=consolidateXml)
