
Test of command-line pbvalidate tool

  $ DATA=$TESTDIR/../data

  $ pbvalidate
  usage: pbvalidate [-h] [-v] [-c] [--verbose] [-q] [--debug] [--quick]
                    [--max MAX_ERRORS] [--max-records MAX_RECORDS]
                    [--type {BAM,Fasta,AlignmentSet,ConsensusSet,SubreadSet,BarcodeSet,ReferenceSet,HdfSubreadSet}]
                    [--index] [--strict] [-x XUNIT_OUT] [--unaligned]
                    [--unmapped] [--aligned] [--mapped]
                    [--contents {SUBREAD,CCS}] [--reference REFERENCE]
                    [--permissive-headers]
                    file
  pbvalidate: error: too few arguments
  [2]
  $ pbvalidate --help | grep -c PacBio
  1

Good Fasta file

  $ pbvalidate $DATA/test0.fasta

Good BAM file

  $ pbvalidate $DATA/tst_1_subreads.bam

Awful BAM file

  $ pbvalidate --index $DATA/tst_2_subreads.bam
  22 spec violation(s) detected:
    This file has not been sorted by position, or the header has not been updated.
    Missing corresponding .pbi index file
    Missing platform (PL) for read group 2f48aec3
    Mismatch between specified and expected read group ID: 2f48aec3 in file, but computed as 3f58e5b8
    The chemistry information for read group 2f48aec3 is either missing or cannot be interpreted: ('null', 'null', '9.8')
    The basecaller version number '9.8' in read group '2f48aec3' is not one of the allowed values; it is probably being misreported by an upstream program such as baz2bam or bax2bam.
    The read group 2f48aec3 declares the tag 'Ipd' for pulse features, but the encoding scheme is not specified or not recognized
    Range specified in QNAME movie1/54130/0_10 conflicts with QS and QE tags
      [1 similar error(s) not displayed]
    The length of the sequence of movie1/54130/0_10 and the length indicated by qStart/qEnd disagree (10 versus 0)
      [2 similar error(s) not displayed]
    Value '2001.0' of tag 'rq' for qname movie1/54130/0_10 is not one of the expected values or within the required range.
      [1 similar error(s) not displayed]
    The HQRegionSNR field ('sn') for the alignment of movie1/54130/0_10 is uninitialized
      [1 similar error(s) not displayed]
    The CIGAR string for the alignment of movie1/54130/0_10 contains the ambiguous 'M' character to specify either a match or mismatch
      [1 similar error(s) not displayed]
    Adjacent insertion and deletion in CIGAR string for QNAME movie1/54130/10_20
    Query name 'movie1_54130_10-20' is not in the expected format!
    The alignment of movie1_54130_10-20 is not mapped to the reference sequence
    QNAME movie1/54130/0_10 occurs more than once (tStart = 1, 1, 1)
  [1]

Now with Xunit output (exit code 0):

  $ pbvalidate --quiet --index --xunit-out tst_2_subreads_xunit.xml $DATA/tst_2_subreads.bam
  $ grep -c 'failures="17"' tst_2_subreads_xunit.xml
  1

Now ignoring some header issues

  $ pbvalidate --quick --index --permissive-headers $DATA/tst_2_subreads.bam
  20 spec violation(s) detected:
    This file has not been sorted by position, or the header has not been updated.
    Missing corresponding .pbi index file
    Missing platform (PL) for read group 2f48aec3
    Mismatch between specified and expected read group ID: 2f48aec3 in file, but computed as 3f58e5b8
    The read group 2f48aec3 declares the tag 'Ipd' for pulse features, but the encoding scheme is not specified or not recognized
    Range specified in QNAME movie1/54130/0_10 conflicts with QS and QE tags
      [1 similar error(s) not displayed]
    The length of the sequence of movie1/54130/0_10 and the length indicated by qStart/qEnd disagree (10 versus 0)
      [2 similar error(s) not displayed]
    Value '2001.0' of tag 'rq' for qname movie1/54130/0_10 is not one of the expected values or within the required range.
      [1 similar error(s) not displayed]
    The HQRegionSNR field ('sn') for the alignment of movie1/54130/0_10 is uninitialized
      [1 similar error(s) not displayed]
    The CIGAR string for the alignment of movie1/54130/0_10 contains the ambiguous 'M' character to specify either a match or mismatch
      [1 similar error(s) not displayed]
    Adjacent insertion and deletion in CIGAR string for QNAME movie1/54130/10_20
    Query name 'movie1_54130_10-20' is not in the expected format!
    The alignment of movie1_54130_10-20 is not mapped to the reference sequence
    QNAME movie1/54130/0_10 occurs more than once (tStart = 1, 1, 1)
  [1]

  $ pbvalidate --quick --aligned --content=CCS $DATA/tst_1_subreads.bam
  1 spec violation(s) detected:
    File was expected to contain only CCS reads, but this file contains standard reads
  [1]

And now a dataset XML file

  $ pbvalidate $DATA/tst_1.alignmentset.xml
  $ pbvalidate --type=ReferenceSet $DATA/tst_1.alignmentset.xml
  1 spec violation(s) detected:
    Unexpected error reading dataset: [Errno 5] Cannot create PacBio.DataSet.ReferenceSet from PacBio.DataSet.AlignmentSet: '*/tst_1.alignmentset.xml'.  This prevents any further validation functions from being run. (glob)
  [1]
  $ pbvalidate $DATA/tst_1.subreadset.xml

Missing XML declaration

  $ pbvalidate $DATA/tst_1b.subreadset.xml
  1 spec violation(s) detected:
    This XML document is either missing the header or the encoding type is missing or wrong; all DataSet XMLs should explicitly specify UTF-8 encoding.
  [1]

  $ pbvalidate $DATA/tst_1c.subreadset.xml
  1 spec violation(s) detected:
    This XML document is either missing the header or the encoding type is missing or wrong; all DataSet XMLs should explicitly specify UTF-8 encoding.
  [1]

Mysterious PyXB error

  $ pbvalidate --strict $DATA/tst_1d.subreadset.xml | grep -c "XML schema error"
  1

Now a dataset with the awful .bam file.  We can't actually do much with this
because there's no .pbi file (pbindex doesn't like this .bam).

  $ pbvalidate $DATA/tst_2_subreads.xml
  1 spec violation(s) detected:
    Unexpected error reading dataset: IndexedBamReader requires bam.pbi index file to read *.  This prevents any further validation functions from being run. (glob)
  [1]

#Now a copy of tst_1.alignmentset.xml, with misleading file name
#FIXME disabled for now pending further discussion
#  $ pbvalidate $DATA/tst_1.ccs.xml
#  1 spec violation(s) detected:
#    The dataset file is named misleadingly - 'tst_1.ccs.xml' suggests ConsensusReadSet, but the actual type is AlignmentSet
#  [1]

Now a ReferenceSet XML file (i.e. Fasta sequences)

  $ pbvalidate $DATA/reference.xml
  1 spec violation(s) detected:
    Sequence 'ecoliK12_pbi_March2013_2955000_to_2980000' does not have line wrapping
  [1]

  $ pbvalidate --strict $DATA/reference.xml
  2 spec violation(s) detected:
    The dataset file reference.xml is named incorrectly - datasets of type 'ReferenceSet' should have the extension '.referenceset.xml'.
    Sequence 'ecoliK12_pbi_March2013_2955000_to_2980000' does not have line wrapping
  [1]

  $ pbvalidate --type=AlignmentSet $DATA/reference.xml
  1 spec violation(s) detected:
    Unexpected error reading dataset: [Errno 5] Cannot create PacBio.DataSet.AlignmentSet from PacBio.DataSet.ReferenceSet: '*/reference.xml'.  This prevents any further validation functions from being run. (glob)
  [1]


bamSieve:

  $ bamSieve  --show-zmws $DATA/tst_1_subreads.bam
  54130
