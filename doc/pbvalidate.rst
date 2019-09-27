pbvalidate - a tool for validating spec compliance
=====================================================

``pbvalidate`` is a tool used to validate that files produced by PacBio software
are compliant with our own internal specifications.

Supported input formats:

 - BAM

 - FASTA

 - DataSet XML

See `here <http://pacbiofileformats.readthedocs.org/en/3.0/>`_
for further information about each format's requirements.

Command-line usage
------------------

Output of ``pbvalidate --help``::

  usage: pbvalidate [-h] [--version] [--log-file LOG_FILE]
                    [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL} | --debug | --quiet | -v]
                    [-c] [--quick] [--max MAX_ERRORS]
                    [--max-records MAX_RECORDS]
                    [--type {BAM,Fasta,AlignmentSet,ConsensusSet,ConsensusAlignmentSet,SubreadSet,BarcodeSet,ContigSet,ReferenceSet,GmapReferenceSet}]
                    [--index] [--strict] [-x XUNIT_OUT] [--unaligned]
                    [--unmapped] [--aligned] [--mapped]
                    [--contents {SUBREAD,CCS}] [--reference REFERENCE]
                    [--permissive-headers]
                    file

  Utility for validating files produced by PacBio software against our own
  internal specifications.

  positional arguments:
    file                  BAM, FASTA, or DataSet XML file

  optional arguments:
    -h, --help            show this help message and exit
    --version             show program's version number and exit
    --log-file LOG_FILE   Write the log to file. Default(None) will write to
                          stdout. (default: None)
    --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                          Set log level (default: CRITICAL)
    --debug               Alias for setting log level to DEBUG (default: False)
    --quiet               Alias for setting log level to CRITICAL to suppress
                          output. (default: False)
    -v, --verbose         Set the verbosity level. (default: None)
    -c
    --quick               Limits validation to the first 100 records (plus file
                          header); equivalent to --max-records=100 (default:
                          False)
    --max MAX_ERRORS      Exit after MAX_ERRORS have been recorded (DEFAULT:
                          check entire file) (default: None)
    --max-records MAX_RECORDS
                          Exit after MAX_RECORDS have been inspected (DEFAULT:
                          check entire file) (default: None)
    --type {BAM,Fasta,AlignmentSet,ConsensusSet,ConsensusAlignmentSet,SubreadSet,BarcodeSet,ContigSet,ReferenceSet,GmapReferenceSet}
                          Use the specified file type instead of guessing
                          (default: None)
    --index               Require index files (.fai or .pbi) (default: False)
    --strict              Turn on additional validation, primarily for DataSet
                          XML (default: False)
    -x XUNIT_OUT, --xunit-out XUNIT_OUT
                          Xunit test results for Jenkins (default: None)

  bam:
    BAM options

    --unaligned           Specify that the file should contain only unmapped
                          alignments (DEFAULT: no requirement) (default: None)
    --unmapped            Alias for --unaligned (default: None)
    --aligned             Specify that the file should contain only mapped
                          alignments (DEFAULT: no requirement) (default: None)
    --mapped              Alias for --aligned (default: None)
    --contents {SUBREAD,CCS}
                          Enforce read type (default: None)
    --reference REFERENCE
                          Path to optional reference FASTA file, used for
                          additional validation of mapped BAM records (default:
                          None)
    --permissive-headers  Don't check chemistry/basecaller versions (default:
                          False)

  fasta:
    Fasta options

Examples
--------

Validating a BAM file::

  $ pbvalidate in.subreads.bam

Validating a FASTA file::

  $ pbvalidate in.fasta

Validating a DataSet XML file::

  $ pbvalidate in.subreadset.xml

Validating a BAM file and its index (.pbi)::

  $ pbvalidate --index in.subreads.bam

Validating a BAM file, exiting after 10 errors detected::

  $ pbvalidate --max 10 in.subreads.bam

Validating up to 100 records in a BAM file::

  $ pbvalidate --max-records 100 in.subreads.bam

Same as above (equivalent to --max-records=100)::

  $ pbvalidate --quick in.subreads.bam

Validating a BAM file, at a desired log level::

  $ pbvalidate --log-level=INFO in.subreads.bam

Validating a BAM file, writing log messages to a file rather than stdout::

  $ pbvalidate --log-file validation_results.log in.subreads.bam

Validating a BAM file, generating XUnit-formatted results::

  $ pbvalidate --xunit-out validation_results.xml in.subreads.bam
