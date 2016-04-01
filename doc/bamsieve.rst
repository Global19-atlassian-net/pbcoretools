bamSieve - a tool for dataset reduction
=======================================

``bamSieve`` is a utility for the generation of sub-datasets by filtering reads
on a per-ZMW basis (keeping all subreads within a read together), inspired on
the ``baxSieve`` program for RSII datasets.
Although it is BAM-centric it has some support for
dataset XML and will propagate metadata, as well as scraps BAM files in
the special case of SubreadSets.  We use this extensively
internally for testing and debugging purposes, e.g. generating minimal test
datasets containing a handful of reads.

The program operates in two modes: **whitelist/blacklist**
where the ZMWs to keep or discard are explicitly specified, or
**percentage/count** mode, where a fraction of ZMWs is randomly selected.
ZMWs may be whitelisted or blacklisted in one of several ways:

  - as a comma-separated list on the command line
  - as a flat text file, one ZMW per line
  - as another PacBio BAM or dataset of any type


Command-line usage
------------------

Output of ``bamSieve --help``::

  usage: bamSieve [-h] [--version] [--log-file LOG_FILE]
                  [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL} | --debug | --quiet | -v]
                  [--show-zmws] [--whitelist WHITELIST] [--blacklist BLACKLIST]
                  [--percentage PERCENTAGE] [-n COUNT] [-s SEED]
                  [--ignore-metadata] [--anonymize] [--barcodes]
                  input_bam [output_bam]
  
  Tool for subsetting a BAM or PacBio DataSet file based on either a whitelist
  of hole numbers or a percentage of reads to be randomly selected.
  
  positional arguments:
    input_bam             Input BAM or DataSet from which reads will be read
    output_bam            Output BAM or DataSet to which filtered reads will be
                          written (default: None)
  
  optional arguments:
    -h, --help            show this help message and exit
    --version             show program's version number and exit
    --log-file LOG_FILE   Write the log to file. Default(None) will write to
                          stdout. (default: None)
    --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                          Set log level (default: WARN)
    --debug               Alias for setting log level to DEBUG (default: False)
    --quiet               Alias for setting log level to CRITICAL to suppress
                          output. (default: False)
    -v, --verbose         Set the verbosity level. (default: None)
    --show-zmws           Print a list of ZMWs and exit (default: False)
    --whitelist WHITELIST
                          Comma-separated list of ZMWs, or file containing
                          whitelist of one hole number per line, or BAM/DataSet
                          file from which to extract ZMWs (default: None)
    --blacklist BLACKLIST
                          Opposite of --whitelist, specifies ZMWs to discard
                          (default: None)
    --percentage PERCENTAGE
                          If you prefer to recover a percentage of a SMRTcell
                          rather than a specific list of reads specify that
                          percentage (range 0-100) here (default: None)
    -n COUNT, --count COUNT
                          Recover a specific number of ZMWs picked at random
                          (default: None)
    -s SEED, --seed SEED  Random seed for selecting a percentage of reads
                          (default: None)
    --ignore-metadata     Discard input DataSet metadata (default: False)
    --anonymize           Randomize sequences for privacy (default: False)
    --barcodes            Indicates that the whitelist or blacklist contains
                          barcode indices instead of ZMW numbers (default:
                          False)


Examples
--------

Pulling out two ZMWs from a BAM file::

  $ bamSieve --whitelist 111111,222222 full.subreads.bam sample.subreads.bam

The dataset equivalent::

  $ bamSieve --whitelist 111111,222222 full.subreadset.xml sample.subreadset.xml

Using a text whitelist::

  $ bamSieve --whitelist zmws.txt full.subreads.bam sample.subreads.bam

Using another BAM or dataset as a whitelist::

  $ bamSieve --whitelist mapped.alignmentset.xml full.subreads.bam mappable.subreads.bam

Generate a whitelist from a dataset::

  $ bamSieve --show-zmws mapped.alignmentset.xml > mapped_zmws.txt

Anonymizing a dataset::

  $ bamSieve --whitelist zmws.txt --ignore-metadata --anonymize full.subreadset.xml anonymous_sample.subreadset.xml

Removing a read::

  $ bamSieve --blacklist 111111 full.subreadset.xml filtered.subreadset.xml

Selecting 0.1% of reads::

  $ bamSieve --percentage 0.1 full.subreads.bam random_sample.subreads.bam

Selecting a different 0.1% of reads::

  $ bamSieve --percentage 0.1 --seed 98765 full.subreads.bam random_sample.subreads.bam

Selecting just two ZMWs/reads at random::

  $ bamSieve --count 2 full.subreads.bam two_reads.subreads.bam

Selecting by barcode::

  $ bamSieve --barcodes --whitelist 4,7 full.subreads.bam two_barcodes.subreads.bam

Generating a tiny BAM file that contains only mappable reads::

  $ bamSieve --whitelist mapped.subreads.bam full.subreads.bam mappable.subreads.bam
  $ bamSieve --count 4 mappable.subreads.bam tiny.subreads.bam

Splitting a dataset into two halves::

  $ bamSieve --percentage 50 full.subreadset.xml split.1of2.subreadset.xml
  $ bamSieve --blacklist split.1of2.subreadset.xml full.subreadset.xml split.2of2.subreadset.xml
