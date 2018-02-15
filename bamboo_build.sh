#!/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
set -ex

NX3PBASEURL=http://nexus/repository/unsupported/pitchfork/gcc-6.4.0
export PATH=$PWD/bin:/mnt/software/a/anaconda2/4.2.0/bin:$PATH
export PYTHONUSERBASE=$PWD
export CFLAGS="-I/mnt/software/a/anaconda2/4.2.0/include"
PIP="pip --cache-dir=$bamboo_build_working_directory/.pip"
module load gcc

rm -rf bin lib include share
mkdir  bin lib include share

# HACK to put binaries on path
if [ ! -z "$PB_TOOLS_BIN" ]; then
  if [ ! -d "$PB_TOOLS_BIN" ]; then
    echo "ERROR: $PB_TOOLS_BIN is not a valid directory"
    exit 1
  fi
  echo "Symlinking to executables in smrttools installation..."
  ln -sfn $PB_TOOLS_BIN/pbindex  bin/
  ln -sfn $PB_TOOLS_BIN/pbmerge  bin/
  ln -sfn $PB_TOOLS_BIN/bax2bam  bin/
  ln -sfn $PB_TOOLS_BIN/bam2fasta  bin/
  ln -sfn $PB_TOOLS_BIN/bam2fastq  bin/
  ln -sfn $PB_TOOLS_BIN/samtools bin/
else
  echo "WARNING: smrttools not available, some tests will be skipped"
fi

$PIP install --user \
  iso8601
$PIP install --user \
  $NX3PBASEURL/pythonpkgs/xmlbuilder-1.0-cp27-none-any.whl \
  $NX3PBASEURL/pythonpkgs/tabulate-0.7.5-cp27-none-any.whl \
  $NX3PBASEURL/pythonpkgs/pysam-0.13-cp27-cp27mu-linux_x86_64.whl \
  $NX3PBASEURL/pythonpkgs/avro-1.7.7-cp27-none-any.whl
ln -sfn ../data repos/PacBioTestData/pbtestdata/data
$PIP install --user --upgrade pylint
$PIP install --user -e repos/PacBioTestData
$PIP install --user -e repos/pbcommand
$PIP install --user -e repos/pbcore
$PIP install --user -r requirements-ci.txt
$PIP install --user -r requirements-dev.txt

$PIP install --user -e ./

make pylint
make test
