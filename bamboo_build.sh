#!/bin/bash
set -euox pipefail
nproc

rm -rf build
mkdir -p build/{bin,lib,include,share}
export PYTHONUSERBASE=$PWD/build
export PATH=$PYTHONUSERBASE/bin:$PATH

if [[ -z ${bamboo_repository_branch_name+x} ]]; then
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/private/otherbins/internalall/bin
elif [[ ${bamboo_repository_branch_name} == develop ]]; then
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/private/otherbins/internalall/bin
elif [[ ${bamboo_repository_branch_name} == master ]]; then
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_master_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_master_smrttools-incremental_installdir/private/otherbins/internalall/bin
else
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/private/otherbins/internalall/bin
fi

# HACK to put binaries on path
if [ ! -z "$PB_TOOLS_BIN" ]; then
  if [ ! -d "$PB_TOOLS_BIN" ]; then
    echo "ERROR: $PB_TOOLS_BIN is not a valid directory"
    exit 1
  fi
  echo "Symlinking to executables in smrttools installation..."
  ln -sfn $PB_TOOLS_BIN/pbindex   build/bin/
  ln -sfn $PB_TOOLS_BIN/pbmerge   build/bin/
  ln -sfn $PB_TOOLS_BIN/bax2bam   build/bin/
  ln -sfn $PB_TOOLS_BIN/bam2fasta build/bin/
  ln -sfn $PB_TOOLS_BIN/bam2fastq build/bin/
  ln -sfn $PB_TOOLS_BIN/samtools  build/bin/
  ln -sfn $PB_TOOLS_INTERNAL_BIN/pbmerge  build/bin/
else
  echo "WARNING: smrttools not available, some tests will be skipped"
fi

set +ve
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load python/2
set -ve

if [[ -z ${bamboo_repository_branch_name+x} ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
elif [[ ${bamboo_repository_branch_name} == develop ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
elif [[ ${bamboo_repository_branch_name} == master ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/master
else
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
fi

PIP="pip --cache-dir=${bamboo_build_working_directory:-$PWD}/.pip"

#$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE -e repos/PacBioTestData
#$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE -e repos/pbcommand
#$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE -e repos/pbcore
$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE -r requirements-ci.txt
$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE -r requirements-dev.txt
#$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE "pylint<2.0.0"
$PIP install --user --no-compile --no-index --find-link $WHEELHOUSE pbtestdata pbcommand pbcore

$PIP install --user -e ./

make test
make pylint
