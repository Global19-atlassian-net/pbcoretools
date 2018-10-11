#!/bin/bash
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load python/2.7.9-mobs-pbcoretools
set -euo pipefail
set -x

if [[ -z ${bamboo_repository_branch_name+x} ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/private/otherbins/internalall/bin
elif [[ ${bamboo_repository_branch_name} == develop ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/private/otherbins/internalall/bin
elif [[ ${bamboo_repository_branch_name} == master ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/master
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_master_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_master_smrttools-incremental_installdir/private/otherbins/internalall/bin
  git -C repos/pbcore checkout master
  git -C repos/pbcommand checkout master
else
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
  PB_TOOLS_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/smrtcmds/bin
  PB_TOOLS_INTERNAL_BIN=/pbi/dept/secondary/builds/links/current_develop_smrttools-incremental_installdir/private/otherbins/internalall/bin
fi

rm -rf build
mkdir -p build/{bin,lib,include,share}

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

PIP="pip --cache-dir=${bamboo_build_working_directory:-$PWD}/.pip"
export PATH=$PWD/build/bin:$PATH
export PYTHONUSERBASE=$PWD/build

$PIP install --no-compile --user --find-link $WHEELHOUSE -e repos/PacBioTestData
$PIP install --no-compile --user --find-link $WHEELHOUSE -e repos/pbcommand
$PIP install --no-compile --user --find-link $WHEELHOUSE -e repos/pbcore
$PIP install --no-compile --user --find-link $WHEELHOUSE -r requirements-ci.txt
$PIP install --no-compile --user --find-link $WHEELHOUSE -r requirements-dev.txt
$PIP install --no-compile --user --find-link $WHEELHOUSE "pylint<2.0.0"

$PIP install --user -e ./

make pylint
make test
