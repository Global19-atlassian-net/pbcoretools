#!/bin/bash -ex

source /mnt/software/Modules/current/init/bash
module load hdf5-tools/1.8.16
module load python/2.7.9

mkdir -p tmp/bin
python /mnt/software/v/virtualenv/13.0.1/virtualenv.py tmp/venv
source tmp/venv/bin/activate
# HACK to put binaries on path
if [ ! -z "$PB_TOOLS_BIN" ]; then
  echo "Symlinking to executables in smrttools installation..."
  ln -s $PB_TOOLS_BIN/pbindex tmp/bin
  ln -s $PB_TOOLS_BIN/pbmerge tmp/bin/
  ln -s $PB_TOOLS_BIN/bax2bam tmp/bin
  ln -s $PB_TOOLS_BIN/bam2bam tmp/bin
  ln -s $PB_TOOLS_BIN/samtools tmp/bin
  export PATH=$PATH:$PWD/tmp/bin
else
  echo "WARNING: smrttools not available, some tests will be skipped"
fi

(cd repos/PacBioTestData && make python)
(cd repos/pbcommand && make install)
pip install --no-cache-dir h5py
pip install -r requirements-ci.txt
(cd repos/pbcore && make install)
pip install -r requirements-dev.txt

python setup.py install
make test
