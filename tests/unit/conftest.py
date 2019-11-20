import os
import pytest

from pbcore.io.dataset.utils import which


def _check_constools():
    return which('dataset') and which('samtools') and which('pbmerge') and which('pbindex')


def _check_bam2fastx():
    return which('bam2fasta') and which('bam2fastq')


def _check_xmllint():
    return which('xmllint')


def _internal_data():
    return os.path.exists("/pbi/dept/secondary/siv/testdata")


def pytest_runtest_setup(item):
    for mark in item.iter_markers():
        if mark.name == 'internal_data':
            if not _internal_data():
                pytest.skip(
                    "need access to '/pbi/dept/secondary/siv/testdata'")
        elif mark.name == 'constools':
            if not _check_constools():
                pytest.skip("need 'pbindex'/'samtools'/'pbmerge'")
        elif mark.name == 'bam2fastx':
            if not _check_bam2fastx():
                pytest.skip("need 'bam2fastx'")
        elif mark.name == 'xmllint':
            if not _check_xmllint():
                pytest.skip("need 'xmllint'")
        else:
            raise LookupError("Unknown pytest mark: '{}'".format(mark.name))
