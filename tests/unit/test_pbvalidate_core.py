import logging
import tempfile
import os
import sys

from pbcore.io import ReaderBase

from pbcoretools.pbvalidate.core import (run_validators,
                                         ValidatorErrorContext,
                                         ValidatorContextFirstError,
                                         ValidatorContextMaxErrors,
                                         ValidateFile,
                                         ValidateRecord)

log = logging.getLogger(__name__)


class ValidateTxtFile(ValidateFile):

    def validate(self, path):
        return path.endswith(".txt")


class ValidateTxtCatRecord(ValidateRecord):

    def validate(self, record):
        return "cat" in record


class ValidateTxtDogRecord(ValidateRecord):

    def validate(self, record):
        return "dog" in record


class ValidateBad (ValidateRecord):

    def validate(self, record):
        assert False


class TextFileReader(ReaderBase):

    def __iter__(self):
        for i in self.file:
            yield i.rstrip()


def _get_tmp(fn):
    return tempfile.NamedTemporaryFile(suffix="."+fn).name


def _write_example_file(lines, path):
    with open(path, 'w+') as w:
        w.write("\n".join(lines))


def _to_max_errors(max_errors):
    # A little more difficult to initialize
    def _wrapper(errors, metrics):
        return ValidatorContextMaxErrors(errors, metrics, max_errors)
    return _wrapper


def _to_max_records(max_records):
    def _wrapper(errors, metrics):
        return ValidatorContextMaxRecords(errors, metrics, max_records)
    return _wrapper


class TestValidatorsContext:
    FILE_PATH = _get_tmp("tmp_file.txt")
    DATA = ["cat dog bird", "dog"]

    def setup_method(self, method):
        with open(self.FILE_PATH, 'w+') as f:
            f.write("\n".join(self.DATA))

    def teardown_method(self, method):
        os.remove(self.FILE_PATH)

    def test_1(self):
        """This is a valid file format"""
        validators = [ValidateTxtCatRecord(), ValidateTxtDogRecord(),
                      ValidateTxtFile()]
        file_path = _get_tmp("file.txt")
        contents = ["cat dog ", "cat dog bird", "cat dog tree"]
        _write_example_file(contents, file_path)
        errors, metrics = run_validators(ValidatorContextFirstError, file_path,
                                         TextFileReader, validators)
        os.remove(file_path)
        assert len(errors) == 0

    def test_2(self):
        """This file has many errors"""
        validators = [
            ValidateTxtCatRecord(), ValidateTxtDogRecord(), ValidateTxtFile()]
        file_path = _get_tmp("file.doc")
        contents = ["cat dog ", "cat dog bird", "cat dog tree"]
        # records with errors
        contents.extend(["fish"] * 5)
        _write_example_file(contents, file_path)
        errors, metrics = run_validators(_to_max_errors(2), file_path,
                                         TextFileReader, validators)
        os.remove(file_path)
        assert len(errors) == 2
        errors, metrics = run_validators(_to_max_records(1), file_path,
                                         TextFileReader, validators)
        assert len(errors) == 1

    def test_3(self):
        """Test for consistent behavior when a validator is broken"""
        validators = [ValidateTxtCatRecord(), ValidateTxtDogRecord(),
                      ValidateTxtFile(), ValidateBad()]
        file_path = _get_tmp("file.txt")
        contents = ["cat dog ", "cat dog bird", "cat dog tree"]
        contents.extend(["fish"] * 1)
        _write_example_file(contents, file_path)
        errors, metrics = run_validators(ValidatorErrorContext, file_path,
                                         TextFileReader, validators)
        assert len(errors) == 6
        os.remove(file_path)
        # TODO figure out how to check for log error output showing traceback
        # of exception raised by failed validator
