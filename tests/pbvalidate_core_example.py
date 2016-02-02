
import functools
import logging
import sys

from test_pbvalidate_core import *
from test_pbvalidate_core import _write_example_file
from pbvalidate.core import (ValidatorContextMaxErrors,
                             ValidatorContextFirstError,
                             ValidatorErrorContext)
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)


def example_01():
    """This is a valid file format"""

    validators = [
        ValidateTxtCatRecord(), ValidateTxtDogRecord(), ValidateTxtFile()]

    file_path = "file.txt"
    contents = ["cat dog ", "cat dog bird", "cat dog tree"]
    _write_example_file(contents, file_path)
    errors, metrics = run_validators(
        ValidatorContextFirstError, file_path, TextFileReader, validators)

    return errors, metrics


def example_02(context_klass):
    """This file has many errors"""

    validators = [
        ValidateTxtCatRecord(), ValidateTxtDogRecord(), ValidateTxtFile()]
    file_path = "file.txt"
    contents = ["cat dog ", "cat dog bird", "cat dog tree"]
    # records with errors
    contents.extend(["fish"] * 5)
    _write_example_file(contents, file_path)
    errors, metrics = run_validators(
        context_klass, file_path, TextFileReader, validators)
    return errors, metrics


def _to_max_errors(max_errors):
    # A little more difficult to initialize
    def _wrapper(errors, metrics):
        return ValidatorContextMaxErrors(errors, metrics, max_errors)
    return _wrapper


def main():

    def _log_results(name, e, m):
        log.info("Test {x} with {n} errors {e} and metrics {m}".format(
            x=name, e=errors, m=m, n=len(errors)))

    setup_log(log, level=logging.DEBUG)

    quick_ex_02 = functools.partial(example_02, ValidatorContextFirstError)
    full_ex_02 = functools.partial(example_02, ValidatorErrorContext)
    # the MaxError context is little bit more difficult to initialize
    max_errors_ex_02 = functools.partial(example_02, _to_max_errors(2))

    examples = [example_01, quick_ex_02, full_ex_02, max_errors_ex_02]

    for example in examples:
        log.info("running example {e}".format(e=example))
        errors, metrics = example()
        _log_results(example_01.__name__, errors, metrics)
        log.info("")

    log.info("exiting main.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
