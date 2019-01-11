#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from setuptools import setup, Extension, find_packages
import os.path
import sys

REQUIREMENTS_TXT = "requirements.txt"

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    print("pbcoretools requires Python 2.7")
    sys.exit(-1)

globals = {}
execfile("pbcoretools/__init__.py", globals)
__VERSION__ = globals["__VERSION__"]

def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        reqs = [line for line in f if not line.startswith("#")]
    return reqs


def _get_local_requirements(file_name):
    return _get_requirements(_get_local_file(file_name))


setup(
    name = 'pbcoretools',
    version=__VERSION__,
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    description="Python CLI tools and add-ons for reading and writing PacBio® data files",
    license=open('LICENSES.txt').read(),
    packages = find_packages('.'),
    package_dir = {'':'.'},
    zip_safe = False,
    entry_points = {"console_scripts": [
        "dataset = pbcoretools.dataset:main",
        'pbvalidate = pbcoretools.pbvalidate.main:main',
        'bamsieve = pbcoretools.bamsieve:main',
        'pbtools-gather = pbcoretools.tasks2.gather:main'
    ]},
    install_requires=_get_local_requirements(REQUIREMENTS_TXT),
    test_requires=("pbtestdata",))
