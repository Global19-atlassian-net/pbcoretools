#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import sys

test_deps = [
    'coverage',
    'jinja2',
    'jsonschema',
    'nose',
    'pbtestdata',
    'pytest',
    'sphinx',
    'xmlbuilder',
]

setup(
    name='pbcoretools',
    version='0.7.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    description='Python CLI tools and add-ons for reading and writing PacBio® data files',
    license='BSD-3-Clause-Clear',
    packages=find_packages('.'),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'dataset = pbcoretools.dataset:main',
            'pbvalidate = pbcoretools.pbvalidate.main:main',
            'bamsieve = pbcoretools.bamsieve:main',
            'pbtools-gather = pbcoretools.tasks.gather:main',
        ]},
    install_requires=[
        'numpy >= 1.7.1',
        'pysam >= 0.15.1',
        'pbcore >= 2.0.0',
        'pbcommand >= 2.0.0',
    ],
    test_requires=test_deps,
    extras_require={
        'test': test_deps},
    python_requires='>=3.7',
)
