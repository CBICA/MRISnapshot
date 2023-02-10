#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
import setuptools

#this_directory = Path(__file__).parent
#long_description = (this_directory / 'README.md').read_text()
long_description = 'TODO'

setup(name='mri_snapshot',
    version='0.0.1-dev',
    description='QC tool for verification of datasets with MRI images and derived maps',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/gurayerus/MRISnapshot',
    author='Guray Erus',
    author_email='guray.erus@pennmedicine.upenn.edu',
    license='MIT',
    zip_safe=False,
    install_requires=[
        'numpy',
        'pandas',
        'nibabel',
        'pillow',
        'scipy',
        'matplotlib'
    ],
    scripts=[
        'create_report', 'prep_dataset'
    ],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Operating System :: Unix',
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
)
