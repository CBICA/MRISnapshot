#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
import setuptools
from glob import glob
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / 'README.rst').read_text()

setup(name='mrisnapshot',
    version='0.0.4-dev',
    description='QC tool for verification of datasets with MRI images and derived maps',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/CBICA/MRISnapshot',
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
    entry_points={
        'console_scripts': [
            'mrisnapshot_prep_data = MRISnapshot.prep_data:main',
            'mrisnapshot_create_report = MRISnapshot.create_report:main',
        ]        
    },    
    data_files=[
        ('', glob('MRISnapshot/js_templates/*.js'))
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
