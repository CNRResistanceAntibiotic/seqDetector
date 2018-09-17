#!/usr/bin/python3

# setup.py file install libraries wanted for seqassembler
# launch it with pip install -e .

# Install setuptools if not already present.
from setuptools import setup, find_packages
import glob

setup(
    name='seqdetector',
    version='1.0.1',
    description='seqdetector: pipeline CNR Resistance for detection',
    packages=find_packages(),
    author='Richard Bonnet',
    author_email='rbonnet@chu-clermontferrand.fr',
    url='https://github.com/CNRResistanceAntibiotic/seqDetector',
    scripts=glob.glob('scripts/*'),
    install_requires=['', 'numpy', 'pandas'],
    license='GPLv3',
    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],

)
