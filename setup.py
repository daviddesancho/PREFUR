#!/usr/bin/env python

# Setup script for bestmsm package

import os
from setuptools import setup,find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
		name='PREFUR',
		version='0.2dev',
		description='',
		url='http://github.com/daviddesancho/PREFUR',
		author='David De Sancho',
		author_email='daviddesancho.at.gmail.com',
		license='GPLv3.0',
		packages=find_packages(),
		keywords= "protein folding kinetics",
		long_description=read('README.md'),
		classifiers = ["""\
				Development Status :: 1 - Planning
				Operating System :: POSIX :: Linux
				Operating System :: MacOS
				Programming Language :: Python :: 2.7
				Topic :: Scientific/Engineering :: Bio-Informatics
				Topic :: Scientific/Engineering :: Chemistry
				"""]
		)
