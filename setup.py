#!/usr/bin/env python

from distutils.core import setup

setup(name='Distutils',
      version='1.0',
      description='Resolving complex tandem repeats with long reads',
      author='Ajay Ummat, Ali Bashir',
      author_email='ali.bashir@mssm.edu',
      url='http://research.mssm.edu/bashir/',
      packages=['gmmCluster', 'pacIO', 'pacModel', 'pHmm', 'psCount'],
      scripts=['pacMonStr_V1.py'],
      
     )
