from setuptools import setup, find_packages
from kmerfeatures import __version__
import glob
import os

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

pkgs = find_packages(exclude=('test'))  # also Util?

# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.

setup(name='kmerfeatures',
      version=__version__,
      # packages=['kmerfeatures'],
      # test_suite='pytest.collector',
      # tests_require=['pytest'],
      description=('Kmer pipeline for generating features and'
                   'running SIEVE models on input sequences'),
      url='http://github.com/biodataganache/KmerPipeline/',
      author='@christinehc, @biodataganache',
      author_email='christine.chang@pnnl.gov',
      license='MIT',
      packages=pkgs,
      entry_points={
          'console_scripts': ['kmerfeatures = kmerfeatures.cli:main']
          },
      package_data={'': ['Snakefile']},
      install_requires=required,
      include_package_data=True,
      keywords=[],
      zip_safe=False)
