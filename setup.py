from setuptools import setup, find_packages
import glob
import os

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from smk import __version__, _program

setup(name='kmer',
      version=__version__,
      packages=['smk'],
      test_suite='pytest.collector',
      tests_require=['pytest'],
      description=('Kmer pipeline for generating features and'
                   'running SIEVE models on input sequences)',
      url='http://github.com/biodataganache/KmerPipeline/',
      author='@christinehc, @biodataganache',
      author_email='christine.chang@pnnl.gov',
      license='MIT',
      entry_points="""
      [console_scripts]
      {program} = snakemake.cli:main
      """.format(program = _program),
      install_requires=required,
      include_package_data=True,
      keywords=[],
      zip_safe=False)
