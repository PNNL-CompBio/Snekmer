from setuptools import setup, find_packages
import glob
import os

version = {}
with open("snekmer/version.py") as fp:
    exec(fp.read(), version)

with open("README.md") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

# with open('requirements.txt') as f:
#     required = [x for x in f.read().splitlines() if not x.startswith("#")]

pkgs = find_packages(exclude=("test"))  # also Util?

# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.

setup(
    name="snekmer",
    version=__version__,
    # packages=['kmerfeatures'],
    # test_suite='pytest.collector',
    # tests_require=['pytest'],
    description=(
        "Kmer pipeline for generating features and"
        "running SIEVE models on input sequences"
    ),
    long_description=readme,
    url="http://github.com/PNNL-CompBio/Snekmer/",
    author="@christinehc, @biodataganache",
    author_email="christine.chang@pnnl.gov",
    license=license,
    packages=pkgs,
    entry_points={"console_scripts": ["snekmer = snekmer.cli:main"]},
    package_data={"": ["rules/*.smk", "scripts/*.py", "templates/*.html"]},
    # install_requires=required,
    include_package_data=True,
    keywords=[],
    zip_safe=False,
)
