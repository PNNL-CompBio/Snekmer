"""process.smk: Module for input file processing.

author: @christinehc

"""

# imports
import gzip
from datetime import datetime
from os import makedirs, remove
from os.path import dirname, exists, join
from shutil import copy

import snekmer as skm
from pandas import DataFrame


# define rules
rule unzip:
    wildcard_constraints:
        uz=r".*\.(?:fa|fna|faa|fasta)$",
    input:
        join("input", "{uz}.gz"),
    output:
        unzipped=join("input", "{uz}"),
        zipped=join("input", "zipped", "{uz}.gz"),
    run:
        os.makedirs(dirname(output.zipped), exist_ok=True)
        copy(input[0], output.zipped)
        # unzip and save file contents
        with gzip.open(input[0], "rb") as openf, open(output.unzipped, "wb") as savef:
            file_content = openf.readlines()
            for line in file_content:
                savef.write(line)
        remove(input[0])
