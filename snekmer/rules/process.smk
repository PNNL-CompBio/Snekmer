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
    input:
        join("input", "{uz}.gz"),
    output:
        unzipped=join("input", "{uz}"),
        zipped=join("input", "zipped", "{uz}.gz"),
    run:
        # preserve zipped file
        copy(input[0], output.zipped)

        # unzip and save file contents
        with gzip.open(input[0], "rb") as openf, open(output.unzipped, "wb") as savef:
            file_content = openf.readlines()
            for line in file_content:
                savef.write(line)

        remove(input[0])


        # rule perform_kmer_walk:
        #     input:
        #         fasta=get_fasta_files
        #     output:
        #         # need to fix code to properly generate an output...
        #     run:
        #         skm.walk.kmer_walk(input.fasta)

