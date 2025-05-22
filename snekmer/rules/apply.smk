# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version

min_version("6.0")


# load snakemake modules
module process:
    snakefile:
        "process.smk"
    config:
        config


module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


import copy
import csv

# built-in imports
import gzip
import itertools
import json
import pickle
import struct
import sys
import time
from datetime import datetime
from glob import glob
from itertools import islice, product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, split, splitext
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from Bio import SeqIO
from scipy import spatial
from scipy.stats import rankdata

import snekmer as skm
from pkg_resources import resource_filename



# Note:
# Pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = glob(join(input_dir, "*"))
# base_file = glob(join(input_dir,"base" "*"))xf
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
annot_files = glob(join("annotations", "*.ann"))
compare_file = glob(join("counts", "*.csv"))
confidence_file = glob(join("confidence", "*.csv"))
decoy_stats_file = glob(join("stats", "*.csv"))
# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
# seq-annotation-scores
# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]
# isolate basenames for all files
FAS = list(FA_MAP.keys())
FAV = list(FA_MAP.values())
# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

threshold_type = config["learnapp"]["threshold"]
selection_type = config["learnapp"]["selection"]


wildcard_constraints:
    dataset=FAS,
    FAS=FAS,


rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        *[
            (
                expand(
                    join(out_dir, "apply", "seq-annotation-scores-{nb}.csv"), nb=FAS
                )
                if config["learnapp"]["save_apply_associations"]
                else []
            )
        ],
        expand(join(out_dir, "apply", "kmer-summary-{nb}.csv"), nb=FAS),


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(
            input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"
        ),
    output:
        data=join(out_dir, "vector", "{nb}.npz"),
        kmerobj=join(out_dir, "kmerize", "{nb}.kmers"),
    log:
        join(out_dir, "kmerize", "log", "{nb}.log"),


rule apply:
    input:
        data=join(out_dir, "vector", "{nb}.npz"),
        annotation=expand("{an}", an=annot_files),
        compare_associations=expand("{comp}", comp=compare_file),
        confidence_associations=expand("{conf}", conf=confidence_file),
        decoy_stats=expand("{decoy}", decoy=decoy_stats_file),
    params:
        selection_type=config["learnapp"]["selection"],
        threshold_type=config["learnapp"]["threshold"],
    output:
        seq_ann=(
            expand(join(out_dir, "apply", "seq-annotation-scores-{nb}.csv"), nb=FAS)
            if config["learnapp"]["save_apply_associations"]
            else []
        ),
        kmer_summary=join(out_dir, "apply", "kmer-summary-{nb}.csv"),
    log:
        join(out_dir, "apply", "log", "{nb}.log"),

    script:
        resource_filename("snekmer", join("scripts/apply.py"))