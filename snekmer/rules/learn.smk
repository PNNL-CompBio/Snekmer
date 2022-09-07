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


# built-in imports
import gzip
import json
import pickle
import struct
from datetime import datetime
from glob import glob
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, splitext, split

import numpy as np
import pandas as pd
import seaborn as sns
import snekmer as skm
from Bio import SeqIO

# collect all fasta-like files, unzipped filenames, and basenames
input_dir = "input" if (("input_dir" not in config) or (str(config["input_dir"]) == "None")) else config["input_dir"]
input_files = glob(join(input_dir, "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]

annot_files = glob(join("annotations", "*.ann"))

# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]

# isolate basenames for all files
FAS = list(FA_MAP.keys())

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

# define output files to be created by snekmer
rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),  # require unzipping
        join(out_dir, "learn", "kmer-associations.csv"),  # require learning


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join(input_dir, "{uz}"),
        zipped=join(input_dir, "zipped", "{uz}.gz"),


# build kmer count vectors for each basis set
use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),

# WORKFLOW to learn kmer associations
# collect all seq files and generate mega-cluster
rule learn:
    input:
        kmerobj=expand(join("output", "kmerize", "{fa}.kmers"), fa=NON_BGS),
        data=expand(join("output", "vector", "{fa}.npz"), fa=NON_BGS),
        annotation=expand("{an}", an=annot_files)
    output:
        table=join(out_dir, "learn", "kmer-associations.csv")
    log:
        join(out_dir, "learn", "log", "learn.log"),
    run:
        # log script start time
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        # load in the kmerbasis and vectors
        for f in input.data:
            df, kmerlist = skm.io.load_npz(f)
            seqids = df["sequence_id"]
            vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)

        # vecs is now an numpy array with columns corresponding to kmerlist
        #      each row is a protein with names given by seqids
        # seqids: >tr|A0A0C2DJP9|A0A0C2DJP9_9STAP Lipid kinase OS=Salinicoccus roseus OX=45670 GN=SN16_09640 PE=3 SV=1

        # read in the annotation files which should be in the format
        #    column 1: sequence id
        #    column 2: annotation
        # annotations: A0A0C2DJP9	IPR016064
        annots = list()
        for f in input.annotation:
            annots.append(pd.read_table(f))
        annotations = pd.concat(annots)

        # The plan is to build a matrix of kmer-annotation associations
        #    by calculating the probability that a kmer is associated with
        #    each annotation (label). Options should give the ability to
        #    filter out low information associations from the output
        # And then output an output/learn/kmer-associations.csv file (or some
        #    format file)
        # Code to go here:


        
        skm.utils.log_runtime(log[0], start_time)
