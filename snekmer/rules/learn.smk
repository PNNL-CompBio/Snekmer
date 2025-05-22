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
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, split, splitext
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from Bio import SeqIO
from scipy.interpolate import interp1d
from scipy.stats import rankdata
import random
import snekmer as skm
from scipy.ndimage import gaussian_filter1d
import sklearn.metrics.pairwise
import itertools
import os
import shutil
from pkg_resources import resource_filename


# Note:
# Pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
inputDir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
inputFiles = glob(join(inputDir, "*"))
# base_file = glob(join(inputDir,"base" "*"))
zipped = [fa for fa in inputFiles if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(inputFiles, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
annotFiles = glob(join("annotations", "*.ann"))
baseCounts = glob(join("base", "counts", "*.csv"))
baseConfidence = glob(join("base", "confidence", "*.csv"))
baseFamilyCheckpoint = glob(join("base", "thresholds", "*.csv"))

uzMap = {skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped}
faMap = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in uzMap.items()]
# isolate basenames for all files
FAS = list(faMap.keys())
# parse any background files
backgroundFiles = glob(join(inputDir, "background", "*"))
if len(backgroundFiles) > 0:
    backgroundFiles = [
        skm.utils.split_file_ext(basename(f))[0] for f in backgroundFiles
    ]
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
outDir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

output_prefixes = (
    ["vector"] if not config["learnapp"]["fragmentation"] else ["vector", "vector_frag"]
)

if (
    config["learnapp"]["selection"] != "top_hit"
    and config["learnapp"]["threshold"] == "None"
):
    raise Exception(
        "The only selection method that allows for None is `top_hit`. Other methods inherently use a threshold"
    )

# define output files to be created by snekmer
rule all:
    input:
        # Vector files
        expand(join(outDir, "vector", "{nb}.npz"), nb=FAS),
        # Kmer counts files for learning
        expand(join(outDir, "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        join(outDir, "learn", "kmer-counts-total.csv"),
        # Fragmentation outputs (only if enabled)
        expand(join(outDir, "fragmented", "{nb}.fasta"), nb=FAS)
        if config["learnapp"]["fragmentation"]
        else [],
        expand(join(outDir, "vector_frag", "{nb}.npz"), nb=FAS)
        if config["learnapp"]["fragmentation"]
        else [],
        # Evaluation application sequences (or fragments if fragmentation is enabled)
        expand(
            join(
                outDir,
                (
                    "eval_apply_sequences"
                    if not config["learnapp"]["fragmentation"]
                    else "eval_apply_frag"
                ),
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        # Evaluation configuration outputs
        join("output", "eval_conf", "family_summary_stats.csv"),
        join("output", "eval_conf", "global-confidence-scores.csv"),
        # Duplicated some outputs for each of use on the user side.
        join(outDir, "apply_inputs", "counts", "kmer-counts-total.csv"),
        join(outDir, "apply_inputs", "stats", "family_summary_stats.csv"),
        join(outDir, "apply_inputs", "confidence", "global-confidence-scores.csv"),


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join(inputDir, "{uz}"),
        zipped=join(inputDir, "zipped", "{uz}.gz"),


if config["learnapp"]["fragmentation"]:

    rule fragmentation:
        input:
            fasta=lambda wildcards: join(
                inputDir, f"{wildcards.nb}.{faMap[wildcards.nb]}"
            ),
        output:
            fasta_out=join(outDir, "fragmented", "{nb}.fasta"),
        params:
            version=config["learnapp"]["version"],
            frag_length=config["learnapp"]["frag_length"],
            location=config["learnapp"]["location"],
            min_length=config["learnapp"]["min_length"],
            seed=config["learnapp"]["seed"],
        script:
            resource_filename("snekmer", join("scripts/learn_fragment.py"))


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(
            outDir if wildcards.prefix == "vector_frag" else inputDir,
            "fragmented" if wildcards.prefix == "vector_frag" else "",
            f"{wildcards.nb}.{faMap[wildcards.nb]}",
        ),
    output:
        data=join(outDir, "{prefix}", "{nb}.npz"),
        kmerobj=join(outDir, "kmerize_{prefix}", "{nb}.kmers"),
    log:
        join(outDir, "{prefix}_kmerize", "log", "{nb}.log"),


# WORKFLOW to learn kmer associations
rule learn:
    input:
        data=join(outDir, "vector", "{nb}.npz"),
        annotation=expand("{an}", an=annotFiles),
    output:
        counts=join(outDir, "learn", "kmer-counts-{nb}.csv"),
    log:
        join(outDir, "learn", "log", "learn-{nb}.log"),
    script:
        resource_filename("snekmer", join("scripts/learn_learn.py"))

rule merge:
    input:
        counts=expand(join(outDir, "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        baseCounts=expand("{bf}", bf=baseCounts),
    output:
        totals=join(outDir, "learn", "kmer-counts-total.csv"),
    log:
        join(outDir, "learn", "log", "merge.log"),
    script:
        resource_filename("snekmer", join("scripts/learn_merge.py"))


rule eval_apply_reverse_seqs:
    input:
        data=join(
            outDir,
            (
                "vector"
                if config["learnapp"]["fragmentation"] == False
                else "vector_frag"
            ),
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annotFiles),
        compareAssociations=join(outDir, "learn", "kmer-counts-total.csv"),
    output:
        apply=join(
            outDir,
            (
                "eval_apply_reversed"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag"
            ),
            "seq-annotation-scores-{nb}.csv.gz",
        ),
    log:
        join(outDir, "eval_apply_reversed", "log", "{nb}.log"),
    script:
        resource_filename("snekmer", join("scripts/learn_eval_apply_reverse_seqs.py"))


rule reverseDecoy_evaluations:
    input:
        evalApplyData=expand(
            join(
                outDir,
                (
                    "eval_apply_reversed"
                    if config["learnapp"]["fragmentation"] == False
                    else "eval_apply_frag"
                ),
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        baseFamilyCheckpoint=expand("{bc}", bc=baseFamilyCheckpoint),
    output:
        familyStats=join(outDir, "eval_conf", "family_summary_stats.csv"),
        checkpoint=join(outDir, "eval_conf", "family_stats_checkpoint.csv"),
    script:
        resource_filename("snekmer", join("scripts/learn_reverse_decoy_evaluations.py"))
    

rule eval_apply_sequences:
    input:
        data=join(
            outDir,
            (
                "vector"
                if config["learnapp"]["fragmentation"] == False
                else "vector_frag"
            ),
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annotFiles),
        compareAssociations=join(outDir, "learn", "kmer-counts-total.csv"),
    output:
        apply=join(
            outDir,
            (
                "eval_apply_sequences"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag"
            ),
            "seq-annotation-scores-{nb}.csv.gz",
        ),
    log:
        join(outDir, "eval_apply_sequences", "log", "{nb}.log"),
    script:
        resource_filename("snekmer", join("scripts/learn_eval_apply_sequences.py"))


rule evaluate:
    input:
        evalApplyData=expand(
            join(
                outDir,
                (
                    "eval_apply_sequences"
                    if config["learnapp"]["fragmentation"] == False
                    else "eval_apply_frag"
                ),
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        baseConfidence=expand("{bc}", bc=baseConfidence),
        reverseDecoyStats=join(outDir, "eval_conf", "family_summary_stats.csv"),
    output:
        evalGlob=join(outDir, "eval_conf", "global-confidence-scores.csv"),
    params:
        modifier=config["learnapp"]["conf_weight_modifier"],
    log:
        join(outDir, "eval_conf", "log", "conf.log"),
    script:
        resource_filename("snekmer", join("scripts/learn_evaluate_sequences.py"))


rule copy_results_for_apply:
    input:
        kmer_counts_total=join(outDir, "learn", "kmer-counts-total.csv"),
        family_stats=join(outDir, "eval_conf", "family_summary_stats.csv"),
        global_conf_scores=join(outDir, "eval_conf", "global-confidence-scores.csv"),
    output:
        kmer_counts_total=join(
            outDir, "apply_inputs", "counts", "kmer-counts-total.csv"
        ),
        family_stats=join(outDir, "apply_inputs", "stats", "family_summary_stats.csv"),
        global_conf_scores=join(
            outDir, "apply_inputs", "confidence", "global-confidence-scores.csv"
        ),
    run:
        target_dir = os.path.join(outDir, "apply_inputs")
        os.makedirs(target_dir, exist_ok=True)
        shutil.copy(input.kmer_counts_total, output.kmer_counts_total)
        shutil.copy(input.family_stats, output.family_stats)
        shutil.copy(input.global_conf_scores, output.global_conf_scores)
