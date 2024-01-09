# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version

min_version("6.0")


module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


# built-in imports
from glob import glob
from itertools import product
from os.path import basename, dirname, join
from pkg_resources import resource_filename

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import snekmer as skm

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
gz_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}.gz"))
seq_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}"))

# check input file size
for f, e in zip(seq_input.filename, seq_input.ext):
    skm.utils.check_n_seqs(
        join("input", f"{f}.{e}"), config["model"]["cv"], show_warning=False
    )

# add unzipped gz files to total input list
for f, e in zip(gz_input.filename, gz_input.ext):
    getattr(seq_input, "filename").append(f)
    getattr(seq_input, "ext").append(e)

# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)


ruleorder: unzip > vectorize > cluster > cluster_report


# define output files to be created by snekmer
rule all:
    input:
        expand(
            join(out_dir, "kmerize", "{f}.{e}.kmers"),
            f=seq_input.filename,
            e=seq_input.ext,
        ),
        join(out_dir, "Snekmer_Cluster_Report.html"),


# if any files are gzip zipped, unzip them
rule unzip:
    input:
        join("input", "{f}.{e}.gz"),
    output:
        unzipped=join("input", "{f}.{e}"),
        zipped=join("input", "zipped", "{f}.{e}.gz"),
    wildcard_constraints:
        f="\w+",
        e="fasta|fna|faa|fa",
    script:
        resource_filename("snekmer", join("scripts", "unzip.py"))


# build kmer count vectors for each basis set
use rule vectorize from kmerize with:
    input:
        fasta=join("input", "{f}.{e}"),
    output:
        data=join(out_dir, "vector", "{f}.{e}.npz"),
        kmerobj=join(out_dir, "kmerize", "{f}.{e}.kmers"),
    wildcard_constraints:
        f="\w+",
        e="fasta|fna|faa|fa",
    log:
        join(out_dir, "kmerize", "log", "{f}.{e}.log"),


# [in-progress] kmer walk
# if config['walk']:
# use rule perform_kmer_walk from process with:
# output:


# UNSUPERVISED WORKFLOW
# collect all seq files and generate mega-cluster
rule cluster:
    input:
        kmerobj=expand(
            join(out_dir, "kmerize", "{f}.{e}.kmers"),
            f=seq_input.filename,
            e=seq_input.ext,
        ),
        data=expand(
            join(out_dir, "vector", "{f}.{e}.npz"),
            f=seq_input.filename,
            e=seq_input.ext,
        ),
    output:
        figs=directory(join(out_dir, "cluster", "figures")),
        table=join(out_dir, "cluster", "snekmer.csv"),
    log:
        join(out_dir, "cluster", "log", "cluster.log"),
    script:
        resource_filename("snekmer", join("scripts", "cluster.py"))


rule cluster_report:
    input:
        figs=rules.cluster.output.figs,
        table=rules.cluster.output.table,
    output:
        join(out_dir, "Snekmer_Cluster_Report.html"),
    script:
        resource_filename("snekmer", join("scripts", "cluster_report.py"))
