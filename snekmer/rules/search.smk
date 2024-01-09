"""search.smk: Module for supervised modeling from kmer vectors.

author: @christinehc, @biodataganache

"""
# snakemake config
from snakemake.utils import min_version

min_version("6.0")  # force snakemake v6.0+ (required for modules)


ruleorder: vectorize > search


module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


# built-in imports
import gzip
import json
import pickle
from ast import literal_eval
from datetime import datetime
from glob import glob
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, splitext
from pkg_resources import resource_filename

# import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier

# external libraries
import snekmer as skm

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
gz_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}.gz"))
seq_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}"))

model_input = glob_wildcards(
    join(config["model_dir"], "{family,\w+}.{ext,fasta|fna|faa|fa}.model")
)

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


# define output files to be created by snekmer
rule all:
    input:
        join(config["basis_dir"], "search_kmers.txt"),  # require common basis
        expand(
            join(out_dir, "vector", "{f}.{e}.npz"),
            f=seq_input.filename,
            e=seq_input.ext,
        ),
        expand(
            join(out_dir, "search", "{m}.{me}", "{f}.{e}.csv"),
            m=model_input.family,
            me=model_input.ext,
            f=seq_input.filename,
            e=seq_input.ext,
        ),
        # require search
        join(out_dir, "Snekmer_Search_Report.html"),


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


rule common_basis:  # build kmer count vectors for each basis set
    input:
        kmerobjs=expand(
            join(config["basis_dir"], "{m}.{me}.kmers"),
            m=model_input.family,
            me=model_input.ext,
        ),
    output:
        kmerbasis=join(config["basis_dir"], "search_kmers.txt"),
    log:
        join(config["basis_dir"], "log", "common_basis.log"),
    run:
        common_basis = list()
        for kobj in input.kmerobjs:
            kmers = skm.io.load_pickle(kobj)

            # consolidate overly long lists of duplicate kmers
            if len(common_basis) > 1e10:
                common_basis = list(set(common_basis))
            common_basis.extend(list(kmers.kmer_set.kmers))

            # capture common basis set -- is faster than np.unique
        common_basis = set(common_basis)
        common_basis = sorted(list(common_basis))

        # output type: plaintext (no pandas) would likely be more compact
        with open(output.kmerbasis, "w") as f:
            for kmer in common_basis:
                f.write(f"{kmer}\n")
                # df = pd.DataFrame({'common': common_basis})
                # df.to_csv(output.kmerbasis, index=False)



use rule vectorize from kmerize with:
    input:
        fasta=join("input", "{f}.{e}"),
        kmerbasis=rules.common_basis.output.kmerbasis,
    output:
        data=join(out_dir, "vector", "{f}.{e}.npz"),
        kmerobj=join(out_dir, "kmerize", "{f}.{e}.kmers"),
    wildcard_constraints:
        f="\w+",
        e="fasta|fna|faa|fa",
    log:
        join(out_dir, "kmerize", "log", "{f}.{e}.log"),


rule search:
    input:
        vecs=join(out_dir, "vector", "{f}.{e}.npz"),  # change to data=join(out_dir, "vector", "{nb}.npz")
        model=join(config["model_dir"], "{m}.{me}.model"),
        kmerobj=join(config["basis_dir"], "{m}.{me}.kmers"),
        scorer=join(config["score_dir"], "{m}.{me}.scorer"),
    output:
        results=join(out_dir, "search", "{m}.{me}", "{f}.{e}.csv"),
    script:
        resource_filename("snekmer", join("scripts", "search.py"))


rule search_report:
    input:
        files=expand(
            join(out_dir, "search", "{m}.{me}", "{f}.{e}.csv"),
            m=model_input.family,
            me=model_input.ext,
            f=seq_input.filename,
            e=seq_input.ext,
        ),
    output:
        join(out_dir, "Snekmer_Search_Report.html"),
    run:
        file_dir = dirname(dirname(input.files[0]))

        # search
        search_vars = dict(
            page_title="Snekmer Search Report",
            title="Snekmer Search Results",
            text="See the below links to access Snekmer Search results.",
        )

        skm.report.create_report_many_csvs(file_dir, search_vars, "search", output[0])
