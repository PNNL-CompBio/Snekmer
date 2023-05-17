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
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = glob(join(input_dir, "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]

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
        expand(join("input", "{uz}"), uz=UZS),  # require unzipping
        join(out_dir, "Snekmer_Cluster_Report.html"),


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join("input", "{uz}"),
        zipped=join("input", "zipped", "{uz}.gz"),


# build kmer count vectors for each basis set
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


# [in-progress] kmer walk
# if config['walk']:
# use rule perform_kmer_walk from process with:
# output:


# UNSUPERVISED WORKFLOW
# collect all seq files and generate mega-cluster
rule cluster:
    input:
        kmerobj=expand(join("output", "kmerize", "{fa}.kmers"), fa=NON_BGS),
        data=expand(join("output", "vector", "{fa}.npz"), fa=NON_BGS),
        bg=BGS,
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
    run:
        # check for figures
        if str(config["cluster"]["cluster_plots"]) == "True":
            fig_params = {
                "image1_name": "PCA Explained Variance",
                        "image1_path": skm.report.correct_rel_path(
                            join(input.figs, "pca_explained_variance_curve.png")
                        ),
                        "image2_name": "Clusters (UMAP)",
                        "image2_path": skm.report.correct_rel_path(
                            join(input.figs, "umap.png")
                        ),
                        "image3_name": "Clusters (t-SNE)",
                        "image3_path": skm.report.correct_rel_path(
                            join(input.figs, "tsne.png")
                        ),
                    }
                else:
                    fig_params = {
                        "image1_name": "",
                        "image1_path": None,
                        "image2_name": "",
                        "image2_path": None,
                        "image3_name": "",
                        "image3_path": None,
                    }

                    # cluster
                cluster_vars = {
                    "page_title": "Snekmer Cluster Report",
                    "title": "Snekmer Cluster Results",
                    "text": (
                "Snekmer clustering results are linked below. "
                "If `cluster_plots` are enabled in the config, "
                "they will be shown below."
            ),
            "dir": dirname(skm.report.correct_rel_path(input.table)),
            "table": skm.report.correct_rel_path(input.table),
            **fig_params,
        }

        skm.report.create_report(cluster_vars, "cluster", output[0])
