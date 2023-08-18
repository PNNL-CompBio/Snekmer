# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version

min_version("6.0")


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

# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])

# collect all fasta-like files, unzipped filenames, and basenames
gz_input = glob_wildcards(join("input", "[!zipped/_]{filename}.{ext}.gz"))
all_input = glob_wildcards(
    join("input", f"{{filename}}.{{ext,({'|'.join(config['input_file_exts'])})}}")
)

# check input file size
for f, e in zip(all_input.filename, all_input.ext):
    skm.utils.check_n_seqs(
        join("input", f"{f}.{e}"), config["model"]["cv"], show_warning=False
    )

# add unzipped gz files to total input list
for f, e in zip(gz_input.filename, gz_input.ext):
    getattr(all_input, "filename").append(f)
    getattr(all_input, "ext").append(e)

# define unique reference ids by kmer file and get seq counts for each file
NSEQS = {
    f: skm.utils.count_n_seqs(join("input", f"{f}.{e}"))
    for f, e in zip(all_input.filename, all_input.ext)
}
FA_REFIDS = skm.utils.get_ref_ids(
    [join("input", f"{f}.{e}") for f, e in zip(all_input.filename, all_input.ext)]
)

# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

# show warnings if files excluded
onstart:
    [
        skm.utils.check_n_seqs(
            join("input", f"{f}.{e}"), config["model"]["cv"], show_warning=True
        )
        for f, e in zip(all_input.filename, all_input.ext)
    ]


# define output files to be created by snekmer
rule all:
    input:
        expand(
            join("input", "zipped", "{gzf}.{gze}.gz"),
            zip,
            gzf=gz_input.filename,
            gze=gz_input.ext,
        ),
        join(out_dir, "Snekmer_Cluster_Report.html"),


rule unzip:
    input:
        join("input", "{gzf}.{gze}.gz"),
    output:
        unzipped=join("input", "{gzf}.{gze}"),
        zipped=join("input", "zipped", "{gzf}.{gze}.gz"),
    script:
        resource_filename("snekmer", join("scripts", "unzip.py"))


# build kmer count vectors for each basis set
rule kmerize:
    input:
        fasta=join("input", "{f}.{e}"),
    output:
        # data=join(out_dir, "vector", "{nb}.npz"),
        kmertable=join(out_dir, "kmerize", "{f}.{e}.ktable"),
    params:
        ref_ids=FA_REFIDS,
    log:
        join(out_dir, "kmerize", "log", "{f}.{e}.log"),
    script:
        resource_filename("snekmer", join("scripts", "kmerize.py"))


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
