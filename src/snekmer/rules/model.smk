"""model.smk: Module for supervised kmer-based annotation models.

author: @christinehc

"""
from collections import namedtuple

# snakemake config
from snakemake.utils import min_version

min_version("6.0")  # force snakemake v6.0+ (required for modules)


# imports
from glob import glob
from itertools import product
from os.path import basename, dirname, join, splitext
from pkg_resources import resource_filename

import matplotlib.pyplot as plt
import numpy as np
import snekmer as skm
from Bio import SeqIO
from sklearn.linear_model import LogisticRegression


# # change matplotlib backend to non-interactive
# plt.switch_backend("Agg")


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
        expand(
            join(out_dir, "scoring", "sequences", "{f}.{e}.csv.gz"),
            zip,
            f=all_input.filename,
            e=all_input.ext,
        ),
        join(out_dir, "kmerize", "merged.ktable"),
        # expand(join(out_dir, "model", "{nb}.model"), nb=NON_BGS),  # require model-building
        # join(out_dir, "Snekmer_Model_Report.html"),


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


# build family score basis
rule merge_vectors:
    input:
        all_seqs=expand(
            join("input", "{af}.{ae}"), zip, af=all_input.filename, ae=all_input.ext
        ),
        all_tables=expand(
            join(out_dir, "kmerize", "{af}.{ae}.ktable"),
            zip,
            af=all_input.filename,
            ae=all_input.ext,
        ),
    output:
        labels=join(out_dir, "kmerize", "labels.npz"),
        table=join(out_dir, "kmerize", "merged.ktable"),
    params:
        nseqs=NSEQS,
    script:
        resource_filename("snekmer", join("scripts", "merge_vectors.py"))


rule score:
    input:
        family_table=join(out_dir, "kmerize", "{f}.{e}.ktable"),
        all_seqs=expand(
            join("input", "{af}.{ae}"), zip, af=all_input.filename, ae=all_input.ext
        ),
        labels=join(out_dir, "kmerize", "labels.npz"),
        merged_table=join(out_dir, "kmerize", "merged.ktable"),
        # bg=BGS,
    output:
        scores=join(out_dir, "scoring", "sequences", "{f}.{e}.csv.gz"),
        weights=join(out_dir, "scoring", "weights", "{f}.{e}.csv.gz"),
        scorer=join(out_dir, "scoring", "{f}.{e}.scorer"),
        # matrix=join(out_dir, "scoring", "{nb}.matrix"),
    log:
        join(out_dir, "score", "log", "{f}.{e}.log"),
    params:
        nseqs=NSEQS,
    script:
        resource_filename("snekmer", join("scripts", "score.py"))


rule model:
    input:
        raw=rules.score.input.merged_table,
        # data=rules.score.output.data,
        weights=rules.score.output.weights,
        kmertable=rules.score.input.family_table,
        scores=rules.score.output.scores,
    output:
        model=join(out_dir, "model", "{f}.{e}.model"),
        results=join(out_dir, "model", "results", "{f}.{e}.csv"),
        figs=directory(join(out_dir, "model", "figures", "{f}.{e}")),
    script:
        resource_filename("snekmer", join("scripts", "model.py"))


rule model_report:
    input:
        results=expand(
            join(out_dir, "model", "results", "{af}.{ae}.csv"),
            zip,
            af=all_input.filename,
            ae=all_input.ext,
        ),
        figs=expand(
            join(out_dir, "model", "figures", "{af}.{ae}"),
            zip,
            af=all_input.filename,
            ae=all_input.ext,
        ),
    output:
        join(out_dir, "Snekmer_Model_Report.html"),
    run:
        fig_dir = dirname(input.figs[0])
        tab_dir = dirname(input.results[0])

        model_vars = dict(
            page_title="Snekmer Model Report",
            title="Snekmer Model Results",
            text=(
                "Classifier model results "
                f"({config['model']['cv']}-Fold Cross-Validation) "
                f"are below."
            ),
            dir=skm.report.correct_rel_path(tab_dir),
        )

        skm.report.create_report_many_images(
            fig_dir, tab_dir, model_vars, "model", output[0]
        )
