"""model.smk: Module for supervised kmer-based annotation models.

author: @christinehc

"""

# snakemake config
from snakemake.utils import min_version

min_version("6.0")  # force snakemake v6.0+ (required for modules)


# load modules
module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


# imports
from glob import glob
from itertools import product
from os.path import basename, dirname, isdir, join, splitext
from pkg_resources import resource_filename

import matplotlib.pyplot as plt
import numpy as np
import snekmer as skm
from Bio import SeqIO
from sklearn.linear_model import LogisticRegression


# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])


# collect all fasta-like files, unzipped filenames, and basenames
gz_input = glob_wildcards(
    join("input", "zipped", "{filename,\w+}.{ext,fasta|fna|faa|fa}.gz")
)
seq_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}"))
bg_input = glob_wildcards(
    join("input", "background", "{filename,\w+}.{ext,fasta|fna|faa|fa}")
)
print(seq_input, bg_input, gz_input)
# gz_exts = [e.split(".")[0] for e in gz_input.raw_ext]

if len(bg_input.filename) > 0:

    ruleorder: unzip > vectorize > vectorize_background > score_background > score_with_background > model > model_report

else:

    ruleorder: unzip > vectorize > score > model > model_report


# check input file size
for f, e in zip(seq_input.filename, seq_input.ext):
    skm.utils.check_n_seqs(
        join("input", f"{f}.{e}"), config["model"]["cv"], show_warning=False
    )

# add unzipped gz files to total input list
for f, e in zip(gz_input.filename, gz_input.ext):
    getattr(seq_input, "filename").append(f)
    getattr(seq_input, "ext").append(e)
print(seq_input)

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
        for f, e in zip(seq_input.filename, seq_input.ext)
    ]


# define output files to be created by snekmer
output = [
    expand(
        join("input", "{gzf}.{gze}"),
        zip,
        gzf=gz_input.filename,
        gze=gz_input.ext,
    ),
    expand(
        join(out_dir, "kmerize", "{f}.{e}.kmers"),
        f=seq_input.filename,
        e=seq_input.ext,
    ),
    expand(
        join(out_dir, "model", "{sf}.model"), sf=seq_input.filename
    ),  # require model-building for non-background seqs
    join(out_dir, "Snekmer_Model_Report.html"),
]
if len(bg_input.filename) > 0:
    output.append(
        expand(
            join(out_dir, "background", "{bf}.{be}.npz"),
            bf=bg_input.filename,
            be=bg_input.ext,
        )
    )
    output.append(join(out_dir, "background", "combined_background.npz"))


rule all:
    input:
        output,
        # expand(join(out_dir, "vector", "background", "{b}.npz", b=BGS)),


rule unzip:
    input:
        join("input", "zipped", "{gzf}.{gze}.gz"),
    output:
        unzipped=join("input", "{gzf}.{gze}"),
        # zipped=join("input", "zipped", "{gzf}.{gze}.gz"),
    wildcard_constraints:
        gzf="\w+",
        gze="fasta|fna|faa|fa",
    script:
        resource_filename("snekmer", join("scripts", "unzip.py"))


# build kmer count vectors for each basis set
use rule vectorize from kmerize with:
    input:
        fasta=join("input", "{f}.{e}"),
    output:
        data=join(out_dir, "vector", "{f}.{e}.npz"),
        kmerobj=join(out_dir, "kmerize", "{f}.{e}.kmers"),
    log:
        join(out_dir, "kmerize", "log", "{f}.{e}.log"),


if len(bg_input.filename) > 0:

    use rule vectorize from kmerize as vectorize_background with:
        input:
            fasta=join("input", "background", "{bf}.{be}"),
        output:
            data=join(out_dir, "background", "vector", "{bf}.{be}.npz"),
            kmerobj=join(out_dir, "background", "kmerize", "{bf}.{be}.kmers"),
        log:
            join(out_dir, "background", "kmerize", "log", "{bf}.{be}.log"),

    rule score_background:
        input:
            data=join(out_dir, "background", "vector", "{bf}.{be}.npz"),
            kmerobj=join(out_dir, "background", "kmerize", "{bf}.{be}.kmers"),
        output:
            scores=join(out_dir, "background", "{bf}.{be}.npz"),
        script:
            resource_filename(
                "snekmer", join("scripts", "model_vectorize_background.py")
            )

    rule combine_background:
        input:
            scores=expand(
                join(out_dir, "background", "{bf_a}.{be_a}.npz"),
                bf_a=bg_input.filename,
                be_a=bg_input.ext,
            ),
        output:
            join(out_dir, "background", "combined_background.npz"),
        script:
            resource_filename("snekmer", join("scripts", "model_combine_background.py"))

    rule score_with_background:
        input:
            data=expand(join(out_dir, "vector", "{sf}.npz"), sf=seq_input.filename),
            kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
            bg=join(out_dir, "background", "combined_background.npz"),
        output:
            data=join(out_dir, "scoring", "sequences", "{f}.csv.gz"),
            weights=join(out_dir, "scoring", "weights", "{f}.csv.gz"),
            scorer=join(out_dir, "scoring", "{f}.scorer"),
            matrix=join(out_dir, "scoring", "{f}.matrix"),
        log:
            join(out_dir, "scoring", "log", "{f}.log"),
        script:
            resource_filename("snekmer", join("scripts/model_score_background.py"))

else:

    # build family score basis
    rule score:
        input:
            kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
            data=expand(join(out_dir, "vector", "{sf}.npz"), sf=seq_input.filename),
        output:
            data=join(out_dir, "scoring", "sequences", "{f}.csv.gz"),
            weights=join(out_dir, "scoring", "weights", "{f}.csv.gz"),
            scorer=join(out_dir, "scoring", "{f}.scorer"),
            matrix=join(out_dir, "scoring", "{f}.matrix"),
        log:
            join(out_dir, "scoring", "log", "{f}.log"),
        script:
            resource_filename("snekmer", join("scripts/model_score.py"))


rule model:
    input:
        kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
        raw=join(out_dir, "vector", "{f}.npz"),
        data=join(out_dir, "scoring", "sequences", "{f}.csv.gz"),
        weights=join(out_dir, "scoring", "weights", "{f}.csv.gz"),
        matrix=join(out_dir, "scoring", "{f}.matrix"),
    output:
        model=join(out_dir, "model", "{f}.model"),
        results=join(out_dir, "model", "results", "{f}.csv"),
        figs=directory(join(out_dir, "model", "figures", "{f}")),
    script:
        resource_filename("snekmer", join("scripts/model_model.py"))


rule model_report:
    input:
        results=expand(
            join(out_dir, "model", "results", "{sf}.csv"), sf=seq_input.filename
        ),
        figs=expand(join(out_dir, "model", "figures", "{sf}"), sf=seq_input.filename),
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
