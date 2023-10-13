"""model.smk: Module for supervised kmer-based annotation models.

author: @christinehc

"""
# snakemake config
from snakemake.utils import min_version

min_version("6.0")  # force snakemake v6.0+ (required for modules)


# load modules
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


# # change matplotlib backend to non-interactive
# plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = [f for f in glob(join(input_dir, "*")) if not isdir(f)]
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
    and skm.utils.check_n_seqs(fa, config["model"]["cv"], show_warning=False)
]


# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# final file map: checks that files are large enough for model building
# FA_MAP = {
#     k: v for k, v in f_map.items() if skm.utils.check_n_seqs(k, config["model"]["cv"])
# }

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


# show warnings if files excluded
onstart:
    [
        skm.utils.check_n_seqs(fa, config["model"]["cv"], show_warning=True)
        for fa in input_files
    ]


# define output files to be created by snekmer
rule all:
    input:
        expand(join("input", "{uz}"), uz=UZS),  # require unzipping
        expand(join(out_dir, "model", "{f}.model"), f=FAS),  # require model-building
        join(out_dir, "Snekmer_Model_Report.html"),


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join("input", "{uz}"),
        zipped=join("input", "zipped", "{uz}.gz"),


# build kmer count vectors for each basis set
use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.f}.{FA_MAP[wildcards.f]}"),
    output:
        data=join(out_dir, "vector", "{f}.npz"),
        kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
    log:
        join(out_dir, "kmerize", "log", "{f}.log"),


# build family score basis
rule score:
    input:
        kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
        data=expand(join(out_dir, "vector", "{nb}.npz"), nb=NON_BGS),
    output:
        data=join(out_dir, "scoring", "sequences", "{f}.csv.gz"),
        weights=join(out_dir, "scoring", "weights", "{f}.csv.gz"),
        scorer=join(out_dir, "scoring", "{f}.scorer"),
        matrix=join(out_dir, "scoring", "{f}.matrix"),
    params:
        bg=BGS,
    log:
        join(out_dir, "scoring", "log", "{f}.log"),
    script:
        resource_filename("snekmer", join("scripts/model_score.py"))


rule model:
    input:
        raw=rules.score.input.data,
        data=rules.score.output.data,
        weights=rules.score.output.weights,
        kmerobj=rules.score.input.kmerobj,
        matrix=rules.score.output.matrix,
    output:
        model=join(out_dir, "model", "{f}.model"),
        results=join(out_dir, "model", "results", "{f}.csv"),
        figs=directory(join(out_dir, "model", "figures", "{f}")),
    script:
        resource_filename("snekmer", join("scripts/model_model.py"))


rule model_report:
    input:
        results=expand(join(out_dir, "model", "results", "{f}.csv"), f=FAS),
        figs=expand(join(out_dir, "model", "figures", "{f}"), f=FAS),
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
