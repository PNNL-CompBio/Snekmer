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
from os.path import basename, dirname, join, splitext
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
input_files = glob(join(input_dir, "*.*"))
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

# get seq counts for each file
NSEQS = {skm.utils.split_file_ext(f)[0]: skm.utils.count_n_seqs(f) for f in unzipped}

# final file map: checks that files are large enough for model building
# FA_MAP = {
#     k: v for k, v in f_map.items() if skm.utils.check_n_seqs(k, config["model"]["cv"])
# }

# get unzipped filenames
UZS = [
    f"{f}.{ext}"
    for f, ext in UZ_MAP.items()
    if skm.utils.check_n_seqs(fa, config["model"]["cv"], show_warning=False)
]

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
        expand(join(out_dir, "scoring", "sequences", "{nb}.csv.gz"), nb=NON_BGS),
        # expand(join(out_dir, "model", "{nb}.model"), nb=NON_BGS),  # require model-building
        # join(out_dir, "Snekmer_Model_Report.html"),


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join("input", "{uz}"),
        zipped=join("input", "zipped", "{uz}.gz"),


# build kmer count vectors for each basis set
use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        # data=join(out_dir, "vector", "{nb}.npz"),
        kmertable=join(out_dir, "kmerize", "{nb}.ktable"),
    log:
        join(out_dir, "kmerize", "log", "{nb}.log"),


# build family score basis
rule merge_vectors:
    input:
        all_seqs=expand(join("input", "{uz}"), uz=UZS),
        all_tables=expand(join(out_dir, "kmerize", "{fa}.ktable"), fa=NON_BGS),
    output:
        labels=join(out_dir, "kmerize", "labels.npz"),
        table=join(out_dir, "kmerize", "merged.ktable"),
    params:
        nseqs=NSEQS,
    script:
        resource_filename("snekmer", join("scripts", "merge_vectors.py"))


rule score:
    input:
        family_table=join(out_dir, "kmerize", "{nb}.ktable"),
        all_seqs=expand(join("input", "{uz}"), uz=UZS),
        labels=join(out_dir, "kmerize", "labels.npz"),
        merged_table=join(out_dir, "kmerize", "merged.ktable"),
        bg=BGS,
    output:
        scores=join(out_dir, "scoring", "sequences", "{nb}.csv.gz"),
        weights=join(out_dir, "scoring", "weights", "{nb}.csv.gz"),
        scorer=join(out_dir, "scoring", "{nb}.scorer"),
        # matrix=join(out_dir, "scoring", "{nb}.matrix"),
    log:
        join(out_dir, "score", "log", "{nb}.log"),
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
        model=join(out_dir, "model", "{nb}.model"),
        results=join(out_dir, "model", "results", "{nb}.csv"),
        figs=directory(join(out_dir, "model", "figures", "{nb}")),
    script:
        resource_filename("snekmer", join("scripts", "model.py"))


rule model_report:
    input:
        results=expand(join(out_dir, "model", "results", "{nb}.csv"), nb=NON_BGS),
        figs=expand(join(out_dir, "model", "figures", "{nb}"), nb=NON_BGS),
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
