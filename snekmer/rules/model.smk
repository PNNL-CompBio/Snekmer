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
gz_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}.gz"))
seq_input = glob_wildcards(join("input", "{filename,\w+}.{ext,fasta|fna|faa|fa}"))
bg_input = glob_wildcards(
    join("input", "background", "{filename,\w+}.{ext,fasta|fna|faa|fa}")
)
bgz_input = glob_wildcards(
    join("input", "background", "{filename,\w+}.{ext,fasta|fna|faa|fa}.gz")
)

if len(bg_input.filename) > 0:

    ruleorder: unzip > vectorize > vectorize_background > score > model > model_report

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
for f, e in zip(bgz_input.filename, bgz_input.ext):
    getattr(bg_input, "filename").append(f)
    getattr(bg_input, "ext").append(e)

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

    # raise error if no background files supplied in a bg mode
    # if (len(bg_input.filename) == 0) and (
    #     config["score"]["method"].lower() != "family"
    # ):
    #     raise FileNotFoundError(
    #         f"Score method `{config['score']['method']}`"
    #         " requires background files (none found)."
    #     )



# define output files to be created by snekmer
output = [
    # expand(
    #     join("input", "{gzf}.{gze}"),
    #     zip,
    #     gzf=gz_input.filename,
    #     gze=gz_input.ext,
    # ),
    expand(
        join(out_dir, "kmerize", "{f}.{e}.kmers"),
        zip,
        f=seq_input.filename,
        e=seq_input.ext,
    ),
    expand(
        join(out_dir, "model", "{sf}.{se}.model"),
        zip,
        sf=seq_input.filename,
        se=seq_input.ext,
    ),  # require model-building for non-background seqs
    join(out_dir, "Snekmer_Model_Report.html"),
]
if len(bg_input.filename) > 0:
    output.append(
        expand(
            join(out_dir, "scoring", "{f}.{e}.matrix"),
            zip,
            f=seq_input.filename,
            e=seq_input.ext,
        )
    )
    # output.append(
    #     expand(
    #         join(out_dir, "background", "{bf}.{be}.npz"),
    #         zip,
    #         bf=bg_input.filename,
    #         be=bg_input.ext,
    #     )
    # )
    # output.append(
    #     expand(
    #         join(out_dir, "background", "{f}.{e}.npz"),
    #         zip,
    #         f=seq_input.filename,
    #         e=seq_input.ext,
    #     )
    # )


rule all:
    input:
        output,
        # expand(join(out_dir, "vector", "background", "{b}.npz", b=BGS)),


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


if len(bg_input.filename) > 0:

    rule unzip_background:
        input:
            join("input", "background", "{bf}.{be}.gz"),
        output:
            unzipped=join("input", "background", "{bf}.{be}"),
            zipped=join("input", "background", "zipped", "{bf}.{be}.gz"),
        wildcard_constraints:
            bf="\w+",
            be="fasta|fna|faa|fa",
        script:
            resource_filename("snekmer", join("scripts", "unzip.py"))

    use rule vectorize from kmerize as vectorize_background with:
        input:
            fasta=join("input", "background", "{bf}.{be}"),
        output:
            data=join(out_dir, "background", "vector", "{bf}.{be}.npz"),
            kmerobj=join(out_dir, "background", "kmerize", "{bf}.{be}.kmers"),
        wildcard_constraints:
            bf="\w+",
            be="fasta|fna|faa|fa",
        log:
            join(out_dir, "background", "kmerize", "log", "{bf}.{be}.log"),

    rule combine_background:
        input:
            kmerobj=join(out_dir, "kmerize", "{f}.{e}.kmers"),
            data=expand(
                join(out_dir, "background", "vector", "{bf_a}.{be_a}.npz"),
                zip,
                bf_a=bg_input.filename,
                be_a=bg_input.ext,
            ),
        output:
            join(out_dir, "background", "{f}.{e}.npz"),
        script:
            resource_filename("snekmer", join("scripts", "combine_background.py"))


def score_input(wildcards):
    input_files = {
        "kmerobj": join(out_dir, "kmerize", f"{wildcards.f}.{wildcards.e}.kmers"),
        "data": expand(
            join(out_dir, "vector", "{sf}.{se}.npz"),
            zip,
            sf=seq_input.filename,
            se=seq_input.ext,
        ),
    }
    if config["score"]["method"] in ["background", "bg", "combined"]:
        input_files["bg"] = join(
            out_dir, "background", f"{wildcards.f}.{wildcards.e}.npz"
        )
    return input_files


# build family score basis
rule score:
    input:
        unpack(score_input),
    output:
        data=join(out_dir, "scoring", "sequences", "{f}.{e}.csv.gz"),
        weights=join(out_dir, "scoring", "weights", "{f}.{e}.csv.gz"),
        scorer=join(out_dir, "scoring", "{f}.{e}.scorer"),
        matrix=join(out_dir, "scoring", "{f}.{e}.matrix"),
    log:
        join(out_dir, "scoring", "log", "{f}.{e}.log"),
    script:
        resource_filename("snekmer", join("scripts", "score_with_background.py"))


def model_input(wildcards):
    input_files = {
        "kmerobj": join(out_dir, "kmerize", f"{wildcards.f}.{wildcards.e}.kmers"),
        "raw": join(out_dir, "vector", f"{wildcards.f}.{wildcards.e}.npz"),
        "data": join(
            out_dir, "scoring", "sequences", f"{wildcards.f}.{wildcards.e}.csv.gz"
        ),
        "weights": join(
            out_dir, "scoring", "weights", f"{wildcards.f}.{wildcards.e}.csv.gz"
        ),
        "matrix": join(out_dir, "scoring", f"{wildcards.f}.{wildcards.e}.matrix"),
    }
    if config["score"]["method"] in ["background", "bg", "combined"]:
        input_files["bg"] = join(
            out_dir, "background", f"{wildcards.f}.{wildcards.e}.npz"
        )
    return input_files


rule model:
    input:
        unpack(model_input),
    output:
        model=join(out_dir, "model", "{f}.{e}.model"),
        results=join(out_dir, "model", "results", "{f}.{e}.csv"),
        figs=directory(join(out_dir, "model", "figures", "{f}.{e}")),
    script:
        resource_filename("snekmer", join("scripts/model_model.py"))


rule model_report:
    input:
        results=expand(
            join(out_dir, "model", "results", "{sf}.{se}.csv"),
            zip,
            sf=seq_input.filename,
            se=seq_input.ext,
        ),
        figs=expand(
            join(out_dir, "model", "figures", "{sf}.{se}"),
            zip,
            sf=seq_input.filename,
            se=seq_input.ext,
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
