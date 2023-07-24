"""model.smk: Module for supervised kmer-based annotation models.

author: @tnitka

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
input_files = glob(join(input_dir, "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
#    and skm.utils.check_n_seqs(fa, config["model"]["cv"], show_warning=False)
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
# onstart:
#    [
#        skm.utils.check_n_seqs(fa, config["model"]["cv"], show_warning=True)
#        for fa in input_files
#    ]

# set number of permutations to test
n_iter = (config["motif"]["n"])


# define output files to be created by snekmer
rule all:
    input:
        expand(join("input", "{uz}"), uz=UZS),  # require unzipping
#    	join(config["basis_dir"], "search_kmers.txt"), # require common basis
#        expand(join(out_dir, "model", "{nb}.model"), nb=NON_BGS),  # require model-building
    	expand(join(out_dir, "scoring", "weights", "{nb}.csv.gz"), nb=NON_BGS), # require scoring
    	expand(join(out_dir, "motif", "p_values", "{nb}.csv.gz"), nb=NON_BGS), # require motif identification

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
        data=join(out_dir, "vector", "{nb}.npz"),
        kmerobj=join(out_dir, "kmerize", "{nb}.kmers"),
    log:
        join(out_dir, "kmerize", "log", "{nb}.log"),

#rule common_basis:  # build kmer count vectors for each basis set
#    input:
#        kmerobjs=expand(join(config["basis_dir"], "{nb}.kmers"), nb=NON_BGS),
#    output:
#        kmerbasis=join(config["basis_dir"], "search_kmers.txt"),
#    log:
#        join(config["basis_dir"], "log", "common_basis.log"),
#    run:
#        common_basis = list()
#        for kobj in input.kmerobjs:
#            kmers = skm.io.load_pickle(kobj)

            # consolidate overly long lists of duplicate kmers
#            if len(common_basis) > 1e10:
#                common_basis = list(set(common_basis))
#            common_basis.extend(list(kmers.kmer_set.kmers))

        # capture common basis set -- is faster than np.unique
#        common_basis = set(common_basis)
#        common_basis = sorted(list(common_basis))

        # output type: plaintext (no pandas) would likely be more compact
#        with open(output.kmerbasis, "w") as f:
#            for kmer in common_basis:
#                f.write(f"{kmer}\n")


# build family score basis
rule score:
    input:
        kmerobj=join(out_dir, "kmerize", "{nb}.kmers"),
        data=expand(join(out_dir, "vector", "{fa}.npz"), fa=NON_BGS),
        bg=BGS,
    output:
        data=join(out_dir, "scoring", "sequences", "{nb}.csv.gz"),
        weights=join(out_dir, "scoring", "weights", "{nb}.csv.gz"),
        scorer=join(out_dir, "scoring", "{nb}.scorer"),
        matrix=join(out_dir, "scoring", "{nb}.matrix"),
    log:
        join(out_dir, "scoring", "log", "{nb}.log"),
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
        model=join(out_dir, "model", "{nb}.model"),
        results=join(out_dir, "model", "results", "{nb}.csv"),
        figs=directory(join(out_dir, "model", "figures", "{nb}")),
    script:
        resource_filename("snekmer", join("scripts/model_model.py"))
        
# rule permute:
#    input:
#        raw=rules.score.input.data,
#        data=rules.score.output.data,
#        weights=rules.score.output.weights,
#        kmers=rules.common_basis.output.kmerbasis,
#        kmerobj=rules.score.input.kmerobj,
#        matrix=rules.score.output.matrix,
#    output:
#        data=temp(join(out_dir, "motif", "perm_data", "{nb}_{n}.csv.gz")),
#    script:
#        resource_filename("snekmer", join("scripts/motif_permute.py"))
        
rule rescore:
    input:
        raw=rules.score.input.data,
        data=rules.score.output.data,
        weights=rules.score.output.weights,
#        kmers=rules.common_basis.output.kmerbasis,
        kmerobj=rules.score.input.kmerobj,
        matrix=rules.score.output.matrix,
    output:
        data=temp(join(out_dir, "motif", "scores", "{nb}_{n}.csv.gz")),
    script:
        resource_filename("snekmer", join("scripts/motif_rescore.py"))
        
rule motif:
    input:
        raw=rules.score.input.data,
        data=rules.score.output.data,
        weights=rules.score.output.weights,
        perm_scores=expand(join(out_dir, "motif", "scores", "{{nb}}_{n}.csv.gz"), n=range(n_iter)),
#        kmers=rules.common_basis.output.kmerbasis,
        kmerobj=rules.score.input.kmerobj,
        matrix=rules.score.output.matrix,
    output:
        data=join(out_dir, "motif", "scores", "{nb}.csv.gz"),
        p_values=join(out_dir, "motif", "p_values", "{nb}.csv.gz"),
    script:
        resource_filename("snekmer", join("scripts/motif_motif.py"))
    