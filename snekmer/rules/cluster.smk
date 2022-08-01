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
import gzip
import json
import pickle
import struct
from datetime import datetime
from glob import glob
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, splitext, split

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import snekmer as skm
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import (StratifiedKFold, cross_val_score,
                                     train_test_split)
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier
# from umap import UMAP

# make this conditional on the library being installed
import bsf

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
input_files = glob(join("input", "*"))
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
bg_files = glob(join("input", "background", "*"))
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
        join(out_dir, "cluster", "cluster.clust"),  # require cluster-building


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
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),


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
    output:
        clusters=join(out_dir, "cluster", "cluster.clust"),
        #bsf=directory(join(out_dir, "cluster", "bsf")),
        #bsfmat=join(out_dir,"cluster","bsf_matrix.csv"),
        #distmat=join(out_dir,"cluster","distance_matrix.csv"),
        figs=directory(join(out_dir, "cluster", "figures")),
        table=join(out_dir, "cluster", "clusters.csv")
    log:
        join(out_dir, "cluster", "log", "cluster.log"),
    run:
        # log script start time
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        # parse all data and label background files
        label = config["score"]["lname"]

        # assuming all kmer basis sets are identical, grab the first
        kmer = skm.io.load_pickle(input.kmerobj[0])

        # tabulate vectorized seq data
        data = list()
        for f in input.data:
            data.append(skm.io.load_npz(f))

        data = pd.concat(data, ignore_index=True)
        data["background"] = [f in BGS for f in data["filename"]]

        # log conversion step runtime
        skm.utils.log_runtime(log[0], start_time, step="load_npz")

        # define feature matrix of kmer vectors not from background set
        bg, non_bg = data[data["background"]], data[~data["background"]]
        full_feature_matrix = skm.utils.to_feature_matrix(data["sequence_vector"].values)
        feature_matrix = skm.utils.to_feature_matrix(non_bg["sequence_vector"].values)

        # currently not used
        if len(bg) > 0:
            bg_feature_matrix = skm.utils.to_feature_matrix(bg["sequence_vector"].values)

        # for bsf we need to make the input matrix a distance matrix
        if config["cluster"]["method"] in ["bsf","hdbsf", "aggbsf"]:
            output_bsf = join(split(output.clusters)[0], "bsf")
            if not exists(output_bsf):
                makedirs(output_bsf)

            model = skm.cluster.KmerClustering(config["cluster"]["method"], config["cluster"]["params"])

            # the uint64 type is required by bsf
            full_feature_matrix = np.array(full_feature_matrix, dtype=np.uint64)

            # the filename that bsf creates you can't specify exactly
            # but we can reconstruct it here - this will likely break
            # because I think bsf breaks large input matrices into chunks
            # and will output multiple files
            nsign, vlen = full_feature_matrix.shape

            bsfname = join(output_bsf, "bin_%d_0_0_%d_%d_bsf.txt.bin" % (nsign,nsign,nsign))

            # maybe we've already generated this matrix?
            # THIS DOESN'T WORK SINCE SNAKEMAKE REMOVES ALL THESE FILES
            #     BEFORE STARTING
            if not exists(bsfname):
                # We can break this job in to chunks by setting the nsign to smaller
                #    than the row size - but we'll need to figure out how to read
                #    in the bin files and use them
                #bsf.analysis_with_chunk(full_feature_matrix, 10000, "bsf.txt", "%s/" % output.bsf)
                bsf.analysis_with_chunk(full_feature_matrix, nsign, "bsf.txt", "%s/" % output_bsf)
            else:
                print("Using existing BSF similarity matrix...")

            with open(bsfname, "rb") as f:
                bsf_mat = np.frombuffer(f.read(), dtype=np.uint32).reshape((nsign, nsign))
                #bsf_mat = np.frombuffer(f.read(), dtype=np.double).reshape((nsign, nsign))

            # this is great for diagnostics - but takes a lot of time/space for large
            # comparisons
            #np.savetxt(output.bsfmat, bsf_mat, fmt="%.3f", delimiter=",")

            # bsf only outputs the upper triangle so we'll copy that to the lower triangle
            bsf_mat = bsf_mat + bsf_mat.T - np.diag(np.diag(bsf_mat))

            # and fill the diagonals with max which is skipped by bsf
            np.fill_diagonal(bsf_mat, 100)

            # now transform the similarity matrix output by bsf to a
            #     distance matrix - using inverse Jaccard similarity

            # this gives Jaccard similarity (since all vectors are the same length)
            # and subtracting it from 1 gives a distance
            #bsf_mat = 1.0-(bsf_mat/100)
            bsf_mat = 100-bsf_mat

            if "dist_thresh" in config["cluster"]:
                dist_thresh = config["cluster"]["dist_thresh"]
                bsf_mat[bsf_mat > dist_thresh] = 100

            # this is great for diagnostics - but takes a lot of time/space for large
            # comparisons
            #np.savetxt(output.distmat, bsf_mat, fmt="%.3f", delimiter=",")

            data["cluster"] = model.predict(bsf_mat)
            with open(output.clusters, "wb") as f:
                pickle.dump(model, f)


            data = data.drop(columns=["sequence_vector"]).to_csv(output.table, index=False)
            skm.utils.log_runtime(log[0], start_time, step="clustering")

        else:
            # fit and save fitted clusters
            # bsf here
            model = skm.cluster.KmerClustering(config["cluster"]["method"], config["cluster"]["params"])

            model.fit(full_feature_matrix)
            # save output
            #data["cluster"] = model.predict(data["sequence_vector"].values)
            data["cluster"] = model.predict(full_feature_matrix)

            with open(output.clusters, "wb") as f:
                pickle.dump(model, f)


            data = data.drop(columns=["sequence_vector"]).to_csv(output.table, index=False)

            # log time to compute clusters
            skm.utils.log_runtime(log[0], start_time, step="clustering")

        if not exists(output.figs):
            makedirs(output.figs)

        if "cluster_plots" in config["cluster"] and config["cluster"]["cluster_plots"] == "True":
        # insert plots here?
            fig, ax = skm.plot.explained_variance_curve(full_feature_matrix)
            fig.savefig(join(output.figs, "pca_explained_variance_curve.png"))
            plt.close("all")

            # plot tsne
            fig, ax = skm.plot.cluster_tsne(full_feature_matrix, model.model.labels_)
            fig.savefig(join(output.figs, "tsne_clusters.png"))
            plt.close("all")

        # plot umap
        # model_embedding = model.predict(full_feature_matrix)
        # umap_embedding = UMAP(metric="jaccard", n_components=2).fit_transform(full_feature_matrix)
        # fig, ax = plt.subplots(dpi=150)
        # sns.scatterplot(
        #     x=umap_embedding[:, 0],
        #     y=umap_embedding[:, 1],
        #     hue=model.model.labels_,
        #     alpha=0.2,
        #     ax=ax,
        # )

        # fig.savefig(join(output.figs, "umap_clusters.png"))
        # plt.close("all")

        # record script endtime
        skm.utils.log_runtime(log[0], start_time)
