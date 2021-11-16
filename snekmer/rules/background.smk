# imports
# built-in imports
import json
import pickle
from datetime import datetime
from glob import glob
from itertools import (product, repeat)
from multiprocessing import Pool
from os import makedirs
from os.path import (basename, dirname, exists, join, splitext)

# external libraries
import snekmer as skm
import numpy as np
import matplotlib.pyplot as plt
from pandas import (DataFrame, read_csv, read_json)
from Bio import SeqIO

# change matplotlib backend to non-interactive
plt.switch_backend('Agg')

# do i even need this? would prolly specify files in the main snakefile
# define background files
bg_files = config['input']
if isinstance(bg_files, str):
    bg_files = [bg_files]
BGS = [
    skm.utils.split_file_ext(basename(f))[0] for flist in [
        glob(join("input", bg)) for bg in bg_files
    ]
    for f in flist
]
# NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files

# main rule
rule all:
    input:
        expand(join("{bg}.pkl"), bg=BGS)

# load main workflow and use to process background files
module main_workflow:
    snakefile: "../Snakefile"
    config: config

# unzip any zipped bg files
use rule unzip from main_workflow as main_unzip with:
    input:
        join("input", "{bg}.fasta")
        # lambda wildcards: join("input", f"{wildcards.uz}.{uz_map[wildcards.uz]}.gz")
    output:
        join("input", "{uz}.{uzext}")
    params:
        outdir=join("input", 'zipped')


# input: background files
# output: background score for a given family given the input background files

# step 1: process bg files using the same params as in the snakefile
# this would prolly occur via the main workflow?
rule kmerize:
    input:
        fasta=join("input", "{bg}.fasta")
        params="kmer params",
    output:
        join("output", "features", "{bg}.txt"),  # kmer vectors
        join("output", "features", "{bg}.csv")  # kmer scores
    run:

# step 1.9 ?: compute full kmer vectors?

# step 2: standardize bg files with the same basis set as family of interest
rule standardize:
    input:
    output:

# step 3: score background
rule score_background:
    input:
        join("{bg}")
    output:
        join("{bg} score")
    run:
        # get kmers for this particular set of sequences
        kmers = skm.io.read_output_kmers(input.kmers)

        # parse all data and label background files
        # modify this to only parse a single background file
        label = config['lname']
        data = skm.io.vecfiles_to_df(
            input.files, labels=config['score']['labels'], label_name=label
        )
        data['background'] = [skm.utils.split_file_ext(f)[0] in BGS for f in data['filename']]

        # parse family names and only add if some are valid
        families = [
            skm.utils.get_family(fn, regex=config['input']['regex'])
            for fn in data['filename']
            ]
        if any(families):
            label = 'family'
            data[label] = families

        # define feature matrix of kmer vectors not from background set
        bg, non_bg = data[data['background']], data[~data['background']]
        full_feature_matrix = skm.score.to_feature_matrix(data['vector'].values)
        feature_matrix = skm.score.to_feature_matrix(non_bg['vector'].values)
        bg_feature_matrix = skm.score.to_feature_matrix(bg['vector'].values)
        # print(feature_matrix.T.shape, bg_feature_matrix.T.shape, np.array(kmers).shape)

        # compute class probabilities
        labels = non_bg[label].values
        class_probabilities = skm.score.feature_class_probabilities(
            feature_matrix.T, labels, kmers=kmers
        )  #.drop(columns=['count'])

        # compute background sequence probabilities
        bg_labels = bg[label].values
        # print(labels)
        # print(bg_labels)
        bg_probabilities = skm.score.feature_class_probabilities(
            bg_feature_matrix.T, bg_labels, kmers=kmers
        )

        # background family probability scores
        for fam in bg[label].unique():
            bg_scores = bg_probabilities[
                bg_probabilities['label'] == fam
            ]['score'].values

            # normalize by max bg score
            # bg_norm = np.max(bg_scores)

            # get background scores for sequences
            bg_only_scores = skm.score.apply_feature_probabilities(
                full_feature_matrix, bg_scores, scaler=config['score']['scaler'],
                **config['score']['scaler_kwargs']
            )

            # normalize by max bg score attained by a sequence
            bg_norm = np.max(bg_only_scores)
            print("bg_norm", bg_norm)

            data[f"{fam}_background_score"] = bg_only_scores / bg_norm

        # assign family probability scores
        for fam in non_bg[label].unique():
            scores = class_probabilities[
                class_probabilities['label'] == fam
            ]['score'].values

            # normalize by sum of all positive scores
            # norm = np.sum([s for s in scores if s > 0])

            # include background sequences for score generation
            total_scores = skm.score.apply_feature_probabilities(
                full_feature_matrix, scores, scaler=config['score']['scaler'],
                **config['score']['scaler_kwargs']
            )

            # normalize by max score
            norm = np.max(total_scores)
            if norm == 0:
                norm = 1.0
            print("norm", norm)

            # assign percent score based on max positive score
            data[f"{fam}_score"] = total_scores / norm

            # weight family score by (1 - normalized bg score)
            if fam in bg[label].unique():
                data[f"{fam}_score_background_weighted"] = [
                    total * (1 - bg) for total, bg in zip(
                        data[f"{fam}_score"], data[f"{fam}_background_score"]
                    )
                ]

                # old scoring method
                data[f"{fam}_background_subtracted_score"] = [
                    total - bg for total, bg in zip(
                        data[f"{fam}_score"], data[f"{fam}_background_score"]
                    )
                ]

        # [IN PROGRESS] compute clusters
        clusters = skm.score.cluster_feature_matrix(full_feature_matrix)
        data['cluster'] = clusters

        # save all files to respective outputs
        delete_cols = ['vec', 'vector']
        for col in delete_cols:
            if col in data.columns:
                data = data.drop(columns=col)
            if col in class_probabilities.columns:
                class_probabilities = class_probabilities.drop(columns=col)
        data.to_csv(output.df, index=False)
        class_probabilities.to_csv(output.scores, index=False)


# use rule preprocess from main_workflow as main_preprocess with:
# use rule generate_kmers from main_workflow as main_generate_kmers with:
