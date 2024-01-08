# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

from datetime import datetime
from os import makedirs
from os.path import dirname, exists, join

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import snekmer as skm
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from umap import UMAP

# ---------------------------------------------------------
# (optional) BSF handling
# ---------------------------------------------------------

try:
    # make this conditional on the library being installed
    import bsf

    BSF_PRESENT = True

except ImportError:
    BSF_PRESENT = False

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

# log script start time
start_time = datetime.now()
with open(snakemake.log[0], "a") as f:
    f.write(f"start time:\t{start_time}\n")

# parse all data and label background files
label = config["score"]["lname"]

# assuming all kmer basis sets are identical, grab the first
kmer = skm.io.load_pickle(snakemake.input.kmerobj[0])

# tabulate vectorized seq data
data, kmers = list(), list()
for dfile, kfile in zip(sorted(snakemake.input.data), sorted(snakemake.input.kmerobj)):
    kmerlist, df = skm.io.load_npz(dfile)
    data.append(df)

    kobj = skm.io.load_pickle(kfile)
    klist = kobj.kmer_set._kmerlist
    kmers.append(klist)
# basis.extend(klist)

# make a superset of kmers
kmerbasis = np.unique(np.hstack(kmers))

basis = skm.vectorize.KmerVec(config["alphabet"], config["k"])
basis.set_kmer_set(kmerbasis)

for i in range(len(data)):
    df, kmerlist = data[i], kmers[i]
    vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)

    converted = basis.harmonize(vecs, kmerlist)
    df["sequence_vector"] = converted.tolist()

data = pd.concat(data, ignore_index=True)
# data["background"] = [f in snakemake.input.bg for f in data["filename"]]

# log conversion step runtime
skm.utils.log_runtime(snakemake.log[0], start_time, step="load_npz")

# define feature matrix of kmer vectors not from background set
# bg, non_bg = data[data["background"]], data[~data["background"]]
full_feature_matrix = skm.utils.to_feature_matrix(data["sequence_vector"].values)
# feature_matrix = skm.utils.to_feature_matrix(non_bg["sequence_vector"].values)

# filter by kmer occurrence
# print("Filtering?")
if str(config["cluster"]["min_rep"]) != "None":
    min_rep = max_rep = config["cluster"]["min_rep"]
    print(f"Filtering by minimum {min_rep}")
    ksums = full_feature_matrix.sum(axis=0)
    full_feature_matrix = full_feature_matrix[:, ksums > min_rep]
    print(full_feature_matrix.shape)

if str(config["cluster"]["max_rep"]) != "None":
    max_rep = config["cluster"]["max_rep"]
    print(f"Filtering by maximum {max_rep}")
    ksums = full_feature_matrix.sum(axis=0)
    full_feature_matrix = full_feature_matrix[:, ksums < max_rep]
    print(full_feature_matrix.shape)

# currently not used
# if len(bg) > 0:
#     bg_feature_matrix = skm.utils.to_feature_matrix(bg["sequence_vector"].values)
# kludge to check on options for agglomerative Clustering
if (
    "n_clusters" in config["cluster"]["params"]
    and config["cluster"]["params"]["n_clusters"] == "None"
):
    config["cluster"]["params"]["n_clusters"] = None
    config["cluster"]["params"]["compute_full_tree"] = True

# initialize clustering model
model = skm.cluster.KmerClustering(
    config["cluster"]["method"], config["cluster"]["params"]
)

# for bsf, make the input matrix a distance matrix
if config["cluster"]["method"] in [
    "density-jaccard",
    "hdensity-jaccard",
    "agglomerative-jaccard",
]:
    # print BSF status
    output_bsf = join(dirname(snakemake.output.table), "bsf")
    if not exists(output_bsf):
        makedirs(output_bsf)

    # if available, use BSF to create clustering distance matrix
    if BSF_PRESENT:
        print("Using BSF to compute distance matrix...")
        full_feature_matrix = np.array(
            full_feature_matrix, dtype=np.uint64
        )  # bsf requires uint64

        # the filename that bsf creates you can't specify exactly
        # but we can reconstruct it here - this will likely break
        # because I think bsf breaks large input matrices into chunks
        # and will output multiple files
        nsign, vlen = full_feature_matrix.shape

        bsfname = join(
            output_bsf, "bin_%d_0_0_%d_%d_bsf.txt.bin" % (nsign, nsign, nsign)
        )

        if not exists(bsfname):
            # We can break this job in to chunks by setting the nsign to smaller
            #    than the row size - but we'll need to figure out how to read
            #    in the bin files and use them
            # bsf.analysis_with_chunk(full_feature_matrix, 10000, "bsf.txt", "%s/" % output.bsf)
            bsf.analysis_with_chunk(
                full_feature_matrix, nsign, "bsf.txt", "%s/" % output_bsf
            )
        else:
            print("Using existing BSF similarity matrix...")

        with open(bsfname, "rb") as f:
            bsf_mat = np.frombuffer(f.read(), dtype=np.uint32).reshape((nsign, nsign))
        # bsf_mat = np.frombuffer(f.read(), dtype=np.double).reshape((nsign, nsign))

        # this is great for diagnostics - but takes a lot of time/space for large
        # comparisons
        # np.savetxt(output.bsfmat, bsf_mat, fmt="%.3f", delimiter=",")

        # bsf only outputs the upper triangle so we'll copy that to the lower triangle
        bsf_mat = bsf_mat + bsf_mat.T - np.diag(np.diag(bsf_mat))

        # and fill the diagonals with max which is skipped by bsf
        np.fill_diagonal(bsf_mat, 100)

        # now transform the similarity matrix output by bsf to a
        #     distance matrix - using inverse Jaccard similarity

        # this gives Jaccard similarity (since all vectors are the same length)
        # and subtracting it from 1 gives a distance
        # bsf_mat = 1.0 - (bsf_mat / 100)
        bsf_mat = (100 - bsf_mat) / 100.0

    else:
        # If bsf isn't available we'll use scipy to
        # calculate jaccard distance in a similar way
        # but it will be slower and more memory intensive
        # especially for larger clustering jobs
        # calculate Jaccard distance
        res = pdist(full_feature_matrix, "jaccard")
        bsf_mat = squareform(res)

    # dist_thresh = config["cluster"]["dist_thresh"]
    # bsf_mat[bsf_mat > dist_thresh] = 100

    # output matrix for diagnostics - time/space heavy for large comparisons
    if config["cluster"]["save_matrix"] is True:
        np.savetxt(
            join(output_bsf, "bsf_matrix.csv"), bsf_mat, fmt="%.3f", delimiter=",",
        )

    model.fit(bsf_mat)
else:
    # fit and save fitted clusters
    # bsf here
    model.fit(full_feature_matrix)

# save output
data["cluster"] = model.labels_
data = data.drop(columns=["sequence", "sequence_vector"])
data.to_csv(snakemake.output.table, index=False)

# with open(output.clusters, "wb") as f:
#     pickle.dump(model, f)

# always create output figure directory
if not exists(snakemake.output.figs):
    makedirs(snakemake.output.figs)

# log time to compute clusters
skm.utils.log_runtime(snakemake.log[0], start_time, step="clustering")

# optionally generate plots
if str(config["cluster"]["cluster_plots"]) == "True":
    # plot explained variance curve
    fig, ax = skm.plot.explained_variance_curve(full_feature_matrix)
    fig.savefig(join(snakemake.output.figs, "pca_explained_variance_curve.png"))
    plt.close("all")

    # plot tsne
    fig, ax = skm.plot.cluster_tsne(full_feature_matrix, model.model.labels_)
    fig.savefig(join(snakemake.output.figs, "tsne.png"))
    plt.close("all")

    # plot umap
    umap_embedding = UMAP(metric="jaccard", n_components=2).fit_transform(
        full_feature_matrix
    )
    fig, ax = plt.subplots(dpi=150)
    sns.scatterplot(
        x=umap_embedding[:, 0],
        y=umap_embedding[:, 1],
        hue=model.model.labels_,
        alpha=0.2,
        ax=ax,
    )

    fig.savefig(join(snakemake.output.figs, "umap.png"))
    plt.close("all")

# record script endtime
skm.utils.log_runtime(snakemake.log[0], start_time)
