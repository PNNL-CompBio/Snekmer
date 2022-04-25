"""snekmer.smk: v1.0.0 code

author: @christinehc

"""
# snakemake config
from snakemake.utils import min_version

min_version("6.0")  # force snakemake v6.0+ (required for modules)

# imports
import pickle
from datetime import datetime
from glob import glob
from itertools import product
from os import makedirs
from os.path import basename, dirname, exists, join, splitext

import matplotlib.pyplot as plt
import numpy as np
import snekmer as skm
from Bio import SeqIO
from pandas import DataFrame
from sklearn.model_selection import StratifiedKFold


# load modules
module process_input:
    snakefile:
        "process_input.smk"
    config:
        config


# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
input_files = glob(join("input", "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input"]["file_extensions"])
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
    config["alphabet"], config["k"], nested=config["output"]["nested_dir"]
)


# define output files to be created by snekmer
rule all:
    input:
        expand(join("input", "{uz}"), uz=UZS),  # require unzipping
        # expand(join("output", "kmerize", "{nb}", "{fa}.json.gz"), nb=NON_BGS, fa=FAS),  # correctly build features
        expand(join("output", "scoring", "{nb}.pkl"), nb=NON_BGS),  # require model-building


# if any files are gzip zipped, unzip them
use rule unzip from process_input with:
    output:
        join("input", "{uz}"),


# define rules
# read and process parameters from config
rule vectorize:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        data=join("output", "vector", "{nb}.npz"),
        # vecs=join("output", "vector", "{nb}.npy"),
        # ids=join("output", "seq_id", "{nb}.npy"),
        kmerobj=join("output", "kmerize", "{nb}.pkl"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),
    run:
        fasta = SeqIO.parse(input.fasta, "fasta")

        # initialize kmerization object
        kmer = skm.vectorize.KmerVec(alphabet=config["alphabet"], k=config["k"])

        vecs, seqs, ids = list(), list(), list()
        for f in fasta:
            vecs.append(kmer.reduce_vectorize(f.seq))
            # seqs.append(
            #     reduce(f.seq, alphabet=config["alphabet"], mapping=skm.alphabet.FULL_ALPHABETS)
            # )
            ids.append(f.id)

        # save seqIO output and transformed vecs
        # np.save(output.vecs, vecs)
        # np.save(output.seqs, seqs)
        # np.save(output.ids, ids)
        np.savez_compressed(output.data, ids=ids, seqs=seqs, vecs=vecs)

        with open(output.kmerobj, "wb") as f:
            pickle.dump(kmer, f)


# def get_ids_vecs(wildcards):
#     """Retrieve ID file and vec file in matching order"""
#     ids = expand(join("output", "seq_id", "{fa}.npy"), fa=FILES)
#     basenames = [basename(f) for f in ids]
#     vecs = [join("output", "vector", f) for f in basenames]
#     return ids, vecs


rule score:
    input:
        kmerobj=join("output", "kmerize", "{nb}.pkl"),
        data=join("output", "vector", "{nb}.npz"),
        # vecs=expand(join("output", "vector", "{fa}.npy"), fa=NON_BGS),
        # ids=expand(join("output", "seq_id", "{fa}.npy"), fa=NON_BGS),
    output:
        df=join("output", "scoring", "sequences", "{nb}.csv.gz"),
        scores=join("output", "scoring", "weights", "{nb}.csv.gz"),
        scorer=join("output", "scoring", "{nb}.pkl"),
    log:
        join("output", "scoring", "log", "{nb}.log"),
    run:
        # log script start time
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        # get kmers for this particular set of sequences
        with open(input.kmerobj, "rb") as f:
            kmer = pickle.load(f)

        # parse all data and label background files
        label = config["score"]["lname"]

        loaded = np.load(input.data)

        for ids, vecs in zip(loaded["ids"], loaded["vecs"]):
            # ids = np.load(idf)
            # vecs = np.load(vecf)
            data = DataFrame(
                {
                    "filename": splitext(basename(idf))[0],
                    "sequence_id": list(ids),
                    "sequence_vector": list(vecs),
                }
            )

        data["background"] = [f in BGS for f in data["filename"]]

        # log conversion step runtime
        skm.utils.log_runtime(log[0], start_time, step="files_to_df")

        # parse family names and only add if some are valid
        families = [
            skm.utils.get_family(
                skm.utils.split_file_ext(fn)[0], regex=config["input"]["regex"]
            )
            for fn in data["filename"]
        ]
        if any(families):
            data[label] = families

        # binary T/F for classification into family
        family = skm.utils.get_family(wildcards.nb)
        binary_labels = [True if value == family else False for value in data[label]]

        # define k-fold split indices
        if config["model"]["cv"] > 1:
            cv = StratifiedKFold(n_splits=config["model"]["cv"], shuffle=True)

            # stratify splits by [0,1] family assignment
            for n, (i_train, _) in enumerate(
                cv.split(data["sequence_vector"], binary_labels)
            ):
                data[f"train_cv-{n + 1:02d}"] = [idx in i_train for idx in data.index]

        elif config["model"]["cv"] in [0, 1]:
            i_train, _ = train_test_split(data.index, stratify=binary_labels)
            data["train"] = [idx in i_train for idx in data.index]

        # generate family scores and object
        scorer = skm.model.KmerScorer()
        scorer.fit(
            list(kmer.kmer_set.kmers),
            data,
            skm.utils.get_family(wildcards.nb, regex=config["input"]["regex"]),
            label_col=label,
            vec_col="sequence_vector",
            **config["score"]["scaler_kwargs"],
        )

        # append scored sequences to dataframe
        data = data.merge(
            DataFrame(scorer.scores["sample"]), left_index=True, right_index=True
        )
        if data.empty:
            raise ValueError("Blank df")

        # save score loadings
        class_probabilities = (
            DataFrame(scorer.probabilities, index=scorer.kmers.basis)
            .reset_index()
            .rename(columns={"index": "kmer"})
        )

        # log time to compute class probabilities
        skm.utils.log_runtime(log[0], start_time, step="class_probabilities")

        # save all files to respective outputs
        delete_cols = ["vec", "sequence_vector"]
        for col in delete_cols:
            # if col in data.columns:
            #     data = data.drop(columns=col)
            if col in class_probabilities.columns:
                class_probabilities = class_probabilities.drop(columns=col)
        data.to_csv(output.df, index=False, compression="gzip")
        class_probabilities.to_csv(output.scores, index=False, compression="gzip")
        with open(output.scorer, "wb") as f:
            pickle.dump(scorer, f)

        # record script endtime
        skm.utils.log_runtime(log[0], start_time)
