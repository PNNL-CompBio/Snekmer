from datetime import datetime

import pandas as pd
import snekmer as skm
from sklearn.model_selection import StratifiedKFold

rule score:
    input:
        kmerobj=join("output", "kmerize", "{nb}.pkl"),
        data=expand(join("output", "vector", "{fa}.npz"), fa=NON_BGS),
    output:
        data=join("output", "scoring", "sequences", "{nb}.csv.gz"),
        weights=join("output", "scoring", "weights", "{nb}.csv.gz"),
        scorer=join("output", "scoring", "{nb}.pkl"),
    log:
        join("output", "scoring", "log", "{nb}.log"),
    run:
        # log script start time
        start_time = datetime.now()
        label = (
            config["score"]["lname"]
            if str(config["score"]["lname"]) != "None"
            else "label"
        )  # e.g. "family"
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        # get kmers for this particular set of sequences
        kmer = skm.io.load_pickle(input.kmerobj)

        # tabulate vectorized seq data
        data = list()
        for f in input.data:
            data.append(skm.io.load_npz(f))

        data = pd.concat(data, ignore_index=True)
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
            pd.DataFrame(scorer.scores["sample"]), left_index=True, right_index=True
        )
        if data.empty:
            raise ValueError("Blank df")

        # save score loadings
        class_probabilities = (
            pd.DataFrame(scorer.probabilities, index=scorer.kmers.basis)
            .reset_index()
            .rename(columns={"index": "kmer"})
        )

        # log time to compute class probabilities
        skm.utils.log_runtime(log[0], start_time, step="class_probabilities")

        # save all files to respective outputs
        delete_cols = ["vec", "sequence_vector"]
        for col in delete_cols:
            if col in class_probabilities.columns:
                class_probabilities = class_probabilities.drop(columns=col)
        data.drop(columns="sequence_vector").to_csv(
        output.data, index=False, compression="gzip"
        )
        class_probabilities.to_csv(output.weights, index=False, compression="gzip")
        with open(output.scorer, "wb") as f:
            pickle.dump(scorer, f)

        # record script endtime
        skm.utils.log_runtime(log[0], start_time)
