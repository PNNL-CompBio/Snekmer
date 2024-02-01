# imports
from os.path import basename
import numpy as np
import snekmer as skm

# get kmers for this particular set of sequences
family = snakemake.wildcards.m
print(f"starting {family}")
kmer = skm.io.load_pickle(snakemake.input.kmerobj)
model = skm.io.load_pickle(snakemake.input.model)
scorer = skm.io.load_pickle(snakemake.input.scorer)
print(f"loaded model {family}")

# load vectorized sequences, score, and predict scores
kmerlist, df = skm.io.load_npz(snakemake.input.vecs)

print(f"making feature matrix {family}")
vecs = np.vstack(df["sequence_vector"].values)
df = df.drop(columns=["sequence_vector", "sequence"])  # clear memory

print(f"getting scores {family}")
scores = scorer.predict(vecs, kmerlist[0])
print(f"making predictions {family}")
predictions = model.predict(scores.reshape(-1, 1))
print(f"getting probabilities {family}")
predicted_probas = model.model.predict_proba(scores.reshape(-1, 1))

# display results (score, family assignment, and probability)
df["score"] = scores  # scorer output
df["in_family"] = [p == 1 for p in predictions]
df["probability"] = [p[1] for p in predicted_probas]
df["filename"] = f"{snakemake.wildcards.f}.{snakemake.wildcards.e}"
df["model"] = basename(snakemake.input.model)

# df = df.drop(columns=["sequence_vector", "sequence"])
df.to_csv(snakemake.output.results, index=False)
