# imports
from os.path import basename
import snekmer as skm

# simplify variable name
family = snakemake.wildcards.m

# get kmers for this particular set of sequences
# print(f"starting {family}")
kmer = skm.io.load_pickle(snakemake.input.kmerobj)
model = skm.io.load_pickle(snakemake.input.model)
scorer = skm.io.load_pickle(snakemake.input.scorer)
# print(f"loaded model {family}")

# load vectorized sequences, score, and predict scores
kmerlist, df = skm.io.load_npz(snakemake.input.vecs)
filename = snakemake.wildcards.f

# print(f"making feature matrix {family}")
vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)

# print(f"getting scores {family}")
scores = scorer.predict(vecs, kmerlist[0])
# print(f"making predictions {family}")
predictions = model.predict(scores.reshape(-1, 1))
# print(f"getting probabilities {family}")
predicted_probas = model.model.predict_proba(scores.reshape(-1, 1))

# display results (score, family assignment, and probability)
df["score"] = scores  # scorer output
df["in_family"] = [True if p == 1 else False for p in predictions]
df["probability"] = [p[1] for p in predicted_probas]
df["filename"] = f"{filename}.{FILE_MAP[filename]}"
df["model"] = basename(snakemake.input.model)

df = df.drop(columns=["sequence_vector", "sequence"])
df.to_csv(snakemake.output.results, index=False)