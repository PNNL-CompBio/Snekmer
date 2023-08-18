# imports
import pickle
from Bio import SeqIO
import snekmer as skm


# initialize kmerization object
kmerizer = skm.vectorize.KmerVecs(
    k=snakemake.config["k"],
    alphabet=snakemake.config["alphabet"],
    ref_id=snakemake.params.ref_ids[snakemake.wildcards.f],
)

# TODO: read kmerbasis if present

min_filter = 0
# if hasattr(input, "kmerbasis") and exists(input.kmerbasis):
#     kmerbasis = skm.io.read_kmers(input.kmerbasis)

#     # quick way to get the number of proteins in
#     #     the fasta file so we can set up an array
#     #     ahead of time
#     nprot = len([1 for line in open(input.fasta) if line.startswith(">")])

# else:
#     # we make our own kmerbasis and filter for minimum
#     #    number of occurrences, etc.
#     # we will only allow filtering by number of kmers
#     # if we're not using a basis set as input -
#     if "min_filter" in config:
#         min_filter = config["min_filter"]

#         # make basis
#     kmerbasis = {}
#     fasta = SeqIO.parse(input.fasta, "fasta")

#     nprot = 0
#     for f in fasta:
#         nprot += 1
#         these = kmer.reduce_vectorize(f.seq)
#         for key in these:
#             if key in kmerbasis:
#                 kmerbasis[key] += 1
#             else:
#                 kmerbasis[key] = 1

#     kmerbasis = np.array(list(kmerbasis.keys()))[
#         np.array(list(kmerbasis.values())) > min_filter
#     ]

# kmer.set_kmer_set(kmerbasis)

# (re)read fasta using bioconda obj
fasta = SeqIO.parse(snakemake.input.fasta, "fasta")
seqs = [str(f.seq) for f in fasta]
kmerizer.kmerize(seqs)

# unpack table
# codes = kmerizer.kmer_codes
# table = kmerizer.table
# matrix = skm.vectorize.unpack_table(table, len(seqs), kmers=codes)

# pre-allocate an array to keep results
# vecs = np.zeros((nprot, len(kmerbasis)))

# I question whether we need to keep the reduced seqs here
# seqs, ids, lengths = list(), list(), list()
# n = 0
# for f in fasta:
#     addvec = kmer.reduce_vectorize(f.seq)
#     vecs[n][np.isin(kmerbasis, addvec)] = 1
#     n += 1
#     seqs.append(
#         skm.vectorize.reduce(
#             f.seq,
#             alphabet=config["alphabet"],
#             mapping=skm.alphabet.FULL_ALPHABETS,
#         )
#     )
#     ids.append(f.id)
#     lengths.append(len(f.seq))

# save seqIO output and transformed vecs
# np.savez_compressed(
#     output.data,
#     kmerlist=kmerbasis,
#     ids=ids,
#     seqs=seqs,
#     vecs=vecs,
#     lengths=lengths,
# )

skm.io.save_pickle(kmerizer, snakemake.output.kmertable)
