# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import random

from typing import List

from Bio import SeqIO

import snekmer as skm
from snekmer.learn import fragment


# -----------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

random.seed(snakemake.params.seed)  # setting the seed for randomness


with open(snakemake.input.fasta, "r") as f:
    fasta_sequences = SeqIO.parse(f, "fasta")

    with open(snakemake.output.fasta_out, "w") as file:
        for fasta in fasta_sequences:
            title_line, sequence = fasta.description, str(fasta.seq)

            fragments = fragment(
                sequence,
                snakemake.params.version,
                snakemake.params.frag_length,
                snakemake.params.location,
                snakemake.params.min_length
            )

            for i, frag in enumerate(fragments):
                file.write(
                    ">" + title_line + " Fragment=" + str(i) + "\n"
                )
                file.write(frag + "\n")