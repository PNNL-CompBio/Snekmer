# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import random
from typing import List

from Bio import SeqIO

import snekmer as skm
from snekmer.learn import fragment


# ---------------------------------------------------------
# Core logic
# ---------------------------------------------------------

def run_fragmentation(
    fasta_in: str,
    fasta_out: str,
    seed: int,
    version: str,
    frag_length: int,
    location: str,
    min_length: int,
) -> None:
    """
    Fragment sequences from a FASTA file and write fragments to a new FASTA.

    Parameters
    ----------
    fasta_in : str
        Path to input FASTA file.
    fasta_out : str
        Path to output FASTA file where fragments will be written.
    seed : int
        Random seed for reproducible fragmentation.
    version : str
        Fragmentation version (passed through to `snekmer.learn.fragment`).
    frag_length : int
        Target fragment length.
    location : str
        Fragmentation location mode (passed to `fragment`).
    min_length : int
        Minimum fragment length to retain.
    """
    random.seed(seed)

    with open(fasta_in, "r") as f:
        fasta_sequences = SeqIO.parse(f, "fasta")

        with open(fasta_out, "w") as file:
            for fasta in fasta_sequences:
                title_line, sequence = fasta.description, str(fasta.seq)

                fragments: List[str] = fragment(
                    sequence,
                    version,
                    frag_length,
                    location,
                    min_length,
                )

                for i, frag in enumerate(fragments):
                    file.write(f">{title_line} Fragment={i}\n")
                    file.write(frag + "\n")


# ---------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------

if "snakemake" in globals():
    run_fragmentation(
        fasta_in=snakemake.input.fasta,
        fasta_out=snakemake.output.fasta_out,
        seed=snakemake.params.seed,
        version=snakemake.params.version,
        frag_length=snakemake.params.frag_length,
        location=snakemake.params.location,
        min_length=snakemake.params.min_length,
    )
