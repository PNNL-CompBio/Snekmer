# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import random
from typing import List
from Bio import SeqIO
import snekmer as skm


# -----------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

random.seed(params.seed)  # setting the seed for randomness


def fragment(
    sequence: str,
    version: str,
    length: int,
    location: str,
    minLength: int
    ) -> List[str]:
    """
    Fragment a given sequence.

    Args:
        sequence (str): the sequence to fragment.
        version (str): fragmentation method, either "absolute" or "percent".
        length (int): length for fragmentation. If version is "percent", this is treated as a percentage.
        location (str): where the fragmentation happens - "start", "end", or "random".
        min_length (int): minimum length a fragment should be to be retained.

    Returns:
        List[str]: the retained fragment(s).
    """
    frags = []
    if version == "absolute":
        for i in range(0, len(sequence) - length + 1):
            frags.append(sequence[i : i + length])

    elif version == "percent":
        actualLength = int(len(sequence) * (length / 100))
        for i in range(0, len(sequence) - actualLength + 1):
            frags.append(sequence[i : i + actualLength])

    # Filter fragments based on minLength
    frags = [frag for frag in frags if len(frag) >= minLength]

    # Retention logic based on the location parameter
    if location == "start":
        if frags:
            frags = [frags[0]]
    elif location == "end":
        if frags:
            frags = [frags[-1]]
    elif location == "random":
        if frags:
            chosenIndex = random.randint(0, len(frags) - 1)
            frags = [frags[chosenIndex]]

    return frags

with open(input.fasta, "r") as f:
    fastaSequences = SeqIO.parse(f, "fasta")

    with open(output.fasta_out, "w") as file:
        for fasta in fastaSequences:
            titleLine, sequence = fasta.description, str(fasta.seq)

            fragments = fragment(
                sequence,
                params.version,
                params.fragLength,
                params.location,
                params.minLength,
            )

            for i, frag in enumerate(fragments):
                file.write(
                    ">" + titleLine + " Fragment=" + str(i) + "\n"
                )
                file.write(frag + "\n")