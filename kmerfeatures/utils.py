"""utils: Utility functions for the kmerfeatures pipeline.

author: @christinehc
"""
# imports
from Bio import SeqIO


# functions
def read_fasta(fasta):
    """Read fasta file and obtain sequences, IDs, and mappings.

    Parameters
    ----------
    fasta : str
        Filename or /path/to/fasta/file.

    Returns
    -------
    (list, dict)
        Tuple containing (id_list, id2seq), where:
        id_list : list
            List of IDs from fasta file.
        id2seq: dict
            Dictionary mapping IDs to sequences (id2seq[id] = seq).

    """
    # read in sequences from the fasta file
    id_list, id2seq = [], {}
    with open(fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            id_list.append(record.id)
            id2seq[record.id] = str(record.seq)
    return id_list, id2seq


def read_example_index(example_indexfile):
    """Read example index file and parse into dictionary.

    Parameters
    ----------
    example_indexfile : str
        /path/to/example_indexfile

    Returns
    -------
    dict
        Dictionary of example indices.

    """
    example_index = {}
    with open(example_indexfile, "r") as f:
        for line in f.readlines():
            seq_id = line.split()[0]
            example_index[seq_id] = 1.0
    return example_index
