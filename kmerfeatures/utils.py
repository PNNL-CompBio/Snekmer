"""utils: Utility functions for the kmerfeatures pipeline.

author: @christinehc
"""
# imports
from Bio import SeqIO


# functions
def read_fasta(fasta, include_map=True):
    """Read fasta file and obtain sequences, IDs, and mappings.

    Parameters
    ----------
    fasta : str
        Filename or /path/to/fasta/file.
    include_map : bool
        If True, outputs the dictionary mapping of IDs to sequences.
        If False, only outputs the sequence list and ID list.
        (default: True)

    Returns
    -------
    (list, list) or (list, list, dict)
        Tuple containing (seq_list, id_list, id2seq), where:
        seq_list : list
            List of sequences (as strings) from fasta file.
        id_list : list
            List of IDs from fasta file.
        id2seq: dict
            Dictionary mapping IDs to sequences (id2seq[id] = seq).
            Only included as an output if include_map=True.

    """
    # read in sequences from the fasta file
    seq_list, id_list, id2seq = list(), list(), dict()
    with open(fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_list.append(str(record.seq))
            id_list.append(record.id)
            id2seq[record.id] = str(record.seq)
    if include_map:
        return seq_list, id_list, id2seq
    return seq_list, id_list


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
