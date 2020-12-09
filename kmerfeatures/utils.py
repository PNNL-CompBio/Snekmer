"""utils: Utility functions for the kmerfeatures pipeline.

author: @christinehc
"""
# imports
import h5py

from Bio import SeqIO

# global
DEFAULT_HDF = {'map_function': 'reduced_alphabet_0',
               ''}


# functions
def read_fasta(fasta):
    """Read fasta file and obtain sequences, IDs, and mappings.

    Parameters
    ----------
    fasta : str
        Filename or /path/to/fasta/file.

    Returns
    -------
    (list, list)
        Tuple containing (seq_list, id_list, id2seq), where:
        seq_list : list
            List of sequences (as strings) from fasta file.
        id_list : list
            List of IDs from fasta file.

    """
    # read in sequences from the fasta file
    seq_list, id_list, = list(), list()
    with open(fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_list.append(str(record.seq))
            id_list.append(record.id)
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
        Returns empty dictionary if no index file given.

    """
    example_index = {}
    if example_indexfile:
        with open(example_indexfile, "r") as f:
            for line in f.readlines():
                seq_id = line.split()[0]
                example_index[seq_id] = 1.0
    return example_index


def save_hdf(filename, files=None, types=None):
    with h5py.File(filename, 'w') as f:
        if files is None:
            raise AttributeError('No data to save!')

        # save map function
        for key in files.keys():
            # if no dtype, use var length string dtype for arrays
            if key not in types.keys() and isinstance(files[key], list):
                dt = h5py.special_dtype(vlen=str)
            # or use native dtype for all other unspecified dtypes
            elif key not in types.keys():
                dt = type(files[key])
            else:
                dt = types[key]

            # save each item in a dictionary as a separate dataset
            if isinstance(files[key], dict):
                for subkey in files[key].keys():
                    f.create_dataset(f"{key}/{subkey}",
                                     data=files[key][subkey], dtype=dt)
            else:
                f.create_dataset(key, data=files[key], dtype=dt)

        # save ids, seqs, and residues in separate group
        f.create_dataset("filtered/sequences", data=filt_seq_list, dtype=dt)
        f.create_dataset("filtered/ids", data=filt_sid_list, dtype=dt)
        f.create_dataset("filtered/residues", data=filt_residues, dtype=bool)

        for
        f.create_dataset("filtered/example_index", data=filt_example_index)



def load_hdf(filename):
    with h5py.File(filename, 'r') as f:
