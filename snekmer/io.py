"""io: Input-output handling in Snekmer.

author: @christinehc

"""
# imports
import gzip
import json
import pickle
import re
from os.path import basename, join, splitext
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from .alphabet import ALPHABET_ORDER


# functions
def load_pickle(filename: str, mode: str = "rb") -> Any:
    """Load a pickled object (wrapper for `pickle.load`).

    Parameters
    ----------
    filename : str
        Description of parameter `filename`.
    mode : str
        Description of parameter `mode` (the default is "rb").

    Returns
    -------
    Any
        Description of returned object.

    Raises
    ------
    ExceptionName
        Why the exception is raised.

    """
    with open(filename, mode) as f:
        return pickle.load(f)


def load_npz(
    filename: str,
    columns: Dict[str, str] = {
        "ids": "sequence_id",
        "seqs": "sequence",
        "vecs": "sequence_vector",
    },
) -> pd.DataFrame:
    """Compile .npz results into dataframe.

    Parameters
    ----------
    filename : str
        /path/to/filename for input .npz file.
    columns : Dict[str, str]
        Mapping for output data column names (the default is
            {
                "ids": "sequence_id",
                "seqs": "sequence",
                "vecs": "sequence_vector",
            }
        ).

    Returns
    -------
    pd.DataFrame
        Tabulated .npz data.

    """
    data = np.load(filename)

    # fill in df based on desired output col names
    df = {"filename": splitext(basename(filename))[0]}

    for in_col, out_col in columns.items():
        df.update({out_col: list(data[in_col])})

        # get seq column for sequence lengths
        if "seq" in in_col:
            df.update({f"{out_col}_length": [len(s) for s in data[in_col]]})

    return pd.DataFrame(df)


def read_kmers(filename: str) -> List[str]:
    """Extract kmer identifiers from Snekmer output file.

    Parameters
    ----------
    filename : str
        /path/to/output/file.txt

    Returns
    -------
    list
        Array of strings (kmer identifiers)

    """
    kmers = list()
    with open(filename) as f:
        for line in f:
            kmers.append(line.strip())  # .split("\t")
            # parse kmer outputs if detected
            # if re.findall(r"^KMER", line_data[0]):
            # kmers += line_data
    # prefix = re.search(r"^KMER-[\d]+-[A-Za-z]+-", kmers[0]).group()
    # return [s.strip().replace(prefix, "") for s in kmers]
    return kmers


def vecfiles_to_df(
    files: List[str], labels: Optional[List[str]] = None, label_name: str = "label"
) -> pd.DataFrame:
    """Load multiple sequence files and parse into common dataframe.

    Parameters
    ----------
    files : list of str
        List of (/path/to/file1.ext1, /path/to/file2.ext2, ..., etc.)
    labels : list of str or None (default: None)
        Labels for specified files; if None, excludes labels.
    label_name : str (default: 'label')
        Desired name for data column containing labels.

    Returns
    -------
    pandas.DataFrame
        Tabulated data of the form:

        | filename | seq_id | vector | vec_shape |
        |----------|--------|--------|-----------|
        | file1    | seq1   | array1 | shape1    |
        ...
        ... etc.

    """
    data = {"filename": [], "seq_id": [], "vector": [], "vec_shape": []}
    for afile in files:
        with gzip.open(afile, "rt") as f:
            # with open(afile, 'r') as f:
            tmp = json.load(f)
            data["filename"] += [basename(afile)] * len(tmp["seq_id"])
            data["seq_id"] += tmp["seq_id"]
            data["vector"] += tmp["vector"]
            data["vec_shape"] += [np.array(arr).shape for arr in tmp["vector"]]
    if str(labels) != "None":
        if len(labels) != len(data):
            raise ValueError("Number of labels must equal number of sequences.")
        data[label_name] = labels
    return pd.DataFrame(data)


def define_output_dir(alphabet: Union[str, int], k: int, nested: bool = False) -> str:
    """Create output directory name using AAR parameters.

    Parameters
    ----------
    alphabet : str or int
        Name or numerical identifier of alphabet.
        See `snekmer.alphabet.ALPHABETS` for descriptions.
    k : int
        Kmer length.
    nested : bool (default: False)
        If True, creates nested directories with alphabet, k.
        If False, returns simple output directory.

    Returns
    -------
    str
        Name of output directory, given nested directory parameters.

    """
    if not nested:
        return "output"
    if not isinstance(alphabet, str):
        alphabet = ALPHABET_ORDER[alphabet]
    return join("output", alphabet, f"k-{k:02}")
