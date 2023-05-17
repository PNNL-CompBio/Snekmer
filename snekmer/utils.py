"""utils: Snekmer utility functions.

author: @christinehc

"""
# imports
import collections.abc
import datetime
import re

from os.path import basename, splitext
from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray


# functions
def _format_timedelta(timedelta: datetime.timedelta) -> str:
    """Format datetime.timedelta object as str with units.

    Parameters
    ----------
    timedelta : datetime.timedelta object
        A datetime.timedelta object.

    Returns
    -------
    str
        Formatted string with the following format:
            %%h %%m %%.%%s
            (e.g. '01h 41m 20.01s')

    """
    # get milliseconds and round
    milliseconds = round(timedelta.microseconds / 1000)

    # get H, M, S
    minutes, seconds = divmod(timedelta.seconds, 60)
    hours, minutes = divmod(minutes, 60)

    return f"{hours}h {minutes:02}m {seconds:02}.{milliseconds}s"


def log_runtime(
    filename: str, start_time: datetime.datetime, step: Optional[str] = None
) -> None:
    """Append runtimes to specified log file.

    Parameters
    ----------
    filename : str
        /path/to/log/file.
    start_time : datetime.datetime object
        Previous `datetime.now()` result.
    step : str or None (default: None)
        Optional name for step

    Returns
    -------
    None
        Appends runtime information to file and returns no output.

    """
    end_time = datetime.datetime.now()
    if step is None:
        step = "end"

    with open(filename, "a") as f:
        f.write(f"start time:\t{start_time}\n")
        f.write(f"{step} time:\t{end_time}\n")
        f.write(f"total time:\t{_format_timedelta(end_time - start_time)}")
    return


def check_list(array: Any) -> bool:
    """Check whether input is a sequence, list, or array-like object.

    Parameters
    ----------
    array : type
        Input array (or any parameter type).

    Returns
    -------
    bool
        Returns True if input is a list/sequence, array, or Series.

    """
    if not isinstance(array, (collections.abc.Sequence, np.ndarray, pd.Series)):
        return False
    return True


def get_family(
    filename: str, regex: str = r"[a-z]{3}[A-Z]{1}", return_first: bool = True
) -> Union[str, List[str]]:
    """Extract family from filename given regular expression format.

    Note: If no family matching the given regular expression is
    found, the original filename is returned, stripped of any file
    extension (i.e. the file's basename).

    Parameters
    ----------
    filename : str
        path/to/filename.ext
    regex : str or r-string or None
        Regular expression for matching a family name
        (default: "[a-z]{3}[A-Z]{1}").
        The default expression is 3 lowercase letters followed
            by one uppercase letter. To write a custom regular
            expression, see https://docs.python.org/3/library/re.html
            for more details on using the built-in re library.
        If None, returns the full file basename.
    return_first : bool
        If True, returns only the first occurrence (default: True).
        If False, returns a list of all occurring regex matches.

    Returns
    -------
    str or list
        Family name or names


    """
    # extract and simplify file basename
    filename = basename(filename)

    # account for directories
    if "." not in filename:  # and filename[-1] == "/"
        filename = f"{filename}.dir"
    s = "_".join(filename.split(".")[:-1]).replace("-", "_").replace(" ", "_")

    # explicitly define regex as an r-string
    if regex is not None:
        regex = r"{}".format(regex)
    else:
        return filename
    search = re.search(regex, s)

    # return list output
    if not return_first:
        re.findall(regex, s)

    # return string output
    if search is not None:
        return search.group()
    return filename.split(".")[0]


def split_file_ext(filename: str) -> Tuple[str, str]:
    """Split file.ext into (file, ext).

    Ignores ".gz" for gzipped files; e.g. "file.ext.gz" returns
    ("file", "ext") rather than ("file.ext", "gz").

    Parameters
    ----------
    filename : str
        /path/to/file.ext.

    Returns
    -------
    (str, str)
        Tuple containing (file, ext) for file.ext.

    """
    filename = basename(filename)

    # for compressed files, returns (filename, ext) without .gz
    if splitext(filename)[1] == ".gz":
        return (
            splitext(splitext(filename)[0])[0],
            splitext(splitext(filename)[0])[1].lstrip("."),
        )

    # otherwise, returns (filename, ext)
    return splitext(filename)[0], splitext(filename)[1].lstrip(".")


def to_feature_matrix(
    array: List, length_array: Optional[Union[List, NDArray, pd.Series]] = None
) -> NDArray:
    """Create properly shaped feature matrix for kmer scoring.

    Parameters
    ----------
    array : numpy.ndarray or list or array-like
        2-D array-like collection of kmer vectors.
        The assumed format is rows = sequences, cols = kmers.

    Returns
    -------
    numpy.ndarray of numpy.ndarrays
        2D array version of the 2D array-like input.

    """
    if length_array is None:
        length_array = np.ones(len(array))
    array = [np.array(a) / length for a, length in zip(array, length_array)]
    return np.asarray(array)


def count_n_seqs(filename: str) -> int:
    """Count number of sequences in a file.

    Parameters
    ----------
    filename : str
        /path/to/sequence/file.fasta
        Other text formatted sequence files (.faa, etc.) also work.

    Returns
    -------
    int
        Number of sequences contained within the input file.

    """
    return len([1 for line in open(filename) if line.startswith(">")])


def check_n_seqs(filename: str, k: int, show_warning: bool = True) -> bool:
    """Check that a file contains at least k sequences.

    Parameters
    ----------
    filename : str
        /path/to/sequence/file.fasta
        Other text formatted sequence files (.faa, etc.) also work.
    k : int
        Minimum threshold.
    show_warning : bool, optional
        When True, if len(file) < k, a warning is displayed;
        by default True.

    Returns
    -------
    bool
        True if len(file) < k; False otherwise.

    """
    n_seqs = len([1 for line in open(filename) if line.startswith(">")])
    if (n_seqs < k) and (show_warning):
        print(
            f"\nWARNING: {filename} contains an insufficient"
            " number of sequences for model cross-validation"
            " and will thus be excluded from Snekmer modeling."
            f" ({k} folds specified in config; {n_seqs}"
            " sequence(s) detected.)\n"
        )
    return n_seqs >= k
