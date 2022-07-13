"""features: Output functions for kmer features.

author: @christinehc

"""
# imports
from typing import Any, Dict, IO, List, Optional


# functions
def output_features(
    save_to: str,
    fformat: str,
    feature_sets: Optional[List[List[Any]]] = None,
    labels: Optional[List] = None,
    mode: str = "w",
    **kwargs,
) -> None:
    """Generate output features based on input fasta file.

    Parameters
    ----------
    save_to : str
        File path to output file (default: None).
    fformat : str
        File format; select one option from the following list:
            ('gist', 'sieve', 'matrix', 'both')
    feature_sets : type
        Description of parameter `feature_sets`.
    labels : list
        Description of parameter `labels`.
    mode : str
        File write mode (default: 'w').
            For details, see documentation for Python built-in
            function open().

    Returns
    -------
    None
        Feature output is written to the `save_to` file path.
        No output is returned.

    """
    # # ensure valid format is specified
    # if fformat not in ['gist', 'sieve', 'matrix', 'both']:
    #     raise ValueError(
    #         ('File format must be one of the following:'
    #          ' "gist", "sieve", "both", or "matrix".')
    #         )

    # update any gist files
    if fformat in ["gist", "both"]:
        _output_gist(
            save_to, labels=labels, feature_sets=feature_sets, mode=mode, **kwargs
        )

    # update any sieve files
    if fformat in ["sieve", "both"]:
        _output_sieve(save_to, feature_sets=feature_sets, mode=mode, **kwargs)

    # update matrix files
    if fformat == "matrix":
        _output_matrix(save_to, labels=labels, feature_sets=feature_sets, mode=mode)


def _output_gist(
    filename: str,
    labels: List[str] = None,
    feature_sets: Optional[List[str]] = None,
    mode: str = "w",
) -> None:
    """Output features in gist format.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    labels : list
        (default: None)
    feature_sets : type
        Description of parameter `feature_sets`.
    mode : str
        (default: 'w')

    Returns
    -------
    None
        Feature output is written to files; no output is returned.

    """
    train_out = f"{filename.split('.')[0]}.train"
    class_out = f"{filename.split('.')[0]}.class"

    # update train
    with open(train_out, mode) as tf, open(class_out, mode) as cf:
        tf.write("corner")
        cf.write("corner\tclass\n")

        if labels:
            for label in labels:
                tf.write("\t%s" % label)
        else:
            for i in range(len(feature_sets[0]) - 1):
                tf.write("\tlabel%d" % i)
        for features in feature_sets:
            _output_gist_features(tf, features, mode)
            _output_gist_class(tf, features, mode)


def _output_sieve(
    filename: str,
    feature_sets: List[str],
    labels: Optional[List[str]] = None,
    mode: str = "w",
    **kwargs,
) -> None:
    """Output features in sieve format.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    feature_sets : type
        Description of parameter `feature_sets`.
    labels : list
        (default: None)
    mode : str
        (default: 'w')

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """

    def output_sieve_features(
        features: List[str],
        filename: str,
        example_index: Optional[Dict[str, float]] = None,
    ) -> None:
        """Write features to SIEVE output file.

        Parameters
        ----------
        features : type
            List of specified features.
        filename : str
            File handle for naming of output files.
        example_index : dict
            Description of parameter `example_index` (default: None).

        Returns
        -------
        None
            Feature output is written to file; no output is returned.

        """
        # parse first item in feature list as feature ID
        fid = features[0]
        if example_index is not None:
            value = example_index.get(fid, 0.0)
        else:
            value = 0.0

        with open(filename, "a") as f:
            f.write("pattern\t%s\t%d\n" % (fid, len(features) - 1))
            f.write("\tinput\t%s" % fid)
            for ft in features[1:]:
                f.write("\t%s" % ft)
            f.write("\n")
            f.write("\toutput\t%s\t%d\n" % (fid, value))
            f.flush()

    # pattern_out = f"{filename}.pattern"
    for features in feature_sets:
        output_sieve_features(features, filename, **kwargs)


def _output_matrix(
    filename: str,
    labels: Optional[List[str]] = None,
    feature_sets: List = None,
    mode: str = "w",
) -> None:
    """Output features in matrix format.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    labels : list
        (default: None)
    feature_sets : type
        Description of parameter `feature_sets`.
    mode : str
        (default: 'w')

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """
    with open(filename, mode) as f:
        if labels is not None:
            f.write("%s" % labels[0])
            for label in labels[1:]:
                f.write("\t%s" % label)
            f.write("\n")
            f.flush()

    if feature_sets is not None:
        for features in feature_sets:
            _output_gist_features(filename, features, mode)


def _output_gist_features(filename: str, features: List[str], mode: str = "w") -> None:
    """Write features to gist output file.

    Parameters
    ----------
    filename : type
        Description of parameter `file`.
    features : type
        Description of parameter `features`.

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """
    with open(filename, mode) as f:
        f.write("%s" % features[0])
        for ft in features[1:]:
            f.write("\t%s" % ft)
        f.write("\n")
    # filename.flush()


def _output_gist_class(
    filename: str,
    features: List[str],
    example_index: Optional[Dict[str, float]] = None,
    mode: str = "w",
) -> None:
    """Write gist class to specified output file.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    features : list
        List of specified features.
    example_index : dict
        Description of parameter `example_index` (default: None).

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """
    pattern = f"{filename.split('.')[0]}.pattern"
    with open(pattern, mode) as f:
        fid = features[0]
        if example_index is not None:
            value = example_index.get(fid, -1)
        else:
            value = -1
        f.write("%s\t%d\n" % (features[0], value))
        # filename.flush()

