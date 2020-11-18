"""features: Output feature processing for the Kmer pipeline."""

# imports


# functions
def output_features(format, feature_sets=None, output_filename=None, labels=None, mode="w", **kw):
    """Generate output features based on input fasta file.

    Generate

    Parameters
    ----------
    format : str
        File format; select one option from the following list:
            ('gist', 'sieve', 'matrix', 'both')
    feature_sets : type
        Description of parameter `feature_sets`.
    output_filename : type
        Description of parameter `output_filename`.
    labels : list
        Description of parameter `labels`.
    mode : str
        File write mode (default: 'w').
            For details, see documentation for Python built-in
            function open().
    **kw : type
        Description of parameter `**kw`.

    Returns
    -------
    type
        Description of returned object.

    """
    # ensure valid format is specified
    if format not in ['gist', 'sieve', 'matrix', 'both']:
        raise ValueError(
            ('File format must be one of the following:'
             ' "gist", "sieve", "both", or "matrix".')
            )

    # update any gist files
    if format in ["gist", "both"]:
        train_out = f"{output_filename}.train"
        class_out = f"%{output_filename}.class"

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
                _output_gist_features(features=features, file=tf, **kw)
                _output_gist_class(features=features, file=cf, **kw)

    # update any sieve files
    if format in ["sieve", "both"]:
        pattern_out = f"{output_filename}.pattern"
        with open(pattern_out, mode) as pf:
            for features in feature_sets:
                _output_sieve_features(features=features, file=pf, **kw)

    # update matrix files
    if format == "matrix":
        file_out = f"{output_filename}.txt"
        with open(file_out, mode) as of:
            if labels:
                of.write("%s" % labels[0])
                for label in labels[1:]:
                    of.write("\t%s" % label)
                of.write("\n")
            if feature_sets:
                for features in feature_sets:
                    _output_gist_features(features=features, file=of, **kw)


def _output_gist_features(features=None, file=None, **kw):
    """Write features to gist output file.

    Parameters
    ----------
    features : type
        Description of parameter `features`.
    file : type
        Description of parameter `filehandle`.
    **kw : type
        Description of parameter `**kw`.

    Returns
    -------
    type
        Description of returned object.

    """
    file.write("%s" % features[0])
    for f in features[1:]:
        file.write("\t%s" % f)
    file.write("\n")
    file.flush()


def _output_sieve_features(features=None, file=None, example_index={}, **kw):
    """Write features to SIEVE output file.

    Parameters
    ----------
    features : type
        List of specified features.
    file : str
        File handle for naming of output files.
    example_index : dict
        Description of parameter `example_index`.
    **kw : type
        Description of parameter `**kw`.

    Returns
    -------
    type
        Description of returned object.

    """
    # parse first item in feature list as feature ID
    fid = features[0]
    value = example_index.get(fid, 0.0)
    file.write("pattern\t%s\t%d\n" % (fid, len(features) - 1))
    file.write("\tinput\t%s" % fid)
    for f in features[1:]:
        file.write("\t%s" % f)
    file.write("\n")
    file.write("\toutput\t%s\t%d\n" % (fid, value))
    file.flush()


def _output_gist_class(features=None, file=None, example_index={}, **kw):
    """Write gist class to specified output file.

    Parameters
    ----------
    features : list
        List of specified features.
    file : str
        File handle for naming of output files.
    example_index : dict
        Description of parameter `example_index`.
    **kw : type
        Description of parameter `**kw`.

    Returns
    -------
    type
        Description of returned object.

    """
    fid = features[0]
    value = example_index.get(fid, -1)
    file.write("%s\t%d\n" % (features[0], value))
    file.flush()
