"""transform: Transformation module for Kmer pipeline.

author: @christinehc
"""

# imports
import random
import re
from multiprocessing import Pool

from Bio import SeqIO
from .alphabet import (ALPHABETS, ALPHABET_ORDER, StandardAlphabet, get_alphabets)


# functions
def identity(character, mapping=None):
    """Return self.

    Used as a placeholder map function when applicable.

    Parameters
    ----------
    character : object
        Object to return.
    mapping : NoneType
        For consistency with reduce_alphabet

    Returns
    -------
    object
        Returned object.

    """
    return character


def baseconvert(n, k, digits="ACDEFGHIKLMNPQRSTVWY"):
    """Generate a kmer sequence from input integer representation.

    Parameters
    ----------
    n : int
        Integer representation of a kmer sequence.
    k : int
        Length of the kmer (i.e. protein sequence length).
    digits : str
        Digits to use in determining the base for conversion.
            (default: "ACDEFGHIKLMNPQRSTVWY")

    Returns
    -------
    str
        Kmer sequence, in string format.
        Returns empty string if inputs are invalid.

    """
    assert len(digits) > 0, "digits argument cannot be empty string."
    # digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    base = len(digits)

    try:
        n = int(n)
    except (TypeError, ValueError):
        return ""

    if n < 0 or base < 2 or base > 36:
        return ""

    # parse integer by digits base to populate sequence backwards
    s = ""
    while n != 0:
        r = int(n % base)
        s = digits[r] + s
        n = n / base

    # fill in any remaining empty slots with first character
    if len(s) < k:
        for i in range(len(s), k):
            s = "%s%s" % (digits[0], s)

    return s


def reduce_alphabet(character, mapping):
    """Reduce alphabet according to pre-defined alphabet mapping.

    Parameters
    ----------
    character : type
        Description of parameter `character`.
    mapping : str
        Name of map function.

    Returns
    -------
    type
        Description of returned object.

    """
    for key in mapping.keys():
        if character in key:
            return mapping[key]
    return None


def parse_mapping(alphabet):
    """Parse alphabet input into mapping parameters.

    Mappings are created as follows:
        If the input map function is a list, this format is assumed:
            (residues, map_name, mapping)
        If the input map function is a string, the mapping is defined
        from pre-existing alphabets.

    Parameters
    ----------
    alphabet : None or list or str
        Alphabet for sequences. Can be one of the following:
            None
                No sequence (use generic alphabets)
            list
                (residues, map_name, mapping)
                    str, str, dict
                Specifications for a random alphabet to use
            str : e.g. "reduced_alphabet_N"
                Use a reduced alphabet (N = 0, 1, 2, 3, 4, or 5)


    Returns
    -------
    type
        residues, map_name, alphabet

    """
    success = 0

    # do not reduce alphabet
    if (alphabet == "None") or (alphabet is None):
        mapping = None
        residues = StandardAlphabet
        map_name = str(alphabet)
        alphabet_map = identity
        success += 1

    # use pre-existing alphabets
    if (alphabet in ALPHABETS) or (alphabet in range(len(ALPHABETS))):
        if (isinstance(alphabet, int)):
            alphabet = ALPHABET_ORDER[alphabet]
        try:
            mapping = get_alphabets()[alphabet]
            residues = get_alphabets()[alphabet]['_keys']
            map_name = alphabet.upper()
            alphabet_map = reduce_alphabet
            success += 1
        except AttributeError as e:
            raise ValueError(
                "alphabet string must be one of the following:"
                f" {list(ALPHABETS.keys())}."
                ) from e

    # use a random alphabet to apply to many sequences
    elif isinstance(alphabet, list):
        residues = alphabet[0]
        map_name = alphabet[1]
        mapping = alphabet[2]
        alphabet_map = reduce_alphabet
        success += 1

    # if all else fails
    if success == 0:
        print("alphabet_map:\t", alphabet_map)
        raise ValueError("Invalid map function specified.")

    return residues, map_name, alphabet_map, mapping


def map_characters(k_map, alphabet_map, mapping):
    """Apply character mapping as specified.

    Parameters
    ----------
    k_map : str
        String of mapped characters.
    alphabet_map : type
        Description of parameter `alphabet_map`.
    mapping : type
        Description of parameter `mapping`.

    Returns
    -------
    type
        Description of returned object.

    """
    k_string = ""
    for char in k_map:
        map_char = alphabet_map(char, mapping)

        # omits unrecognized characters (may be undesireable in some cases)
        if map_char is None:
            continue
        k_string += map_char
    return k_string


def generate_labels(k, alphabet, filter_list=None):
    """Generate labels for mapped residues.

    Parameters
    ----------
    k : int
        Sequence length k of the kmer.
    alphabet : str
        Name of alphabet (default: None). If None, no mapping is
        applied.
    residues : str or iterable
        Residues alphabet; the default standard alphabet is below.
        (default: "AILMVFYWSTQNCHDEKRGP")
    filter_list : list
        List of filter sequences (default: None).

    Returns
    -------
    list
        List of labels in the format 'KMER-{k}-{map_name}-{mapped residues}'.

    """
    # if residues or alphabet not specified, set generically
    # alphabet = alphabet or identity

    residues, map_name, _, _ = parse_mapping(alphabet)

    # we should provide the ability to include a bunch of strings
    # that can be used to vectorize
    # if pow(len(residues), k) <= 2e4:
    #     raise RuntimeError(
    #         "Given parameters will generate >20k inputs."
    #     )

    labels = []

    # if there is a filter list, return labels for listed filters
    if filter_list:
        for filt in filter_list:
            label = "KMER-%d-%s-%s" % (k, map_name, filt)
            labels.append(label)
        return labels

    for bit in range(pow(len(residues), k)):
        label = "KMER-%d-%s-%s" % (k, map_name, baseconvert(bit, k,
                                                            digits=residues))
        labels.append(label)
    return labels


def set_sequence_endpoints(sequence, k, start, end):
    """Set logical start and end indices for a given sequence.

    Parameters
    ----------
    sequence : str or iterable
        Sequence for which to determine endpoints.
    k : int
        Kmer length (k units).
    start : int
        Start index of the sequence.
    end : int
        End index of the sequence.

    Returns
    -------
    type
        Description of returned object.

    """
    start = start or 0
    end = end or len(sequence) - k
    # set logical start and end indices for sequence
    if start < 0:
        start = len(sequence) + start
        if start < 0:
            start = 0
    return start, end


def exclude_from_string(k_string, exclude):
    """Exclude string if at least one given substring is detected.

    Parameters
    ----------
    k_string : str
        String to check.
    exclude : list
        List of substrings to search for in k_string.

    Returns
    -------
    bool
        Returns True if a specified substring is found in k_string.
        Returns False if none of the specified substrings are found
            in k_string, or if no exclusion list is given.

    """
    if exclude:
        if not isinstance(exclude, list):
            exclude = [exclude]
        for bit in exclude:
            if k_string.find(bit) > 0:
                return True
    return False


def vectorize_string(sequence, k, alphabet, start=0, end=False,
                     feature_dict=None, filter_list=None, exclude=None,
                     return_dict=False, verbose=False, log_file=False):
    """Transform a protein sequence into a vector representation.

    Parameters
    ----------
    sequence : str
        Protein sequence.
    k : int
        Sequence length k of the kmer.
    start : int
        Start index of the sequence (for sequence slicing).
    end : int
        End index of the sequence (for sequence slicing).
    alphabet : str
        Name of the map function (e.g. "reduced_alphabet_0")
        (default: "None").
    feature_dict : dict
        Dictionary of {feature ID: feature} (default: None).
    filter_list : list
        List of filter sequences (default: None).
    exclude : list
        List of sequence strings for exclusion (default: None).
        This is from the kmer_walk approach and should be shorter
        sequences (though they may also be the same length).
    return_dict : bool
        If True, returns {kmer: count} for kmers (default: False).
    verbose : bool
        If True, prints verbose output during run (default: False).
    log_file : str
        /path/to/log.file for verbose outputs (default: False)
        If False, pipes verbose outputs to console instead.

    Returns
    -------
    list
        List containing kmer counts.

    """
    if alphabet == "None":
        alphabet = None

    # if residues or alphabet not specified, set generically
    # alphabet = alphabet or identity
    _, _, alphabet_map, mapping = parse_mapping(alphabet)

    # we should provide the ability to include a bunch of strings
    # that can be used to vectorize
    # if pow(len(residues), k) <= 2e4:
    #     raise RuntimeError(
    #         "Given parameters will generate >20k inputs."
    #     )

    # set logical start and end indices for sequence
    start, end = set_sequence_endpoints(sequence, k, start, end)

    results = feature_dict or {}
    if filter_list:
        results = {key: 0 for key in filter_list}

    for i in range(start, end):
        k_map = sequence[i: i + k]

        # perform mapping to a reduced alphabet
        k_string = map_characters(k_map, alphabet_map, mapping)

        if verbose and log_file:
            with open(log_file, 'a') as f:
                f.write(f"{i}\t{k_map}\t{k_string}\t1\n")
        elif verbose and not log_file:
            print(f"{i}\t{k_map}\t{k_string}\t1", )

        # filter unrecognized characters or filter from list
        if (len(k_string) < k) or (filter_list and k_string
                                   not in filter_list):
            continue

        # FILTER HERE

        # exclude specified substrings from results
        if exclude_from_string(k_string, exclude):
            continue

        # initialize value in dictionary if missing
        if k_string not in results:
            results[k_string] = 0
        results[k_string] += 1

    if filter_list:
        results = {item: results[item] for item in filter_list}
    if return_dict:
        return results
    return list(results.values())


def scramble_sequence(sequence_id, sequence, n=1, id_modifier=False,
                      first_residue=1, example_index=None):
    """Scramble sequences given an identifier and its sequence.

    Given a sequence and an identifier, returns a list of n scrambled
    sequences and numbered identifiers.

    Parameters
    ----------
    sequence_id : type
        Identifier.
    sequence : str
        Sequence to be scrambled.
    n : int
        Number of scrambled sequences to produce (default: 1).
    id_modifier : bool
        If True, adds modifier to sequence ID to specify shuffle
        parameters (detault: False).
    first_residue : int
        Represents the first index at which to scramble the sequence.
        e.g. when first_residue=1, protects the first residue from
        being scrambled, i.e. protecting an N-terminal M (default: 1).
        When first_residue=0, scrambles the full sequence.
    example_index : type
        (default: None)

    Returns
    -------
    type
        Description of returned object.

    """
    if id_modifier:
        id_list = [f"{sequence_id}_shuffle_{i}" for i in range(n)]
    else:
        id_list = [sequence_id for i in range(n)]

    start_residue = sequence[:first_residue]
    seq = [char for char in sequence[first_residue:]]

    scrambled_seqs = []
    for sid in id_list:
        random.shuffle(seq)
        shuffled = [char for char in "".join([start_residue] + seq)]
        scrambled_seqs.append(shuffled)

        # change indices to -1 for shuffled
        example_index[sid] = -1.0  # this confuses me-- mistake?

    return id_list, scrambled_seqs, example_index


def make_n_terminal_fusions(sequence_id, filename):
    """Fuse sequence with sequences in fasta file.

    Parameters
    ----------
    sequence_id : type
        Description of parameter `sequence_id`.
    filename : str
        Filename (e.g. "/path/to/file.fasta").

    Returns
    -------
    (list, list)
        Tuple containing (id_list, sequence_list), where:
            id_list : List of fused sequence IDS
            sequence_list : List of sequences in fasta file

    """
    sequence_list, id_list = [], []
    with open(filename, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            id_string = f"{record.id}-{sequence_id}"
            id_list.append(id_string)
            sequence_list.append(str(record.seq))

    return id_list, sequence_list


def randomize_alphabet(alphabet):
    """Create a randomized residue alphabet.

    Parameters
    ----------
    alphabet : str
        Name of the map function (e.g. "reduced_alphabet_0").

    Returns
    -------
    dict
        Dictionary containing randomized residue mappings.

    """
    alpha = get_alphabets()[alphabet]

    rand_alphabet = {}
    # this gives us a non-repeating string for use as new keys
    rand_str = ''.join(random.sample("ACDEFGHIKLMNPQRSTVWY", 20))
    for key_a in alpha.keys():
        if key_a == "_key":
            rand_alphabet["_key"] = alpha["_key"]
            continue
        key_r = rand_str[0:len(key_a)]

        # trim off used sequence bits
        rand_str = rand_str[len(key_a):]
        rand_alphabet[key_r] = alpha[key_a]

    return rand_alphabet
