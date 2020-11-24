"""transform: Transformation module for Kmer pipeline.

author: @christinehc
"""

# imports
import re

from Util.SIEVEInit import StandardAlphabet
from Util.SIEVEInit import get_alphabets


# global variables and dictionary mappings
RESIDUES = {0: "SV", 1: "APFNDKC", 2: "CAP", 3: "LHC", 4: "LHB"}
MAPFN2RESIDUE = {f"reduced_alphabet_{n}": RESIDUES[n] for n in range(5)}
MAPFN2NAME = {f"reduced_alphabet_{n}": f"RED{n}" for n in range(5)}


# functions
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
        r = n % base
        s = digits[r] + s
        n = n / base

    # fill in any remaining empty slots with first character
    if len(s) < k:
        for i in range(len(s), k):
            s = "%s%s" % (digits[0], s)

    return s


def parse_map_function(map_function, mapping=None):
    """Short summary.

    Parameters
    ----------
    map_function : None or list or str
        Map function for sequences. Can be one of the following:
            None
                No sequence (use generic alphabets)
            list
                Specifications for a random alphabet to use
            str : e.g. "reduced_alphabet_N"
                Use a reduced alphabet (N = 0, 1, 2, 3, or 4)

    mapping : type
        Mapping specification for sequence (?)

    Returns
    -------
    type
        residues, map_name, map_function, kwargs['mapping']

    """
    def reduce_alphabet(character, mapping=None, **kwargs):
        for key in mapping.keys():
            if character in key:
                return mapping[key]

    # specifically for when we create a random alphabet
    #   that we want to apply to a lot of sequences
    if isinstance(map_function, list):
        mapping = map_function[2]
        residues = map_function[0]
        map_name = map_function[1]
        map_function = reduce_alphabet

    elif isinstance(map_function, str):
        try:
            mapfn = map_function
            i_mapfn = int(re.search("(?<=reduced_alphabet_)[0-4]",
                                    mapfn).group())
            map_function = reduce_alphabet
            mapping = get_alphabets()[mapfn]
            residues = RESIDUES[i_mapfn]
            map_name = MAPFN2NAME[i_mapfn]
        except AttributeError as e:
            raise ValueError(
                ('map_function string must be in this format:'
                 ' "reduced_alphabet_n", with n = 0,1,2,3,4')
                ) from e

    return residues, map_name, map_function, mapping


def string_vectorize(sequence=None, kmer=3, start=None,
                     end=None, map_function=None, return_labels=False,
                     feature_dict=None, filter_list=None, exclusion_list=None,
                     return_dict=None, kmer_output=None,
                     residues=StandardAlphabet, **kwargs):
    """Short summary.

    Parameters
    ----------
    sequence : type
        Description of parameter `sequence`.
    kmer : int
        Sequence length k of the kmer.
    start : type
        Description of parameter `start`.
    end : type
        Description of parameter `end`.
    map_function : type
        Description of parameter `map_function`.
    return_labels : type
        Description of parameter `return_labels`.
    feature_dict : type
        Description of parameter `feature_dict`.
    filter_list : type
        Description of parameter `filter_list`.
    exclusion_list : type
        Description of parameter `exclusion_list`.
    return_dict : type
        Description of parameter `return_dict`.
    kmer_output : type
        Description of parameter `kmer_output`.
    residues : str
        (default: "AILMVFYWSTQNCHDEKRGP")

    Returns
    -------
    type
        why is this returning so many things ????

    """
    # appears extraneous
    # def identity(character, **kwargs):
    #     return character
    # if residues or map_function not specified, set generically
    # map_function = map_function #  or identity
    # residues = residues or StandardAlphabet
    mapping = None
    map_name = "NAT"

    residues, map_name, map_function, mapping = parse_map_function(map_function)

    # we should provide the ability to include a bunch of strings
    # that can be used to vectorize
    error_message = "Error: given parameters will generate >20k inputs."
    assert pow(len(residues), kmer) <= 2e4, error_message

    # return only labels
    if return_labels:
        labels = []
        if filter_list:
            for filter in filter_list:
                label = "KMER-%d-%s-%s" % (kmer, map_name, filter)
                labels.append(label)
            return labels

        else:
            for bit in range(pow(len(residues), kmer)):
                label = "KMER-%d-%s-%s" % (kmer, map_name, baseconvert(bit, k=kmer, digits=residues))
                labels.append(label)
            return labels

    tstart = start or 0
    if tstart < 0:
        tstart = len(sequence) + tstart
        if tstart < 0:
            tstart = 0

    tend = end or len(sequence)-kmer

    results = {}
    if feature_dict:
        results = feature_dict

    elif filter_list:

        # there is a more elegant way of doing this- and probably clever too
        for item in filter_list:
            results[item] = 0

    for i in range(tstart, tend):
        kmap = sequence[i:i+kmer]

        # do mapping to a reduced alphabet, e.g.
        kstring = ""
        for char in kmap:
            cc = map_function(char, **kwargs)

            # this has the effect of omitting unrecognized characters, which may be undesireable in some cases
            if cc == None:
                continue
            kstring += cc

        if len(kstring) < kmer:
            # this happens when there are unrecognized characters
            # and we need to not include these
            continue

        # if we don't find the kstring in the filter_list
        #   then we skip.
        #if filter_dict:
        #    print(filter_dict.keys())
        #    print(kstring)

        if kmer_output:
          print(i, "\t", kmap, "\t", kstring, "\t1", )

        if filter_list and not kstring in filter_list:
            #print("hoopla")
            continue

        #FILTER HERE

        # a list of strings to exclude from consideration
        #   this is from the kmer_walk approach and should
        #   be shorter sequences (though could also be the
        #   same length)
        if exclusion_list:
            breaker = 0
            for bit in exclusion_list:
                if kstring.find(bit) > 0:
                    breaker = 1
                    break
            if breaker: next

        if not kstring in results:
            results[kstring] = 0
        results[kstring] += 1

    if return_dict:
        return results

    if filter_list:
        results_ordered = []
        for item in filter_list:
            results_ordered.append(results[item])
        return results_ordered

    return results.values()

    # old stuff below
    results_vector = []
    for j in range(0, pow(len(residues), kmer)):
        results_vector.append(0.0)

    for key in results.keys():
        # base conversion- though we're using base 20 (A-Y, minus B, J, O, U) and not base 25 (A-Y)
        # for the standard aa alphabet. We can also use various reduced alphabets
        x = 0

        if len(key) < kmer:
            continue

        for k in range(kmer):
            # provide a unique index for all possible strings
            x += residues.find(key[k]) * pow(len(residues), k)
        results_vector[x] = results[key]

    return results_vector


def scramble_sequence(sequence_id, sequence, n=1, id_modifier=False,
                      first_residue=1, example_index=None):
    """Return scrambled sequences given a sequence and identifier.

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

    start_residue = sequence[0:first_residue]
    scramble = [char for char in sequence[first_residue:]]

    seq_list = []
    for i in range(len(id_list)):
        random.shuffle(scramble)
        shuffled = [char for char in "".join([start_residue] + scramble)]
        seq_list.append(shuffled)

        example_index[sid] = -1  # this confuses me-- is this mistakenly written?

    return id_list, seq_list, example_index


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
            id_string = "%s-%s" % (record.id, sequence_id)
            id_list.append(id_string)
            sequence_list.append(str(record.seq))

    return id_list, sequence_list

def randomize_alphabet(map_function):
    """Short summary.

    Parameters
    ----------
    map_function : type
        Description of parameter `map_function`.

    Returns
    -------
    type
        Description of returned object.

    """
    alpha = get_alphabets()[map_function]

    residues = alpha["_keys"]
    map_name = f"RND{map_function[-1]}"

    rand_alphabet = {}

    # this gives us a non-repeating string for use as new keys
    rand_str = ''.join(random.sample("ACDEFGHIKLMNPQRSTVWY", 20))
    for key_a in alpha.keys():
        if key_a == "_key":
            rand_alphabet["_key"] = alpha["_key"]
            continue
        key_r = rand_str[0:len(key_a)]

        # trim off this sequence
        rand_str = rand_str[len(key_a):]
        rand_alphabet[key_r] = alpha[key_a]  # ?? what about rand_str

    return rand_alphabet
