"""transform: Transformation module for Kmer pipeline."""

# imports
import re

from Utils.SIEVEInit import StandardAlphabet, get_alphabets

# global variables and dictionary mappings
RESIDUES = {0: "SV", 1: "APFNDKC", 2: "CAP", 3: "LHC", 4: "LHB"}
MAPFN2RESIDUE = {f"reduced_alphabet_{n}": RESIDUES[n] for n in range(5)}
MAPFN2NAME = {f"reduced_alphabet_{n}": f"RED{n}" for n in range(5)}


# functions
def baseconvert(n, k=None, digits=None, **kwargs):
    """Convert positive decimal integer n to equivalent in another.

    Parameters
    ----------
    n : type
        Description of parameter `n`.
    k : type
        Description of parameter `k`.
    digits : type
        Description of parameter `digits`.
    **kwargs : type
        Optional keyword arguments.

    Returns
    -------
    type
        Description of returned object.

    """
    #digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    digits = digits or "ACDEFGHIKLMNPQRSTVWY"
    base = len(digits)

    try:
        n = int(n)
        base = int(base)
    except (TypeError, ValueError, SyntaxError):
        return ""

    if n < 0 or base < 2 or base > 36:
        return ""

    s = ""

    while 1:
        r = n % base

        s = digits[r] + s
        n = n / base
        n = int(n)
        if n == 0:
            break

    if len(s) < k:
        for i in range(len(s), k):
            s = "%s%s" % (digits[0], s)

    return s


def _parse_map_function(map_function, **kwargs):
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

    **kwargs : type
        Optional keyword arguments for specification.
            e.g. mapping

    Returns
    -------
    type
        residues, map_name, map_function, kwargs['mapping']

    """
    def reduce_alphabet(character, mapping=None, **kw):
        for key in mapping.keys():
            if character in key:
                return mapping[key]

    # specifically for when we create a random alphabet
    #   that we want to apply to a lot of sequences
    if isinstance(map_function, list):
        kwargs["mapping"] = map_function[2]
        residues = map_function[0]
        map_name = map_function[1]
        map_function = reduce_alphabet

    elif isinstance(map_function, str):
        try:
            mapfn = map_function
            i_mapfn = int(re.search("(?<=reduced_alphabet_)[0-4]",
                                    mapfn).group())
            map_function = reduce_alphabet
            kwargs['mapping'] = get_alphabets()[mapfn]
            residues = RESIDUES[i_mapfn]
            map_name = MAPFN2NAME[i_mapfn]
        except AttributeError as e:
            raise ValueError(
                ('map_function string must be in this format:'
                 ' "reduced_alphabet_n", with n = 0,1,2,3,4')
                ) from e

    return residues, map_name, map_function, kwargs['mapping']


def string_vectorize(residues=None, sequence=None, kmer=3, start=None,
                     end=None, map_function=None, return_labels=False,
                     feature_dict=None, filter_list=None, exclusion_list=None,
                     return_dict=None, kmer_output=None, **kw):
    """Short summary.

    Parameters
    ----------
    residues : type
        Description of parameter `residues`.
    sequence : type
        Description of parameter `sequence`.
    kmer : type
        Description of parameter `kmer`.
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
    **kw : type
        Description of parameter `**kw`.

    Returns
    -------
    type
        Description of returned object.

    """
    def identity(character, **kw):
        return character

    # if residues or map_function not specified, set generically
    residues = residues or StandardAlphabet
    map_function = map_function or identity
    mapping = None
    map_name = "NAT"

    residues, map_name, map_function, kwargs['mapping'] = _parse_map_function(map_function)

    # we should provide the ability to include a bunch of strings
    # that can be used to vectorize
    # if pow(len(residues), kmer) > 20000:
    #    sys.stderr.write("Error: this will generate more than 20k inputs\n")
    #    return

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
        for c in kmap:
            cc = map_function(c, **kw)

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
