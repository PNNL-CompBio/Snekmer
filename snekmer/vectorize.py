"""vectorize: Create kmer vectors.
author: @christinehc

"""
import itertools
from collections import Counter
from typing import Set, Union

import numpy as np

from snekmer.alphabets import ALPHABET, FULL_ALPHABETS, get_alphabet, get_alphabet_keys


# generate all possible kmer combinations
def _generate(alphabet: Union[str, int], k: int):
    for c in itertools.product(alphabet, repeat=k):
        yield "".join(c)


# iterator object for kmer basis set given alphabet and k
class KmerSet:
    def __init__(self, alphabet: Union[str, int], k: int):
        self.alphabet = alphabet
        self.k = k
        self._kmerlist = list(_generate(get_alphabet_keys(alphabet), k))

    @property
    def kmers(self):
        return iter(self._kmerlist)


# manually reduce alphabet
def reduce(
    string: str, alphabet: Union[str, int], mapping: dict = FULL_ALPHABETS
) -> str:
    """Short summary.

    Parameters
    ----------
    string : str
        Description of parameter `string`.
    alphabet : Union[str, int]
        Description of parameter `alphabet`.
    mapping : dict
        Description of parameter `mapping` (the default is FULL_ALPHABETS).

    Returns
    -------
    str
        Description of returned object.

    Raises
    ------
    ExceptionName
        Why the exception is raised.

    """
    alphabet_map: dict = get_alphabet(alphabet, mapping=mapping)
    return string.translate(string.maketrans(alphabet_map))


class KmerVec:
    def __init__(self, alphabet: Union[str, int], k: int):
        self.alphabet = alphabet
        self.k = k
        self.kmer_gen = None
        self.char_set = get_alphabet_keys(alphabet)
        self.vector = None
        self.kmer_set = KmerSet(alphabet, k)

    # iteratively get all kmers in a string
    def _kmer_gen_str_limit(self, sequence):
        """Generator object for segment in string format"""
        i = 0
        n = len(sequence) - self.k + 1

        # iterate character-by-character
        while i < n:
            kmer = sequence[i : i + self.k]
            if set(kmer) <= self.char_set:
                yield kmer
            i += 1

    # not used
    @staticmethod
    def _kmer_gen_str(sequence, k):
        """Generator object for segment in string format"""
        for n in range(0, len(sequence) - k + 1):
            yield sequence[n : n + k]

    # apply alphabet reduction
    @staticmethod
    def _reduce(sequence: str, alphabet_map: dict) -> str:
        return sequence.translate(sequence.maketrans(alphabet_map))

    # generate kmer vectors with bag-of-words approach
    def vectorize(self, sequence: str) -> np.ndarray:
        N = len(self.char_set) ** self.k

        alphabet_map = get_alphabet(self.alphabet, mapping=FULL_ALPHABETS)
        sequence = self._reduce(sequence, alphabet_map=alphabet_map)
        kmers = list(self._kmer_gen_str_limit(sequence))
        kmer2count = Counter(kmers)

        # Convert to vector of counts
        vector = np.zeros(N)
        for i, word in enumerate(self.kmer_set.kmers):
            vector[i] += kmer2count[word]

        # Convert to frequencies
        # vector /= sum(kmer2count.values())

        return vector
