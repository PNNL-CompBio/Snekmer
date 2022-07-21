"""vectorize: Create kmer vectors.

author: @christinehc

"""
import itertools
from collections import Counter
from typing import Dict, Generator, Set, Union

import numpy as np
from numpy.typing import NDArray
from .alphabet import FULL_ALPHABETS, get_alphabet, get_alphabet_keys
from .utils import check_list

# store kmer basis set and transform new vectors into fitted basis
class KmerBasis:
    """Store kmer basis set and perform basis set transforms.

    Attributes
    ----------
    basis : list
        Description of attribute `basis`.
    basis_order : type
        Description of attribute `basis_order`.

    """

    def __init__(self):
        self.basis = []
        self.basis_order = {}

    def set_basis(self, basis):
        """Specify kmer basis set.

        Parameters
        ----------
        basis : list or array-like of str
            Ordered array of kmer strings.

        Returns
        -------
        self: object
            Fitted KmerBasis object.

        """
        if not check_list(basis):
            raise TypeError("`basis` input must be list or array-like.")

        self.basis = basis
        self.basis_order = {i: k for i, k in enumerate(basis)}

    def transform(self, vector, vector_basis):
        """Apply basis to new vector with separate basis.

        e.g. If the basis is set to a list of p kmers, the input
            vector array of size (m, n) -> (m, p).

        Note: Order is preserved from kmer arrays. Kmer count
        vectors are assumed to follow the order from corresponding
        kmer basis sets.

        Parameters
        ----------
        vector : list or array-like
            Array of size (m, n).
            Contains m vectors built from n kmers in the basis.
        vector_basis : list or array-like of str
            Array of size (n,) of ordered kmers.

        Returns
        -------
        list or array-like
            Transformed array of size (m, p).

        """
        if not check_list(vector_basis):
            raise TypeError("`vector_basis` input must be list or array-like.")

        if not isinstance(vector, np.ndarray):
            vector = np.asarray(vector)

        # make sure input vector matches shape of vector basis
        try:
            vector_size = vector.shape[1]
        except IndexError:
            vector_size = len(vector)

        if vector_size != len(vector_basis):
            raise ValueError(
                "Vector and supplied basis shapes"
                " must match (vector shape ="
                f" {vector.shape}"
                " and len(vector_basis) ="
                f" {len(vector_basis)})."
            )

        # get index order of kmers in the vector basis set
        vector_basis_order = {
            k: i if k in self.basis else np.nan for i, k in enumerate(vector_basis)
        }

        # convert vector basis into set basis
        i_convert = list()
        for i in range(len(self.basis)):
            kmer = self.basis_order[i]  # get basis set kmer in correct order
            idx = vector_basis_order[kmer]  # locate kmer in the new vector
            i_convert.append(idx)

        # correctly index into ND or 1D array
        try:
            return vector[:, i_convert]
        except IndexError:
            return vector[i_convert]


# generate all possible kmer combinations
def _generate(alphabet: Set[str], k: int):
    for c in itertools.product(alphabet, repeat=k):
        yield "".join(c)


# iterator object for kmer basis set given alphabet and k
class KmerSet:
    """Given alphabet and k, creates iterator for kmer basis set.
    """

    def __init__(self, alphabet: Union[str, int], k: int):
        self.alphabet = alphabet
        self.k = k
        self._kmerlist = list(_generate(get_alphabet_keys(alphabet), k))

    @property
    def kmers(self):
        return iter(self._kmerlist)


# manually reduce alphabet
def reduce(
    sequence: str, alphabet: Union[str, int], mapping: dict = FULL_ALPHABETS
) -> str:
    """Reduce sequence into character space of specified alphabet.

    Parameters
    ----------
    sequence : str
        Input sequence.
    alphabet : Union[str, int]
        Alphabet name or number (see `snekmer.alphabet`).
    mapping : dict
        Defined mapping for alphabet (the default is FULL_ALPHABETS).

    Returns
    -------
    str
        Transformed sequence.

    """
    sequence = str(sequence).rstrip("*")
    alphabet_map: Dict[str, str] = get_alphabet(alphabet, mapping=mapping)
    return sequence.translate(sequence.maketrans(alphabet_map))


class KmerVec:
    def __init__(self, alphabet: Union[str, int], k: int):
        self.alphabet = alphabet
        self.k = k
        self.char_set = get_alphabet_keys(alphabet)
        self.vector = None
        self.kmer_set = KmerSet(alphabet, k)

    # iteratively get all kmers in a string
    def _kmer_gen(self, sequence: str) -> Generator[str, None, None]:
        """Generator object for segment in string format"""
        i = 0
        n = len(sequence) - self.k + 1

        # iterate thru sequence in blocks of length k
        while i < n:
            kmer = sequence[i : i + self.k]
            if set(kmer) <= self.char_set:
                yield kmer
            i += 1

    # not used: iterate using range() vs. while loop
    @staticmethod
    def _kmer_gen_str(sequence: str, k: int) -> Generator[str, None, None]:
        """Generator object for segment in string format"""
        for n in range(0, len(sequence) - k + 1):
            yield sequence[n : n + k]

    # generate kmer vectors with bag-of-words approach
    def vectorize(self, sequence: str) -> NDArray:
        """Transform sequence into representative kmer vector.

        Parameters
        ----------
        sequence : str
            Input sequence.

        Returns
        -------
        NDArray
            Vector representation of sequence as kmer counts vector.

        """
        N = len(self.char_set) ** self.k

        kmers = list(self._kmer_gen(sequence))
        kmer2count = Counter(kmers)

        # Convert to vector of counts
        vector = np.zeros(N)
        for i, word in enumerate(self.kmer_set.kmers):
            vector[i] += kmer2count[word]

        # Convert to frequencies
        # vector /= sum(kmer2count.values())

        return vector

    def reduce_vectorize(self, sequence: str) -> NDArray:
        """Simplify and vectorize sequence into reduced kmer vector.

        Parameters
        ----------
        sequence : str
            Input sequence.

        Returns
        -------
        NDArray
            Vector representation of sequence as reduced kmer vector.

        """
        N = len(self.char_set) ** self.k

        reduced = reduce(sequence, alphabet=self.alphabet, mapping=FULL_ALPHABETS)
        kmers = list(self._kmer_gen(reduced))
        kmer2count = Counter(kmers)

        # Convert to vector of counts
        vector = np.zeros(N)
        for i, word in enumerate(self.kmer_set.kmers):
            vector[i] += kmer2count[word]

        # Convert to frequencies
        # vector /= sum(kmer2count.values())

        return vector
