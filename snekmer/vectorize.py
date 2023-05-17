"""vectorize: Create kmer vectors.

author: @christinehc

"""
import itertools
import random
import string
from collections import Counter
from typing import Dict, Generator, Optional, Set, Union

import numpy as np
from biotite.sequence import GeneralSequence, LetterAlphabet, ProteinSequence
from biotite.sequence.align import KmerAlphabet, KmerTable
from numpy.typing import NDArray
from ._version import __version__
from .alphabet import (
    ALPHABETS,
    FULL_ALPHABETS,
    check_valid,
    get_alphabet,
    get_alphabet_name,
)
from .utils import check_list

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


# ---------
# NEW CODE
# ---------


def unpack_table(
    table, ref_id: Union[int, tuple[int, int]], kmers: Optional[NDArray] = None
) -> NDArray:
    """_summary_

    Parameters
    ----------
    table : _type_
        _description_
    ref_id : Union[int, tuple[int, int]]
        _description_
    kmers : Optional[NDArray], optional
        _description_, by default None

    Returns
    -------
    NDArray
        _description_
    """
    if kmers is None:
        kmers = table.get_kmers()

    if isinstance(ref_id, int):
        ref_id = 0, ref_id

    i_min, i_max = ref_id
    matrix = np.zeros((i_max - i_min, len(kmers)))

    for j, kmer in enumerate(kmers):
        subset = Counter(table[kmer][:, 0])
        for i, count in subset.items():
            matrix[i - i_min][j] = count
    return matrix


def _get_alphabet_chars(
    alphabet: Union[str, int], mapping: dict = ALPHABETS
) -> Dict[str, str]:
    """Short summary.

    Parameters
    ----------
    alphabet : Union[str, int]
        Alphabet name (as str) or alphabet id (as int).
        Must be one of the following:
            0: "hydro",
            1: "standard",
            2: "solvacc",
            3: "hydrocharge",
            4: "hydrostruct",
            5: "miqs"
        or None.
    mapping : dict
        All alphabet maps (the default is ALPHABETS).

    Returns
    -------
    dict
        Dictionary map of amino acids to alphabet character.

    Raises
    ------
    ValueError
        Raised if alphabet not in pre-defined list.

    """
    check_valid(alphabet)

    # enforce alphabet name
    alphabet_name = get_alphabet_name(alphabet)
    if alphabet_name == "None":
        return "".join(sorted(mapping[alphabet].values()))
    return sorted(mapping[alphabet]["_keys"])


def _generate_kmers(codes, alphabet, as_mapping=False):
    for c in codes:
        if not as_mapping:
            yield "".join(alphabet.decode(c))
        else:
            yield c, "".join(alphabet.decode(c))


def _generate_vector(kmerlist, table):
    for kmer in kmerlist:
        yield len(table[kmer])


def _generate_id(size=2, chars=string.ascii_uppercase + string.digits):
    return "".join(random.choice(chars) for _ in range(size))


class KmerVecs:
    def __init__(self, k, alphabet, family=None, ref_id=None):
        self.k = k
        self.alphabet_name = get_alphabet_name(alphabet)
        self.alphabet = LetterAlphabet(_get_alphabet_chars(self.alphabet_name))
        self.ref_id = ref_id

    def kmerize(self, seqs):
        if self.ref_id is None:
            self.ref_id = (0, len(seqs))
        table = KmerTable.from_sequences(
            k=self.k,
            sequences=[
                GeneralSequence(
                    sequence=reduce(s, self.alphabet_name), alphabet=self.alphabet
                )
                for s in seqs
            ],
            ref_ids=np.arange(self.ref_id[0], self.ref_id[1]),
            alphabet=self.alphabet,
        )
        self.table = table

    @property
    def kmer_alphabet(self):
        return self.table.kmer_alphabet

    @property
    def kmer_codes(self):
        return self.table.get_kmers()

    @property
    def kmers(self):
        return _generate_kmers(self.table.get_kmers(), self.kmer_alphabet)

    # ["".join(self.kmer_alphabet.decode(kmer)) for kmer in self.table.get_kmers()]

    @property
    def kmer_map(self):
        return _generate_kmers(
            self.table.get_kmers(), self.kmer_alphabet, as_mapping=True
        )

    # return counts vec
    @property
    def vector(self):
        """Get vector as generator (direct access via dict(self.vector))"""
        # if not hasattr(self, "_vector"):
        self._vector = _generate_vector(self.kmer_codes, self.table)
        # {
        #     "".join(self.kmer_alphabet.decode(kmer)): len(self.table[kmer])
        #     for kmer in self.kmer_codes
        # }
        return self._vector
