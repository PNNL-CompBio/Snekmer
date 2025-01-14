"""alphabet: Built-in alphabet mappings for Snekmer.

2021.08.11 - Alphabet names have been reworked per @wichne
             (@christinehc)

author(s): @biodataganache, @wichne

"""
# imports
from typing import Dict, Mapping, Set, Union

# define standard amino acid alphabet
StandardAlphabet = "AILMVFYWSTQNCHDEKRGP"
AA_SELF_MAPPING = {a: a for a in StandardAlphabet}

# post-translational modification
PTM_CHARS = "-_!^#$@.%&"
PTM_SELF_MAPPING = {c: c for c in PTM_CHARS}

# define alphabet names and ordering
ALPHABET_ORDER = {
    0: "hydro",
    1: "standard",
    2: "solvacc",
    3: "hydrocharge",
    4: "hydrostruct",
    5: "miqs",
}

# define alphabets (from Utils.SIEVEInit)
ALPHABETS = {
    # 2-value hydrophobicity alphabet taken
    # from Arnold, et al. 2009. PLoS Pathogens 5(4): e1000376
    "hydro": {"SFTNKYEQCWPHDR": "S", "VMLAIG": "V", "_keys": "SV"},
    # 'standard' reduction alphabet taken
    # from Arnold, et al. 2009. PLoS Pathogens 5(4): e1000376
    "standard": {
        "AGILMV": "A",  # hydrophobic
        "PH": "P",  # hydrophilic
        "FWY": "F",  # aromatic
        "NQST": "N",  # polar
        "DE": "D",  # acidic
        "KR": "K",  # alkaline
        "C": "C",  # ionizable
        "_keys": "APFNDKC",
    },
    # Solvent accessibility alphabet
    # from Bacardit, et al. 2009. BMC Bioinformatics 10:6
    "solvacc": {"CILMVFWY": "C", "AGHST": "A", "PDEKNQR": "P", "_keys": "CAP"},
    # 2-value hydrophobicity with charged residues as a third
    # category. Made by @biodataganache.
    "hydrocharge": {
        "SFTNYQCWPH": "L",  # hydrophilic (L-ove)
        "VMLAIG": "H",  # hydrophobic (H-ate)
        "KNDER": "C",  # charged (C-harged)
        "_keys": "LHC",
    },
    # 2-value hydrophobicity with structural-breakers as a third category
    # Made by @biodataganache
    "hydrostruct": {"SFTNKYEQCWHDR": "L", "VMLAI": "H", "PG": "B", "_keys": "LHB"},
    # MIQS alphabet
    # Made by @wichne
    "miqs": {
        "A": "A",  # Alanine
        "C": "C",  # Cysteine
        "DEN": "D",  # acidicish
        "FWY": "F",  # aromatic
        "G": "G",  # glycine
        "H": "H",  # histidine
        "ILMQV": "I",  # hydrophobicish
        "KR": "K",  # alkaline
        "P": "P",  # proline
        "ST": "S",  # hydroxyl
        "_keys": "ACDFGHIKPS",
    },
    # # identity
    # "None": {**AA_SELF_MAPPING, "_keys": ALL_AMINO_ACIDS}  # OU
    # post-translational modification alphabet
    "ptm": {
        **AA_SELF_MAPPING,
        **PTM_SELF_MAPPING,
        "_keys": StandardAlphabet + PTM_CHARS,
    },
    "None": AA_SELF_MAPPING,
}

# reconfigure alphabet dict into "long-form"
FULL_ALPHABETS: Dict[str, dict] = {a: {} for a in ALPHABETS.keys()}
for alphabet, mapping in ALPHABETS.items():
    for k, v in mapping.items():
        if k == "_keys":
            continue
        elif len(k) > 1:
            FULL_ALPHABETS[alphabet].update({k[i]: v for i in range(len(k))})
        else:
            FULL_ALPHABETS[alphabet].update({k: v})

# create generic alphabet identifiers
ALPHABET_ID = {
    f"RED{n}": {v: k for k, v in ALPHABETS[ALPHABET_ORDER[n]].items()}
    for n in range(len(ALPHABET_ORDER))
}

# alphabet to alphabet code
ALPHABET2ID = {ALPHABET_ORDER[n]: f"RED{n}" for n in range(len(ALPHABET_ORDER))}


# SIEVEInit alphabet grabbing function
def get_alphabets():
    """Return all alphabet mappings.

    Returns
    -------
    dict of dicts
        {'alphabet': {'key': 'mapping', ...}} for each alphabet

    """
    return ALPHABETS


def check_valid(alphabet: Union[str, int]) -> None:
    """Check validity of input alphabet vs. defined list.

    Parameters
    ----------
    alphabet : Union[str, int]
        Alphabet name or identifier.

    Raises
    ------
    ValueError
        If `alphabet` does not match a pre-defined alphabet.

    """
    if (
        # raise error if alphabet not included in pre-defined set
        (
            # doesn't match integer designations
            (alphabet not in range(len(ALPHABETS)))
            # doesn't match str name designations
            and (alphabet not in ALPHABETS)
        )
        # or doesn't match None (no mapping)
        and (str(alphabet) != "None")
    ):  # and config['alphabet'] != 'custom':
        raise ValueError(
            "Invalid alphabet specified; alphabet must be a"
            " string (see snekmer.alphabet) or integer"
            " n between"
            f" {min(list(ALPHABET_ORDER.keys()))}"
            " and"
            f" {max(list(ALPHABET_ORDER.keys()))}"
            "."
        )
    return


def get_alphabet(
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

    # enforce string for dict key
    if alphabet is None:
        alphabet = str(alphabet)

    if isinstance(alphabet, int):
        alphabet = ALPHABET_ORDER[alphabet]
    return mapping[alphabet]


def get_alphabet_name(
    alphabet: Union[str, int], mapping: dict = ALPHABETS
) -> Dict[str, str]:
    """Get alphabet name given any input type.

    Parameters
    ----------
    alphabet : Union[str, int]
        Alphabet name (as str) or alphabet id (as int).
        Must be one of the follwing:
            0: "hydro",
            1: "standard",
            2: "solvacc",
            3: "hydrocharge",
            4: "hydrostruct",
            5: "miqs"
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

    # enforce string for dict key
    if alphabet is None:
        alphabet = str(alphabet)

    if isinstance(alphabet, int):
        alphabet = ALPHABET_ORDER[alphabet]
    return alphabet


def get_alphabet_keys(
    alphabet: Union[str, int], mapping: Dict[str, dict] = FULL_ALPHABETS
) -> Set[str]:
    """Retrieve keys for specified alphabet.

    Parameters
    ----------
    alphabet : Union[str, int]
        Description of parameter `alphabet`.
    mapping : Mapping[dict]
        Description of parameter `mapping` (the default is FULL_ALPHABETS).

    Returns
    -------
    dict
        Description of returned object.

    """
    # enforce string for dict key
    if alphabet is None:
        alphabet = str(alphabet)

    alphabet_map = get_alphabet(alphabet, mapping)
    if "_keys" in alphabet_map.keys():
        alphabet_map.pop("_keys")
    return set(alphabet_map.values())


# def add_alphabet(alphabet_name, mapping):
#     return
