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
        "KNDR": "C",  # charged (C-harged)
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
}

# reconfigure alphabet dict into "long-form"
FULL_ALPHABETS: dict = {a: {} for a in ALPHABETS.keys()}
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


def check_valid(alphabet):
    """Check validity of input alphabet vs. defined list.

    Parameters
    ----------
    alphabet : str
        Alphabet name.

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

    if isinstance(alphabet, int):
        alphabet = ALPHABET_ORDER[alphabet]
    return mapping[alphabet]


def get_alphabet_keys(
    alphabet: Union[str, int], mapping: dict = FULL_ALPHABETS
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
    alphabet_map = get_alphabet(alphabet, mapping)
    if "_keys" in alphabet_map.keys():
        alphabet_map.pop("_keys")
    return set(alphabet_map.values())


# def add_alphabet(alphabet_name, mapping):
#     return


ORGANISM_ORDER = [
    "agrobacterium_tumefaciens_c58_uwash",
    "agrobacterium_tumefaciens_c58_cereon",
    "agrobacterium_tumefaciens_a348",
    "bacillus_anthracis",
    "bacillus_subtilis",
    "bordetella_pertussis",
    "brucella_melitensis",
    "brucella_suis",
    "campylobacter_jejuni",
    "chlamydia_trachomatis",
    "clostridium_perfringens",
    "escherichia_coli",
    "halobacterium_sp",
    "helicobacter_pylori",
    "listeria_monocytogenes",
    "methanococcus_jannaschii",
    "methanococcus_maripaludis",
    "mycobacterium_bovis",
    "mycobacterium_tuberculosis",
    "neisseria_meningitidis",
    "pseudomonas_aeruginosa",
    "pyrococcus_abyssi",
    "pyrococcus_furiosus",
    "rhodopseudomonas_palustris",
    "rickettsia_conorii",
    "rickettsia_prowazekii",
    "salmonella_typhimurium",
    "shewanella_oneidensis",
    "shigella_flexneri",
    "staphylococcus_aureus",
    "thermotoga_maritima",
    "vibrio_cholerae",
    "vibrio_parahaemolyticus",
    "vibrio_vulnificus",
    "yersinia_pestis",
    "yersinia_pestis_bgi_91001",
    "yersinia_pestis_bgi_kim",
    "yersinia_pestis_bgi_co92",
    "arabidopsis_thaliana",
    "caenorhabditis_elegans",
    "canis_familiaris",
    "drosophila_melanogaster",
    "encephalitozoon_cuniculi",
    "homo_sapiens",
    "magnaporthe_grisea",
    "mus_musculus",
    "mus_musculus_bgi",
    "oryza_sativa_japonica_fl",
    "oryza_sativa_indica_9311",
    "oryza_sativa_japonica_syngenta",
    "pan_troglodytes",
    "plasmodium_falciparum",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
]

ORGANISM_DICT = {
    "agrobacterium_tumefaciens_a348": 0.0,
    "agrobacterium_tumefaciens_c58_cereon": 0.0,
    "agrobacterium_tumefaciens_c58_uwash": 0.0,
    "arabidopsis_thaliana": 0.0,
    "bacillus_anthracis": 0.0,
    "bacillus_subtilis": 0.0,
    "bordetella_pertussis": 0.0,
    "brucella_melitensis": 0.0,
    "brucella_suis": 0.0,
    "caenorhabditis_elegans": 0.0,
    "campylobacter_jejuni": 0.0,
    "canis_familiaris": 0.0,
    "chlamydia_trachomatis": 0.0,
    "clostridium_perfringens": 0.0,
    "drosophila_melanogaster": 0.0,
    "encephalitozoon_cuniculi": 0.0,
    "escherichia_coli": 0.0,
    "halobacterium_sp": 0.0,
    "helicobacter_pylori": 0.0,
    "homo_sapiens": 0.0,
    "listeria_monocytogenes": 0.0,
    "magnaporthe_grisea": 0.0,
    "methanococcus_jannaschii": 0.0,
    "methanococcus_maripaludis": 0.0,
    "mus_musculus": 0.0,
    "mus_musculus_bgi": 0.0,
    "mycobacterium_bovis": 0.0,
    "mycobacterium_tuberculosis": 0.0,
    "neisseria_meningitidis": 0.0,
    "oryza_sativa_indica_9311": 0.0,
    "oryza_sativa_japonica_fl": 0.0,
    "oryza_sativa_japonica_syngenta": 0.0,
    "pan_troglodytes": 0.0,
    "plasmodium_falciparum": 0.0,
    "pseudomonas_aeruginosa": 0.0,
    "pyrococcus_abyssi": 0.0,
    "pyrococcus_furiosus": 0.0,
    "rattus_norvegicus": 0.0,
    "rhodopseudomonas_palustris": 0.0,
    "rickettsia_conorii": 0.0,
    "rickettsia_prowazekii": 0.0,
    "saccharomyces_cerevisiae": 0.0,
    "salmonella_typhimurium": 0.0,
    "shewanella_oneidensis": 0.0,
    "shigella_flexneri": 0.0,
    "staphylococcus_aureus": 0.0,
    "thermotoga_maritima": 0.0,
    "vibrio_cholerae": 0.0,
    "vibrio_parahaemolyticus": 0.0,
    "vibrio_vulnificus": 0.0,
    "yersinia_pestis": 0.0,
    "yersinia_pestis_bgi_91001": 0.0,
    "yersinia_pestis_bgi_co92": 0.0,
    "yersinia_pestis_bgi_kim": 0.0,
}
