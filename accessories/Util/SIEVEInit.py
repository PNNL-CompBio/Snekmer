ALPHABETS = {
    # 2-value hydrophobicity alphabet taken
    # from Arnold, et al. 2009. PLoS Pathogens 5(4): e1000376
    "reduced_alphabet_0":{"SFTNKYEQCWPHDR":"S",
                          "VMLAIG":"V",
                          "_keys":"SV"},

    # 'standard' reduction alphabet taken
    # from Arnold, et al. 2009. PLoS Pathogens 5(4): e1000376
    "reduced_alphabet_1":{"AGILMV":"A",   # hydrophobic
                          "PH":"P",       # hydrophilic
                          "FWY":"F",      # aromatic
                          "NQST":"N",     # polar
                          "DE":"D",       # acidic
                          "KR":"K",       # alkaline
                          "CY":"C",       # ionizable
                          "_keys":"APFNDKC"},      

    # Solvent accessibility alphabet
    # from Bacardit, et al. 2009. BMC Bioinformatics 10:6
    "reduced_alphabet_2":{"CILMVFWY":"C",
                          "AGHST":"A",
                          "PDEKNQR":"P",
                          "_keys":"CAP"}}

ALPHABETCODE = {"RED0":{v:k for k, v in ALPHABETS["reduced_alphabet_0"].items()},
                "RED1":{v:k for k, v in ALPHABETS["reduced_alphabet_1"].items()},
                "RED2":{v:k for k, v in ALPHABETS["reduced_alphabet_2"].items()}}

StandardAlphabet = "AILMVFYWSTQNCHDEKRGP"

def get_alphabets():
    return ALPHABETS

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
                "saccharomyces_cerevisiae"]

ORGANISM_DICT = {
                "agrobacterium_tumefaciens_a348":0.0,
                "agrobacterium_tumefaciens_c58_cereon":0.0,
                "agrobacterium_tumefaciens_c58_uwash":0.0,
                "arabidopsis_thaliana":0.0,
                "bacillus_anthracis":0.0,
                "bacillus_subtilis":0.0,
                "bordetella_pertussis":0.0,
                "brucella_melitensis":0.0,
                "brucella_suis":0.0,
                "caenorhabditis_elegans":0.0,
                "campylobacter_jejuni":0.0,
                "canis_familiaris":0.0,
                "chlamydia_trachomatis":0.0,
                "clostridium_perfringens":0.0,
                "drosophila_melanogaster":0.0,
                "encephalitozoon_cuniculi":0.0,
                "escherichia_coli":0.0,
                "halobacterium_sp":0.0,
                "helicobacter_pylori":0.0,
                "homo_sapiens":0.0,
                "listeria_monocytogenes":0.0,
                "magnaporthe_grisea":0.0,
                "methanococcus_jannaschii":0.0,
                "methanococcus_maripaludis":0.0,
                "mus_musculus":0.0,
                "mus_musculus_bgi":0.0,
                "mycobacterium_bovis":0.0,
                "mycobacterium_tuberculosis":0.0,
                "neisseria_meningitidis":0.0,
                "oryza_sativa_indica_9311":0.0,
                "oryza_sativa_japonica_fl":0.0,
                "oryza_sativa_japonica_syngenta":0.0,
                "pan_troglodytes":0.0,
                "plasmodium_falciparum":0.0,
                "pseudomonas_aeruginosa":0.0,
                "pyrococcus_abyssi":0.0,
                "pyrococcus_furiosus":0.0,
                "rattus_norvegicus":0.0,
                "rhodopseudomonas_palustris":0.0,
                "rickettsia_conorii":0.0,
                "rickettsia_prowazekii":0.0,
                "saccharomyces_cerevisiae":0.0,
                "salmonella_typhimurium":0.0,
                "shewanella_oneidensis":0.0,
                "shigella_flexneri":0.0,
                "staphylococcus_aureus":0.0,
                "thermotoga_maritima":0.0,
                "vibrio_cholerae":0.0,
                "vibrio_parahaemolyticus":0.0,
                "vibrio_vulnificus":0.0,
                "yersinia_pestis":0.0,
                "yersinia_pestis_bgi_91001":0.0,
                "yersinia_pestis_bgi_co92":0.0,
                "yersinia_pestis_bgi_kim":0.0
             }
