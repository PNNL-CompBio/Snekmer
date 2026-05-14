# require snakemake 9.0+ (required for modules)
from snakemake.utils import min_version

min_version("9.0")


# load snakemake modules
module process:
    snakefile:
        "process.smk"
    config:
        config


module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


from glob import glob
from itertools import product
from os.path import basename, dirname, join

import os
import shutil

import snekmer as skm
from importlib.resources import files


# collect all fasta-like files, unzipped filenames, and basenames
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = glob(join(input_dir, "*"))
# base_file = glob(join(input_dir,"base" "*"))xf
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.removesuffix(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.removesuffix(".gz").endswith(f".{ext}")
]
compare_file = glob(join("counts", "*.csv"))
confidence_file = glob(join("confidence", "*.csv"))
decoy_stats_file = glob(join("stats", "*.csv"))
# map extensions to basename (basename.ext.gz -> {basename: ext})
uz_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
fa_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
# seq-annotation-scores
# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in uz_map.items()]
# isolate basenames for all files
FAS = list(fa_map.keys())
FAV = list(fa_map.values())
# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
non_bgs, bgs = [f for f in FAS if f not in bg_files], bg_files
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)
out_dir = str(out_dir)

threshold_type = config["learn_apply"]["threshold"]
selection_type = config["learn_apply"]["selection"]

apply_cfg = config.get("learn_apply", {})
concat_results = apply_cfg.get("apply_output", "snekmer_results.csv")
extra_all = [concat_results] if concat_results is not None else []


onstart:
    import sys as _sys
    from Bio import SeqIO as _SeqIO

    _errors = []
    _warns = []

    if not FAS:
        _errors.append(
            f"No FASTA files found in '{input_dir}/'. "
            f"Expected extensions: {config['input_file_exts']}"
        )

    for _label, _files, _subdir in [
        ("k-mer counts",      compare_file,     "counts/"),
        ("confidence scores", confidence_file,  "confidence/"),
        ("family stats",      decoy_stats_file, "stats/"),
    ]:
        if not _files:
            _errors.append(
                f"No {_label} CSV found in '{_subdir}'. "
                f"Run 'snekmer learn' first, then copy apply_inputs/ files here."
            )

    for _fa in unzipped:
        try:
            _recs = list(_SeqIO.parse(_fa, "fasta"))
        except Exception as _exc:
            _errors.append(f"Cannot parse FASTA '{_fa}': {_exc}")
            continue
        if not _recs:
            _errors.append(
                f"FASTA file has no sequences: '{_fa}'\n"
                f"  Verify the file uses FASTA format (records start with '>')."
            )
            continue
        _short = [r for r in _recs if len(r.seq) < config["k"]]
        if _short:
            _warns.append(
                f"{len(_short)}/{len(_recs)} sequence(s) in '{_fa}' are shorter "
                f"than k={config['k']} and will yield no k-mers."
            )
        _sample = "".join(str(r.seq).upper() for r in _recs[:10])
        if _sample and sum(1 for c in _sample if c in "ATCGN") / len(_sample) > 0.9:
            _warns.append(
                f"'{_fa}' may contain nucleotide sequences. "
                f"Snekmer expects amino acid (protein) input."
            )

    # ── Parameter checks ─────────────────────────────────────────────────────
    _VALID_SEL = {"top_hit", "greatest_distance", "combined_distance"}

    if config["k"] < 1:
        _errors.append(f"k must be ≥ 1, got k={config['k']}.")

    try:
        from snekmer.alphabet import get_alphabet_keys as _gak
        _n_sym = len(_gak(config["alphabet"]))
        _vocab = _n_sym ** config["k"]
        if _vocab > 50_000_000:
            _errors.append(
                f"Vocabulary too large: alphabet '{config['alphabet']}' has {_n_sym} "
                f"symbols, k={config['k']} → {_vocab:,} k-mers. This will exhaust memory.\n"
                f"  Use a smaller k or a coarser alphabet (e.g., alphabet=0, k=6)."
            )
        elif _vocab > 1_000_000:
            _warns.append(
                f"Large vocabulary: alphabet '{config['alphabet']}' has {_n_sym} symbols, "
                f"k={config['k']} → {_vocab:,} k-mers. Runtime and memory may be high."
            )
    except Exception:
        pass

    _sel = config.get("learn_apply", {}).get("selection", "top_hit")
    if _sel not in _VALID_SEL:
        _errors.append(
            f"Unknown selection method '{_sel}'. "
            f"Valid options: {', '.join(sorted(_VALID_SEL))}"
        )

    if _sel == "combined_distance":
        _la = config.get("learn_apply", {})
        _wt, _wd = _la.get("weight_top", 0.7), _la.get("weight_distance", 0.3)
        if abs(_wt + _wd - 1.0) > 0.01:
            _warns.append(
                f"weight_top ({_wt}) + weight_distance ({_wd}) = {_wt + _wd:.2f} "
                f"(expected 1.0). Combined-distance scores may be poorly calibrated."
            )

    for _w in _warns:
        print(f"\nWARNING: {_w}", file=_sys.stderr)
    if _errors:
        for _e in _errors:
            print(f"\nINPUT ERROR: {_e}", file=_sys.stderr)
        raise Exception("\nInput validation failed — fix the errors above and re-run.")


def resource_path(package: str, *parts) -> str:
    """
    Re-create pkg_resources.resource_filename()
    but using importlib.resources.
    """
    return str(files(package).joinpath(*parts))


rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        *[
            (
                expand(
                    join(out_dir, "apply", "seq_annotation_scores_{nb}.csv"), nb=FAS
                )
                if config["learn_apply"]["save_apply_associations"]
                else []
            )
        ],
        *extra_all,
        expand(join(out_dir, "apply", "kmer_summary_{nb}.csv"), nb=FAS),
        join(out_dir, "Snekmer_Apply_Report.html"),


use rule unzip from process with:
    wildcard_constraints:
        uz=r".*\.(?:fa|fna|faa|fasta)$",
    output:
        unzipped=join(input_dir, "{uz}"),
        zipped=join(input_dir, "zipped", "{uz}.gz"),


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(
            input_dir, f"{wildcards.nb}.{fa_map[wildcards.nb]}"
        ),
    output:
        data=join(out_dir, "vector", "{nb}.npz"),
        kmerobj=join(out_dir, "kmerize", "{nb}.kmers"),
    message:
        "Kmerizing and re-encoding Amino acids in {input.fasta}. Output written to {output.data}."


rule apply:
    input:
        data=join(out_dir, "vector", "{nb}.npz"),
        compare_associations=expand("{comp}", comp=compare_file),
        confidence_associations=expand("{conf}", conf=confidence_file),
        decoy_stats=expand("{decoy}", decoy=decoy_stats_file),
    params:
        selection_type=config["learn_apply"]["selection"],
        threshold_type=config["learn_apply"]["threshold"],
    output:
        seq_ann=(
            join(out_dir, "apply", "seq_annotation_scores_{nb}.csv")
            if config["learn_apply"]["save_apply_associations"]
            else []
        ),
        kmer_summary=join(out_dir, "apply", "kmer_summary_{nb}.csv"),
    message:
        "Running Snekmer Apply on {input.data}. Output written to {output.kmer_summary}."
    script:
        resource_path("snekmer", "scripts", "apply.py")


if concat_results is not None:

    rule concat_kmer_summary:
        """Concatenate all sample kmer summaries into one CSV"""
        input:
            expand(join(out_dir, "apply", "kmer_summary_{nb}.csv"), nb=FAS),
        output:
            concat_results,
        message:
            "Writing consolidated k-mer summary to {output}"
        script:
            resource_path("snekmer", "scripts", "apply_concat.py")


rule apply_report:
    input:
        seq_scores=(
            expand(join(out_dir, "apply", "seq_annotation_scores_{nb}.csv"), nb=FAS)
            if config["learn_apply"]["save_apply_associations"]
            else []
        ),
        kmer_sum=expand(join(out_dir, "apply", "kmer_summary_{nb}.csv"), nb=FAS),
        concat=(
            [config["learn_apply"]["apply_output"]]
            if config.get("learn_apply", {}).get("apply_output")
            else []
        ),
    output:
        report=join(out_dir, "Snekmer_Apply_Report.html"),
    message:
        "Generating full Snekmer Apply Report at {output.report}"
    script:
        resource_path("snekmer", "scripts", "apply_report.py")
