# require snakemake v9 (required for modules)
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
from os.path import basename, join
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
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.removesuffix(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.removesuffix(".gz").endswith(f".{ext}")
]
annot_files = glob(join("annotations", "*.ann"))
base_counts = glob(join("base", "counts", "*.csv"))
base_confidence = glob(join("base", "confidence", "*.csv"))
base_family_checkpoint = glob(join("base", "thresholds", "*.csv"))

uz_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
fa_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in uz_map.items()]

# isolate basenames for all files
FAS = list(fa_map.keys())
# parse any background files
background_files = glob(join(input_dir, "background", "*"))
if len(background_files) > 0:
    background_files = [
        skm.utils.split_file_ext(basename(f))[0] for f in background_files
    ]

# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)
out_dir = str(out_dir)

if (
    config["learn_apply"]["selection"] != "top_hit"
    and config["learn_apply"]["threshold"] == "None"
):
    raise Exception(
        "The only selection method that allows for None is `top_hit`. Other methods inherently use a threshold"
    )


def resource_path(package: str, *parts) -> str:
    """
    Re-create pkg_resources.resource_filename()
    but using importlib.resources.
    """
    return str(files(package).joinpath(*parts))


rule all:
    input:
        # Vector files
        expand(join(out_dir, "vector", "vector", "{nb}.npz"), nb=FAS),
        # Kmer counts for learning
        expand(join(out_dir, "learn", "kmer_counts_{nb}.csv"), nb=FAS),
        join(out_dir, "learn", "kmer_counts_total.csv"),
        # Fragmentation outputs (only if enabled)
        expand(join(out_dir, "fragmented", "{nb}.fasta"), nb=FAS)
        if config["learn_apply"]["fragmentation"]
        else [],
        expand(join(out_dir, "vector", "vector_frag", "{nb}.npz"), nb=FAS)
        if config["learn_apply"]["fragmentation"]
        else [],
        # Forward evaluation scores
        expand(
            join(
                out_dir,
                "evaluate",
                (
                "eval_apply_sequences"
                    if not config["learn_apply"]["fragmentation"]
                    else "eval_apply_frag"
                ),
                "seq_annotation_scores_{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        # Reverse evaluation scores
        expand(
            join(
                out_dir,
                "evaluate",
                (
                "eval_apply_reversed"
                    if not config["learn_apply"]["fragmentation"]
                    else "eval_apply_reversed_frag"
                ),
                "seq_annotation_scores_{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        # Evaluation-level summaries & confidence
        join(out_dir, "eval_conf", "family_summary_stats.csv"),
        join(out_dir, "eval_conf", "global_confidence_scores.csv"),
        # Apply Inputs Copy
        join("apply_inputs", "counts", "kmer_counts_total.csv"),
        join("apply_inputs", "stats", "family_summary_stats.csv"),
        join("apply_inputs", "confidence", "global_confidence_scores.csv"),
        # Report
        join(out_dir, "Snekmer_Learn_Report.html"),


use rule unzip from process with:
    wildcard_constraints:
        uz = r".*\.(?:fa|fna|faa|fasta)$"
    output:
        unzipped=join(input_dir, "{uz}"),
        zipped=join(input_dir, "zipped", "{uz}.gz"),


if config["learn_apply"]["fragmentation"]:

    rule fragmentation:
        input:
            fasta=lambda wildcards: join(
                input_dir, f"{wildcards.nb}.{fa_map[wildcards.nb]}"
            ),
        output:
            fasta_out=join(out_dir, "fragmented", "{nb}.fasta"),
        message:
            "Fragmenting sequences in {input.fasta}. Output written to {output.fasta_out}."
        params:
            version=config["learn_apply"]["version"],
            frag_length=config["learn_apply"]["frag_length"],
            location=config["learn_apply"]["location"],
            min_length=config["learn_apply"]["min_length"],
            seed=config["learn_apply"]["seed"],
        script:
            resource_path("snekmer", "scripts", "learn_fragment.py")


prefix = "vector"
use rule vectorize from kmerize as vectorize_vector with:
    input:
        fasta=lambda wc: join(input_dir, f"{wc.nb}.{fa_map[wc.nb]}"),
    output:
        data=join(out_dir, "vector", "vector", "{nb}.npz"),
        kmerobj=join(out_dir, "kmerize", "kmer", "{nb}.kmers"),


prefix = "vector_frag"
use rule vectorize from kmerize as vectorize_frag with:
    input:
        fasta=lambda wc: join(out_dir, "fragmented", f"{wc.nb}.fasta"),
    output:
        data=join(out_dir, "vector", "vector_frag", "{nb}.npz"),
        kmerobj=join(out_dir, "kmerize", "kmer_frag", "{nb}.kmers"),


rule learn:
    input:
        data=join(out_dir, "vector", "vector", "{nb}.npz"),
        annotation=expand("{an}", an=annot_files),
    output:
        counts=join(out_dir, "learn", "kmer_counts_{nb}.csv"),
    message:
        "Building kmer-association matrix from {input.data}. Output written to {output.counts}."
    script:
        resource_path("snekmer", "scripts", "learn_learn.py")


rule merge:
    input:
        counts=expand(join(out_dir, "learn", "kmer_counts_{nb}.csv"), nb=FAS),
        base_counts=expand("{bf}", bf=base_counts),
    output:
        totals=join(out_dir, "learn", "kmer_counts_total.csv"),
    message:
        "Merging individual k-mer association matrix files into consolidated {output.totals}."
    script:
        resource_path("snekmer", "scripts", "learn_merge.py")


rule eval_apply_reverse_seqs:
    input:
        data=join(
            out_dir,
            "vector",
            (
            "vector"
                if not config["learn_apply"]["fragmentation"]
                else "vector_frag"
            ),
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annot_files),
        compare_associations=join(out_dir, "learn", "kmer_counts_total.csv"),
    output:
        apply=join(
            out_dir,
            "evaluate",
            (
            "eval_apply_reversed"
                if not config["learn_apply"]["fragmentation"]
                else "eval_apply_reversed_frag"
            ),
            "seq_annotation_scores_{nb}.csv.gz",
        ),
    message:
        "Using Apply to test reversed (decoy) sequences in {input.data}. Output written to {output.apply}."
    script:
        resource_path("snekmer", "scripts", "learn_eval_apply_reverse_seqs.py")


rule reverse_decoy_evaluations:
    input:
        eval_apply_data=expand(
            join(
                out_dir,
                "evaluate",
                (
                "eval_apply_reversed"
                    if not config["learn_apply"]["fragmentation"]
                    else "eval_apply_reversed_frag"
                ),
                "seq_annotation_scores_{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        base_family_checkpoint=expand("{bc}", bc=base_family_checkpoint),
    output:
        family_stats=join(out_dir, "eval_conf", "family_summary_stats.csv"),
        checkpoint=join(out_dir, "eval_conf", "family_stats_checkpoint.csv"),
    message:
        "Evaluating reverse decoy sequences and writing family stats to {output.family_stats}."
    script:
        resource_path("snekmer", "scripts", "learn_reverse_decoy_evaluations.py")


rule eval_apply_sequences:
    input:
        data=join(
            out_dir,
            "vector",
            (
            "vector"
                if not config["learn_apply"]["fragmentation"]
                else "vector_frag"
            ),
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annot_files),
        compare_associations=join(out_dir, "learn", "kmer_counts_total.csv"),
    output:
        apply=join(
            out_dir,
            "evaluate",
            (
            "eval_apply_sequences"
                if not config["learn_apply"]["fragmentation"]
                else "eval_apply_frag"
            ),
            "seq_annotation_scores_{nb}.csv.gz",
        ),
    message:
        "Using Apply to test normal sequences in {input.data}. Output written to {output.apply}."
    script:
        resource_path("snekmer", "scripts", "learn_eval_apply_sequences.py")


rule evaluate:
    input:
        eval_apply_data=expand(
            join(
                out_dir,
                "evaluate",
                (
                "eval_apply_sequences"
                    if not config["learn_apply"]["fragmentation"]
                    else "eval_apply_frag"
                ),
                "seq_annotation_scores_{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        base_confidence=expand("{bc}", bc=base_confidence),
        reverse_decoy_stats=join(out_dir, "eval_conf", "family_summary_stats.csv"),
    output:
        eval_glob=join(out_dir, "eval_conf", "global_confidence_scores.csv"),
    message:
        "Calculating global confidence scores based on Apply results. Output written to {output.eval_glob}."
    params:
        modifier=config["learn_apply"]["conf_weight_modifier"],
    script:
        resource_path("snekmer", "scripts", "learn_evaluate_sequences.py")


rule copy_results_for_apply:
    input:
        kmer_counts_total=join(out_dir, "learn", "kmer_counts_total.csv"),
        family_stats=join(out_dir, "eval_conf", "family_summary_stats.csv"),
        global_conf_scores=join(out_dir, "eval_conf", "global_confidence_scores.csv"),
    output:
        kmer_counts_total=join("apply_inputs", "counts", "kmer_counts_total.csv"),
        family_stats=join("apply_inputs", "stats", "family_summary_stats.csv"),
        global_conf_scores=join(
            "apply_inputs", "confidence", "global_confidence_scores.csv"
        ),
    message:
        "Copying files needed for downstream apply workflow to local apply_inputs directory."
    run:
        target_dir = os.path.join("apply_inputs")
        os.makedirs(target_dir, exist_ok=True)
        shutil.copy(str(input.kmer_counts_total), str(output.kmer_counts_total))
        shutil.copy(str(input.family_stats), str(output.family_stats))
        shutil.copy(str(input.global_conf_scores), str(output.global_conf_scores))


rule learn_report:
    input:
        per_sample=expand(join(out_dir, "learn", "kmer_counts_{nb}.csv"), nb=FAS),
        total_counts=join(out_dir, "learn", "kmer_counts_total.csv"),
        family_stats=join(out_dir, "eval_conf", "family_summary_stats.csv"),
        family_checkpoint=join(out_dir, "eval_conf", "family_stats_checkpoint.csv"),
        global_conf=join(out_dir, "eval_conf", "global_confidence_scores.csv"),
    output:
        report=join(out_dir, "Snekmer_Learn_Report.html"),
    message:
        "Generating full Snekmer Learn Report at {output.report}"
    params:
        fragmentation=config["learn_apply"]["fragmentation"],
    script:
        resource_path("snekmer", "scripts", "learn_report.py")

