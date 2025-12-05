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
