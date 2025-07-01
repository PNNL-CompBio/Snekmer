# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version

min_version("6.0")


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
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
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

threshold_type = config["learnapp"]["threshold"]
selection_type = config["learnapp"]["selection"]

apply_cfg = config.get("learnapp", {})
concat_results = apply_cfg.get("apply_output", "snekmer_results.csv")
extra_all = [concat_results] if concat_results is not None else []


def resource_path(package: str, *parts) -> str:
    """
    Re-create pkg_resources.resource_filename()
    but using importlib.resources.
    """
    return str(files(package).joinpath(*parts))


wildcard_constraints:
    dataset=FAS,
    FAS=FAS,


rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        *[
            (
                expand(
                    join(out_dir, "apply", "seq-annotation-scores-{nb}.csv"), nb=FAS
                )
                if config["learnapp"]["save_apply_associations"]
                else []
            )
        ],
        *extra_all,
        expand(join(out_dir, "apply", "kmer-summary-{nb}.csv"), nb=FAS),
        join(out_dir, "Snekmer_Apply_Report.html"),


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
        selection_type=config["learnapp"]["selection"],
        threshold_type=config["learnapp"]["threshold"],
    output:
        seq_ann=(
            expand(join(out_dir, "apply", "seq-annotation-scores-{nb}.csv"), nb=FAS)
            if config["learnapp"]["save_apply_associations"]
            else []
        ),
        kmer_summary=join(out_dir, "apply", "kmer-summary-{nb}.csv"),
    message:
        "Running Snekmer Apply on {input.data}. Output written to {output.kmer_summary}."
    script:
        resource_path("snekmer", "scripts", "apply.py")


if concat_results is not None:

    rule concat_kmer_summary:
        """NEW: concatenate all per-sample k-mer summaries into one CSV"""
        input:
            expand(join(out_dir, "apply", "kmer-summary-{nb}.csv"), nb=FAS),
        output:
            concat_results,
        message:
            "Writing consolidated k-mer summary to {output}"
        shell:
            """                                       
            # write header from the first file
            head -n 1 {input[0]} > {output}
            # append rows from every remaining file (skip header lines)
            for f in {input}; do
                tail -n +2 "$f" >> {output}
            done
            """


rule apply_report:
    input:
        files=expand(join(out_dir, "apply", "kmer-summary-{f}.csv"), f=FAS),
    output:
        join(out_dir, "Snekmer_Apply_Report.html"),
    message:
        "Basic HTML report written to {output}"
    run:
        file_dir = dirname(dirname(input.files[0]))

        # apply
        apply_vars = dict(
            page_title="Snekmer Apply Report",
            title="Snekmer Apply Results",
            text="See the below links to access Snekmer Apply results.",
        )

        skm.report.create_report_many_csvs(file_dir, apply_vars, "apply", output[0])
