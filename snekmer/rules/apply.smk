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


wildcard_constraints:
    dataset=FAS,
    FAS=FAS,


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
        # seq_ann=(
        #     expand(join(out_dir, "apply", "seq_annotation_scores_{nb}.csv"), nb=FAS)
        #     if config["learn_apply"]["save_apply_associations"]
        #     else []
        # ),
        seq_ann=join(out_dir, "apply", "seq_annotation_scores_{nb}.csv")
        if config["learn_apply"]["save_apply_associations"]
        else [],
        kmer_summary=join(out_dir, "apply", "kmer_summary_{nb}.csv"),
    message:
        "Running Snekmer Apply on {input.data}. Output written to {output.kmer_summary}."
    script:
        resource_path("snekmer", "scripts", "apply.py")


if concat_results is not None:

    rule concat_kmer_summary:
        """NEW: concatenate all sample k-mer summaries into one CSV"""
        input:
            expand(join(out_dir, "apply", "kmer_summary_{nb}.csv"), nb=FAS),
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
        seq_scores=(
            expand(join(out_dir, "apply", "seq_annotation_scores_{nb}.csv"), nb=FAS)
            if config["learn_apply"]["save_apply_associations"]
            else []
        ),
        kmer_sum=expand(join(out_dir, "apply", "kmer_summary_{nb}.csv"), nb=FAS),
        concat=config.get("learn_apply", {}).get("apply_output", None),
    output:
        report=join(out_dir, "Snekmer_Apply_Report.html"),
    message:
        "Generating full Snekmer Apply Report at {output.report}"
    script:
        resource_path("snekmer", "scripts", "apply_report.py")


    # run:
        # import os, glob, pandas as pd
        # from os.path import dirname, relpath, join
        # from datetime import datetime

        # report_src = dirname(dirname(input.kmer_sum[0]))
        # seq_scores_rel = (
        #     sorted(
        #         [
        #             relpath(p, report_src).replace(os.sep, "/")
        #             for p in glob.glob(
        #                 join(report_src, "apply", "seq_annotation_scores_*.csv")
        #             )
        #         ]
        #     )
        #     if input.seq_scores
        #     else []
        # )

        # kmer_sum_rel = sorted(
        #     [
        #         relpath(p, report_src).replace(os.sep, "/")
        #         for p in glob.glob(join(report_src, "apply", "kmer_summary_*.csv"))
        #     ]
        # )

        # concat_rel = None
        # if input.concat:
        #     concat_rel = relpath(input.concat, report_src).replace(os.sep, "/")

        # vector_rel = sorted(
        #     [
        #         relpath(p, report_src).replace(os.sep, "/")
        #         for p in glob.glob(join(report_src, "vector", "*.npz"))
        #     ]
        # )

        # kmerize_rel = sorted(
        #     [
        #         relpath(p, report_src).replace(os.sep, "/")
        #         for p in glob.glob(join(report_src, "kmerize", "*.kmers"))
        #     ]
        # )

        # overview = (
        #     "The Apply pipeline computes cosine similarities between the learned k-mer association matrix and each newly provided sequence, "
        #     f"applies the <strong>{config['learn_apply']['selection']}</strong> selection method using family-specific thresholds "
        #     "from <code>family_summary_stats.csv</code> and the global confidence mapping from <code>global-confidence-scores.csv</code>, "
        #     f"and then writes <strong>{len(kmer_sum_rel)}</strong> summary CSV files under <code>apply/</code>."
        #     + (
        #         f" A consolidated results table is also available at <code>{concat_rel}</code>."
        #         if concat_rel
        #         else ""
        #     )
        # )

        # desc = {
        #     "kmerize": "<p><strong>K-mer Extraction:</strong> Each fasta file was parsed into a `.kmers` object.</p>",
        #     "vector": "<p><strong>Vectorization:</strong> Each sequence encoded as a binary k-mer `.npz` vector.</p>",
        #     "scores": (
        #         "<p><strong>Annotation Scores:</strong> "
        #         "For each fasta, cosine similarities against the learned k-mer matrix "
        #         "were computed.  These CSVs include all raw score comparisons.</p>"
        #     ),
        #     "concat": "<p><strong>Consolidated Summary:</strong> All kmer summaries merged into one CSV.</p>",
        # }

        # selection_method = config["learn_apply"]["selection"]
        # threshold_type = config["learn_apply"]["threshold"]
        # weight_top = config["learn_apply"].get("weight_top")
        # weight_distance = config["learn_apply"].get("weight_distance")

        # threshold_blurb = (
        #     f"<p><strong>Threshold:</strong> using <strong>{threshold_type}</strong> values from "
        #     "<code>family_summary_stats.csv</code> as noise cutoffs.</p>"
        # )

        # weight_blurb = ""
        # if selection_method == "combined_distance":
        #     weight_blurb = f"<p><strong>Weights:</strong> cosine similarity score x <strong>{weight_top}</strong> + distance from threshold x <strong>{weight_distance}</strong></p>"

        # sel_descriptions = {
        #     "top_hit": (
        #         "<p><strong>Selection Method: Top Hit</strong> selects the family with the highest cosine similarity. "
        #         "</p>"
        #     ),
        #     "greatest_distance": (
        #         "<p><strong>Selection Method: Greatest Distance</strong> computes (score - threshold) per family "
        #         "and picks the one with the largest positive difference above its threshold, "
        #         "favoring truly above-noise hits.</p>"
        #     ),
        #     "combined_distance": (
        #         "<p><strong>Selection Method: Combined Distance</strong> computes a weighted sum of raw score and "
        #         "(score - threshold), then picks the family maximizing that metric.</p>"
        #     ),
        # }

        # method_blurb = sel_descriptions.get(
        #     selection_method,
        #     f"<p><strong>Selection: {selection_method}</strong> - custom method.</p>",
        # )

        # desc["summary"] = (
        #     method_blurb + threshold_blurb + weight_blurb + "\n<ul>\n"
        #     "  <li><strong>Prediction:</strong> the family with the top cosine-similarity score.</li>\n"
        #     "  <li><strong>Score Delta (Δ):</strong> the difference between the top-rank and second-rank scores "
        #     "(Δ = top₁ - top₂), indicating how clear the best match was.</li>\n"
        #     "  <li><strong>Confidence:</strong> a global confidence value for that Δ, "
        #     "computed by comparing true vs. decoy T/F rates across all samples and interpolating over Δ.</li>\n"
        #     "</ul>\n"
        # )

        # file_info = {}
        # for group in [
        #     kmerize_rel,
        #     vector_rel,
        #     seq_scores_rel,
        #     kmer_sum_rel,
        #     [concat_rel] if concat_rel else [],
        # ]:
        #     for f in group:
        #         full = join(report_src, f)
        #         size = round(os.path.getsize(full) / 1024, 1)
        #         mtime = datetime.fromtimestamp(os.path.getmtime(full)).strftime(
        #             "%Y-%m-%d %H:%M"
        #         )
        #         file_info[f] = {"size": size, "mtime": mtime}

        # apply_vars = dict(
        #     page_title="Snekmer Apply Report",
        #     title="Snekmer Apply Pipeline Results",
        #     overview_text=overview,
        #     section_desc=desc,
        #     kmerize_rel=kmerize_rel,
        #     vector_rel=vector_rel,
        #     seq_scores_rel=seq_scores_rel,
        #     kmer_sum_rel=kmer_sum_rel,
        #     concat_rel=concat_rel,
        #     file_info=file_info,
        # )

        # skm.report.create_report_many_csvs(
        #     report_src, apply_vars, "apply", output.report
        # )
