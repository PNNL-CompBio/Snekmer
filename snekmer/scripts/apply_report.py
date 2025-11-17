# snekmer/scripts/apply_report.py

import os
import glob
from os.path import dirname, relpath, join
from datetime import datetime

import pandas as pd
import snekmer as skm  # for skm.report.create_report_many_csvs

# Use config
config = snakemake.config


def main():
    # --------------------------------------------
    # Derive report root from kmer_sum location
    # --------------------------------------------
    kmer_sum_inputs = list(snakemake.input.kmer_sum)
    if not kmer_sum_inputs:
        raise RuntimeError("apply_report: no kmer_sum inputs found.")
    report_src = dirname(dirname(kmer_sum_inputs[0]))

    # --------------------------------------------
    # Collect file lists as paths relative to report_src
    # --------------------------------------------
    # seq_scores are optional depending on config["learn_apply"]["save_apply_associations"]
    seq_scores_input = list(snakemake.input.seq_scores)
    if seq_scores_input:
        seq_scores_rel = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob(
                    join(report_src, "apply", "seq_annotation_scores_*.csv")
                )
            ]
        )
    else:
        seq_scores_rel = []

    kmer_sum_rel = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob.glob(join(report_src, "apply", "kmer_summary_*.csv"))
        ]
    )

    concat_rel = None
    concat_input = getattr(snakemake.input, "concat", None)
    if concat_input:
        concat_rel = relpath(concat_input, report_src).replace(os.sep, "/")

    vector_rel = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob.glob(join(report_src, "vector", "*.npz"))
        ]
    )

    kmerize_rel = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob.glob(join(report_src, "kmerize", "*.kmers"))
        ]
    )

    # --------------------------------------------
    # Overview text
    # --------------------------------------------
    selection_method = config["learn_apply"]["selection"]
    overview = (
        "The Apply pipeline computes cosine similarities between the learned k-mer association matrix and each newly provided sequence, "
        f"applies the <strong>{selection_method}</strong> selection method using family-specific thresholds "
        "from <code>family_summary_stats.csv</code> and the global confidence mapping from <code>global-confidence-scores.csv</code>, "
        f"and then writes <strong>{len(kmer_sum_rel)}</strong> summary CSV files under <code>apply/</code>."
        + (
            f" A consolidated results table is also available at <code>{concat_rel}</code>."
            if concat_rel
            else ""
        )
    )

    # --------------------------------------------
    # Method descriptions (unchanged from your rule)
    # --------------------------------------------
    desc = {
        "kmerize": "<p><strong>K-mer Extraction:</strong> Each fasta file was parsed into a `.kmers` object.</p>",
        "vector": "<p><strong>Vectorization:</strong> Each sequence encoded as a binary k-mer `.npz` vector.</p>",
        "scores": (
            "<p><strong>Annotation Scores:</strong> "
            "For each fasta, cosine similarities against the learned k-mer matrix "
            "were computed.  These CSVs include all raw score comparisons.</p>"
        ),
        "concat": "<p><strong>Consolidated Summary:</strong> All kmer summaries merged into one CSV.</p>",
    }

    threshold_type = config["learn_apply"]["threshold"]
    weight_top = config["learn_apply"].get("weight_top")
    weight_distance = config["learn_apply"].get("weight_distance")

    threshold_blurb = (
        f"<p><strong>Threshold:</strong> using <strong>{threshold_type}</strong> values from "
        "<code>family_summary_stats.csv</code> as noise cutoffs.</p>"
    )

    weight_blurb = ""
    if selection_method == "combined_distance":
        weight_blurb = (
            f"<p><strong>Weights:</strong> cosine similarity score x <strong>{weight_top}</strong> "
            f"+ distance from threshold x <strong>{weight_distance}</strong></p>"
        )

    sel_descriptions = {
        "top_hit": (
            "<p><strong>Selection Method: Top Hit</strong> selects the family with the highest cosine similarity. "
            "</p>"
        ),
        "greatest_distance": (
            "<p><strong>Selection Method: Greatest Distance</strong> computes (score - threshold) per family "
            "and picks the one with the largest positive difference above its threshold, "
            "favoring truly above-noise hits.</p>"
        ),
        "combined_distance": (
            "<p><strong>Selection Method: Combined Distance</strong> computes a weighted sum of raw score and "
            "(score - threshold), then picks the family maximizing that metric.</p>"
        ),
    }

    method_blurb = sel_descriptions.get(
        selection_method,
        f"<p><strong>Selection: {selection_method}</strong> - custom method.</p>",
    )

    desc["summary"] = (
        method_blurb + threshold_blurb + weight_blurb + "\n<ul>\n"
        "  <li><strong>Prediction:</strong> the family with the top cosine-similarity score.</li>\n"
        "  <li><strong>Score Delta (Δ):</strong> the difference between the top-rank and second-rank scores "
        "(Δ = top₁ - top₂), indicating how clear the best match was.</li>\n"
        "  <li><strong>Confidence:</strong> a global confidence value for the corresponding Δ. "
        "Computed by comparing true vs. decoy T/F rates across all samples and interpolating over Δ.</li>\n"
        "</ul>\n"
    )

    # --------------------------------------------
    # File metadata (size + mtime)
    # --------------------------------------------
    file_info = {}
    for group in [
        kmerize_rel,
        vector_rel,
        seq_scores_rel,
        kmer_sum_rel,
        [concat_rel] if concat_rel else [],
    ]:
        for f in group:
            if not f:
                continue
            full = join(report_src, f)
            size = round(os.path.getsize(full) / 1024, 1)
            mtime = datetime.fromtimestamp(os.path.getmtime(full)).strftime(
                "%Y-%m-%d %H:%M"
            )
            file_info[f] = {"size": size, "mtime": mtime}

    # --------------------------------------------
    # Variables passed to the HTML template
    # --------------------------------------------
    apply_vars = dict(
        page_title="Snekmer Apply Report",
        title="Snekmer Apply Pipeline Results",
        overview_text=overview,
        section_desc=desc,
        kmerize_rel=kmerize_rel,
        vector_rel=vector_rel,
        seq_scores_rel=seq_scores_rel,
        kmer_sum_rel=kmer_sum_rel,
        concat_rel=concat_rel,
        file_info=file_info,
    )

    # --------------------------------------------
    # Render report
    # --------------------------------------------
    skm.report.create_report_many_csvs(
        report_src, apply_vars, "apply", snakemake.output.report
    )


if __name__ == "__main__":
    main()
