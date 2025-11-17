# snekmer/scripts/learn_report.py

import os
from glob import glob
from os.path import dirname, relpath, join
from datetime import datetime

import pandas as pd
import snekmer as skm  # for skm.report.create_report_many_csvs


# Snakemake provides a global `snakemake` object here.
# Inputs / outputs / params are accessed as snakemake.input, snakemake.output, snakemake.params

def main():
    frag_enabled = snakemake.params.fragmentation

    # total_counts is something like: <out_dir>/learn/kmer_counts_total.csv
    # report_src is the top-level <out_dir> directory
    total_counts = snakemake.input.total_counts
    report_src = dirname(dirname(total_counts))

    # ------------------------------------------------------------------
    # Collect all relevant files (paths relative to report_src)
    # ------------------------------------------------------------------
    vector_files = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob(join(report_src, "vector", "vector", "*.npz"))
        ]
    )
    learn_files = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob(join(report_src, "learn", "*.csv"))
        ]
    )
    eval_seq = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob(join(report_src, "evaluate", "eval_apply_sequences", "*.csv.gz"))
        ]
    )
    eval_rev = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob(join(report_src, "evaluate", "eval_apply_reversed", "*.csv.gz"))
        ]
    )

    fam_stats = relpath(snakemake.input.family_stats, report_src).replace(os.sep, "/")
    glob_conf = relpath(snakemake.input.global_conf, report_src).replace(os.sep, "/")
    fam_checkpoint = relpath(snakemake.input.family_checkpoint, report_src).replace(
        os.sep, "/"
    )

    kmer_obj = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob(join(report_src, "kmerize", "kmer", "*.kmers"))
        ]
    )

    # apply_inputs is rooted at the working dir, not report_src
    apply_files = sorted(
        [
            relpath(p, report_src).replace(os.sep, "/")
            for p in glob("apply_inputs/**/*.csv", recursive=True)
        ]
    )

    frag_fastas = []
    vector_frag = []
    kmer_frag = []
    if frag_enabled:
        frag_fastas = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob(join(report_src, "fragmented", "*.fasta"))
            ]
        )
        vector_frag = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob(join(report_src, "vector", "vector_frag", "*.npz"))
            ]
        )
        kmer_frag = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob(join(report_src, "kmerize", "kmer_frag", "*.kmers"))
            ]
        )

    # ------------------------------------------------------------------
    # File metadata (size + mtime)
    # ------------------------------------------------------------------
    all_file_lists = [
        vector_files,
        learn_files,
        eval_seq,
        eval_rev,
        [fam_stats, fam_checkpoint, glob_conf],
        kmer_obj,
        apply_files,
    ]
    if frag_enabled:
        all_file_lists += [frag_fastas, vector_frag, kmer_frag]

    file_info = {}
    for lst in all_file_lists:
        for relp in lst:
            if relp.startswith("apply_inputs/"):
                full = relp  # already relative to CWD
            else:
                full = join(report_src, relp)
            size = round(os.path.getsize(full) / 1024, 1)
            mtime = datetime.fromtimestamp(os.path.getmtime(full)).strftime(
                "%Y-%m-%d %H:%M"
            )
            file_info[relp] = {"size": size, "mtime": mtime}

    # ------------------------------------------------------------------
    # Overview text (same as in your run: block)
    # ------------------------------------------------------------------
    per_sample_files = list(snakemake.input.per_sample)

    total_df = pd.read_csv(snakemake.input.total_counts)
    fam_df = pd.read_csv(snakemake.input.family_stats)

    overview = (
        f"The Learn pipeline processed {len(per_sample_files)} FASTA files. "
        f"K-mer counts files were written under the learn directory, and then merged into "
        f"`learn/kmer_counts_total.csv` ({total_df.shape[0]} x "
        f"{total_df.shape[1]} matrix). "
        f"We identified {fam_df.shape[0]} families to generate "
        "family-specific thresholds, and produced a global confidence score mapping based on cosine similarity scores. "
        "The required inputs for the Apply pipeline have been copied to the apply_inputs directory for downstream workflows."
    )

    # ------------------------------------------------------------------
    # Section descriptions (unchanged)
    # ------------------------------------------------------------------
    desc = {
        "kmerize": (
            "<p><strong>K-mer Extraction:</strong> "
            "Parse each FASTA to enumerate all unique k-mers and save them as .kmers objects.</p>"
            "<ul>"
            "<li><strong>K-mer objects:</strong> <code>kmerize/kmer/*.kmers</code> store the k-mer files.</li>"
            "</ul>"
        ),
        "vector": (
            "<p><strong>Vectorization:</strong> "
            "Convert each sequence into a binary k-mer presence/absence vector (.npz) for downstream processing.</p>"
            "<ul>"
            "<li><strong>Sequence vectors (.npz):</strong> Binary k-mer presence/absence for each FASTA.</li>"
            "</ul>"
        ),
        "fragmentation": (
            "<p><strong>Fragmentation:</strong> "
            "When enabled, each input FASTA is broken into fragments based on config params including: "
            "<code>version</code>, <code>frag_length</code>, <code>location</code>, and <code>seed</code> settings.</p>"
            "<ul>"
            "<li><strong>Fragments:</strong> <code>fragmented/{nb}.fasta</code> files.</li>"
            "</ul>"
        ),
        "vector_frag": (
            "<p><strong>Vectorization of Fragments:</strong> "
            "Convert each fragment into a binary k-mer presence/absence vector (.npz).</p>"
            "<ul>"
            "<li><strong>Fragment Vectors:</strong> <code>vector/vector_frag/{nb}.npz</code></li>"
            "<li><strong>Fragment KMERS:</strong> <code>kmerize/kmer_frag/{nb}.kmers</code></li>"
            "</ul>"
        ),
        "learn": (
            "<p><strong>Learn (Count & Merge):</strong> "
            "Count k-mer occurrences per sequence, then consolidate into a global k-mer association matrix.</p>"
            "<ul>"
            "<li><strong>K-mer Counts:</strong> One CSV per fasta file produced by <code>generate_kmer_counts()</code>.</li>"
            "<li><strong>Merged Matrix:</strong> <code>learn/kmer_counts_total.csv</code> sums all counts into one table.</li>"
            "</ul>"
        ),
        "evaluate": (
            "<p><strong>Evaluation:</strong> "
            "Compute cosine similarities between sequence k-mer counts and the learned matrix for both original and reversed inputs.</p>"
            "<ul>"
            "<li><strong>Apply on decoys:</strong> <code>evaluate/eval_apply_reversed/</code> holds reversed (decoy) scores.</li>"
            "<li><strong>Apply on original:</strong> <code>evaluate/eval_apply_sequences/</code> holds real sequence scores.</li>"
            "</ul>"
        ),
        "eval_conf": (
            "<p><strong>Thresholding & Confidence:</strong> "
            "Use reversed (decoy) scores to derive per-family noise thresholds, then compute a global confidence mapping based on cosine similarity scores.</p>"
            "<ul>"
            "<li><strong>Family Thresholds:</strong> <code>eval_conf/family_summary_stats.csv</code> with count, sum, sumSqr, min/max, and percentiles (reservoir-sampled).</li>"
            "<li><strong>Checkpoint:</strong> <code>eval_conf/family_stats_checkpoint.csv</code> when additional data is added in Learn and merging stats.</li>"
            "<li><strong>Global-Confidence Map:</strong> <code>eval_conf/global_confidence_scores.csv</code> showing confidence vs. cosine score difference.</li>"
            "</ul>"
        ),
        "apply_inputs": (
            "<p><strong>Staged Outputs:</strong> "
            "Key files are copied here for the downstream “apply” pipeline.</p>"
            "<ul>"
            "<li><strong>Counts:</strong> <code>apply_inputs/counts/kmer_counts_total.csv</code></li>"
            "<li><strong>Stats:</strong> <code>apply_inputs/stats/family_summary_stats.csv</code></li>"
            "<li><strong>Confidence:</strong> <code>apply_inputs/confidence/global_confidence_scores.csv</code></li>"
            "</ul>"
        ),
    }

    learn_counts = [f for f in learn_files if "total" not in f]
    learn_total = [f for f in learn_files if "total" in f]

    learn_vars = dict(
        page_title="Snekmer Learn Report",
        title="Snekmer Learn Pipeline Results",
        overview_text=overview,
        section_desc=desc,
        vector_files=vector_files,
        frag_fastas=frag_fastas,
        vector_frag=vector_frag,
        kmer_frag=kmer_frag,
        learn_counts=learn_counts,
        learn_total=learn_total,
        eval_seq=eval_seq,
        eval_rev=eval_rev,
        fam_stats=fam_stats,
        fam_checkpoint=fam_checkpoint,
        glob_conf=glob_conf,
        kmer_obj=kmer_obj,
        apply_inputs=apply_files,
        file_info=file_info,
    )

    # Render the report
    skm.report.create_report_many_csvs(
        report_src,
        learn_vars,
        "learn",
        snakemake.output.report,
    )


if __name__ == "__main__":
    main()
