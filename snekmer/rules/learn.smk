# force snakemake v6.0(required for modules)
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
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
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
        shutil.copy(input.kmer_counts_total, output.kmer_counts_total)
        shutil.copy(input.family_stats, output.family_stats)
        shutil.copy(input.global_conf_scores, output.global_conf_scores)


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
    run:
        import os, glob, pandas as pd
        from os.path import dirname, relpath, join
        from datetime import datetime

        frag_enabled = params.fragmentation
        report_src = dirname(dirname(input.total_counts))

        vector_files = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob(join(report_src, "vector", "vector", "*.npz"))
            ]
        )
        learn_files = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob(join(report_src, "learn", "*.csv"))
            ]
        )
        eval_seq = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob(
                    join(report_src, "evaluate", "eval_apply_sequences", "*.csv.gz")
                )
            ]
        )
        eval_rev = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob(
                    join(report_src, "evaluate", "eval_apply_reversed", "*.csv.gz")
                )
            ]
        )

        fam_stats = relpath(input.family_stats, report_src).replace(os.sep, "/")
        glob_conf = relpath(input.global_conf, report_src).replace(os.sep, "/")
        fam_checkpoint = relpath(input.family_checkpoint, report_src).replace(
            os.sep, "/"
        )

        kmer_obj = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob(join(report_src, "kmerize", "kmer", "*.kmers"))
            ]
        )

        apply_files = sorted(
            [
                relpath(p, report_src).replace(os.sep, "/")
                for p in glob.glob("apply_inputs/**/*.csv", recursive=True)
            ]
        )

        frag_fastas = []
        vector_frag = []
        kmer_frag = []
        if frag_enabled:
            frag_fastas = sorted(
                [
                    relpath(p, report_src).replace(os.sep, "/")
                    for p in glob.glob(join(report_src, "fragmented", "*.fasta"))
                ]
            )
            vector_frag = sorted(
                [
                    relpath(p, report_src).replace(os.sep, "/")
                    for p in glob.glob(
                        join(report_src, "vector", "vector_frag", "*.npz")
                    )
                ]
            )
            kmer_frag = sorted(
                [
                    relpath(p, report_src).replace(os.sep, "/")
                    for p in glob.glob(
                        join(report_src, "kmerize", "kmer_frag", "*.kmers")
                    )
                ]
            )

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
                    full = relp
                else:
                    full = join(report_src, relp)
                size = round(os.path.getsize(full) / 1024, 1)
                mtime = datetime.fromtimestamp(os.path.getmtime(full)).strftime(
                    "%Y-%m-%d %H:%M"
                )
                file_info[relp] = {"size": size, "mtime": mtime}

        overview = (
            f"The Learn pipeline processed {len(input.per_sample)} FASTA files. "
            f"K-mer counts files were written under the learn directory, and then merged into "
            f"`learn/kmer_counts_total.csv` ({pd.read_csv(input.total_counts).shape[0]} x "
            f"{pd.read_csv(input.total_counts).shape[1]} matrix). "
            f"We identified {pd.read_csv(input.family_stats).shape[0]} families to generate "
            "family-specific thresholds, and produced a global confidence score mapping based on cosine similarity scores. "
            "The required inputs for the Apply pipeline have been copied to the apply_inputs directory for downstream workflows."
        )
        desc = {
            "kmerize": (
            "<p><strong>K-mer Extraction:</strong> "
                "Parse each FASTA to enumerate all unique k-mers and save them as .kmers objects.</p>"
                "<ul>"
                "<li><strong>K-mer objects:</strong> <code>kmerize/kmer/*.kmers</code> store the k-mer definitions.</li>"
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
                "When enabled, each input FASTA is broken into fragments per your "
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

        skm.report.create_report_many_csvs(
            report_src, learn_vars, "learn", output.report
        )
