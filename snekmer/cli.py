"""cli: Command line interface for Snekmer.

author: @christinehc, @biodataganache, @snakemake

"""
# imports
import argparse
import os
import sys
import json
import shlex
import subprocess
import textwrap
import tempfile

from multiprocessing import cpu_count
from importlib.resources import files, as_file
from snekmer import __version__


# Actual alphabets (id -> (name, size, description))
ALPHABETS = {
    0: ("hydro", 2, "2-value hydrophobicity alphabet"),
    1: ("standard", 7, "“Standard” reduction alphabet"),
    2: ("solvacc", 3, "Solvent accessibility alphabet"),
    3: (
        "hydrocharge",
        3,
        "2-value hydrophobicity with charged residues as a third category",
    ),
    4: (
        "hydrostruct",
        3,
        "2-value hydrophobicity with structural-breakers as a third category",
    ),
    5: ("miqs", 10, "MIQS alphabet3"),
    None: ("None", 20, "No reduced alphabet"),
}


# Mode registry — single source of truth for names and descriptions
_MODES = {
    "cluster":          "Unsupervised clustering workflow.",
    "model":            "Train supervised models + cross-validation reports.",
    "search":           "Score sequences against trained models.",
    "learn":            "Build annotation-associated k-mer distributions + confidence evaluation.",
    "apply":            "Predict annotations using outputs from learn.",
    "easy-learn-apply": "Guided front-end that runs learn then apply end-to-end.",
}


def _mode_table() -> str:
    """Return a formatted table of modes and their descriptions."""
    width = max(len(m) for m in _MODES)
    lines = []
    for mode, desc in _MODES.items():
        lines.append(f"  {mode:<{width}}  {desc}")
    return "\n".join(lines)


class _SnekmerMainParser(argparse.ArgumentParser):
    """Main parser with concise, organized help and targeted error messages."""

    def print_help(self, file=None):
        if file is None:
            file = sys.stdout
        print(
            f"Snekmer {__version__} — protein sequence fingerprinting via amino acid reduction.\n"
            "\n"
            "Usage:\n"
            "  snekmer <mode> [options]\n"
            "  snekmer <mode> --help\n"
            "\n"
            "Modes:\n"
            f"{_mode_table()}\n"
            "\n"
            "Global options (accepted by all modes):\n"
            "  --k N           K-mer length (default: 8).\n"
            "  --alphabet N    Reduced alphabet encoding 0-5 or name (default: 2 = solvacc).\n"
            "  --cores N       CPU cores to use (default: all).\n"
            "  --dry-run       Show what would be done without executing.\n"
            "  --configfile    Path(s) to YAML/JSON config file(s).\n"
            "  -v, --version   Print version and exit.\n"
            "  -h, --help      Show this help message and exit.\n"
            "\n"
            "Run 'snekmer <mode> --help' for full options for a specific mode.",
            file=file,
        )

    def error(self, message):
        # Strip the default argparse preamble and show a focused message
        clean = message
        # "argument mode: invalid choice: 'X' (choose from ...)" → "unknown mode: 'X'"
        if "argument mode: invalid choice:" in message:
            bad = message.split("invalid choice:")[1].split("(")[0].strip()
            clean = f"unknown mode: {bad}"
        print(f"snekmer: error: {clean}\n", file=sys.stderr)
        print("Valid modes:", file=sys.stderr)
        print(_mode_table(), file=sys.stderr)
        print(
            "\nRun 'snekmer <mode> --help' for options, or 'snekmer --help' for an overview.",
            file=sys.stderr,
        )
        sys.exit(2)


# ------------------------- argparse helpers -------------------------
class SnekmerHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


class _SnekmerSubParser(argparse.ArgumentParser):
    """Base class for all Snekmer mode subparsers.

    - Defaults to SnekmerHelpFormatter so argument groups show descriptions
      and defaults in --help output.
    - Collapses the usage line to 'usage: snekmer <mode> [options]' so it
      does not expand into a wall of flags.
    - Overrides error() to suppress the usage line and print a concise
      message with a pointer to --help.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", SnekmerHelpFormatter)
        kwargs.setdefault("usage", "%(prog)s [options]")
        super().__init__(*args, **kwargs)

    def error(self, message):
        print(f"{self.prog}: error: {message}\n", file=sys.stderr)
        print(f"Run '{self.prog} --help' for usage.", file=sys.stderr)
        sys.exit(2)


class StoreSeenAction(argparse.Action):
    """
    Like store/store_const, but also records whether the user explicitly provided the option.
    This lets us:
      - show defaults in help,
      - avoid overriding config.yaml unless the user actually passed a flag,
      - still apply defaults when running without any config file.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        # NOTE: argparse passes [] (not None) for nargs=0 flags.
        if (values is None or values == []) and getattr(self, "const", None) is not None:
            val = self.const
        else:
            val = values
        setattr(namespace, self.dest, val)
        setattr(namespace, f"{self.dest}_seen", True)


def _alphabet_options_block() -> str:
    lines = ["Alphabets (k-mer recoding):"]
    for k in [0, 1, 2, 3, 4, 5]:
        name, size, desc = ALPHABETS[k]
        lines.append(f"  {k}: {name} (size {size}) — {desc}")
    name, size, desc = ALPHABETS[None]
    lines.append(f"  None: {name} (size {size}) — {desc}")
    lines.append("")
    lines.append(
        "You may pass either an integer (0–5) or the alphabet name (e.g. 'hydro'), or 'None'."
    )
    return "\n".join(lines)


def _main_description() -> str:
    return textwrap.dedent(
        """\
Snekmer: A scalable pipeline for protein sequence fingerprinting using amino acid reduction (AAR).

Modes:
  cluster           Unsupervised clustering workflow.
  model             Train supervised models + cross-validation reports.
  search            Score sequences against trained models.
  learn             Build annotation-associated k-mer distributions + confidence evaluation.
  apply             Predict annotations using outputs from learn.
  easy-learn-apply  Guided front-end that runs learn then apply end-to-end.

General usage:
  snekmer <mode> [snakemake arguments] [snekmer parameter overrides]

Important:
  The "Snakemake arguments" below are passed through to Snakemake (they are not Snekmer-specific).
  Snekmer parameters can be provided via config.yaml / --configfile, or overridden via the
  Snekmer parameter flags shown in the section(s) relevant to your mode.

Defaults:
  The defaults shown for Snekmer parameter flags match the template config.yaml defaults.
  These defaults are applied automatically only when no config file is in use, or when a
  flag is explicitly provided.

More help:
  Get help for any subcommand with:
    snekmer <mode> -h

Config precedence:
  1) Default configfile (auto): ./config.yaml (or <DIR>/config.yaml with -d/--directory)
  2) Any explicit --configfile PATH values (in order)
  3) Any Snekmer parameter flags you explicitly provide
  4) Any -C/--config KEY=VALUE overrides (highest)

Running without a config file:
  Use --no-default-configfile (optional) and rely on defaults and/or provide overrides.
"""
    ).strip()


def _main_epilog() -> str:
    return textwrap.dedent(
        f"""\
{_alphabet_options_block()}
"""
    ).rstrip()


def _add_bool_flag(group, dest, name, default, help_true, help_false):
    """
    Add paired --<name> / --no-<name> flags to set a boolean config value.
    Default matches config.yaml template.
    """
    mx = group.add_mutually_exclusive_group()
    mx.add_argument(
        f"--{name}",
        dest=dest,
        action=StoreSeenAction,
        nargs=0,
        const=True,
        default=default,
        help=help_true,
    )
    mx.add_argument(
        f"--no-{name}",
        dest=dest,
        action=StoreSeenAction,
        nargs=0,
        const=False,
        default=default,
        help=help_false,
    )


# ------------------------- Snekmer parameter sections (mode-scoped) -------------------------
def _add_required_args(cfg_required: argparse.ArgumentParser) -> None:
    g = cfg_required.add_argument_group(
        "Snekmer parameters (all modes; defaults match config.yaml)"
    )

    g.add_argument(
        "--k",
        dest="cfg_k",
        type=int,
        default=8,
        action=StoreSeenAction,
        metavar="",
        help="K-mer length.",
    )
    g.add_argument(
        "--alphabet",
        dest="cfg_alphabet",
        type=str,
        default="2",
        action=StoreSeenAction,
        metavar="",
        help="Reduced alphabet encoding (0–5, alphabet name, or 'None'). See alphabets list below.",
    )
    g.add_argument(
        "--input-dir",
        dest="cfg_input_dir",
        type=str,
        default="input",
        action=StoreSeenAction,
        metavar="",
        help="Directory containing input sequence files.",
    )
    g.add_argument(
        "--input-file-exts",
        dest="cfg_input_file_exts",
        nargs="+",
        default=["fasta", "fna", "faa", "fa"],
        action=StoreSeenAction,
        metavar="",
        help="File extensions to consider valid input sequence files (space-separated).",
    )
    g.add_argument(
        "--input-file-regex",
        dest="cfg_input_file_regex",
        type=str,
        default=".*",
        action=StoreSeenAction,
        metavar="",
        help="Regular expression for parsing family/annotation identifiers from filenames.",
    )
    _add_bool_flag(
        g,
        dest="cfg_nested_output",
        name="nested-output",
        default=False,
        help_true="Enable nested output directory structure: {save_dir}/{alphabet}/{k}.",
        help_false="Disable nested output directory structure (flat output layout).",
    )


def _add_score_args(cfg_score: argparse.ArgumentParser) -> None:
    g = cfg_score.add_argument_group("Snekmer Model and Search Parameters: scoring")

    _add_bool_flag(
        g,
        dest="cfg_score_scaler",
        name="score-scaler",
        default=True,
        help_true="Enable k-mer score scaling (applies configured scaler to family scores).",
        help_false="Disable k-mer score scaling.",
    )
    g.add_argument(
        "--score-scaler-n",
        dest="cfg_score_scaler_n",
        type=float,
        default=0.25,
        action=StoreSeenAction,
        metavar="",
        help="Scaler keyword argument 'n' (passed to the k-mer scaler).",
    )
    g.add_argument(
        "--score-labels",
        dest="cfg_score_labels",
        type=str,
        default=None,
        action=StoreSeenAction,
        metavar="",
        help="If None, uses default k-mer label set for scaler. Otherwise uses provided value (string or JSON).",
    )
    g.add_argument(
        "--score-lname",
        dest="cfg_score_lname",
        type=str,
        default=None,
        action=StoreSeenAction,
        metavar="",
        help='Label name (e.g., "family").',
    )


def _add_cluster_args(cfg_cluster: argparse.ArgumentParser) -> None:
    g = cfg_cluster.add_argument_group("Snekmer Cluster Parameters: clustering")

    g.add_argument(
        "--cluster-method",
        dest="cfg_cluster_method",
        type=str,
        default="agglomerative-jaccard",
        action=StoreSeenAction,
        metavar="",
        help='Clustering method (options include "kmeans", "agglomerative", "correlation", "density", "birch", "optics", "hdbscan").',
    )
    g.add_argument(
        "--cluster-n-clusters",
        dest="cfg_cluster_n_clusters",
        type=str,
        default=None,
        action=StoreSeenAction,
        metavar="",
        help="Number of clusters (int) or 'None' (method-dependent).",
    )
    g.add_argument(
        "--cluster-linkage",
        dest="cfg_cluster_linkage",
        type=str,
        default="average",
        action=StoreSeenAction,
        metavar="",
        help='Linkage method for agglomerative clustering (e.g. "average").',
    )
    g.add_argument(
        "--cluster-distance-threshold",
        dest="cfg_cluster_distance_threshold",
        type=float,
        default=0.92,
        action=StoreSeenAction,
        metavar="",
        help="Distance threshold for agglomerative clustering (method-dependent).",
    )
    _add_bool_flag(
        g,
        dest="cfg_cluster_compute_full_tree",
        name="cluster-compute-full-tree",
        default=True,
        help_true="Compute full tree for hierarchical clustering (agglomerative).",
        help_false="Do not compute full tree for hierarchical clustering (agglomerative).",
    )
    _add_bool_flag(
        g,
        dest="cfg_cluster_cluster_plots",
        name="cluster-plots",
        default=False,
        help_true="Generate plots illustrating clustering results.",
        help_false="Do not generate clustering plots.",
    )
    g.add_argument(
        "--cluster-min-rep",
        dest="cfg_cluster_min_rep",
        type=str,
        default=None,
        action=StoreSeenAction,
        metavar="",
        help="Minimum repetition threshold for kmers (int) or 'None'. Kmers below this are discarded.",
    )
    g.add_argument(
        "--cluster-max-rep",
        dest="cfg_cluster_max_rep",
        type=str,
        default=None,
        action=StoreSeenAction,
        metavar="",
        help="Maximum repetition threshold for kmers (int) or 'None'. Kmers above this are discarded.",
    )
    _add_bool_flag(
        g,
        dest="cfg_cluster_save_matrix",
        name="cluster-save-matrix",
        default=False,
        help_true="Save distance matrices (BSF). Not recommended for large datasets.",
        help_false="Do not save distance matrices (BSF).",
    )
    g.add_argument(
        "--cluster-dist-thresh",
        dest="cfg_cluster_dist_thresh",
        type=int,
        default=100,
        action=StoreSeenAction,
        metavar="",
        help="Distance threshold for BSF matrix.",
    )


def _add_model_args(cfg_model: argparse.ArgumentParser) -> None:
    g = cfg_model.add_argument_group("Snekmer Model Parameters: model training")

    g.add_argument(
        "--model-cv",
        dest="cfg_model_cv",
        type=int,
        default=5,
        action=StoreSeenAction,
        metavar="",
        help="Number of cross-validation folds for model evaluation.",
    )
    g.add_argument(
        "--model-random-state",
        dest="cfg_model_random_state",
        type=str,
        default=None,
        action=StoreSeenAction,
        metavar="",
        help="Random state for model evaluation (int) or 'None'.",
    )


def _add_search_args(cfg_search: argparse.ArgumentParser) -> None:
    g = cfg_search.add_argument_group("Snekmer Search Parameters: search inputs")

    g.add_argument(
        "--model-dir",
        dest="cfg_model_dir",
        type=str,
        default="output/model/",
        action=StoreSeenAction,
        metavar="",
        help="Directory containing model object(s) (.model).",
    )
    g.add_argument(
        "--basis-dir",
        dest="cfg_basis_dir",
        type=str,
        default="output/example-model/",
        action=StoreSeenAction,
        metavar="",
        help="Directory containing k-mer basis set(s) (.kmers).",
    )
    g.add_argument(
        "--score-dir",
        dest="cfg_score_dir",
        type=str,
        default="output/scoring/",
        action=StoreSeenAction,
        metavar="",
        help="Directory containing scoring object(s) (.scorer).",
    )


def _add_learn_apply_args(cfg_learn_apply: argparse.ArgumentParser) -> None:
    g = cfg_learn_apply.add_argument_group("Snekmer Learn and Apply Parameters:")

    _add_bool_flag(
        g,
        dest="cfg_la_save_apply_associations",
        name="save-apply-associations",
        default=False,
        help_true="Save large optional outputs containing all cosine similarity scores (increases storage substantially).",
        help_false="Do not save large optional cosine similarity outputs.",
    )
    g.add_argument(
        "--conf-weight-modifier",
        dest="cfg_la_conf_weight_modifier",
        type=int,
        default=20,
        action=StoreSeenAction,
        metavar="",
        help="Weighting modifier for updating confidence when adding data to an existing k-mer count matrix.",
    )

    # ONLY --fragmentation (default False). No --no-fragmentation.
    g.add_argument(
        "--fragmentation",
        dest="cfg_la_fragmentation",
        action=StoreSeenAction,
        nargs=0,
        const=True,
        default=False,
        help="Enable training-data fragmentation (default False).",
    )

    g.add_argument(
        "--fragment-version",
        dest="cfg_la_version",
        type=str,
        default="absolute",
        action=StoreSeenAction,
        metavar="",
        help="Fragment length interpretation: 'absolute' or 'percent'.",
    )
    g.add_argument(
        "--frag-length",
        dest="cfg_la_frag_length",
        type=int,
        default=50,
        action=StoreSeenAction,
        metavar="",
        help="Fragment length (units depend on --fragment-version).",
    )
    g.add_argument(
        "--min-length",
        dest="cfg_la_min_length",
        type=int,
        default=50,
        action=StoreSeenAction,
        metavar="",
        help="Minimum fragment length to retain; shorter fragments are discarded.",
    )
    g.add_argument(
        "--fragment-location",
        dest="cfg_la_location",
        type=str,
        default="random",
        action=StoreSeenAction,
        metavar="",
        help="Fragment location: 'start', 'end', or 'random'.",
    )
    g.add_argument(
        "--seed",
        dest="cfg_la_seed",
        type=int,
        default=999,
        action=StoreSeenAction,
        metavar="",
        help="Random seed for reproducible fragmentation.",
    )

    g.add_argument(
        "--selection",
        dest="cfg_la_selection",
        choices=["top_hit", "greatest_distance", "combined_distance"],
        default="top_hit",
        action=StoreSeenAction,
        metavar="",
        help="Annotation selection method.",
    )
    g.add_argument(
        "--threshold",
        dest="cfg_la_threshold",
        type=str,
        default="Median",
        action=StoreSeenAction,
        metavar="",
        help="Family-specific threshold used for prediction filtering (e.g. 'Median', 'Mean', '90th Percentile', or 'None').",
    )
    g.add_argument(
        "--weight-top",
        dest="cfg_la_weight_top",
        type=float,
        default=0.7,
        action=StoreSeenAction,
        metavar="",
        help="Weight for 'top_hit' when selection method is 'combined_distance'.",
    )
    g.add_argument(
        "--weight-distance",
        dest="cfg_la_weight_distance",
        type=float,
        default=0.3,
        action=StoreSeenAction,
        metavar="",
        help="Weight for 'greatest_distance' when selection method is 'combined_distance'.",
    )
    g.add_argument(
        "--apply-output",
        dest="cfg_la_apply_output",
        type=str,
        default="snekmer_results.csv",
        action=StoreSeenAction,
        metavar="",
        help="Output filename for apply results in single-file format.",
    )


# ------------------------- easy-learn-apply args -------------------------
def _add_easy_la_args(p: argparse.ArgumentParser) -> None:
    """Add all arguments for the easy-learn-apply subcommand."""

    # -- Input/output --------------------------------------------------------
    io = p.add_argument_group("Input / output")
    io.add_argument(
        "--train",
        metavar="PATH",
        default=None,
        help=(
            "Path to training sequences (FASTA file or directory of FASTA files). "
            "If omitted, the wizard will prompt for it."
        ),
    )
    io.add_argument(
        "--query",
        metavar="PATH",
        default=None,
        help=(
            "Path to query sequences to annotate (FASTA file or directory). "
            "If omitted, the wizard will prompt for it."
        ),
    )
    io.add_argument(
        "--output",
        dest="output",
        metavar="DIR",
        default=None,
        help="Output directory for the workspace. If omitted, the wizard will prompt.",
    )

    # -- Annotation ----------------------------------------------------------
    ann = p.add_argument_group("Annotation (choose one)")
    ann_mx = ann.add_mutually_exclusive_group()
    ann_mx.add_argument(
        "--ann",
        dest="ann",
        metavar="PATH",
        default=None,
        help=(
            "Path to an existing annotation file (.ann). "
            "Format: tab-separated with columns 'id' and 'family'."
        ),
    )
    ann_mx.add_argument(
        "--create-ann",
        dest="create_ann",
        action="store_true",
        default=False,
        help=(
            "Generate annotations from training FASTA headers. "
            "Requires headers in the format: >db|FAMILY_LABEL|seqid description "
            "(the field between the first pair of | | becomes the family label)."
        ),
    )

    # -- K-mer parameters ----------------------------------------------------
    kmer = p.add_argument_group("K-mer parameters")
    kmer.add_argument(
        "--k",
        dest="k",
        type=int,
        default=8,
        metavar="N",
        help="K-mer length.",
    )
    kmer.add_argument(
        "--alphabet",
        dest="alphabet",
        type=str,
        default="2",
        metavar="",
        help=(
            "Reduced alphabet encoding (0–5, alphabet name, or 'None'). "
            "2 = solvacc (3-letter). See alphabets list below."
        ),
    )

    # -- Learn/apply options -------------------------------------------------
    la = p.add_argument_group("Learn / apply options")
    la.add_argument(
        "--selection",
        dest="selection",
        choices=["top_hit", "greatest_distance", "combined_distance"],
        default="top_hit",
        metavar="",
        help="Annotation selection method {top_hit, greatest_distance, combined_distance}.",
    )
    la.add_argument(
        "--threshold",
        dest="threshold",
        type=str,
        default="Median",
        metavar="",
        help=(
            "Family-specific score threshold for prediction filtering. "
            "Options: 'Median', 'Mean', '90th Percentile', 'None'."
        ),
    )
    la.add_argument(
        "--apply-output",
        dest="apply_output",
        type=str,
        default="snekmer_results.csv",
        metavar="FILENAME",
        help="Output filename for apply results.",
    )

    # -- Fragmentation (advanced) --------------------------------------------
    frag = p.add_argument_group("Fragmentation (advanced)")
    frag.add_argument(
        "--fragmentation",
        dest="fragmentation",
        action="store_true",
        default=False,
        help="Split sequences into fragments before kmerization.",
    )
    frag.add_argument(
        "--frag-length",
        dest="frag_length",
        type=int,
        default=None,
        metavar="N",
        help="Fragment length in residues (default: 50).",
    )
    frag.add_argument(
        "--min-length",
        dest="min_length",
        type=int,
        default=None,
        metavar="N",
        help="Minimum sequence length to fragment (default: 50).",
    )
    frag.add_argument(
        "--fragment-version",
        dest="fragment_version",
        type=str,
        default=None,
        metavar="",
        help="Fragmentation version (default: absolute).",
    )
    frag.add_argument(
        "--fragment-location",
        dest="fragment_location",
        type=str,
        default=None,
        metavar="",
        help="Fragment location method (default: random).",
    )
    frag.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=None,
        metavar="N",
        help="Random seed for fragmentation (default: 999).",
    )

    # -- Snakemake pass-through ----------------------------------------------
    smk = p.add_argument_group("Snakemake options")
    smk.add_argument(
        "--cores",
        "-c",
        dest="cores",
        type=int,
        default=cpu_count(),
        metavar="N",
        help="CPU cores to use.",
    )
    smk.add_argument(
        "--dry-run",
        "-n",
        dest="dryrun",
        action="store_true",
        help="Show what would be done without executing.",
    )
    smk.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Show additional Snakemake debug output.",
    )
    smk.add_argument(
        "--quiet",
        "-q",
        nargs="*",
        choices=["progress", "rules", "all"],
        default=None,
        help="Reduce Snakemake output.",
    )

    # -- File handling -------------------------------------------------------
    misc = p.add_argument_group("Miscellaneous")
    misc.add_argument(
        "--copy-files",
        dest="copy_files",
        action="store_true",
        help=(
            "Copy input files into the workspace instead of symlinking them "
            "(useful when the workspace will be moved or shared)."
        ),
    )


# ------------------------- parser construction -------------------------
def get_argument_parser():
    parser = {}

    # -------------------------
    # Snakemake arguments (passed through)
    # -------------------------
    parser["smk"] = argparse.ArgumentParser(add_help=False)

    g_smk = parser["smk"].add_argument_group(
        "Snakemake arguments (passed through to Snakemake)"
    )
    g_smk.add_argument(
        "--dry-run",
        "--dryrun",
        "-n",
        dest="dryrun",
        action="store_true",
        help="Do not execute anything, and display what would be done. "
        "If you have a very large workflow, use --dry-run --quiet to just "
        "print a summary of the DAG of jobs.",
    )
    g_smk.add_argument(
        "--configfile",
        nargs="+",
        metavar="PATH",
        help=(
            "Specify or overwrite workflow config file(s). Multiple files overwrite each other "
            "in the given order. Values are available via Snakemake's global config dictionary."
        ),
    )
    g_smk.add_argument(
        "-C",
        "--config",
        nargs="*",
        metavar="KEY=VALUE",
        help="Set or overwrite values in the workflow config object (Snakemake --config KEY=VALUE).",
    )
    g_smk.add_argument("--unlock", action="store_true", help="Unlock the working directory.")
    g_smk.add_argument(
        "--until",
        "-U",
        nargs="+",
        metavar="TARGET",
        help="Run the workflow until it reaches the specified rules or files.",
    )
    g_smk.add_argument(
        "--keepgoing",
        "--keep-going",
        "-k",
        action="store_true",
        help="Continue with independent jobs if a job fails.",
    )
    g_smk.add_argument(
        "--latency",
        "-w",
        "--output-wait",
        "--latency-wait",
        type=int,
        default=30,
        metavar="SECONDS",
        help="Wait given seconds for output files to appear after job completion (filesystem latency).",
    )
    g_smk.add_argument(
        "--touch",
        "-t",
        action="store_true",
        help="Touch output files instead of running commands (mark as up-to-date).",
    )
    g_smk.add_argument(
        "--cores",
        "-c",
        default=cpu_count(),
        type=int,
        metavar="N",
        help="Use at most N CPU cores/jobs in parallel.",
    )
    g_smk.add_argument(
        "--count",
        metavar="N",
        type=int,
        help="Number of files to process (limits DAG size).",
    )
    g_smk.add_argument(
        "--countstart",
        metavar="IDX",
        type=int,
        default=0,
        help="Starting file index (for use with --count).",
    )
    g_smk.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Show additional debug output.",
    )
    g_smk.add_argument(
        "--quiet",
        "-q",
        nargs="*",
        choices=["progress", "rules", "all"],
        default=None,
        help="Reduce Snakemake output (progress/rules/all). If used without args, quiets progress and rules.",
    )
    g_smk.add_argument(
        "--directory",
        "-d",
        metavar="DIR",
        action="store",
        help="Specify working directory (relative paths in the Snakefile use this origin).",
    )
    g_smk.add_argument(
        "--forcerun",
        "-R",
        nargs="*",
        metavar="TARGET",
        help="Force re-execution/creation of the given rules or files.",
    )
    g_smk.add_argument(
        "--list-code-changes",
        "--lc",
        action="store_true",
        help="List output files for which the rule body changed in the Snakefile.",
    )
    g_smk.add_argument(
        "--list-params-changes",
        "--lp",
        action="store_true",
        help="List output files for which defined params changed in the Snakefile.",
    )
    g_smk.add_argument(
    "--scheduler",
    choices=["greedy", "ilp"],
    default=None,
    help="Snakemake scheduler plugin to use.",
    )

    g_smk.add_argument(
        "--scheduler-ilp-solver",
        default=None,
        help="MILP solver to use with the ILP scheduler.",
    )

    g_smk.add_argument(
        "--scheduler-ilp-solver-path",
        metavar="PATH",
        default=None,
        help="PATH to search for ILP scheduler solver binaries.",
    )
    g_smk_cfg = parser["smk"].add_argument_group("Snekmer configfile behavior")
    g_smk_cfg.add_argument(
        "--no-default-configfile",
        action="store_true",
        help="Do not auto-load ./config.yaml (or <DIR>/config.yaml when using -d).",
    )

    g_clust = parser["smk"].add_argument_group(
        "Snakemake cluster execution (passed through)"
    )
    g_clust.add_argument(
        "--clust",
        nargs="+",
        metavar="PATH",
        help="Path to cluster execution YAML configuration file (e.g., for SLURM).",
    )
    g_clust.add_argument(
        "-j",
        "--jobs",
        metavar="N",
        type=int,
        default=1000,
        help="Number of simultaneous jobs to submit to the scheduler.",
    )

    # -------------------------
    # Snekmer parameter parsers (mode-scoped)
    # -------------------------
    parser["cfg_required"] = argparse.ArgumentParser(add_help=False)
    parser["cfg_score"] = argparse.ArgumentParser(add_help=False)
    parser["cfg_cluster"] = argparse.ArgumentParser(add_help=False)
    parser["cfg_model"] = argparse.ArgumentParser(add_help=False)
    parser["cfg_search"] = argparse.ArgumentParser(add_help=False)
    parser["cfg_learn_apply"] = argparse.ArgumentParser(add_help=False)

    _add_required_args(parser["cfg_required"])
    _add_score_args(parser["cfg_score"])
    _add_cluster_args(parser["cfg_cluster"])
    _add_model_args(parser["cfg_model"])
    _add_search_args(parser["cfg_search"])
    _add_learn_apply_args(parser["cfg_learn_apply"])

    # -------------------------
    # main parser
    # -------------------------
    parser["main"] = _SnekmerMainParser(
        description=_main_description(),
        epilog=_main_epilog(),
        formatter_class=SnekmerHelpFormatter,
        parents=[
            parser["smk"],
            parser["cfg_required"],
            parser["cfg_score"],
            parser["cfg_cluster"],
            parser["cfg_model"],
            parser["cfg_search"],
            parser["cfg_learn_apply"],
        ],
    )
    parser["main"].add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="Print version and exit.",
    )

    parser["subparsers"] = parser["main"].add_subparsers(
        title="mode",
        description="Snekmer mode (cluster, model, search, learn, apply, easy-learn-apply).",
        dest="mode",
        parser_class=_SnekmerSubParser,
    )

    # subparsers (each includes only the relevant Snekmer parameter sections)
    parser["cluster"] = parser["subparsers"].add_parser(
        "cluster",
        description="Unsupervised clustering workflow.",
        parents=[parser["smk"], parser["cfg_required"], parser["cfg_cluster"]],
    )
    parser["model"] = parser["subparsers"].add_parser(
        "model",
        description="Train supervised models + cross-validation reports.",
        parents=[
            parser["smk"],
            parser["cfg_required"],
            parser["cfg_score"],
            parser["cfg_model"],
        ],
    )
    parser["search"] = parser["subparsers"].add_parser(
        "search",
        description="Score sequences against trained models.",
        parents=[
            parser["smk"],
            parser["cfg_required"],
            parser["cfg_score"],
            parser["cfg_search"],
        ],
    )
    parser["learn"] = parser["subparsers"].add_parser(
        "learn",
        description="Build annotation-associated k-mer distributions + confidence evaluation.",
        parents=[parser["smk"], parser["cfg_required"], parser["cfg_learn_apply"]],
    )
    parser["apply"] = parser["subparsers"].add_parser(
        "apply",
        description="Predict annotations using outputs from learn.",
        parents=[parser["smk"], parser["cfg_required"], parser["cfg_learn_apply"]],
    )

    parser["easy-learn-apply"] = parser["subparsers"].add_parser(
        "easy-learn-apply",
        description=(
            "Guided front-end that runs learn then apply end-to-end.\n\n"
            "Prompts for training sequences, query sequences, and annotation style,\n"
            "then builds a self-contained workspace and runs both pipeline steps.\n"
            "All prompts can be skipped by supplying the corresponding flags."
        ),
        formatter_class=SnekmerHelpFormatter,
        epilog=_main_epilog(),
    )
    _add_easy_la_args(parser["easy-learn-apply"])

    return parser


def get_main_args():
    parser = get_argument_parser()
    return parser["main"]


# ------------------------- Snakemake v9 bridge helpers -------------------------
def _coerce_scalar(val: str):
    """Best-effort cast of KEY=VALUE overrides (bool/int/float/json/str/None)."""
    low = val.lower()
    if low in {"true", "false"}:
        return low == "true"
    if low in {"none", "null"}:
        return None
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    try:
        return json.loads(val)
    except Exception:
        return val


def _deep_update(dst: dict, src: dict) -> dict:
    """Recursively update dst with src (in-place), return dst."""
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def _cli_config_overrides(args, force_defaults: bool) -> dict:
    """
    Build nested config dict from Snekmer parameter flags.

    - If force_defaults is False, only include values the user explicitly provided.
    - If force_defaults is True (no config files are in use), include defaults too.
    """
    cfg = {}

    def want(dest: str) -> bool:
        return force_defaults or bool(getattr(args, f"{dest}_seen", False))

    def put(path, value):
        d = value
        for key in reversed(path):
            d = {key: d}
        _deep_update(cfg, d)

    # required (all modes)
    if hasattr(args, "cfg_k") and want("cfg_k"):
        put(["k"], args.cfg_k)
    if hasattr(args, "cfg_alphabet") and want("cfg_alphabet"):
        put(["alphabet"], _coerce_scalar(str(args.cfg_alphabet)))
    if hasattr(args, "cfg_input_dir") and want("cfg_input_dir"):
        put(["input_dir"], args.cfg_input_dir)
    if hasattr(args, "cfg_input_file_exts") and want("cfg_input_file_exts"):
        put(["input_file_exts"], list(args.cfg_input_file_exts))
    if hasattr(args, "cfg_input_file_regex") and want("cfg_input_file_regex"):
        put(["input_file_regex"], args.cfg_input_file_regex)
    if hasattr(args, "cfg_nested_output") and want("cfg_nested_output"):
        put(["nested_output"], args.cfg_nested_output)

    # model + search only: score
    if hasattr(args, "cfg_score_scaler") and want("cfg_score_scaler"):
        put(["score", "scaler"], args.cfg_score_scaler)
    if hasattr(args, "cfg_score_scaler_n") and want("cfg_score_scaler_n"):
        put(["score", "scaler_kwargs", "n"], args.cfg_score_scaler_n)
    if hasattr(args, "cfg_score_labels") and want("cfg_score_labels"):
        put(["score", "labels"], args.cfg_score_labels)
    if hasattr(args, "cfg_score_lname") and want("cfg_score_lname"):
        put(["score", "lname"], args.cfg_score_lname)

    # cluster only
    if hasattr(args, "cfg_cluster_method") and want("cfg_cluster_method"):
        put(["cluster", "method"], args.cfg_cluster_method)
    if hasattr(args, "cfg_cluster_n_clusters") and want("cfg_cluster_n_clusters"):
        put(
            ["cluster", "params", "n_clusters"],
            _coerce_scalar(str(args.cfg_cluster_n_clusters)),
        )
    if hasattr(args, "cfg_cluster_linkage") and want("cfg_cluster_linkage"):
        put(["cluster", "params", "linkage"], args.cfg_cluster_linkage)
    if hasattr(args, "cfg_cluster_distance_threshold") and want(
        "cfg_cluster_distance_threshold"
    ):
        put(
            ["cluster", "params", "distance_threshold"],
            args.cfg_cluster_distance_threshold,
        )
    if hasattr(args, "cfg_cluster_compute_full_tree") and want(
        "cfg_cluster_compute_full_tree"
    ):
        put(
            ["cluster", "params", "compute_full_tree"],
            args.cfg_cluster_compute_full_tree,
        )
    if hasattr(args, "cfg_cluster_cluster_plots") and want("cfg_cluster_cluster_plots"):
        put(["cluster", "cluster_plots"], args.cfg_cluster_cluster_plots)
    if hasattr(args, "cfg_cluster_min_rep") and want("cfg_cluster_min_rep"):
        put(["cluster", "min_rep"], _coerce_scalar(str(args.cfg_cluster_min_rep)))
    if hasattr(args, "cfg_cluster_max_rep") and want("cfg_cluster_max_rep"):
        put(["cluster", "max_rep"], _coerce_scalar(str(args.cfg_cluster_max_rep)))
    if hasattr(args, "cfg_cluster_save_matrix") and want("cfg_cluster_save_matrix"):
        put(["cluster", "save_matrix"], args.cfg_cluster_save_matrix)
    if hasattr(args, "cfg_cluster_dist_thresh") and want("cfg_cluster_dist_thresh"):
        put(["cluster", "dist_thresh"], args.cfg_cluster_dist_thresh)

    # model only
    if hasattr(args, "cfg_model_cv") and want("cfg_model_cv"):
        put(["model", "cv"], args.cfg_model_cv)
    if hasattr(args, "cfg_model_random_state") and want("cfg_model_random_state"):
        put(
            ["model", "random_state"],
            _coerce_scalar(str(args.cfg_model_random_state)),
        )

    # search only
    if hasattr(args, "cfg_model_dir") and want("cfg_model_dir"):
        put(["model_dir"], args.cfg_model_dir)
    if hasattr(args, "cfg_basis_dir") and want("cfg_basis_dir"):
        put(["basis_dir"], args.cfg_basis_dir)
    if hasattr(args, "cfg_score_dir") and want("cfg_score_dir"):
        put(["score_dir"], args.cfg_score_dir)

    # learn + apply only
    if hasattr(args, "cfg_la_save_apply_associations") and want(
        "cfg_la_save_apply_associations"
    ):
        put(
            ["learn_apply", "save_apply_associations"],
            args.cfg_la_save_apply_associations,
        )
    if hasattr(args, "cfg_la_conf_weight_modifier") and want(
        "cfg_la_conf_weight_modifier"
    ):
        put(["learn_apply", "conf_weight_modifier"], args.cfg_la_conf_weight_modifier)
    if hasattr(args, "cfg_la_fragmentation") and want("cfg_la_fragmentation"):
        put(["learn_apply", "fragmentation"], args.cfg_la_fragmentation)
    if hasattr(args, "cfg_la_version") and want("cfg_la_version"):
        put(["learn_apply", "version"], args.cfg_la_version)
    if hasattr(args, "cfg_la_frag_length") and want("cfg_la_frag_length"):
        put(["learn_apply", "frag_length"], args.cfg_la_frag_length)
    if hasattr(args, "cfg_la_min_length") and want("cfg_la_min_length"):
        put(["learn_apply", "min_length"], args.cfg_la_min_length)
    if hasattr(args, "cfg_la_location") and want("cfg_la_location"):
        put(["learn_apply", "location"], args.cfg_la_location)
    if hasattr(args, "cfg_la_seed") and want("cfg_la_seed"):
        put(["learn_apply", "seed"], args.cfg_la_seed)
    if hasattr(args, "cfg_la_selection") and want("cfg_la_selection"):
        put(["learn_apply", "selection"], args.cfg_la_selection)
    if hasattr(args, "cfg_la_threshold") and want("cfg_la_threshold"):
        put(["learn_apply", "threshold"], _coerce_scalar(str(args.cfg_la_threshold)))
    if hasattr(args, "cfg_la_weight_top") and want("cfg_la_weight_top"):
        put(["learn_apply", "weight_top"], args.cfg_la_weight_top)
    if hasattr(args, "cfg_la_weight_distance") and want("cfg_la_weight_distance"):
        put(["learn_apply", "weight_distance"], args.cfg_la_weight_distance)
    if hasattr(args, "cfg_la_apply_output") and want("cfg_la_apply_output"):
        put(["learn_apply", "apply_output"], args.cfg_la_apply_output)

    return cfg


def _parse_config_overrides(args, force_defaults: bool):
    """Re-implement minimal parse_config behavior (Snakemake ≥8 removed API)."""
    cfg = _cli_config_overrides(args, force_defaults=force_defaults)

    # -C overrides are always applied (highest precedence)
    if args.config:
        for item in args.config:
            if "=" not in item:
                raise ValueError(f"--config expects KEY=VALUE, got: {item!r}")
            k, v = item.split("=", 1)
            cfg[k] = _coerce_scalar(v)

    # preserve original start/stop semantics
    if args.count is not None:
        cfg = {"start": args.countstart, "stop": args.countstart + args.count, **cfg}

    return cfg


def _quiet_cli(args):
    """Map original quiet behavior to Snakemake CLI flags."""
    if args.quiet is None:
        return []
    if len(args.quiet) == 0:
        return ["--quiet", "progress", "--quiet", "rules"]
    out = []
    for q in args.quiet:
        out.extend(["--quiet", q])
    return out


def _resolve_rulefile(relpath):
    """Resolve packaged rules/<name>.smk using importlib.resources (replaces resource_filename)."""
    res = files("snekmer") / relpath
    return as_file(res)  # context manager yielding filesystem path


def _build_snakemake_cmd(
    rulefile_path, args, configfiles, config=None, keepgoing_override=None
):
    """Translate previous snakemake(...) call into a Snakemake v9 CLI invocation."""
    cmd = ["snakemake", "-s", rulefile_path, "--cores", str(args.cores)]

    if args.dryrun:
        cmd.append("--dry-run")
    if args.touch:
        cmd.append("--touch")
    if args.unlock:
        cmd.append("--unlock")
    if args.list_code_changes:
        cmd.append("--list-code-changes")
    if args.list_params_changes:
        cmd.append("--list-params-changes")
    if args.forcerun:
        cmd.extend(["--forcerun", *args.forcerun])
    if args.until:
        cmd.extend(["--until", *args.until])
    if args.latency is not None:
        cmd.extend(["--latency-wait", str(args.latency)])
    if args.verbose:
        cmd.append("--verbose")
    
    if getattr(args, "scheduler", None):
        cmd.extend(["--scheduler", args.scheduler])

    if getattr(args, "scheduler_ilp_solver", None):
        cmd.extend(["--scheduler-ilp-solver", args.scheduler_ilp_solver])

    if getattr(args, "scheduler_ilp_solver_path", None):
        cmd.extend(
            ["--scheduler-ilp-solver-path", args.scheduler_ilp_solver_path])
        

    cmd.append("--rerun-incomplete")

    if keepgoing_override is True or (keepgoing_override is None and args.keepgoing):
        cmd.append("--keep-going")

    # config files (JSON/YAML)
    if configfiles:
        cmd.extend(["--configfile", *configfiles])

    if args.jobs is not None:
        cmd.extend(["--jobs", str(args.jobs)])

    if args.directory:
        cmd.extend(["--directory", os.path.abspath(args.directory)])

    cmd.extend(_quiet_cli(args))

    if args.clust is not None:
        cluster = (
            "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} "
            "-J {cluster.name} --ntasks-per-node {cluster.ntasks} -p {cluster.partition}"
        )
        cmd.extend(
            [
                "--executor",
                "cluster-generic",
                "--cluster-generic-submit-cmd",
                cluster,
            ]
        )
        for cc in args.clust:
            cmd.extend(["--cluster-config", cc])

    return cmd


def _run(rule_basename, args, keepgoing_override=None):
    """Helper to mirror per-mode Snakemake CLI calls with minimal changes."""
    # resolve configfiles
    if args.configfile is not None:
        configfiles = list(map(os.path.abspath, args.configfile))
    elif getattr(args, "no_default_configfile", False):
        configfiles = []
    else:
        candidate = (
            os.path.join(args.directory, "config.yaml")
            if args.directory is not None
            else os.path.abspath("config.yaml")
        )
        configfiles = [os.path.abspath(candidate)] if os.path.exists(candidate) else []

    force_defaults = len(configfiles) == 0
    config = _parse_config_overrides(args, force_defaults=force_defaults)

    # Write Snekmer-derived config (defaults and/or overrides) to a temporary JSON config file.
    # This avoids Snakemake's CLI --config parsing edge cases.
    tmp_cfg_path = None
    if config:
        tmp = tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".json",
            delete=False,
            prefix="snekmer_cli_config_",
        )
        try:
            json.dump(config, tmp)
            tmp.flush()
            tmp_cfg_path = os.path.abspath(tmp.name)
        finally:
            tmp.close()

        configfiles = [*configfiles, tmp_cfg_path]

    relpath = os.path.join("rules", f"{rule_basename}.smk")
    try:
        with _resolve_rulefile(relpath) as rulefile_path:
            cmd = _build_snakemake_cmd(
                str(rulefile_path),
                args,
                configfiles,
                config=None,
                keepgoing_override=keepgoing_override,
            )
            if args.verbose:
                print(">> Running:", shlex.join(cmd), file=sys.stderr)
            proc = subprocess.run(cmd)
            return proc.returncode
    finally:
        if tmp_cfg_path is not None:
            try:
                os.unlink(tmp_cfg_path)
            except OSError:
                pass


def main():
    parser = get_argument_parser()

    # Route each subcommand through its own subparser so that unrecognized-
    # argument errors show only that subcommand's concise usage, not the full
    # snekmer wall of options.  Unknown or missing modes fall through to the
    # main parser, which now has a clean, organized error handler.
    mode = sys.argv[1] if len(sys.argv) > 1 else None
    if mode in _MODES:
        args = parser[mode].parse_args(sys.argv[2:])
        args.mode = mode
    else:
        args = parser["main"].parse_args()

    if args.quiet is not None and len(args.quiet) == 0:
        args.quiet = ["progress", "rules"]

    if args.mode == "cluster":
        rc = _run("cluster", args, keepgoing_override=args.keepgoing)
    elif args.mode == "model":
        rc = _run("model", args, keepgoing_override=args.keepgoing)
    elif args.mode == "search":
        rc = _run("search", args, keepgoing_override=True)
    elif args.mode == "learn":
        rc = _run("learn", args, keepgoing_override=True)
    elif args.mode == "apply":
        rc = _run("apply", args, keepgoing_override=True)
    elif args.mode == "easy-learn-apply":
        from snekmer.easy import run_easy_learn_apply
        rc = run_easy_learn_apply(args)
    else:
        parser["main"].print_help()
        rc = 2

    sys.exit(rc)


if __name__ == "__main__":
    main()