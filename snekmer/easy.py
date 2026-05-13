"""easy: Guided front-end for the Snekmer learn/apply pipeline.

Provides 'snekmer easy-learn-apply', a wizard that:
  1. Collects training sequences, query sequences, and annotation info from
     the user (interactively or via flags).
  2. Optionally extracts family labels directly from FASTA headers.
  3. Scaffolds a self-contained output directory with learn/ and apply/
     subdirectories, each with the config and inputs the pipeline expects.
  4. Runs snekmer learn, copies the three handoff files, then runs snekmer apply.
"""

import argparse
import gzip
import json
import os
import re
import shutil
import sys
import tempfile
from multiprocessing import cpu_count
from pathlib import Path
from typing import Optional

from Bio import SeqIO


FASTA_EXTS = ("fasta", "fna", "faa", "fa")

# Same regex used internally by learn_learn.py to extract the pipe-enclosed ID
_SEQID_PATTERN = re.compile(r"\|(.*?)\|")


# ---------------------------------------------------------------------------
# Interactive prompt helper
# ---------------------------------------------------------------------------

def _prompt(message: str, default: Optional[str] = None) -> str:
    """Print a prompt and return stripped user input, using default on empty."""
    suffix = f" [{default}]" if default is not None else ""
    try:
        val = input(f"{message}{suffix}: ").strip()
    except (KeyboardInterrupt, EOFError):
        print("\nAborted.", file=sys.stderr)
        sys.exit(1)
    return val if val else (default or "")


# ---------------------------------------------------------------------------
# File collection
# ---------------------------------------------------------------------------

def _collect_fasta_files(path: str) -> list:
    """Return a sorted list of Path objects for FASTA files at path (file or dir)."""
    p = Path(path).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")
    if p.is_file():
        return [p]
    files = []
    for ext in FASTA_EXTS:
        files.extend(p.glob(f"*.{ext}"))
        files.extend(p.glob(f"*.{ext}.gz"))
    if not files:
        raise ValueError(
            f"No FASTA files found in: {path}\n"
            f"  Looked for extensions: {', '.join(f'.{e}' for e in FASTA_EXTS)}"
        )
    return sorted(files)


def _open_fasta(path: Path):
    """Return an open file handle for a FASTA file (plain or gzip-compressed)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


# ---------------------------------------------------------------------------
# Annotation generation from FASTA headers
# ---------------------------------------------------------------------------

def generate_ann_from_fasta(fasta_files: list, ann_path: Path) -> int:
    """Parse FASTA headers and write a tab-separated annotation file (id TAB family).

    Extracts the first pipe-enclosed field from each sequence header as the
    family label:

        >db|FAMILY_LABEL|seqid Description  ->  id=full_seqid, family=FAMILY_LABEL

    Parameters
    ----------
    fasta_files : list of Path
        Training FASTA files to parse.
    ann_path : Path
        Destination path for the generated .ann file.

    Returns
    -------
    int
        Number of successfully annotated sequences.
    """
    seen_families: set = set()
    unparseable = []

    for fasta in fasta_files:
        with _open_fasta(fasta) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                m = _SEQID_PATTERN.search(record.id)
                if m:
                    seen_families.add(m.group(1))
                else:
                    unparseable.append(record.id)

    if unparseable:
        print(
            f"\n  WARNING: {len(unparseable)} sequence(s) have no pipe-enclosed field "
            "in their header and will be excluded from training."
        )
        print("  Examples of unparseable headers:")
        for sid in unparseable[:3]:
            print(f"    {sid}")
        if len(unparseable) > 3:
            print(f"    ... and {len(unparseable) - 3} more.")

    # filter_and_construct extracts the between-pipes field from each sequence ID
    # and looks it up in the annotation "id" column, so id must equal the family label.
    ann_path.parent.mkdir(parents=True, exist_ok=True)
    with open(ann_path, "w", encoding="utf-8") as fh:
        fh.write("id\tfamily\n")
        for family in sorted(seen_families):
            fh.write(f"{family}\t{family}\n")

    return len(seen_families)


# ---------------------------------------------------------------------------
# Workspace scaffolding
# ---------------------------------------------------------------------------

def _coerce_alphabet(alphabet):
    """Return the alphabet value in the form Snakemake's check_valid() expects.

    check_valid() accepts an integer (0-5) or a named string (e.g. 'solvacc').
    A bare numeric string like '2' matches neither, so convert it to int.
    'None' is kept as the string 'None'.
    """
    if isinstance(alphabet, str) and alphabet.lower() != "none":
        try:
            return int(alphabet)
        except ValueError:
            pass
    return alphabet


def _write_config(
    config_path: Path,
    k: int,
    alphabet,
    la_overrides: Optional[dict] = None,
) -> None:
    """Write a config.yaml for the learn or apply working directory.

    JSON is a strict subset of YAML, so writing JSON to config.yaml is
    valid and avoids a PyYAML dependency.
    """
    cfg = {
        "k": k,
        "alphabet": _coerce_alphabet(alphabet),
        "input_dir": "input",
        "input_file_exts": list(FASTA_EXTS),
        "input_file_regex": ".*",
        "nested_output": False,
        "learn_apply": {
            "save_apply_associations": False,
            "fragmentation": False,
            "version": "absolute",
            "frag_length": 50,
            "min_length": 50,
            "location": "random",
            "seed": 999,
            "conf_weight_modifier": 20,
            "selection": "top_hit",
            "threshold": "Median",
            "weight_top": 0.7,
            "weight_distance": 0.3,
            "apply_output": "snekmer_results.csv",
        },
    }
    if la_overrides:
        cfg["learn_apply"].update(la_overrides)
    with open(config_path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2)


def _link_or_copy(src: Path, dst: Path, symlink: bool) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if symlink:
        dst.symlink_to(src)
    else:
        shutil.copy2(src, dst)


def setup_workspace(
    output_dir: Path,
    train_files: list,
    query_files: list,
    ann_file: Path,
    k: int,
    alphabet,
    la_overrides: Optional[dict] = None,
    symlink: bool = True,
):
    """Create learn/ and apply/ subdirectories with inputs and configs.

    Returns
    -------
    (learn_dir, apply_dir) : tuple of Path
    """
    learn_dir = output_dir / "learn"
    apply_dir = output_dir / "apply"

    for d in [
        learn_dir / "input",
        learn_dir / "annotations",
        apply_dir / "input",
        apply_dir / "counts",
        apply_dir / "confidence",
        apply_dir / "stats",
    ]:
        d.mkdir(parents=True, exist_ok=True)

    for src in train_files:
        _link_or_copy(src, learn_dir / "input" / src.name, symlink)

    for src in query_files:
        _link_or_copy(src, apply_dir / "input" / src.name, symlink)

    shutil.copy2(ann_file, learn_dir / "annotations" / "annotations.ann")

    _write_config(learn_dir / "config.yaml", k, alphabet, la_overrides)
    _write_config(apply_dir / "config.yaml", k, alphabet, la_overrides)

    return learn_dir, apply_dir


def handoff_learn_to_apply(learn_dir: Path, apply_dir: Path) -> None:
    """Copy the three required learn outputs into the apply input directories.

    learn produces these via the copy_results_for_apply rule:
      learn/apply_inputs/counts/kmer_counts_total.csv
      learn/apply_inputs/stats/family_summary_stats.csv
      learn/apply_inputs/confidence/global_confidence_scores.csv

    apply expects them at:
      apply/counts/*.csv
      apply/stats/*.csv
      apply/confidence/*.csv
    """
    copies = [
        (
            learn_dir / "apply_inputs" / "counts" / "kmer_counts_total.csv",
            apply_dir / "counts" / "kmer_counts_total.csv",
        ),
        (
            learn_dir / "apply_inputs" / "stats" / "family_summary_stats.csv",
            apply_dir / "stats" / "family_summary_stats.csv",
        ),
        (
            learn_dir / "apply_inputs" / "confidence" / "global_confidence_scores.csv",
            apply_dir / "confidence" / "global_confidence_scores.csv",
        ),
    ]
    for src, dst in copies:
        if not src.exists():
            raise FileNotFoundError(
                f"Expected learn output not found: {src}\n"
                "  The learn step may have failed — check the output above."
            )
        shutil.copy2(src, dst)


# ---------------------------------------------------------------------------
# Build a run-args namespace for cli._run()
# ---------------------------------------------------------------------------

def _make_run_args(base_args, directory: str) -> argparse.Namespace:
    """Construct an argparse.Namespace that cli._run() accepts.

    Pulls cores/verbose/quiet/dryrun from base_args where available and
    sets safe defaults for everything else. The working directory is set
    to `directory` so that Snakemake uses the pre-built workspace.
    """
    ns = argparse.Namespace()

    # Snakemake pass-through
    ns.cores = getattr(base_args, "cores", None) or cpu_count()
    ns.dryrun = getattr(base_args, "dryrun", False)
    ns.touch = False
    ns.unlock = False
    ns.list_code_changes = False
    ns.list_params_changes = False
    ns.forcerun = None
    ns.until = None
    ns.latency = 30
    ns.verbose = getattr(base_args, "verbose", False)
    ns.quiet = getattr(base_args, "quiet", None)
    ns.keepgoing = True
    ns.clust = None
    ns.jobs = 1000
    ns.scheduler = None
    ns.scheduler_ilp_solver = None
    ns.scheduler_ilp_solver_path = None
    ns.count = None
    ns.countstart = 0

    # Config: auto-discover config.yaml written to the working directory.
    # Setting configfile=None and no_default_configfile=False makes _run()
    # look for <directory>/config.yaml automatically.
    ns.configfile = None
    ns.no_default_configfile = False
    ns.config = None

    # Working directory managed by easy-learn-apply
    ns.directory = directory

    return ns


# ---------------------------------------------------------------------------
# Wizard sub-steps (keep run_easy_learn_apply readable)
# ---------------------------------------------------------------------------

def _wizard_inputs(args, interactive: bool):
    """Prompt for training and query paths; return (train_files, query_files)."""
    train_path = getattr(args, "train", None)
    if not train_path:
        train_path = _prompt("Step 1  Training sequences (file or directory path)")

    query_path = getattr(args, "query", None)
    if not query_path:
        query_path = _prompt("Step 2  Query sequences (file or directory path)")

    train_files = _collect_fasta_files(train_path)
    query_files = _collect_fasta_files(query_path)

    if interactive:
        print(
            f"  Found {len(train_files)} training file(s), "
            f"{len(query_files)} query file(s)."
        )
    return train_files, query_files


def _wizard_annotation(args, interactive: bool, train_files: list) -> Path:
    """Resolve annotation source; return path to a ready-to-use .ann temp file.

    Checks --ann (existing file) and --create-ann (generate from headers) flags.
    When neither is supplied, interactively asks the user.
    The caller is responsible for unlinking the returned path when done.
    """
    ann_file_given = getattr(args, "ann", None)
    create_ann = getattr(args, "create_ann", False)

    if not ann_file_given and not create_ann:
        print(
            "\nStep 3  How are your training sequences annotated?\n"
            "\n"
            "  [1] Family labels are embedded in FASTA headers (between | | characters)\n"
            "      Example:  >db|TIGR04183|seqid Description text\n"
            "                       ^^^^^^^^\n"
            "                  this field becomes the family label\n"
            "      (equivalent to passing --create-ann)\n"
            "\n"
            "  [2] I have a separate annotation file (.ann)\n"
            "      Format: tab-separated with columns: id  family\n"
            "      (equivalent to passing --ann <path>)\n"
        )
        choice = _prompt("  Choice", default="1")
        if choice.strip() == "1":
            create_ann = True
        else:
            ann_file_given = _prompt("\n  Path to annotation file (.ann)")

    # Use a system temp file so output_dir is not needed at this stage
    fd, tmp_path = tempfile.mkstemp(suffix=".ann", prefix="snekmer_ann_")
    os.close(fd)
    ann_tmp = Path(tmp_path)

    if create_ann:
        if interactive:
            print("  Extracting family labels from training FASTA headers...")
        n_annotated = generate_ann_from_fasta(train_files, ann_tmp)
        if n_annotated == 0:
            ann_tmp.unlink(missing_ok=True)
            raise ValueError(
                "No sequences with pipe-enclosed family labels were found.\n"
                "  Ensure your training FASTA headers use the format:\n"
                "    >db|FAMILY_LABEL|seqid Description text\n"
                "  where FAMILY_LABEL is the annotation for each sequence."
            )
        if interactive:
            print(f"  Annotated {n_annotated} sequences.")
    else:
        ann_path = Path(ann_file_given).resolve()
        if not ann_path.exists():
            ann_tmp.unlink(missing_ok=True)
            raise FileNotFoundError(f"Annotation file not found: {ann_path}")
        shutil.copy2(ann_path, ann_tmp)

    return ann_tmp


def _wizard_output_dir(args) -> Path:
    """Prompt for output directory if not supplied; return resolved Path."""
    output_str = getattr(args, "output", None)
    if not output_str:
        output_str = _prompt("\nStep 4  Output directory", default="snekmer_easy_output")
    return Path(output_str).resolve()


def _collect_la_overrides(args) -> dict:
    """Extract any learn/apply config overrides from args."""
    overrides = {}
    for attr, key in [
        ("selection", "selection"),
        ("threshold", "threshold"),
        ("save_apply_associations", "save_apply_associations"),
        ("fragmentation", "fragmentation"),
        ("fragment_version", "version"),
        ("frag_length", "frag_length"),
        ("min_length", "min_length"),
        ("fragment_location", "location"),
        ("seed", "seed"),
        ("conf_weight_modifier", "conf_weight_modifier"),
        ("weight_top", "weight_top"),
        ("weight_distance", "weight_distance"),
        ("apply_output", "apply_output"),
    ]:
        val = getattr(args, attr, None)
        if val is not None:
            overrides[key] = val
    return overrides


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_easy_learn_apply(args) -> int:
    """Main entry point for 'snekmer easy-learn-apply'.

    Runs without any prompts when --train, --query, and an annotation source
    (--ann or --create-ann) are all supplied.  Any missing value is prompted
    for interactively.  K-mer parameters always use defaults (k=8, alphabet=2)
    unless --k or --alphabet flags are explicitly provided.
    """
    # Deferred import avoids circular dependency (easy <- cli <- easy)
    from snekmer.cli import _run  # noqa: PLC0415

    has_ann = bool(getattr(args, "ann", None) or getattr(args, "create_ann", False))
    interactive = not (getattr(args, "train", None) and getattr(args, "query", None) and has_ann)

    if interactive:
        print("\n=== Snekmer easy-learn-apply ===\n")

    # ---- Collect inputs (wizard steps 1 → 2 → 3 → 4) ----------------------
    try:
        train_files, query_files = _wizard_inputs(args, interactive)      # 1, 2
        ann_tmp = _wizard_annotation(args, interactive, train_files)       # 3
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    output_dir = _wizard_output_dir(args)                                  # 4

    # k and alphabet use defaults silently; only applied if flags were given
    k = int(getattr(args, "k", None) or 8)
    alphabet = getattr(args, "alphabet", None) or 2
    la_overrides = _collect_la_overrides(args)

    # ---- Scaffold workspace ------------------------------------------------
    if interactive:
        print(f"\nSetting up workspace at: {output_dir}")

    learn_dir, apply_dir = setup_workspace(
        output_dir=output_dir,
        train_files=train_files,
        query_files=query_files,
        ann_file=ann_tmp,
        k=k,
        alphabet=alphabet,
        la_overrides=la_overrides or None,
        symlink=not getattr(args, "copy_files", False),
    )
    ann_tmp.unlink(missing_ok=True)

    if interactive:
        print(f"  learn/ : {learn_dir}")
        print(f"  apply/ : {apply_dir}")

    # ---- Run learn ---------------------------------------------------------
    print("\n--- Running snekmer learn ---")
    rc = _run("learn", _make_run_args(args, str(learn_dir)), keepgoing_override=True)
    if rc != 0:
        print(f"\nERROR: snekmer learn failed (exit code {rc})", file=sys.stderr)
        return rc

    # In dry-run mode the learn DAG was printed but no files were created.
    # Show the apply DAG too and exit cleanly — there is nothing to hand off.
    if getattr(args, "dryrun", False):
        print("\n--- Running snekmer apply (dry run) ---")
        _run("apply", _make_run_args(args, str(apply_dir)), keepgoing_override=True)
        return 0

    # ---- Copy handoff files ------------------------------------------------
    print("\nCopying learn outputs to apply input directories...")
    try:
        handoff_learn_to_apply(learn_dir, apply_dir)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if interactive:
        print("  Done.")

    # ---- Run apply ---------------------------------------------------------
    print("\n--- Running snekmer apply ---")
    rc = _run("apply", _make_run_args(args, str(apply_dir)), keepgoing_override=True)
    if rc != 0:
        print(f"\nERROR: snekmer apply failed (exit code {rc})", file=sys.stderr)
        return rc

    # ---- Report ------------------------------------------------------------
    results_filename = la_overrides.get("apply_output", "snekmer_results.csv")
    results_file = apply_dir / results_filename
    print("\n=== Complete ===")
    print(f"Results: {results_file}")
    return 0
