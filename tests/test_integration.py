"""Integration tests using demo data.

Run with:  pytest -m integration
Skip with: pytest -m "not integration"

Test map
--------
test_baseline              easy with all defaults
test_low_memory_params     k=6, alphabet=0 (hydro) — vocab = 2^6 = 64, memory-safe
test_fragmentation         fragmentation=True via CLI flags
test_create_ann            --create-ann: derive annotations from FASTA headers
test_copy_files            --copy-files: copies rather than symlinks input files
test_dry_run               --dry-run smoke test (no output files expected)
test_selection_method      --selection greatest_distance (non-default)
test_standalone_learn      snekmer learn alone → verify handoff files land correctly
test_standalone_apply      learn then apply separately (manual pipeline)
test_additive_database     two snekmer learn runs merged via base/ directory
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent

_BASE = REPO_ROOT / "resources/demo_sequences/learn_apply_inputs"
DEMO = {
    "train":         _BASE / "learn",
    "train_batch_a": sorted((_BASE / "learn").glob("training_sequences_[1-5].fasta")),
    "train_batch_b": sorted((_BASE / "learn").glob("training_sequences_[6-9].fasta"))
                   + [_BASE / "learn" / "training_sequences_10.fasta"],
    "query":         _BASE / "apply" / "test_sequences_1.fasta",
    "ann":           _BASE / "annotations" / "TIGRFAMs_annotation.ann",
    "train_family_headers": _BASE / "learn_family_headers",
}

TIMEOUT = 600  # seconds — full pipeline can be slow on first run


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _snekmer_exe() -> str:
    """Return the snekmer binary co-located with the active Python interpreter.

    Prefers the script next to sys.executable (same venv) so tests exercise
    the locally installed code, not whatever lands first on PATH.
    """
    candidate = Path(sys.executable).parent / "snekmer"
    if candidate.exists():
        return str(candidate)
    fallback = shutil.which("snekmer")
    if fallback:
        return fallback
    pytest.skip("snekmer not found alongside sys.executable or on PATH")
    return ""  # unreachable; satisfies linter (pytest.skip raises)


def _run(*args, cwd=None):
    """Run a snekmer CLI command and return the CompletedProcess."""
    return subprocess.run(
        [_snekmer_exe(), *args],
        capture_output=True,
        text=True,
        timeout=TIMEOUT,
        cwd=cwd,
        check=False,
    )


def _check_results(results_csv: Path, min_rows: int = 1) -> pd.DataFrame:
    """Assert results CSV has required columns and at least min_rows rows."""
    assert results_csv.exists(), f"Results file not found: {results_csv}"
    df = pd.read_csv(results_csv)
    required = {"Sequence", "Prediction", "Score", "delta", "Confidence"}
    missing = required - set(df.columns)
    assert not missing, f"Results missing columns: {missing}"
    assert len(df) >= min_rows, f"Expected >= {min_rows} rows, got {len(df)}"
    return df


def _setup_learn_dir(dest: Path, fasta_files, ann_file: Path) -> None:
    """Populate a learn workspace: input/ + annotations/."""
    (dest / "input").mkdir(parents=True)
    (dest / "annotations").mkdir()
    for f in fasta_files:
        shutil.copy2(f, dest / "input" / Path(f).name)
    shutil.copy2(ann_file, dest / "annotations" / "annotations.ann")


def _run_learn(learn_dir: Path) -> subprocess.CompletedProcess:
    return _run("learn", "--no-default-configfile", "-d", str(learn_dir), "--cores", "1")


# ---------------------------------------------------------------------------
# Session-scoped guard
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def demo_guard():
    """Skip integration tests if demo data or the snekmer CLI is absent."""
    for path in (DEMO["train"], DEMO["query"], DEMO["ann"], DEMO["train_family_headers"]):
        if not path.exists():
            pytest.skip(f"Demo data missing: {path}")
    if not shutil.which("snekmer"):
        pytest.skip("snekmer not found in PATH — install it first")


# ---------------------------------------------------------------------------
# easy tests
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.usefixtures("demo_guard")
class TestEasyLearnApply:
    """End-to-end tests via snekmer easy."""

    def test_baseline(self, tmp_path):
        """Default parameters (k=8, alphabet=2/solvacc) produce valid predictions."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train"]),
            "--query", str(DEMO["query"]),
            "--ann",   str(DEMO["ann"]),
            "--output", str(out),
            "--cores", "1",
        )
        assert r.returncode == 0, f"easy failed:\n{r.stderr}"
        df = _check_results(out / "apply" / "snekmer_results.csv", min_rows=10)
        assert (df["Score"] > 0).any(), "All scores zero — pipeline may be broken"
        assert (df["Confidence"] > 0).any(), "All confidence values zero"

    def test_low_memory_params(self, tmp_path):
        """k=6 + alphabet=0 (hydro): vocab 2^6=64 k-mers — verifies small-alphabet path."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train"]),
            "--query", str(DEMO["query"]),
            "--ann",   str(DEMO["ann"]),
            "--output", str(out),
            "--k", "6",
            "--alphabet", "0",
            "--cores", "1",
        )
        assert r.returncode == 0, f"easy k=6 alphabet=0 failed:\n{r.stderr}"
        _check_results(out / "apply" / "snekmer_results.csv", min_rows=10)

    def test_fragmentation(self, tmp_path):
        """Fragmentation splits sequences into sub-fragments before kmerization."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train"]),
            "--query", str(DEMO["query"]),
            "--ann",   str(DEMO["ann"]),
            "--output", str(out),
            "--fragmentation",
            "--frag-length", "50",
            "--min-length",  "50",
            "--cores", "1",
        )
        assert r.returncode == 0, f"easy --fragmentation failed:\n{r.stderr}"
        _check_results(out / "apply" / "snekmer_results.csv", min_rows=10)

    def test_create_ann(self, tmp_path):
        """--create-ann derives the annotation file from >db|FAMILY|seqid FASTA headers."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train_family_headers"]),
            "--query", str(DEMO["query"]),
            "--create-ann",
            "--output", str(out),
            "--cores", "1",
        )
        assert r.returncode == 0, f"easy --create-ann failed:\n{r.stderr}"
        _check_results(out / "apply" / "snekmer_results.csv", min_rows=10)

    def test_copy_files(self, tmp_path):
        """--copy-files should copy inputs into workspace rather than symlink them."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train"]),
            "--query", str(DEMO["query"]),
            "--ann",   str(DEMO["ann"]),
            "--output", str(out),
            "--copy-files",
            "--cores", "1",
        )
        assert r.returncode == 0, f"easy --copy-files failed:\n{r.stderr}"
        _check_results(out / "apply" / "snekmer_results.csv", min_rows=10)

        learn_inputs = list((out / "learn" / "input").glob("*.fasta"))
        assert learn_inputs, "No FASTA files in learn/input/"
        symlinks = [f for f in learn_inputs if f.is_symlink()]
        assert not symlinks, f"--copy-files produced symlinks: {symlinks}"

    def test_dry_run(self, tmp_path):
        """--dry-run should exit 0 without writing any output files."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train"]),
            "--query", str(DEMO["query"]),
            "--ann",   str(DEMO["ann"]),
            "--output", str(out),
            "--dry-run",
            "--cores", "1",
        )
        assert r.returncode == 0, f"easy --dry-run failed:\n{r.stderr}"
        assert not (out / "apply" / "snekmer_results.csv").exists(), (
            "dry-run should not produce snekmer_results.csv"
        )

    def test_selection_method(self, tmp_path):
        """greatest_distance selection assigns by largest gap from competitors."""
        out = tmp_path / "output"
        r = _run(
            "easy",
            "--train", str(DEMO["train"]),
            "--query", str(DEMO["query"]),
            "--ann",   str(DEMO["ann"]),
            "--output", str(out),
            "--selection", "greatest_distance",
            "--cores", "1",
        )
        assert r.returncode == 0, (
            f"easy --selection greatest_distance failed:\n{r.stderr}"
        )
        _check_results(out / "apply" / "snekmer_results.csv", min_rows=10)


# ---------------------------------------------------------------------------
# Direct learn and apply
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.usefixtures("demo_guard")
class TestLearnApply:
    """Tests for snekmer learn and snekmer apply run as separate commands."""

    def test_standalone_learn(self, tmp_path):
        """snekmer learn alone should produce the three handoff files in apply_inputs/."""
        learn_dir = tmp_path / "learn"
        _setup_learn_dir(learn_dir, DEMO["train"].glob("*.fasta"), DEMO["ann"])

        r = _run_learn(learn_dir)
        assert r.returncode == 0, f"snekmer learn failed:\n{r.stderr}"

        for path in (
            learn_dir / "apply_inputs" / "counts"     / "kmer_counts_total.csv",
            learn_dir / "apply_inputs" / "confidence" / "global_confidence_scores.csv",
            learn_dir / "apply_inputs" / "stats"      / "family_summary_stats.csv",
        ):
            assert path.exists(), f"Handoff file missing: {path}"

    def test_standalone_apply(self, tmp_path):
        """snekmer apply run after snekmer learn produces a valid results file."""
        learn_dir = tmp_path / "learn"
        apply_dir = tmp_path / "apply"
        _setup_learn_dir(learn_dir, DEMO["train"].glob("*.fasta"), DEMO["ann"])

        r_learn = _run_learn(learn_dir)
        assert r_learn.returncode == 0, f"snekmer learn failed:\n{r_learn.stderr}"

        handoff = learn_dir / "apply_inputs"

        # Build apply workspace from learn handoff files
        (apply_dir / "input").mkdir(parents=True)
        for sub in ("counts", "confidence", "stats"):
            (apply_dir / sub).mkdir()

        shutil.copy2(DEMO["query"], apply_dir / "input" / DEMO["query"].name)
        shutil.copy2(
            handoff / "counts"     / "kmer_counts_total.csv",
            apply_dir / "counts"   / "kmer_counts_total.csv",
        )
        shutil.copy2(
            handoff / "confidence" / "global_confidence_scores.csv",
            apply_dir / "confidence" / "global_confidence_scores.csv",
        )
        shutil.copy2(
            handoff / "stats"      / "family_summary_stats.csv",
            apply_dir / "stats"    / "family_summary_stats.csv",
        )

        r_apply = _run(
            "apply", "--no-default-configfile",
            "-d", str(apply_dir), "--cores", "1",
        )
        assert r_apply.returncode == 0, f"snekmer apply failed:\n{r_apply.stderr}"
        _check_results(apply_dir / "snekmer_results.csv", min_rows=10)


# ---------------------------------------------------------------------------
# Additive database (multiple learns merged via base/ directory)
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.usefixtures("demo_guard")
class TestAdditiveDatabase:
    """Verify that the base/ mechanism actually incorporates prior learn data.

    Three learn runs are compared:
      A     — batch A files only                (reference)
      B     — batch B files only, no base/      (control)
      B+A   — batch B files + base/ from A      (additive)

    The additive run's total sequence count must equal A + B (±5%),
    proving the merge combined both datasets rather than overwriting one.
    """

    def test_additive_learn(self, tmp_path):
        learn_a  = tmp_path / "learn_a"
        learn_b  = tmp_path / "learn_b"
        learn_ba = tmp_path / "learn_b_plus_a"

        # --- Run A: batch A alone ---
        _setup_learn_dir(learn_a, DEMO["train_batch_a"], DEMO["ann"])
        r_a = _run_learn(learn_a)
        assert r_a.returncode == 0, f"Batch A learn failed:\n{r_a.stderr}"
        counts_a = learn_a / "apply_inputs" / "counts" / "kmer_counts_total.csv"
        conf_a   = learn_a / "apply_inputs" / "confidence" / "global_confidence_scores.csv"
        assert counts_a.exists(), "Batch A kmer_counts_total.csv missing"
        seqs_a = pd.read_csv(counts_a)["Sequence count"].sum()

        # --- Run B (control): batch B alone, no base/ ---
        _setup_learn_dir(learn_b, DEMO["train_batch_b"], DEMO["ann"])
        r_b = _run_learn(learn_b)
        assert r_b.returncode == 0, f"Batch B (control) learn failed:\n{r_b.stderr}"
        counts_b = learn_b / "apply_inputs" / "counts" / "kmer_counts_total.csv"
        assert counts_b.exists(), "Batch B kmer_counts_total.csv missing"
        seqs_b = pd.read_csv(counts_b)["Sequence count"].sum()

        # --- Run B+A (additive): batch B with base/ from A ---
        _setup_learn_dir(learn_ba, DEMO["train_batch_b"], DEMO["ann"])
        (learn_ba / "base" / "counts").mkdir(parents=True)
        (learn_ba / "base" / "confidence").mkdir(parents=True)
        shutil.copy2(counts_a, learn_ba / "base" / "counts" / "kmer_counts_total.csv")
        shutil.copy2(conf_a,   learn_ba / "base" / "confidence" / "global_confidence_scores.csv")

        r_ba = _run_learn(learn_ba)
        assert r_ba.returncode == 0, f"Additive learn (B + base A) failed:\n{r_ba.stderr}"
        counts_ba = learn_ba / "apply_inputs" / "counts" / "kmer_counts_total.csv"
        assert counts_ba.exists(), "Additive kmer_counts_total.csv missing"
        seqs_ba = pd.read_csv(counts_ba)["Sequence count"].sum()

        # Additive run must exceed either individual run
        assert seqs_ba > seqs_a, (
            f"Additive run ({seqs_ba} seqs) should exceed batch A alone ({seqs_a})"
        )
        assert seqs_ba > seqs_b, (
            f"Additive run ({seqs_ba} seqs) should exceed batch B alone ({seqs_b})"
        )
        # Merged total should be close to A + B (within 5% tolerance for any overlap)
        expected = seqs_a + seqs_b
        assert seqs_ba >= expected * 0.95, (
            f"Additive run ({seqs_ba} seqs) should be ~A+B ({expected}); "
            f"base/ merge may not have incorporated prior data"
        )
