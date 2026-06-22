"""Tests for snekmer/easy.py."""
import json
from pathlib import Path

import pytest

from snekmer.easy import (
    _coerce_alphabet,
    _collect_fasta_files,
    _write_config,
    generate_ann_from_fasta,
    handoff_learn_to_apply,
    setup_workspace,
    validate_inputs,
    validate_params,
)


# ---------------------------------------------------------------------------
# _collect_fasta_files
# ---------------------------------------------------------------------------

class TestCollectFastaFiles:
    def test_single_file_returned(self, fasta_file):
        result = _collect_fasta_files(str(fasta_file))
        assert result == [fasta_file.resolve()]

    def test_directory_returns_all_fasta_files(self, fasta_dir):
        result = _collect_fasta_files(str(fasta_dir))
        assert len(result) == 3
        extensions = {p.suffix for p in result}
        assert extensions == {".fasta", ".faa", ".fna"}

    def test_directory_returns_sorted_paths(self, fasta_dir):
        result = _collect_fasta_files(str(fasta_dir))
        assert result == sorted(result)

    def test_gz_file_accepted(self, fasta_gz_file):
        result = _collect_fasta_files(str(fasta_gz_file))
        assert len(result) == 1
        assert str(result[0]).endswith(".gz")

    def test_gz_in_directory(self, tmp_path, fasta_gz_file):
        d = tmp_path / "gzdir"
        d.mkdir()
        import shutil
        shutil.copy2(fasta_gz_file, d / fasta_gz_file.name)
        result = _collect_fasta_files(str(d))
        assert len(result) == 1

    def test_missing_path_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            _collect_fasta_files(str(tmp_path / "nonexistent.fasta"))

    def test_empty_directory_raises(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(ValueError, match="No FASTA files found"):
            _collect_fasta_files(str(empty))

    def test_non_fasta_files_ignored(self, tmp_path):
        d = tmp_path / "mixed"
        d.mkdir()
        (d / "data.csv").write_text("a,b\n")
        (d / "readme.txt").write_text("hello\n")
        with pytest.raises(ValueError, match="No FASTA files found"):
            _collect_fasta_files(str(d))


# ---------------------------------------------------------------------------
# generate_ann_from_fasta
# ---------------------------------------------------------------------------

class TestGenerateAnnFromFasta:
    def test_writes_header_and_family_rows(self, fasta_file, tmp_path):
        out = tmp_path / "out.ann"
        n = generate_ann_from_fasta([fasta_file], out)
        lines = out.read_text().splitlines()
        assert lines[0] == "id\tfamily"
        assert n == 2  # FAM001 and FAM002

    def test_id_equals_family(self, fasta_file, tmp_path):
        out = tmp_path / "out.ann"
        generate_ann_from_fasta([fasta_file], out)
        data_lines = out.read_text().splitlines()[1:]
        for line in data_lines:
            id_col, fam_col = line.split("\t")
            assert id_col == fam_col

    def test_deduplication(self, fasta_file, tmp_path):
        """Two sequences from FAM001 must produce only one FAM001 row."""
        out = tmp_path / "out.ann"
        n = generate_ann_from_fasta([fasta_file], out)
        assert n == 2  # FAM001 deduplicated, FAM002 unique
        lines = [l for l in out.read_text().splitlines() if l.strip() and l != "id\tfamily"]
        assert len(lines) == 2

    def test_multiple_files(self, fasta_dir, tmp_path):
        files = _collect_fasta_files(str(fasta_dir))
        out = tmp_path / "out.ann"
        n = generate_ann_from_fasta(files, out)
        assert n == 2  # FAM001 and FAM002

    def test_no_pipes_returns_zero(self, fasta_no_pipes, tmp_path):
        out = tmp_path / "out.ann"
        n = generate_ann_from_fasta([fasta_no_pipes], out)
        assert n == 0
        lines = out.read_text().splitlines()
        assert lines == ["id\tfamily"]

    def test_creates_parent_dirs(self, fasta_file, tmp_path):
        out = tmp_path / "nested" / "deep" / "out.ann"
        generate_ann_from_fasta([fasta_file], out)
        assert out.exists()

    def test_rows_sorted(self, fasta_file, tmp_path):
        out = tmp_path / "out.ann"
        generate_ann_from_fasta([fasta_file], out)
        families = [l.split("\t")[0] for l in out.read_text().splitlines()[1:] if l]
        assert families == sorted(families)


# ---------------------------------------------------------------------------
# _coerce_alphabet
# ---------------------------------------------------------------------------

class TestCoerceAlphabet:
    def test_numeric_string_becomes_int(self):
        assert _coerce_alphabet("2") == 2
        assert _coerce_alphabet("0") == 0

    def test_named_string_unchanged(self):
        assert _coerce_alphabet("solvacc") == "solvacc"
        assert _coerce_alphabet("hydro") == "hydro"

    def test_none_string_unchanged(self):
        assert _coerce_alphabet("None") == "None"

    def test_int_unchanged(self):
        assert _coerce_alphabet(2) == 2

    def test_none_string_case_insensitive(self):
        assert _coerce_alphabet("none") == "none"


# ---------------------------------------------------------------------------
# _write_config
# ---------------------------------------------------------------------------

class TestWriteConfig:
    def test_creates_valid_json(self, tmp_path):
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, k=8, alphabet=2)
        cfg = json.loads(cfg_path.read_text())
        assert cfg["k"] == 8
        assert cfg["alphabet"] == 2

    def test_required_keys_present(self, tmp_path):
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, k=8, alphabet=2)
        cfg = json.loads(cfg_path.read_text())
        for key in ("k", "alphabet", "input_dir", "learn_apply"):
            assert key in cfg

    def test_la_overrides_applied(self, tmp_path):
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, k=8, alphabet=2, la_overrides={"selection": "greatest_distance"})
        cfg = json.loads(cfg_path.read_text())
        assert cfg["learn_apply"]["selection"] == "greatest_distance"

    def test_alphabet_coerced(self, tmp_path):
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, k=8, alphabet="2")
        cfg = json.loads(cfg_path.read_text())
        assert cfg["alphabet"] == 2  # string "2" → int

    def test_default_la_values(self, tmp_path):
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, k=8, alphabet=2)
        la = json.loads(cfg_path.read_text())["learn_apply"]
        assert la["selection"] == "top_hit"
        assert la["threshold"] == "Median"


# ---------------------------------------------------------------------------
# setup_workspace
# ---------------------------------------------------------------------------

class TestSetupWorkspace:
    def test_creates_learn_and_apply_dirs(self, tmp_path, fasta_file, ann_file):
        learn_dir, apply_dir = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
        )
        assert learn_dir.is_dir()
        assert apply_dir.is_dir()

    def test_expected_subdirs_exist(self, tmp_path, fasta_file, ann_file):
        learn_dir, apply_dir = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
        )
        for subdir in ("input", "annotations"):
            assert (learn_dir / subdir).is_dir()
        for subdir in ("input", "counts", "confidence", "stats"):
            assert (apply_dir / subdir).is_dir()

    def test_train_files_linked_to_learn_input(self, tmp_path, fasta_file, ann_file):
        learn_dir, _ = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
        )
        linked = learn_dir / "input" / fasta_file.name
        assert linked.exists()

    def test_query_files_linked_to_apply_input(self, tmp_path, fasta_file, ann_file):
        _, apply_dir = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
        )
        linked = apply_dir / "input" / fasta_file.name
        assert linked.exists()

    def test_ann_file_copied(self, tmp_path, fasta_file, ann_file):
        learn_dir, _ = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
        )
        assert (learn_dir / "annotations" / "annotations.ann").exists()

    def test_configs_written(self, tmp_path, fasta_file, ann_file):
        learn_dir, apply_dir = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
        )
        assert (learn_dir / "config.yaml").exists()
        assert (apply_dir / "config.yaml").exists()

    def test_copy_files_flag(self, tmp_path, fasta_file, ann_file):
        learn_dir, _ = setup_workspace(
            output_dir=tmp_path / "out",
            train_files=[fasta_file],
            query_files=[fasta_file],
            ann_file=ann_file,
            k=8,
            alphabet=2,
            symlink=False,
        )
        linked = learn_dir / "input" / fasta_file.name
        assert linked.exists()
        assert not linked.is_symlink()


# ---------------------------------------------------------------------------
# handoff_learn_to_apply
# ---------------------------------------------------------------------------

class TestHandoffLearnToApply:
    def test_copies_three_files(self, tmp_path, learn_output_dir):
        apply_dir = tmp_path / "apply"
        for d in ("counts", "stats", "confidence"):
            (apply_dir / d).mkdir(parents=True)

        handoff_learn_to_apply(learn_output_dir, apply_dir)

        assert (apply_dir / "counts" / "kmer_counts_total.csv").exists()
        assert (apply_dir / "stats" / "family_summary_stats.csv").exists()
        assert (apply_dir / "confidence" / "global_confidence_scores.csv").exists()

    def test_missing_source_raises(self, tmp_path):
        learn_dir = tmp_path / "learn_empty"
        learn_dir.mkdir()
        apply_dir = tmp_path / "apply"
        apply_dir.mkdir()

        with pytest.raises(FileNotFoundError, match="Expected learn output not found"):
            handoff_learn_to_apply(learn_dir, apply_dir)

    def test_content_preserved(self, tmp_path, learn_output_dir):
        apply_dir = tmp_path / "apply"
        for d in ("counts", "stats", "confidence"):
            (apply_dir / d).mkdir(parents=True)

        handoff_learn_to_apply(learn_output_dir, apply_dir)

        content = (apply_dir / "counts" / "kmer_counts_total.csv").read_text()
        assert "FAM001" in content


# ---------------------------------------------------------------------------
# validate_inputs
# ---------------------------------------------------------------------------

class TestValidateInputs:
    def test_valid_inputs_no_errors(self, fasta_file, ann_file):
        warns, errs = validate_inputs([fasta_file], [fasta_file], ann_file, k=4)
        assert not errs

    def test_empty_fasta_error(self, fasta_empty, fasta_file, ann_file):
        _, errs = validate_inputs([fasta_empty], [fasta_file], ann_file, k=4)
        assert any("no sequences" in e for e in errs)

    def test_dna_fasta_warning(self, fasta_dna, fasta_file, ann_file):
        warns, _ = validate_inputs([fasta_file], [fasta_dna], ann_file, k=4)
        assert any("nucleotide" in w for w in warns)

    def test_sequences_shorter_than_k_warning(self, fasta_file, ann_file):
        # fasta_file sequences are 16–20 aa; k=100 makes them all "short"
        warns, _ = validate_inputs([fasta_file], [fasta_file], ann_file, k=100)
        assert any("shorter than k" in w for w in warns)

    def test_ann_wrong_columns_error(self, fasta_file, ann_wrong_cols):
        _, errs = validate_inputs([fasta_file], [fasta_file], ann_wrong_cols, k=4)
        assert any("family" in e for e in errs)

    def test_ann_comma_separated_error(self, fasta_file, ann_comma):
        _, errs = validate_inputs([fasta_file], [fasta_file], ann_comma, k=4)
        assert any("comma" in e.lower() for e in errs)

    def test_ann_no_data_rows_error(self, fasta_file, ann_no_data):
        _, errs = validate_inputs([fasta_file], [fasta_file], ann_no_data, k=4)
        assert any("no data rows" in e for e in errs)

    def test_no_annotation_id_match_error(self, fasta_no_pipes, ann_file):
        # fasta_no_pipes has no pipe-enclosed fields → zero matches
        _, errs = validate_inputs([fasta_no_pipes], [fasta_no_pipes], ann_file, k=4)
        assert any("No training sequences matched" in e for e in errs)

    def test_single_family_error(self, fasta_file, ann_one_family):
        # ann_one_family has FAM001 as id (matches fasta_file headers) but only 1 family
        _, errs = validate_inputs([fasta_file], [fasta_file], ann_one_family, k=4)
        assert any("family/families" in e for e in errs)

    def test_query_errors_reported_separately(self, fasta_file, fasta_empty, ann_file):
        _, errs = validate_inputs([fasta_file], [fasta_empty], ann_file, k=4)
        assert any("Query" in e for e in errs)
        assert not any("Training" in e for e in errs)


# ---------------------------------------------------------------------------
# validate_params
# ---------------------------------------------------------------------------

class TestValidateParams:
    def test_valid_defaults_no_issues(self):
        warns, errs = validate_params(k=8, alphabet=2, la_overrides={})
        assert not errs
        assert not warns

    def test_k_zero_error(self):
        _, errs = validate_params(k=0, alphabet=2, la_overrides={})
        assert any("k must be" in e for e in errs)

    def test_k_negative_error(self):
        _, errs = validate_params(k=-5, alphabet=2, la_overrides={})
        assert any("k must be" in e for e in errs)

    def test_large_vocab_error(self):
        # miqs (alphabet=5, 10 symbols) + k=8 → 100M k-mers
        _, errs = validate_params(k=8, alphabet=5, la_overrides={})
        assert any("Vocabulary too large" in e for e in errs)

    def test_medium_vocab_warning_not_error(self):
        # standard (alphabet=1, 7 symbols) + k=8 → 5.7M k-mers: warn only
        warns, errs = validate_params(k=8, alphabet=1, la_overrides={})
        assert any("Large vocabulary" in w for w in warns)
        assert not errs

    def test_invalid_selection_error(self):
        _, errs = validate_params(k=8, alphabet=2, la_overrides={"selection": "bogus"})
        assert any("selection method" in e for e in errs)

    def test_valid_selection_methods_accepted(self):
        for sel in ("top_hit", "greatest_distance", "combined_distance"):
            _, errs = validate_params(k=8, alphabet=2, la_overrides={"selection": sel})
            assert not any("selection" in e for e in errs), f"Unexpected error for '{sel}'"

    def test_combined_distance_bad_weights_warning(self):
        warns, _ = validate_params(
            k=8, alphabet=2,
            la_overrides={"selection": "combined_distance", "weight_top": 0.9, "weight_distance": 0.9},
        )
        assert any("poorly calibrated" in w for w in warns)

    def test_combined_distance_good_weights_no_warning(self):
        warns, _ = validate_params(
            k=8, alphabet=2,
            la_overrides={"selection": "combined_distance", "weight_top": 0.7, "weight_distance": 0.3},
        )
        assert not any("weight" in w for w in warns)

    def test_frag_length_less_than_k_error(self):
        _, errs = validate_params(
            k=8, alphabet=2,
            la_overrides={"fragmentation": True, "frag_length": 4},
        )
        assert any("frag" in e.lower() for e in errs)

    def test_frag_length_equal_to_k_ok(self):
        _, errs = validate_params(
            k=8, alphabet=2,
            la_overrides={"fragmentation": True, "frag_length": 8},
        )
        assert not any("frag" in e.lower() for e in errs)

    def test_fragmentation_off_ignores_frag_length(self):
        _, errs = validate_params(
            k=8, alphabet=2,
            la_overrides={"fragmentation": False, "frag_length": 2},
        )
        assert not errs
