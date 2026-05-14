"""Docker integration tests.

Run with:  pytest -m docker
Requires:  docker in PATH.

Tests pull jjacobson95/snekmer:latest from DockerHub and run the pipeline
inside it against local demo data.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent

DEMO = {
    "train": REPO_ROOT / "resources/demo_sequences/learn_apply_inputs/learn",
    "query": REPO_ROOT / "resources/demo_sequences/learn_apply_inputs/apply/test_sequences_1.fasta",
    "ann": REPO_ROOT / "resources/demo_sequences/learn_apply_inputs/annotations/TIGRFAMs_annotation.ann",
}

DOCKERHUB_IMAGE = "jjacobson95/snekmer:latest"
RUN_TIMEOUT = 600  # seconds — pipeline run inside container


# ---------------------------------------------------------------------------
# Session fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def docker_image():
    """Pull jjacobson95/snekmer:latest from DockerHub; skip if docker unavailable."""
    if not shutil.which("docker"):
        pytest.skip("docker not found in PATH")
    for key, path in DEMO.items():
        if not path.exists():
            pytest.skip(f"Demo data missing: {path}")
    result = subprocess.run(
        ["docker", "pull", DOCKERHUB_IMAGE],
        capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        pytest.skip(f"Could not pull {DOCKERHUB_IMAGE}:\n{result.stderr}")
    return DOCKERHUB_IMAGE


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.docker
class TestDockerImage:
    def test_snekmer_version(self, docker_image):
        """snekmer --version should exit 0 inside the container."""
        result = subprocess.run(
            ["docker", "run", "--rm", docker_image, "--version"],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, (
            f"snekmer --version failed inside container:\n{result.stderr}"
        )

    def test_easy_learn_apply(self, tmp_path, docker_image):
        """easy-learn-apply should complete inside the container and produce results."""
        out = tmp_path / "output"
        out.mkdir()

        result = subprocess.run(
            [
                "docker", "run", "--rm",
                "-v", f"{DEMO['train']}:/data/train:ro",
                "-v", f"{DEMO['query'].parent}:/data/apply:ro",
                "-v", f"{DEMO['ann'].parent}:/data/ann:ro",
                "-v", f"{out}:/work/output",
                docker_image,
                "easy-learn-apply",
                "--train",  "/data/train",
                "--query",  f"/data/apply/{DEMO['query'].name}",
                "--ann",    f"/data/ann/{DEMO['ann'].name}",
                "--output", "/work/output",
                "--cores",  "1",
            ],
            capture_output=True,
            text=True,
            timeout=RUN_TIMEOUT,
        )
        assert result.returncode == 0, f"Container run failed:\n{result.stderr}"

        results_csv = out / "apply" / "snekmer_results.csv"
        assert results_csv.exists(), "snekmer_results.csv not produced by container run"

        import pandas as pd
        df = pd.read_csv(results_csv)
        required = {"Sequence", "Prediction", "Score", "delta", "Confidence"}
        assert required <= set(df.columns), (
            f"Missing columns in results: {required - set(df.columns)}"
        )
        assert len(df) >= 10, f"Expected >= 10 result rows, got {len(df)}"
