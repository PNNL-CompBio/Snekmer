#!/usr/bin/env bash
set -euo pipefail

# Activate the venv
source /opt/snekmer_env/bin/activate

# Determine mode: "entrypoint.sh learn ..." or symlink-named "learn"
name=$(basename "$0")
if [ "$name" = "entrypoint.sh" ]; then
  mode=${1:-"--help"}
  shift || true
else
  mode=$name
fi

[ -d /work ] && cd /work || true

# Dispatch to the CLI (e.g., snekmer learn/apply/...)
exec snekmer "$mode" "$@"

