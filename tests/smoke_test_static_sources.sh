#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

while IFS= read -r script; do
  bash -n "$script"
done < <(find "$repo_root/scripts_DoNotTouch" "$repo_root/tests" -type f -name '*.sh' -print)

python3 - "$repo_root" <<'PY'
import ast
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
for path in sorted((root / "scripts_DoNotTouch").rglob("*.py")):
    ast.parse(path.read_text(), filename=str(path))
PY

Rscript -e 'root <- commandArgs(TRUE)[1]; files <- list.files(file.path(root, "scripts_DoNotTouch"), pattern="[.][Rr]$", recursive=TRUE, full.names=TRUE); for (file in files) parse(file=file)' "$repo_root"

echo "Shell, Python, and R source syntax checks passed."
