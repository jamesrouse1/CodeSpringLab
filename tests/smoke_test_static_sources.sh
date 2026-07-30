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

grep -Fq 'pca_differential_peaks.png' "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"
grep -Fq 'dba.plotPCA(db, contrast = 1L' "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"
grep -Fq 'all_differential_peaks.bed' "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"

python3 - "$repo_root" <<'PY'
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
manifest = root / "scripts_DoNotTouch/test/manifest_atac/design_matrix.txt"
fastq_dir = root / "scripts_DoNotTouch/test/fastq_atac"
with manifest.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

required = {"sample", "condition", "replicate", "filename"}
if not rows or not required.issubset(rows[0]):
    raise SystemExit("ATAC example manifest is missing required sample metadata columns")

samples = [row["sample"].strip() for row in rows]
if any(not sample for sample in samples) or len(samples) != len(set(samples)):
    raise SystemExit("ATAC example sample identifiers must be non-empty and unique")

condition_counts = {}
for row in rows:
    condition = row["condition"].strip()
    replicate = row["replicate"].strip()
    if not condition or not replicate:
        raise SystemExit("ATAC example condition and replicate values must be non-empty")
    condition_counts[condition] = condition_counts.get(condition, 0) + 1
    files = [item.strip() for item in row["filename"].split(",") if item.strip()]
    if len(files) != 2 or not all((fastq_dir / item).is_file() for item in files):
        raise SystemExit(f"ATAC example FASTQ pair does not resolve for {row['sample']}")

if len(condition_counts) < 2 or any(count < 2 for count in condition_counts.values()):
    raise SystemExit("ATAC example needs at least two replicates in each comparison condition")
PY

echo "Shell, Python, and R source syntax checks passed."
