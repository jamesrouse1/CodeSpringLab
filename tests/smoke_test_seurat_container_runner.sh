#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
runner="$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline.sh"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring-seurat-container.XXXXXX")"
trap 'rm -rf "$work"' EXIT

mkdir -p "$work/bin" "$work/inputs" "$work/output"
printf 'placeholder container image\n' > "$work/codespring-seurat_1.0.0.sif"
printf 'sample_id\tinput_path\nexample\t%s\n' "$work/inputs/example" > "$work/samples.tsv"
printf 'key\tvalue\nnormalization\tsct\n' > "$work/params.tsv"
mkdir -p "$work/inputs/example"

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -Eeuo pipefail' \
  'printf "%s\n" "$@" > "${CSL_TEST_SINGULARITY_LOG:?}"' > "$work/bin/singularity"
chmod +x "$work/bin/singularity"

PATH="$work/bin:$PATH" CSL_TEST_SINGULARITY_LOG="$work/singularity_args.txt" \
  "$runner" seurat "$work/samples.tsv" "$work/output" "$work/params.tsv" inspect "" "$work/codespring-seurat_1.0.0.sif"

grep -Fxq 'exec' "$work/singularity_args.txt"
grep -Fxq -- '--cleanenv' "$work/singularity_args.txt"
grep -Fxq -- "--bind=$work/output" "$work/singularity_args.txt"
grep -Fxq "$work/codespring-seurat_1.0.0.sif" "$work/singularity_args.txt"
grep -Fxq 'Rscript' "$work/singularity_args.txt"
grep -Fq 'scrna_pipeline_seurat.R' "$work/singularity_args.txt"

echo "Shared Seurat container runner smoke test passed."
