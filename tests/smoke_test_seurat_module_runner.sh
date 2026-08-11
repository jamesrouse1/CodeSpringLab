#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
runner="$repo_root/scripts_DoNotTouch/singleCellRNAseq/scrna_pipeline.sh"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring-seurat-module.XXXXXX")"
trap 'rm -rf "$work"' EXIT

mkdir -p "$work/bin" "$work/inputs" "$work/output"
printf 'sample_id\tinput_path\nexample\t%s\n' "$work/inputs/example" > "$work/samples.tsv"
printf 'key\tvalue\nnormalization\tsct\n' > "$work/params.tsv"
mkdir -p "$work/inputs/example"

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -Eeuo pipefail' \
  'printf "%s\n" "$*" >> "${CSL_TEST_MODULE_LOG:?}"' > "$work/bin/module"
chmod +x "$work/bin/module"

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -Eeuo pipefail' \
  'printf "R_LIBS_USER=%s\nR_ENVIRON_USER=%s\nR_PROFILE_USER=%s\nARGS=%s\n" "${R_LIBS_USER:-}" "${R_ENVIRON_USER:-}" "${R_PROFILE_USER:-}" "$*" >> "${CSL_TEST_RSCRIPT_LOG:?}"' > "$work/bin/Rscript"
chmod +x "$work/bin/Rscript"

PATH="$work/bin:$PATH" \
  CSL_TEST_MODULE_LOG="$work/module_args.txt" \
  CSL_TEST_RSCRIPT_LOG="$work/rscript_args.txt" \
  "$runner" seurat "$work/samples.tsv" "$work/output" "$work/params.tsv" inspect

grep -Fxq 'load EB5Modules' "$work/module_args.txt"
grep -Fxq 'load Seurat/5.4.0-foss-2024a-R-4.4.2' "$work/module_args.txt"
grep -Fxq "R_LIBS_USER=$work/output/.codespring_unused_user_library" "$work/rscript_args.txt"
grep -Fxq 'R_ENVIRON_USER=/dev/null' "$work/rscript_args.txt"
grep -Fxq 'R_PROFILE_USER=/dev/null' "$work/rscript_args.txt"
grep -Fq 'scrna_pipeline_seurat.R' "$work/rscript_args.txt"

echo "Seurat module runner smoke test passed."
