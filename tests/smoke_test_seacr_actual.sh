#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring_seacr_smoke.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT
module() { :; }
export -f module

target="$work/target.bedgraph"
for i in $(seq 0 49); do
  start=$((i * 100))
  signal=$((i % 10 + 1))
  printf 'chr1\t%s\t%s\t%s\n' "$start" "$((start + 50))" "$signal" >> "$target"
done

(
  cd "$work"
  bash "$repo_root/scripts_DoNotTouch/SEACR/seacr_cutrun.sh" \
    "$target" none non stringent "$work/test_sample" smoke none
)

[[ -f "$work/test_sample.stringent.bed" ]] || {
  echo "ASSERTION FAILED: SEACR did not create the expected stringent peak file" >&2
  exit 1
}
[[ -s "$work/test_sample_seacr_summary.txt" ]] || {
  echo "ASSERTION FAILED: SEACR did not create its summary" >&2
  exit 1
}
grep -q $'^peak_bed\t' "$work/test_sample_seacr_summary.txt"
grep -q $'^peak_count\t' "$work/test_sample_seacr_summary.txt"
if find "$work" -maxdepth 1 -type f \( -name '*.auc' -o -name '*.threshold.txt' -o -name '*.fdr.txt' \) | grep -q .; then
  echo "ASSERTION FAILED: SEACR left temporary calculation files behind" >&2
  exit 1
fi
if find "$work" -maxdepth 1 -type d -name '.seacr_tmp.*' | grep -q .; then
  echo "ASSERTION FAILED: SEACR temporary directory was not cleaned" >&2
  exit 1
fi

echo "Bundled SEACR actual-algorithm smoke test passed."
