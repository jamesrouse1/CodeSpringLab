#!/usr/bin/env bash
set -Eeuo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring_cellranger_smoke.XXXXXX")"
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/bin" "$work/fastq" "$work/reference" "$work/output"
printf '{}\n' > "$work/reference/reference.json"
printf '@r1\nAAAA\n+\n####\n' | gzip > "$work/fastq/donor1_S1_L001_R1_001.fastq.gz"
printf '@r1\nACGT\n+\n####\n' | gzip > "$work/fastq/donor1_S1_L001_R2_001.fastq.gz"

cat > "$work/bin/cellranger" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
id=""
for arg in "$@"; do [[ "$arg" == --id=* ]] && id="${arg#--id=}"; done
[[ -n "$id" ]]
mkdir -p "$id/outs/filtered_feature_bc_matrix"
if [[ "${FAKE_CELLRANGER_FAIL_ONCE:-0}" == "1" && ! -f "$id/resume_checkpoint" ]]; then
  touch "$id/resume_checkpoint"
  exit 17
fi
printf 'matrix\n' | gzip > "$id/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
printf 'features\n' | gzip > "$id/outs/filtered_feature_bc_matrix/features.tsv.gz"
printf 'barcodes\n' | gzip > "$id/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
printf 'metric,value\nreads,1\n' > "$id/outs/metrics_summary.csv"
printf '<html>summary</html>\n' > "$id/outs/web_summary.html"
printf '%s\n' "$@" > "$id/arguments.txt"
EOF
chmod +x "$work/bin/cellranger"
export PATH="$work/bin:$PATH"
export CELLRANGER_MIN_FREE_GB=0
module() { :; }
export -f module

bash "$repo/scripts_DoNotTouch/singleCellRNAseq/cellranger_count.sh" \
  donor1 "$work/fastq" donor1 "$work/reference" "$work/output" 5000

test -s "$work/output/donor1/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
test -s "$work/output/donor1/donor1_cellranger_summary.txt"
test -f "$work/output/donor1/_CELLRANGER_COMPLETE"
grep -Fq -- '--sample=donor1' "$work/output/donor1/arguments.txt"
grep -Fq -- '--expect-cells=5000' "$work/output/donor1/arguments.txt"
grep -Fq -- '--nosecondary' "$work/output/donor1/arguments.txt"
grep -Fq -- '--create-bam=false' "$work/output/donor1/arguments.txt"
test -s "$work/output/donor1/outs/metrics_summary.csv"
test -s "$work/output/donor1/outs/web_summary.html"
test ! -e "$work/output/.staging_donor1"

# A failed invocation must preserve a stable staging directory, and an
# identical retry must see the checkpoint and complete instead of restarting.
export FAKE_CELLRANGER_FAIL_ONCE=1
if bash "$repo/scripts_DoNotTouch/singleCellRNAseq/cellranger_count.sh" \
  donor2 "$work/fastq" donor1 "$work/reference" "$work/output" 0; then
  echo "Expected the first resumability fixture invocation to fail." >&2
  exit 1
fi
test -f "$work/output/.staging_donor2/cellranger_run/resume_checkpoint"
bash "$repo/scripts_DoNotTouch/singleCellRNAseq/cellranger_count.sh" \
  donor2 "$work/fastq" donor1 "$work/reference" "$work/output" 0
test -f "$work/output/donor2/_CELLRANGER_COMPLETE"
test ! -e "$work/output/.staging_donor2"
echo "Cell Ranger fake-tool smoke test passed."
