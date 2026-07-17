#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring_peak_smoke.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

real_rscript="$(command -v Rscript)"
"$real_rscript" -e 'invisible(parse(file=commandArgs(TRUE)[1])); invisible(parse(file=commandArgs(TRUE)[2]))' \
  "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R" \
  "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind_chip.R"
if grep -q 'example_dataset' "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R"; then
  echo "ASSERTION FAILED: ATAC DiffBind must not fabricate samples based on an example path" >&2
  exit 1
fi
real_python="$(command -v python3)"
"$real_python" -c 'import ast, pathlib, sys; ast.parse(pathlib.Path(sys.argv[1]).read_text())' \
  "$repo_root/scripts_DoNotTouch/bulkChIPseq.py"
if grep -Eq 'GRCm39_M29|gencode\.vM29|gencode\.v42|hg38_p13|mm10' \
  "$repo_root/scripts_DoNotTouch/bulkChIPseq.py"; then
  echo "ASSERTION FAILED: CodeSpringLab ChIP-seq still references a legacy genome" >&2
  exit 1
fi
grep -Fq '^chr([0-9]+|X|Y)$' "$repo_root/scripts_DoNotTouch/bowtie2/bowtie2_PE.sh"
grep -Fq '^chr([0-9]+|X|Y)$' "$repo_root/scripts_DoNotTouch/bowtie2/bowtie2_chip_SE.sh"
while IFS= read -r shell_script; do
  bash -n "$shell_script"
done < <(find "$repo_root/scripts_DoNotTouch" -type f -name '*.sh' -print)

fake_bin="$work/bin"
mkdir -p "$fake_bin"
cat > "$fake_bin/macs2" <<'FAKE_MACS2'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--version" ]]; then echo "macs2 2.2.9.1-fake"; exit 0; fi
name=""; outdir=""; broad=false
while (($#)); do
  case "$1" in
    -n) name="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --broad) broad=true; shift ;;
    *) shift ;;
  esac
done
[[ -n "$name" && -n "$outdir" ]] || { echo "fake macs2 missing name/outdir" >&2; exit 2; }
mkdir -p "$outdir"
if [[ "${FAKE_MACS2_MODE:-success}" == "fail" ]]; then echo "forced failure" >&2; exit 9; fi
ext="narrowPeak"; $broad && ext="broadPeak"
if [[ "${FAKE_MACS2_MODE:-success}" == "zero" ]]; then
  : > "$outdir/${name}_peaks.${ext}"
else
  printf 'chr1\t10\t30\tpeak1\t100\n' > "$outdir/${name}_peaks.${ext}"
fi
printf 'chr\tstart\tend\nchr1\t10\t30\n' > "$outdir/${name}_peaks.xls"
if [[ "${FAKE_MACS2_MODE:-success}" == "traceback" ]]; then
  printf 'Traceback (most recent call last):\nOSError: No space left on device\n' >&2
fi
FAKE_MACS2
cat > "$fake_bin/computeMatrix" <<'FAKE_MATRIX'
#!/usr/bin/env bash
set -euo pipefail
out=""; regions=""
while (($#)); do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    --outFileSortedRegions) regions="$2"; shift 2 ;;
    *) shift ;;
  esac
done
printf 'matrix\n' > "$out"
printf 'chr1\t1\t2\n' > "$regions"
FAKE_MATRIX
cat > "$fake_bin/plotHeatmap" <<'FAKE_HEATMAP'
#!/usr/bin/env bash
set -euo pipefail
out=""
while (($#)); do if [[ "$1" == "-out" ]]; then out="$2"; shift 2; else shift; fi; done
printf 'png\n' > "$out"
FAKE_HEATMAP
cat > "$fake_bin/annotatePeaks.pl" <<'FAKE_HOMER'
#!/usr/bin/env bash
if [[ "${FAKE_HOMER_MODE:-success}" == "fail" ]]; then echo "forced annotation failure" >&2; exit 7; fi
printf 'PeakID\tAnnotation\npeak1\tPromoter\n'
FAKE_HOMER
cat > "$fake_bin/getDifferentialPeaksReplicates.pl" <<'FAKE_HOMER_DIFF'
#!/usr/bin/env bash
set -euo pipefail
printf 'PeakID\tChr\tStart\tEnd\npeak1\tchr1\t10\t30\n'
FAKE_HOMER_DIFF
cat > "$fake_bin/Rscript" <<'FAKE_RSCRIPT'
#!/usr/bin/env bash
set -euo pipefail
script="$1"
if [[ "$(basename "$script")" == "DiffBind_chip.R" ]]; then
  outdir="$2"; reference="$4"; comparison="$5"
else
  outdir="$2"; reference="$5"; comparison="$6"
fi
mkdir -p "$outdir"
prefix="DifferentialPeaks_${comparison}_vs_${reference}_ref"
printf 'Fold\tFDR\n1\t0.01\n' > "$outdir/${prefix}.txt"
if [[ "${FAKE_DIFFBIND_MODE:-success}" == "zero" ]]; then
  : > "$outdir/${prefix}.with_stats.bed"
else
  printf 'chr1\t10\t30\tpeak1\t1\n' > "$outdir/${prefix}.with_stats.bed"
fi
if [[ "$(basename "$script")" == "DiffBind_chip.R" ]]; then
  printf 'status\tcomplete\n' > "$outdir/_DIFFBIND_COMPLETE"
fi
FAKE_RSCRIPT
chmod +x "$fake_bin"/*
export PATH="$fake_bin:$PATH"
module() { :; }
export -f module

assert_file() { [[ -s "$1" ]] || { echo "ASSERTION FAILED: expected non-empty file $1" >&2; exit 1; }; }
assert_absent() { [[ ! -e "$1" ]] || { echo "ASSERTION FAILED: expected absent path $1" >&2; exit 1; }; }

atac_root="$work/atac"
atac_out="$atac_root/macs2/S1"
mkdir -p "$atac_out" "$atac_root/bowtie2/S1"
printf 'chr1\t10\t30\n' > "$atac_root/S1.bed"
printf 'chr1\t1\t2\n' > "$atac_root/tss.bed"
printf 'bw\n' > "$atac_root/bowtie2/S1/S1Aligned.sortedByCoord_removeDup.out.bw"
bash "$repo_root/scripts_DoNotTouch/MACS2/macs2_PE.sh" \
  S1 "$atac_root/S1.bed" mm unused "$atac_out" "$atac_root/tss.bed" 0.05 mm39 "$atac_root/bowtie2"
assert_file "$atac_out/S1_macs2_complete.txt"
assert_file "$atac_out/S1_macs2_summary.txt"
assert_file "$atac_out/S1_peaks.narrowPeak"
if find "$atac_out" -maxdepth 1 -name '.macs2_tmp_*' | grep -q .; then
  echo "ASSERTION FAILED: ATAC MACS2 temporary directory was not cleaned" >&2
  exit 1
fi

rm -f "$atac_out/S1_macs2_complete.txt"
export FAKE_MACS2_MODE=traceback
if bash "$repo_root/scripts_DoNotTouch/MACS2/macs2_PE.sh" \
  S1 "$atac_root/S1.bed" mm unused "$atac_out" "$atac_root/tss.bed" 0.05 mm39 "$atac_root/bowtie2" \
  > "$work/atac_traceback.out" 2>&1; then
  echo "ASSERTION FAILED: ATAC MACS2 accepted a fake internal traceback" >&2
  exit 1
fi
grep -q "internal peak-calling exception" "$work/atac_traceback.out"
assert_absent "$atac_out/S1_macs2_complete.txt"
unset FAKE_MACS2_MODE

export FAKE_MACS2_MODE=zero
bash "$repo_root/scripts_DoNotTouch/MACS2/macs2_PE.sh" \
  S1 "$atac_root/S1.bed" mm unused "$atac_out" "$atac_root/tss.bed" 0.05 mm39 "$atac_root/bowtie2"
assert_file "$atac_out/S1_macs2_complete.txt"
grep -q $'^peak_count\t0$' "$atac_out/S1_macs2_complete.txt"
unset FAKE_MACS2_MODE

chip_root="$work/chip"
chip_out="$chip_root/macs2/T1"
mkdir -p "$chip_out"
printf 'bam\n' > "$chip_root/target.bam"
printf 'bam\n' > "$chip_root/input.bam"
bash "$repo_root/scripts_DoNotTouch/MACS2/macs2_chip_SE.sh" \
  T1 "$chip_root/target.bam" "$chip_root/input.bam" BAM mm 0.01 narrow "$chip_out"
assert_file "$chip_out/T1_macs2_complete.txt"
assert_file "$chip_out/T1_macs2_summary.txt"
assert_file "$chip_out/T1_peaks.narrowPeak"

rm -f "$chip_out/T1_macs2_complete.txt"
export FAKE_MACS2_MODE=traceback
if bash "$repo_root/scripts_DoNotTouch/MACS2/macs2_chip_SE.sh" \
  T1 "$chip_root/target.bam" "$chip_root/input.bam" BAM mm 0.01 narrow "$chip_out" \
  > "$work/chip_traceback.out" 2>&1; then
  echo "ASSERTION FAILED: ChIP MACS2 accepted a fake internal traceback" >&2
  exit 1
fi
grep -q "internal exception" "$work/chip_traceback.out"
assert_absent "$chip_out/T1_macs2_complete.txt"
unset FAKE_MACS2_MODE

diff_root="$work/diffbind"
mkdir -p "$diff_root/design" "$diff_root/peaks" "$diff_root/bams"
printf 'sample\tcondition\tfilename\n' > "$diff_root/design/design_matrix.txt"
printf 'fake R code\n' > "$diff_root/DiffBind.R"
bash "$repo_root/scripts_DoNotTouch/DiffBind/diffbind.sh" \
  "$diff_root/DiffBind.R" "$diff_root/success" "$diff_root/design" "$diff_root/peaks" A B mouse "$diff_root/bams"
assert_file "$diff_root/success/_COMPLETE"
assert_absent "$diff_root/success/_RUN_STARTED"
grep -q $'^status\tcomplete$' "$diff_root/success/_COMPLETE"

export FAKE_HOMER_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/DiffBind/diffbind.sh" \
  "$diff_root/DiffBind.R" "$diff_root/annotation_fail" "$diff_root/design" "$diff_root/peaks" A B mouse "$diff_root/bams" \
  > "$work/diffbind_annotation_fail.out" 2>&1; then
  echo "ASSERTION FAILED: ATAC DiffBind accepted a failed annotation" >&2
  exit 1
fi
assert_file "$diff_root/annotation_fail/_RUN_STARTED"
assert_absent "$diff_root/annotation_fail/_COMPLETE"
unset FAKE_HOMER_MODE

export FAKE_DIFFBIND_MODE=zero
bash "$repo_root/scripts_DoNotTouch/DiffBind/diffbind.sh" \
  "$diff_root/DiffBind.R" "$diff_root/zero" "$diff_root/design" "$diff_root/peaks" A B mouse "$diff_root/bams"
assert_file "$diff_root/zero/_COMPLETE"
assert_file "$diff_root/zero/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt"
grep -q $'^status\tcomplete$' "$diff_root/zero/_COMPLETE"
unset FAKE_DIFFBIND_MODE

chip_diff_root="$work/chip_diffbind"
mkdir -p "$chip_diff_root"
printf 'fake R code\n' > "$chip_diff_root/DiffBind_chip.R"
printf 'SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller\n' > "$chip_diff_root/samples.csv"
bash "$repo_root/scripts_DoNotTouch/DiffBind/diffbind_chip.sh" \
  "$chip_diff_root/DiffBind_chip.R" "$chip_diff_root/success" "$chip_diff_root/samples.csv" \
  A B mouse none
assert_file "$chip_diff_root/success/_COMPLETE"
assert_absent "$chip_diff_root/success/_RUN_STARTED"
assert_absent "$chip_diff_root/success/_DIFFBIND_COMPLETE"
grep -q $'^status\tcomplete$' "$chip_diff_root/success/_COMPLETE"

export FAKE_HOMER_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/DiffBind/diffbind_chip.sh" \
  "$chip_diff_root/DiffBind_chip.R" "$chip_diff_root/annotation_fail" "$chip_diff_root/samples.csv" \
  A B mouse none > "$work/chip_diffbind_annotation_fail.out" 2>&1; then
  echo "ASSERTION FAILED: ChIP DiffBind accepted a failed annotation" >&2
  exit 1
fi
assert_file "$chip_diff_root/annotation_fail/_RUN_STARTED"
assert_absent "$chip_diff_root/annotation_fail/_COMPLETE"
assert_absent "$chip_diff_root/annotation_fail/_DIFFBIND_COMPLETE"
unset FAKE_HOMER_MODE

homer_out="$work/homer_diff"
bash "$repo_root/scripts_DoNotTouch/Homer/homer_diffpeak.sh" \
  "$homer_out" ref_tags target_tags mm39 control treated
assert_file "$homer_out/DiffPeak_treated_vs_control(ref).txt"
assert_file "$homer_out/DiffPeak_treated_vs_control(ref)_annotated.txt"

rm -f "$homer_out/DiffPeak_treated_vs_control(ref).txt" "$homer_out/DiffPeak_treated_vs_control(ref)_annotated.txt"
export FAKE_HOMER_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/Homer/homer_diffpeak.sh" \
  "$homer_out" ref_tags target_tags mm39 control treated > "$work/homer_annotation_fail.out" 2>&1; then
  echo "ASSERTION FAILED: HOMER differential runner accepted a failed annotation" >&2
  exit 1
fi
assert_absent "$homer_out/DiffPeak_treated_vs_control(ref).txt"
assert_absent "$homer_out/DiffPeak_treated_vs_control(ref)_annotated.txt"
unset FAKE_HOMER_MODE

echo "Peak-runner fake-data smoke tests passed."
