#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring_peak_smoke.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

real_rscript="$(command -v Rscript)"
export CSL_REAL_RSCRIPT="$real_rscript"
"$real_rscript" -e 'invisible(parse(file=commandArgs(TRUE)[1])); invisible(parse(file=commandArgs(TRUE)[2]))' \
  "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R" \
  "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind_chip.R"
"$real_rscript" -e 'invisible(parse(file=commandArgs(TRUE)[1]))' \
  "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.R"
if grep -q 'example_dataset' "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R"; then
  echo "ASSERTION FAILED: ATAC DiffBind must not fabricate samples based on an example path" >&2
  exit 1
fi
grep -Fq 'diffbind_sample_sheet.tsv' "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R"
grep -Fq 'dba.plotPCA(dbobject,DBA_CONDITION,label=DBA_ID)' "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R"
grep -Fq 'dba.plotPCA(dbobject,contrast=1,label=DBA_ID)' "$repo_root/scripts_DoNotTouch/DiffBind/DiffBind.R"
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
printf 'PeakID\tAnnotation\n'
if [[ -f "${1:-}" ]]; then
  awk 'BEGIN {OFS="\t"} NF >= 4 {print $4, "Promoter"}' "$1"
else
  printf 'peak1\tPromoter\n'
fi
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
if [[ "$(basename "$script")" == "expand_peakid_stats.R" ]]; then
  exec "$CSL_REAL_RSCRIPT" "$@"
fi
if [[ "$(basename "$script")" == "cutrun_diffbind.R" ]]; then
  printf '%s\n' "$@" > "${FAKE_CUTRUN_DIFFBIND_ARGS:?}"
  exit 0
fi
if [[ "$(basename "$script")" == "DiffBind_chip.R" ]]; then
  outdir="$2"; reference="$4"; comparison="$5"
else
  outdir="$2"; reference="$5"; comparison="$6"
fi
mkdir -p "$outdir"
prefix="DifferentialPeaks_${comparison}_vs_${reference}_ref"
printf 'Fold\tp.value\tFDR\n1\t0.002\t0.01\n' > "$outdir/${prefix}.txt"
if [[ "${FAKE_DIFFBIND_MODE:-success}" == "zero" ]]; then
  : > "$outdir/${prefix}.with_stats.bed"
else
  printf 'chr1\t10\t30\tpeak1|Fold=1|p.value=0.002|FDR=0.01\t1\n' > "$outdir/${prefix}.with_stats.bed"
fi
if [[ "$(basename "$script")" == "DiffBind_chip.R" ]]; then
  printf 'status\tcomplete\n' > "$outdir/_DIFFBIND_COMPLETE"
fi
FAKE_RSCRIPT
chmod +x "$fake_bin"/*
export PATH="$fake_bin:$PATH"
export FAKE_CUTRUN_DIFFBIND_ARGS="$work/cutrun_diffbind_args.txt"
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
assert_file "$atac_out/S1_peaks_annotated.txt"
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

rm -f "$atac_out/S1_macs2_complete.txt" "$atac_out/S1_peaks_annotated.txt"
export FAKE_HOMER_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/MACS2/macs2_PE.sh" \
  S1 "$atac_root/S1.bed" mm unused "$atac_out" "$atac_root/tss.bed" 0.05 mm39 "$atac_root/bowtie2" \
  > "$work/atac_annotation_fail.out" 2>&1; then
  echo "ASSERTION FAILED: ATAC MACS2 accepted a failed integrated HOMER annotation" >&2
  exit 1
fi
assert_absent "$atac_out/S1_macs2_complete.txt"
assert_absent "$atac_out/S1_peaks_annotated.txt"
unset FAKE_HOMER_MODE

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
head -n 1 "$diff_root/success/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt" | grep -q $'p.value\tFDR'

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

cutrun_diff_root="$work/cutrun_diffbind"
mkdir -p "$cutrun_diff_root"
printf 'SampleID\tCellType\tMark\tCondition\tReplicate\tbamReads\tPeaks\tSpikein\tnormalization_mode\n' > "$cutrun_diff_root/samples.tsv"
printf 'fake R code\n' > "$cutrun_diff_root/cutrun_diffbind.R"
bash "$repo_root/scripts_DoNotTouch/DiffBind/qsub_cutrun_diffbind.sh" \
  "$cutrun_diff_root/cutrun_diffbind.R" "$cutrun_diff_root/samples.tsv" "$cutrun_diff_root/output" \
  Veh 1 mouse none AA AKPS Creb 100 shared_test_overlap \
  "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.sh"
grep -Fxq '100' "$FAKE_CUTRUN_DIFFBIND_ARGS"
grep -Fxq 'shared_test_overlap' "$FAKE_CUTRUN_DIFFBIND_ARGS"
if grep -Fxq "$repo_root/scripts_DoNotTouch/DiffBind/cutrun_diffbind.sh" "$FAKE_CUTRUN_DIFFBIND_ARGS"; then
  echo "ASSERTION FAILED: CUT&RUN DiffBind runner path leaked into R arguments" >&2
  exit 1
fi

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
head -n 1 "$chip_diff_root/success/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt" | grep -q $'p.value\tFDR'

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

annotation_root="$work/peak_annotation_project"
mkdir -p \
  "$annotation_root/macs2/S1" \
  "$annotation_root/diffbind/B_vs_A" \
  "$annotation_root/diffbind/legacy" \
  "$annotation_root/cutrun_diffbind/Creb/APC_AA_vs_Veh"
printf 'chr1\t10\t30\tmacs_peak\t100\t.\t12\t8\t6\t5\n' > \
  "$annotation_root/macs2/S1/S1_peaks.narrowPeak"
printf 'chr1\t20\t50\tdiff_peak|Fold=2|p.value=0.001|FDR=0.01\t2\n' > \
  "$annotation_root/diffbind/B_vs_A/DifferentialPeaks_B_vs_A_ref.with_stats.bed"
printf 'seqnames\tstart\tend\twidth\tstrand\tConc\tConc_A\tConc_B\tFold\tp.value\tFDR\nchr1\t101\t140\t40\t+\t5\t4\t6\t2\t0.001\t0.01\n' > \
  "$annotation_root/diffbind/legacy/DifferentialPeaks_B_vs_A_ref.txt"
printf 'seqnames\tstart\tend\tFold\tFDR\nchr1\t31\t60\t-1.5\t0.02\n' > \
  "$annotation_root/cutrun_diffbind/Creb/APC_AA_vs_Veh/all_differential_peaks.tsv"
printf 'chr1\t70\t90\tcutrun_peak|Fold=3|FDR=0.001\t3\n' > \
  "$annotation_root/cutrun_diffbind/Creb/APC_AA_vs_Veh/significant_differential_peaks.bed"
printf 'chr1\tfake\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n' > "$annotation_root/reference.gtf"
bash "$repo_root/scripts_DoNotTouch/Homer/annotate_peak_results.sh" \
  "$annotation_root" mm39 "$annotation_root/reference.gtf"
assert_file "$annotation_root/peak_annotation/_COMPLETE"
assert_absent "$annotation_root/peak_annotation/_RUN_STARTED"
assert_file "$annotation_root/peak_annotation/peak_annotation_summary.tsv"
assert_file "$annotation_root/macs2/S1/S1_peaks_annotated.txt"
assert_file "$annotation_root/diffbind/B_vs_A/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt"
assert_file "$annotation_root/diffbind/legacy/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt"
assert_file "$annotation_root/cutrun_diffbind/Creb/APC_AA_vs_Veh/all_differential_peaks_annotated_with_stats.txt"
assert_file "$annotation_root/cutrun_diffbind/Creb/APC_AA_vs_Veh/significant_differential_peaks_annotated_with_stats.txt"
head -n 1 "$annotation_root/diffbind/B_vs_A/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt" | grep -q $'p.value\tFDR'
head -n 1 "$annotation_root/diffbind/legacy/DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt" | grep -q $'p.value\tFDR'
grep -q $'^annotated_files\t5$' "$annotation_root/peak_annotation/_COMPLETE"

# A rerun with unchanged inputs must reuse the current annotations instead of
# invoking HOMER again. Force the fake HOMER executable to fail so this test
# proves the reuse path is actually taken.
export FAKE_HOMER_MODE=fail
bash "$repo_root/scripts_DoNotTouch/Homer/annotate_peak_results.sh" \
  "$annotation_root" mm39 "$annotation_root/reference.gtf" > "$work/peak_annotation_reuse.out" 2>&1
unset FAKE_HOMER_MODE
grep -q 'Reusing current annotation' "$work/peak_annotation_reuse.out"
assert_file "$annotation_root/peak_annotation/_COMPLETE"
assert_absent "$annotation_root/peak_annotation/_RUN_STARTED"
grep -q $'^annotated_files\t5$' "$annotation_root/peak_annotation/_COMPLETE"

annotation_fail_root="$work/peak_annotation_failure"
mkdir -p "$annotation_fail_root/macs2/S1"
printf 'chr1\t10\t30\tmacs_peak\t100\t.\t12\t8\t6\t5\n' > \
  "$annotation_fail_root/macs2/S1/S1_peaks.narrowPeak"
printf 'chr1\tfake\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n' > "$annotation_fail_root/reference.gtf"
export FAKE_HOMER_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/Homer/annotate_peak_results.sh" \
  "$annotation_fail_root" mm39 "$annotation_fail_root/reference.gtf" > "$work/peak_annotation_failure.out" 2>&1; then
  echo "ASSERTION FAILED: failed peak annotation returned success" >&2
  exit 1
fi
assert_file "$annotation_fail_root/peak_annotation/_RUN_STARTED"
assert_absent "$annotation_fail_root/peak_annotation/_COMPLETE"
unset FAKE_HOMER_MODE

echo "Peak-runner fake-data smoke tests passed."
