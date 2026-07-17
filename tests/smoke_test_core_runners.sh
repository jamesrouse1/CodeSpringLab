#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring_core_smoke.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

fake_bin="$work/bin"
mkdir -p "$fake_bin"
module() { :; }
export -f module

cat > "$fake_bin/STAR" <<'FAKE_STAR'
#!/usr/bin/env bash
set -euo pipefail
prefix=""
while (($#)); do
  case "$1" in
    --outFileNamePrefix) prefix="$2"; shift 2 ;;
    *) shift ;;
  esac
done
[[ -n "$prefix" ]] || { echo "fake STAR missing output prefix" >&2; exit 2; }
if [[ "${FAKE_STAR_MODE:-success}" == "fail" ]]; then exit 9; fi
printf 'bam\n' > "${prefix}Aligned.sortedByCoord.out.bam"
printf 'transcriptome\n' > "${prefix}Aligned.toTranscriptome.out.bam"
printf 'mapping summary\n' > "${prefix}Log.final.out"
FAKE_STAR

cat > "$fake_bin/samtools" <<'FAKE_SAMTOOLS'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
  quickcheck) exit 0 ;;
  index)
    shift
    [[ "${1:-}" == "-b" ]] && shift
    bam="$1"; output="${2:-${bam}.bai}"
    printf 'index\n' > "$output"
    ;;
  *) echo "unsupported fake samtools command: ${1:-}" >&2; exit 2 ;;
esac
FAKE_SAMTOOLS

cat > "$fake_bin/infer_experiment.py" <<'FAKE_INFER'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${FAKE_INFER_MODE:-success}" == "malformed" ]]; then
  printf 'unable to infer\n'
  exit 0
fi
printf 'line1\nline2\nline3\nline4\nforward 0.80\nreverse 0.10\n'
FAKE_INFER

cat > "$fake_bin/featureCounts" <<'FAKE_FEATURECOUNTS'
#!/usr/bin/env bash
set -euo pipefail
out=""
printf '%s\n' "$*" > "${FAKE_FEATURECOUNTS_ARGS:?}"
while (($#)); do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    *) shift ;;
  esac
done
[[ -n "$out" ]] || { echo "fake featureCounts missing output" >&2; exit 2; }
printf 'Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\ngene1\tchr1\t1\t10\t+\t10\t5\n' > "$out"
printf 'Status\tsample.bam\nAssigned\t5\n' > "${out}.summary"
FAKE_FEATURECOUNTS

cat > "$fake_bin/cutadapt" <<'FAKE_CUTADAPT'
#!/usr/bin/env bash
set -euo pipefail
out1=""; out2=""; inputs=()
while (($#)); do
  case "$1" in
    -o) out1="$2"; shift 2 ;;
    -p) out2="$2"; shift 2 ;;
    -j|-m|-a|-A) shift 2 ;;
    *) inputs+=("$1"); shift ;;
  esac
done
cat "${inputs[0]}" > "$out1"
if [[ -n "$out2" ]]; then cat "${inputs[1]}" > "$out2"; fi
FAKE_CUTADAPT

cat > "$fake_bin/kallisto" <<'FAKE_KALLISTO'
#!/usr/bin/env bash
set -euo pipefail
out=""; inputs=()
printf '%s\n' "$*" > "${FAKE_KALLISTO_ARGS:?}"
while (($#)); do
  case "$1" in
    quant|--single) shift ;;
    -i|-t|-o|-b|-l|-s) [[ "$1" == "-o" ]] && out="$2"; shift 2 ;;
    *) inputs+=("$1"); shift ;;
  esac
done
[[ -n "$out" ]] || { echo "fake kallisto missing output directory" >&2; exit 2; }
for input in "${inputs[@]}"; do cat "$input" >/dev/null; done
if [[ "${FAKE_KALLISTO_MODE:-success}" == "fail" ]]; then exit 9; fi
mkdir -p "$out"
printf 'target_id\tlength\teff_length\test_counts\ttpm\ntx1\t100\t80\t5\t10\n' > "$out/abundance.tsv"
printf '{}\n' > "$out/run_info.json"
printf 'h5\n' > "$out/abundance.h5"
FAKE_KALLISTO

cat > "$fake_bin/rsem-calculate-expression" <<'FAKE_RSEM'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" > "${FAKE_RSEM_ARGS:?}"
argc=$#
((argc >= 3)) || exit 2
prefix="${!argc}"
if [[ "${FAKE_RSEM_MODE:-success}" == "fail" ]]; then exit 9; fi
printf 'gene_id\texpected_count\tTPM\tFPKM\ngene1\t5\t10\t8\n' > "${prefix}.genes.results"
printf 'transcript_id\tgene_id\texpected_count\tTPM\tFPKM\ntx1\tgene1\t5\t10\t8\n' > "${prefix}.isoforms.results"
printf 'temp\n' > "${prefix}.transcript.bam"
mkdir -p "${prefix}.stat"
FAKE_RSEM

cat > "$fake_bin/fastqc" <<'FAKE_FASTQC'
#!/usr/bin/env bash
set -euo pipefail
out=""; read=""
while (($#)); do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    -t) shift 2 ;;
    *) read="$1"; shift ;;
  esac
done
base="$(basename "$read")"; base="${base%.gz}"; base="${base%.fastq}"; base="${base%.fq}"
printf 'html\n' > "$out/${base}_fastqc.html"
printf 'zip\n' > "$out/${base}_fastqc.zip"
FAKE_FASTQC

cat > "$fake_bin/fastq_screen" <<'FAKE_SCREEN'
#!/usr/bin/env bash
set -euo pipefail
out=""; read=""
while (($#)); do
  case "$1" in
    -outdir) out="$2"; shift 2 ;;
    -threads) shift 2 ;;
    *) read="$1"; shift ;;
  esac
done
printf 'screen\n' > "$out/$(basename "$read").screen.txt"
FAKE_SCREEN

cat > "$fake_bin/bowtie2" <<'FAKE_BOWTIE2'
#!/usr/bin/env bash
exit 0
FAKE_BOWTIE2

chmod +x "$fake_bin"/*
export PATH="$fake_bin:$PATH"
export FAKE_FEATURECOUNTS_ARGS="$work/featurecounts_args.txt"
export FAKE_KALLISTO_ARGS="$work/kallisto_args.txt"
export FAKE_RSEM_ARGS="$work/rsem_args.txt"
export FASTQ_SCREEN_BIN="$fake_bin/fastq_screen"

assert_file() { [[ -s "$1" ]] || { echo "ASSERTION FAILED: expected non-empty file $1" >&2; exit 1; }; }
assert_absent() { [[ ! -e "$1" ]] || { echo "ASSERTION FAILED: expected absent path $1" >&2; exit 1; }; }

mkdir -p "$work/genome index" "$work/inputs" "$work/output"
printf '@r1\nACGT\n+\n!!!!\n' > "$work/inputs/lane1_R1.fastq"
printf '@r2\nTGCA\n+\n!!!!\n' > "$work/inputs/lane1_R2.fastq"
printf '@r3\nAAAA\n+\n!!!!\n' > "$work/inputs/lane2_R1.fastq"
printf '@r4\nTTTT\n+\n!!!!\n' > "$work/inputs/lane2_R2.fastq"
gzip -c "$work/inputs/lane1_R1.fastq" > "$work/inputs/lane1_R1.fastq.gz"
gzip -c "$work/inputs/lane1_R2.fastq" > "$work/inputs/lane1_R2.fastq.gz"
printf 'bam\n' > "$work/inputs/input.bam"
printf 'transcriptome\n' > "$work/inputs/transcriptome.bam"
printf 'gtf\n' > "$work/inputs/genes.gtf"
printf 'chr1\t1\t10\tgene1\n' > "$work/inputs/strand.bed"
printf 'kallisto-index\n' > "$work/inputs/transcripts.idx"
printf 'rsem-index\n' > "$work/inputs/rsem.grp"

# Explicit runner arguments must make wrappers independent of the submit directory.
(cd "$work" && bash "$repo_root/scripts_DoNotTouch/STAR/qsub_star_PE.sh" \
  "$work/output/paired sample" "$work/genome index" \
  "$work/inputs/lane1_R1.fastq.gz" "$work/inputs/lane1_R2.fastq.gz" project \
  "$repo_root/scripts_DoNotTouch/STAR/star_PE.sh")
assert_file "$work/output/paired sampleAligned.sortedByCoord.out.bam"
assert_file "$work/output/paired sampleAligned.sortedByCoord.out.bam.bai"

(cd "$work" && bash "$repo_root/scripts_DoNotTouch/STAR/qsub_star_SE.sh" \
  "$work/output/single sample" "$work/genome index" "$work/inputs/lane1_R1.fastq.gz" project \
  "$repo_root/scripts_DoNotTouch/STAR/star_SE.sh")
assert_file "$work/output/single sampleAligned.sortedByCoord.out.bam"

# A failed retry must not leave an old BAM that the App could mistake for success.
printf 'stale\n' > "$work/output/failedAligned.sortedByCoord.out.bam"
export FAKE_STAR_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/STAR/star_SE.sh" \
  "$work/output/failed" "$work/genome index" "$work/inputs/lane1_R1.fastq.gz"; then
  echo "ASSERTION FAILED: forced STAR failure returned success" >&2
  exit 1
fi
assert_absent "$work/output/failedAligned.sortedByCoord.out.bam"
unset FAKE_STAR_MODE

(cd "$work" && bash "$repo_root/scripts_DoNotTouch/featureCounts/qsub_featurecounts_PE.sh" \
  "$work/inputs/input.bam" "$work/inputs/genes.gtf" gene_name "$work/output/paired counts" \
  "$work/inputs/strand.bed" project "$repo_root/scripts_DoNotTouch/featureCounts/featurecounts_PE.sh")
assert_file "$work/output/paired counts_counts.txt"
grep -q -- '-p --countReadPairs' "$FAKE_FEATURECOUNTS_ARGS"

(cd "$work" && bash "$repo_root/scripts_DoNotTouch/featureCounts/qsub_featurecounts_SE.sh" \
  "$work/inputs/input.bam" "$work/inputs/genes.gtf" gene_id "$work/output/single counts" \
  "$work/inputs/strand.bed" project "$repo_root/scripts_DoNotTouch/featureCounts/featurecounts_SE.sh")
assert_file "$work/output/single counts_counts.txt"
if grep -q -- '--countReadPairs' "$FAKE_FEATURECOUNTS_ARGS"; then
  echo "ASSERTION FAILED: single-end featureCounts used paired-end flags" >&2
  exit 1
fi

export FAKE_INFER_MODE=malformed
if bash "$repo_root/scripts_DoNotTouch/featureCounts/featurecounts_SE.sh" \
  "$work/inputs/input.bam" "$work/inputs/genes.gtf" gene_id "$work/output/bad counts" "$work/inputs/strand.bed"; then
  echo "ASSERTION FAILED: malformed strandedness report returned success" >&2
  exit 1
fi
assert_absent "$work/output/bad counts_counts.txt"
unset FAKE_INFER_MODE

# Kallisto must write to the exact App target directory for paired and single-end jobs.
(cd "$work" && bash "$repo_root/scripts_DoNotTouch/Kallisto/qsub_kallisto_PE.sh" \
  "$work/output/kallisto paired" "$work/inputs/transcripts.idx" \
  "$work/inputs/lane1_R1.fastq.gz,$work/inputs/lane2_R1.fastq" \
  "$work/inputs/lane1_R2.fastq.gz,$work/inputs/lane2_R2.fastq" project \
  "$repo_root/scripts_DoNotTouch/Kallisto/kallisto_PE.sh")
assert_file "$work/output/kallisto paired/abundance.tsv"
if find "$work/output/kallisto paired" -maxdepth 1 -type d -name '.kallisto_tmp.*' | grep -q .; then
  echo "ASSERTION FAILED: paired Kallisto FIFO temporary directory was not cleaned" >&2
  exit 1
fi

(cd "$work" && bash "$repo_root/scripts_DoNotTouch/Kallisto/qsub_kallisto_SE.sh" \
  "$work/output/kallisto single" "$work/inputs/transcripts.idx" \
  "$work/inputs/lane1_R1.fastq.gz,$work/inputs/lane2_R1.fastq" project \
  "$repo_root/scripts_DoNotTouch/Kallisto/kallisto_SE.sh")
assert_file "$work/output/kallisto single/abundance.tsv"
assert_absent "$work/output/kallisto single_kallisto"
grep -q -- '--single' "$FAKE_KALLISTO_ARGS"

printf 'stale\n' > "$work/output/kallisto single/abundance.tsv"
export FAKE_KALLISTO_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/Kallisto/kallisto_SE.sh" \
  "$work/output/kallisto single" "$work/inputs/transcripts.idx" "$work/inputs/lane1_R1.fastq.gz"; then
  echo "ASSERTION FAILED: forced Kallisto failure returned success" >&2
  exit 1
fi
assert_absent "$work/output/kallisto single/abundance.tsv"
unset FAKE_KALLISTO_MODE

# RSEM wrappers must be location-independent, distinguish paired/single-end, and clear stale results.
(cd "$work" && bash "$repo_root/scripts_DoNotTouch/RSEM/qsub_RSEM_PE.sh" \
  "$work/inputs/input.bam" "$work/inputs/rsem" gene_id "$work/output/rsem paired" \
  "$work/inputs/strand.bed" "$work/inputs/transcriptome.bam" project \
  "$repo_root/scripts_DoNotTouch/RSEM/RSEM_PE.sh")
assert_file "$work/output/rsem paired.genes.results"
assert_file "$work/output/rsem paired.isoforms.results"
grep -q -- '--paired-end' "$FAKE_RSEM_ARGS"
assert_absent "$work/output/rsem paired.stat"

(cd "$work" && bash "$repo_root/scripts_DoNotTouch/RSEM/qsub_RSEM_SE.sh" \
  "$work/inputs/input.bam" "$work/inputs/rsem" gene_id "$work/output/rsem single" \
  "$work/inputs/strand.bed" "$work/inputs/transcriptome.bam" project \
  "$repo_root/scripts_DoNotTouch/RSEM/RSEM_SE.sh")
assert_file "$work/output/rsem single.genes.results"
if grep -q -- '--paired-end' "$FAKE_RSEM_ARGS"; then
  echo "ASSERTION FAILED: single-end RSEM used --paired-end" >&2
  exit 1
fi

printf 'stale\n' > "$work/output/rsem failed.genes.results"
export FAKE_RSEM_MODE=fail
if bash "$repo_root/scripts_DoNotTouch/RSEM/RSEM_PE.sh" \
  "$work/inputs/input.bam" "$work/inputs/rsem" gene_id "$work/output/rsem failed" \
  "$work/inputs/strand.bed" "$work/inputs/transcriptome.bam"; then
  echo "ASSERTION FAILED: forced RSEM failure returned success" >&2
  exit 1
fi
assert_absent "$work/output/rsem failed.genes.results"
unset FAKE_RSEM_MODE

# Exercise pooled-lane FIFO streaming with real gzip input and a fake trimmer.
bash "$repo_root/scripts_DoNotTouch/cutadapt_PE/cutadapt_PE.sh" 20 A1 A2 \
  "$work/output/trim_R1.fastq" "$work/output/trim_R2.fastq" \
  "$work/inputs/lane1_R1.fastq.gz,$work/inputs/lane2_R1.fastq" \
  "$work/inputs/lane1_R2.fastq.gz,$work/inputs/lane2_R2.fastq"
assert_file "$work/output/trim_R1.fastq"
[[ "$(grep -c '^@r' "$work/output/trim_R1.fastq")" -eq 2 ]] || { echo "ASSERTION FAILED: pooled R1 lanes were not concatenated" >&2; exit 1; }
if find "$work/output" -maxdepth 1 -name '.cutadapt_tmp.*' | grep -q .; then
  echo "ASSERTION FAILED: Cutadapt FIFO temporary directory was not cleaned" >&2
  exit 1
fi

bash "$repo_root/scripts_DoNotTouch/cutadapt_SE/qsub_cutadapt_SE.sh" 20 A1 unused \
  "$work/output/trim_single.fastq" "$work/output/unused.fastq" \
  "$work/inputs/lane1_R1.fastq.gz,$work/inputs/lane2_R1.fastq" \
  "$work/inputs/lane1_R1.fastq.gz,$work/inputs/lane2_R1.fastq" project \
  "$repo_root/scripts_DoNotTouch/cutadapt_SE/cutadapt_SE.sh"
assert_file "$work/output/trim_single.fastq"
[[ "$(grep -c '^@r' "$work/output/trim_single.fastq")" -eq 2 ]] || { echo "ASSERTION FAILED: pooled single-end lanes were not concatenated" >&2; exit 1; }

bash "$repo_root/scripts_DoNotTouch/FastQC/fastqc.sh" \
  "$work/inputs/lane1_R1.fastq.gz,$work/inputs/lane1_R2.fastq.gz" "$work/output/fastqc"
assert_file "$work/output/fastqc/lane1_R1_fastqc.html"
assert_file "$work/output/fastqc/lane1_R2.fastq.gz.screen.txt"

echo "Core RNA runner fake-tool smoke tests passed."
