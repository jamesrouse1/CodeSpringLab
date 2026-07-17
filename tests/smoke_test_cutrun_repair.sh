#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/codespring_cutrun_repair.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT
fake_bin="$work/bin"
mkdir -p "$fake_bin" "$work/picard"
printf 'jar\n' > "$work/picard/picard.jar"
export EBROOTPICARD="$work/picard"
module() { :; }
export -f module

cat > "$fake_bin/samtools" <<'FAKE_SAMTOOLS'
#!/usr/bin/env bash
set -euo pipefail
command="$1"; shift
case "$command" in
  quickcheck) exit 0 ;;
  index)
    [[ "${1:-}" == "-b" ]] && shift
    printf 'index\n' > "${1}.bai"
    ;;
  idxstats) printf 'chr1\t1000\t100\t0\n' ;;
  sort)
    out=""; input=""
    while (($#)); do
      case "$1" in
        -o) out="$2"; shift 2 ;;
        -T|-@) shift 2 ;;
        -n) shift ;;
        *) input="$1"; shift ;;
      esac
    done
    cp "$input" "$out"
    ;;
  view)
    if [[ " $* " == *" -c "* ]]; then printf '100\n'; else printf 'bam\n'; fi
    ;;
  *) echo "unsupported fake samtools command: $command" >&2; exit 2 ;;
esac
FAKE_SAMTOOLS

cat > "$fake_bin/bedtools" <<'FAKE_BEDTOOLS'
#!/usr/bin/env bash
set -euo pipefail
command="$1"; shift
case "$command" in
  bamtobed)
    if [[ " $* " == *" -bedpe "* ]]; then
      printf 'chr1\t10\t20\tchr1\t30\t40\tread1\t60\t+\t-\n'
    else
      printf 'chr1\t10\t40\tread1\t60\t+\n'
    fi
    ;;
  genomecov) printf 'chr1\t10\t40\t1\n' ;;
  *) echo "unsupported fake bedtools command: $command" >&2; exit 2 ;;
esac
FAKE_BEDTOOLS

cat > "$fake_bin/bamCoverage" <<'FAKE_BAMCOVERAGE'
#!/usr/bin/env bash
set -euo pipefail
out=""
while (($#)); do
  case "$1" in
    --outFileName) out="$2"; shift 2 ;;
    *) shift ;;
  esac
done
[[ -n "$out" ]] || exit 2
printf 'bigwig\n' > "$out"
FAKE_BAMCOVERAGE
chmod +x "$fake_bin"/*
export PATH="$fake_bin:$PATH"

sample_dir="$work/sample with spaces"
mkdir -p "$sample_dir"
prefix="$sample_dir/sample1"
printf 'bam\n' > "${prefix}Aligned.sortedByCoord.out.bam"
printf 'bam\n' > "${prefix}Aligned.sortedByCoord_removeDup.out.bam"
printf 'chr1\t1000\n' > "$work/chrom.sizes"

bash "$repo_root/scripts_DoNotTouch/bowtie2/bowtie2_cutrun_PE.sh" \
  "$prefix" unused-index unused-R1 unused-R2 "$work/chrom.sizes" smoke \
  30 1000 keepdup y CPM none spikein 1000 repair

for output in \
  "${prefix}_fragments.bed" \
  "${prefix}_fragments.raw.bedgraph" \
  "${prefix}_fragments.CPM.bedgraph" \
  "${prefix}_fragments.CPM.bw" \
  "${prefix}_alignment_summary.txt" \
  "${prefix}_postprocess_summary.txt"; do
  [[ -s "$output" ]] || { echo "ASSERTION FAILED: missing CUT&RUN repair output $output" >&2; exit 1; }
done
grep -q $'^duplicate_fraction\tNA$' "${prefix}_alignment_summary.txt"
grep -q $'^status\tpost_alignment_repaired$' "${prefix}_postprocess_summary.txt"
if find "$sample_dir" -maxdepth 1 -type d -name 'tmp_cutrun_*' | grep -q .; then
  echo "ASSERTION FAILED: CUT&RUN repair temporary directory was not cleaned" >&2
  exit 1
fi

missing_prefix="$work/missing/sample2"
if bash "$repo_root/scripts_DoNotTouch/bowtie2/bowtie2_cutrun_PE.sh" \
  "$missing_prefix" unused-index unused-R1 unused-R2 "$work/chrom.sizes" smoke \
  30 1000 keepdup y CPM none spikein 1000 repair > "$work/missing.out" 2>&1; then
  echo "ASSERTION FAILED: CUT&RUN repair accepted missing BAMs" >&2
  exit 1
fi
grep -q 'repair requires an existing non-empty BAM' "$work/missing.out"

echo "CUT&RUN post-alignment repair smoke test passed."
