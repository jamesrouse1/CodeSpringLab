#!/usr/bin/env bash
set -Eeuo pipefail

test_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bash "$test_dir/smoke_test_static_sources.sh"
bash "$test_dir/smoke_test_scanpy_container_runner.sh"
bash "$test_dir/smoke_test_cellranger_runner.sh"
Rscript "$test_dir/smoke_test_scrna_processed_object.R" "$(cd "$test_dir/.." && pwd)"
bash "$test_dir/smoke_test_core_runners.sh"
bash "$test_dir/smoke_test_cutrun_repair.sh"
bash "$test_dir/smoke_test_peak_runners.sh"
bash "$test_dir/smoke_test_seacr_actual.sh"
echo "All CodeSpringLab tests passed."
