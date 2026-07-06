#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_shiny
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="${1:-$APP_DIR/shiny_results_config.R}"
HOST_OVERRIDE="${2:-}"
PORT_OVERRIDE="${3:-}"

bash "$APP_DIR/run_rnaseq_results_explorer.sh" "$CONFIG_PATH" "$HOST_OVERRIDE" "$PORT_OVERRIDE"

