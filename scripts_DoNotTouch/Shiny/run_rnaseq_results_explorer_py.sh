#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="${1:-$APP_DIR/shiny_results_config.R}"
HOST_OVERRIDE="${2:-}"
PORT_OVERRIDE="${3:-}"
PYTHON_BIN="${4:-python}"

export RNASEQ_SHINY_CONFIG="$CONFIG_PATH"
export RNASEQ_SHINY_HOST_OVERRIDE="$HOST_OVERRIDE"
export RNASEQ_SHINY_PORT_OVERRIDE="$PORT_OVERRIDE"

exec "$PYTHON_BIN" "$APP_DIR/app_server_py.py"
