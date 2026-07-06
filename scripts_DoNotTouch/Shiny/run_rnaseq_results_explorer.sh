#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="${1:-$APP_DIR/shiny_results_config.R}"
HOST_OVERRIDE="${2:-}"
PORT_OVERRIDE="${3:-}"
CONFIG_PATH="$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$CONFIG_PATH")"

module load EBModules
module load R/4.2.0-foss-2021b

export RNASEQ_SHINY_CONFIG="$CONFIG_PATH"
export RNASEQ_SHINY_HOST_OVERRIDE="$HOST_OVERRIDE"
export RNASEQ_SHINY_PORT_OVERRIDE="$PORT_OVERRIDE"

Rscript -e 'cfg <- new.env(parent = baseenv()); cfg_path <- Sys.getenv("RNASEQ_SHINY_CONFIG"); if (file.exists(cfg_path)) sys.source(cfg_path, envir = cfg); host <- Sys.getenv("RNASEQ_SHINY_HOST_OVERRIDE", unset = ""); port <- Sys.getenv("RNASEQ_SHINY_PORT_OVERRIDE", unset = ""); if (!nzchar(host)) host <- if (exists("host", envir = cfg, inherits = FALSE)) get("host", envir = cfg) else "0.0.0.0"; if (!nzchar(port)) port <- if (exists("port", envir = cfg, inherits = FALSE)) as.character(get("port", envir = cfg)) else "3838"; shiny::runApp(appDir = normalizePath("'"$APP_DIR"'/app_server.R", winslash = "/", mustWork = FALSE), host = host, port = as.integer(port), launch.browser = FALSE)'
