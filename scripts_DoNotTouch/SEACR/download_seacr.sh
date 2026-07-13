#!/usr/bin/env bash
set -euo pipefail

out_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$out_dir"

curl -L -o SEACR_1.3.sh https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.sh
chmod +x SEACR_1.3.sh
echo "Installed SEACR to: ${out_dir}/SEACR_1.3.sh"
