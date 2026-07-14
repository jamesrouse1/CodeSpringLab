#!/usr/bin/env bash
set -euo pipefail

out_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$out_dir"

echo "Refreshing bundled SEACR 1.3 shell and R scripts from FredHutch/SEACR..."
curl -L -o SEACR_1.3.sh https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.sh
curl -L -o SEACR_1.3.R https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.R
chmod +x SEACR_1.3.sh
curl -L -o LICENSE_SEACR https://raw.githubusercontent.com/FredHutch/SEACR/master/LICENSE
echo "SEACR is available at: ${out_dir}/SEACR_1.3.sh and ${out_dir}/SEACR_1.3.R"
