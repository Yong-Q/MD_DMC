#!/bin/bash
set -euo pipefail
source /work/home/acumhdgcrs/dmc_analysis/bin/activate
ROOT="${DMC_ROOT:-/work/home/acumhdgcrs/qy/1/MD_DMC}"
cd "$ROOT"
i="$1"
[ -e "analysis_v4/$i" ] && { echo "SKIP $i"; exit 0; }
python scripts/transport_v4.py --case "runs/$i" --output-dir "analysis_v4/$i" \
  --lag-frac 0.8 --max-lag-ps 12000 --width-min 500 --width-max 6000 --width-step 250 --start-step 100 \
  > "slurm/logs/v4_$i.log" 2>&1 && echo "DONE $i" || echo "FAIL $i"
