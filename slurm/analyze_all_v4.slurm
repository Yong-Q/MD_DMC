#!/bin/bash
#SBATCH --job-name=dmc_post_v4_all
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=13600M
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/%j.out
set -euo pipefail
ROOT="${DMC_ROOT:-${SLURM_SUBMIT_DIR:?Submit from repository root}}"
cd "$ROOT"
LIST="slurm/completed_v4.txt"
: > "$LIST"
for d in runs/*/; do
  c=$(grep -h "^COMPLETED_CHUNK" "$d"/logs/*.log 2>/dev/null | awk '{print $2}' | sort -n | wc -l)
  [ "$c" -eq 15 ] && basename "$d" >> "$LIST"
done
echo "completed cases: $(wc -l < "$LIST")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p analysis_v4
cat "$LIST" | xargs -P 4 -I{} "$ROOT/slurm/run_v4_one.sh" {}
echo "ALL_DONE"
