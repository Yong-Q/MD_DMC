#!/bin/bash
#SBATCH --job-name=dmc_post_v3
#SBATCH --partition=xahcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=13600M
#SBATCH --time=04:00:00
#SBATCH --array=0-167%4
#SBATCH --output=slurm/logs/%A_%a.out
set -euo pipefail
source /work/home/acumhdgcrs/dmc_analysis/bin/activate
ROOT="${DMC_ROOT:-${SLURM_SUBMIT_DIR:?Submit from repository root}}"
cd "$ROOT"
printf -v IDX '%03d' "$SLURM_ARRAY_TASK_ID"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
if [[ "${DMC_KINISI:-0}" == 1 ]]; then
  exec python scripts/transport_v3.py --case "$ROOT/runs/$IDX" --output-dir "$ROOT/analysis_v3/$IDX" --kinisi
else
  exec python scripts/transport_v3.py --case "$ROOT/runs/$IDX" --output-dir "$ROOT/analysis_v3/$IDX"
fi
