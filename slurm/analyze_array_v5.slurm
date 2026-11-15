#!/bin/bash
#SBATCH --job-name=dmc_post_v5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-167%4
set -euo pipefail
ROOT="${DMC_ROOT:-${SLURM_SUBMIT_DIR:?Submit from repository root}}"
cd "$ROOT"
printf -v IDX '%03d' "$SLURM_ARRAY_TASK_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ARGS=()
if [[ "${DMC_KINISI:-0}" == 1 ]]; then ARGS+=(--kinisi); fi
exec python scripts/transport_v5.py --case "$ROOT/runs/$IDX" --output-dir "$ROOT/analysis_v5/$IDX" "${ARGS[@]}"
