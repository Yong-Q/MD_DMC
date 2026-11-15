#!/bin/bash
#SBATCH --job-name=dmc_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-167%4
set -euo pipefail
ROOT="${DMC_ROOT:-${SLURM_SUBMIT_DIR:?Submit from the repository root}}"
cd "$ROOT"
printf -v IDX '%03d' "$SLURM_ARRAY_TASK_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec python scripts/robust_transport.py --case "$ROOT/runs/$IDX" --output-dir "$ROOT/analysis_results/$IDX"
