#!/bin/bash
#SBATCH --job-name=transport_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=08:00:00
set -euo pipefail
ROOT="${DMC_ROOT:-${SLURM_SUBMIT_DIR:?Submit from repository root}}"
cd "$ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ARGS=()
if [[ "${DMC_RESUME:-0}" == 1 ]]; then ARGS+=(--resume); fi
if [[ -n "${DMC_CACHE_ROOT:-}" ]]; then ARGS+=(--cache-root "$DMC_CACHE_ROOT"); fi
exec python scripts/transport_framework.py batch --config config/transport.json --manifest manifest.csv --runs-root "${DMC_RUNS_ROOT:-$ROOT/runs}" --output "${DMC_ANALYSIS_OUTPUT:-$ROOT/transport_analysis}" --workers "${SLURM_CPUS_PER_TASK:-2}" "${ARGS[@]}"
