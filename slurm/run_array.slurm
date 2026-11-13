#!/bin/bash
#SBATCH --job-name=dmcmd
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --array=0-167%8
set -euo pipefail
ROOT="${DMC_ROOT:-${SLURM_SUBMIT_DIR:?Submit from the repository root}}"
cd "$ROOT"
RANKS="${SLURM_NTASKS:-48}"
exec python scripts/run_case.py --index "$SLURM_ARRAY_TASK_ID" --ranks "$RANKS" --launcher srun ${DMC_RESUME:+--resume}
