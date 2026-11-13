# MOF + LiFSI + DMC: 168 systems

Six MOFs × seven DMC/LiFSI ratios (2, 4, 6, 8, 10, 12, 14) × four independently prepared replicas. Each system has 50 LiFSI pairs; the number of DMC molecules is the ratio multiplied by 50. All 168 initial systems passed a 100-step LAMMPS check on the original cluster. `manifest.csv` maps array indices to compositions, replicas and SHA-256 of `system.data`.

These are **initial configurations**, including cases that had already started or finished elsewhere. The repository does not contain trajectories, restart files or prior results. Running a case here starts it from its initial structure; it cannot continue the prior cluster's partial trajectory. `lists/cancelled_before_start_indices.txt` is the list of 115 jobs canceled while pending in the original array, not a live status list. Avoid duplicate runs by checking your own records.

The six materials were selected using historical Li self-diffusion classifications; these labels do not establish conductivity enhancement or suppression. All calculations use the original fixed framework and force-field inputs. Production is 15 chunks of 1 ns at 298 K, after the original relaxation/annealing protocol. Each replica has its own initial structure and velocity seed. Cell size is fixed across concentrations within each MOF.

## 原服务器任务状态（2026-11-13 快照）

168 个初始体系全部已上传。原服务器 Slurm 数组 `121903` 当时的状态如下；这是一份静态快照，之后 6 个运行中的任务可能已结束，请避免重复提交。

| 状态 | 数量 | 处理建议 |
|---|---:|---|
| 排队时取消、未开始生产 MD | 115 | 可在其他机器从头运行 |
| 开始后取消 | 10 | 可在其他机器从头运行；原服务器保留了部分轨迹 |
| 仍在运行 | 6 | 先确认原服务器最终状态，避免重复计算 |
| Slurm 显示完成 | 36 | 优先检查原服务器已有结果 |
| 后处理内存不足 | 1 | 15 ns 轨迹已完成；优先修复后处理，不必默认重跑 MD |

在另一台机器运行全部 **125 个已取消体系**，从仓库根目录提交：

```bash
ids=$(paste -sd, lists/cancelled_all_125_indices.txt)
sbatch --array="${ids}%8" slurm/run_array.slurm
```

- [125 个已取消体系编号](lists/cancelled_all_125_indices.txt)：包含下面两类，可直接用于上面的命令。
- [115 个未开始体系编号](lists/cancelled_before_start_indices.txt)。
- [10 个开始后取消的体系编号](lists/cancelled_after_start_10_indices.txt)：43–52；原服务器已有部分轨迹，但本仓库没有轨迹和重启文件，所以新机器会从初始结构重新运行。
- [168 个体系的原服务器状态明细](lists/server_status_20261113.csv)：包含 MOF、浓度、重复编号、Slurm 状态和当时完成的 1 ns 分段数。`completed_1ns_chunks` 为快照，不是实时状态；Slurm COMPLETED 也不等同于科学结果已收敛。

## Run on another Slurm cluster

1. Clone the repository to a path **without spaces**. Install a LAMMPS build with MPI, MOLECULE, KSPACE and DCD support, and load the appropriate MPI/LAMMPS environment. Use a compatible MPI/Slurm installation. Set `LAMMPS_EXE` to its executable if it is not `lmp`. Adjust the partition, time, memory, rank count and `MPI_TYPE` to your cluster. The provided script defaults to 48 MPI ranks, PMI2 and 72 hours.
2. From the repository root, for a single case: `sbatch --array=0 slurm/run_array.slurm`. For all 125 canceled jobs (including 10 canceled after starting): `ids=$(paste -sd, lists/cancelled_all_125_indices.txt); sbatch --array="${ids}%8" slurm/run_array.slurm`. The `%8` limits concurrent array elements; adjust to your allocation. To select 40 cores, pass `--ntasks=40` to sbatch. Submission starts each case from its original built configuration.
3. Results are written to `runs/NNN/`. The program refuses to overwrite a preexisting run. For a run interrupted after at least one completed chunk, resubmit with `DMC_RESUME=1 sbatch --export=ALL --array=INDEX slurm/run_array.slurm`. Incomplete files and autosaves are preserved under `runs/NNN/interrupted/`; continuation starts from the last fully completed chunk, not an incomplete autosave. If the initial chunk never completed, preserve or rename the whole run directory before starting again.

For a non-Slurm host: `python scripts/run_case.py --index 0 --launcher mpirun --ranks 40 --lammps /path/to/lmp` (or `--launcher none` for a serial executable). Existing runs require `--resume`. Do not run the same index concurrently.

## Analysis

Install the packages in `requirements.txt` into a suitable Python environment. After all 15 chunks finish:

```
python scripts/analyze_diffusion.py --case runs/000
PYTHONPATH="$PWD/src" python scripts/analyze_collective_conductivity.py --case runs/000
```

The first script computes Li and FSI MSD and log–log slopes, with linear MSD fitting rather than Kinisi to reduce memory use. Check the slope and fit window before interpreting diffusion. The second computes collective charge-displacement conductivity (Einstein–Helfand) and several time windows; it requires all chunks. Results belong to each run's `analysis/` directory. Use the four independent replicas to estimate variability; this package does not claim statistically converged transport coefficients. Analysis can still require substantial RAM and time for large trajectories. Run it as a separate appropriately sized cluster job.

The structures are ready to run; original source paths retained in metadata are provenance only and are not consulted for the MD run. Git does not include raw trajectories, restart files, existing outputs, or machine-specific LAMMPS binaries.
