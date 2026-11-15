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

## 双对数区间筛选和电导率（探索版）

`scripts/robust_transport.py` 从一条完整的 15 ns 轨迹计算 Li 自扩散、FSI 氮原子锚点自扩散、以及全部移动原子部分电荷的集体电荷位移。分别为三条曲线寻找窗口，输出双对数图、线性图和区间；没有合格窗口就标记 `no_diffusive_window`，不输出该量的数值。它是筛选与诊断工具，不能凭双对数斜率自动证明扩散或电导率已收敛。

从仓库根目录，加载包含 `requirements.txt` 依赖的 Python 环境后：

```bash
python scripts/robust_transport.py --case runs/000 --output-dir analysis_results/000
# Slurm：自行调整分区、内存和并发。必须从仓库根目录提交。
sbatch --array=0-3%4 slurm/analyze_array.slurm
# 汇总不同独立重复；缺少或未通过筛选的重复会被明确标出。
python scripts/summarize_transport.py --manifest manifest.csv --results-root analysis_results --output analysis_results/summary.csv
```

输出目录必须不存在，以免覆盖已有结果。每个 case 包含 `diffusion_windows.png`、`diffusion_windows.pdf`（每一行左侧双对数图、右侧线性拟合图，阴影为所选区间）、`curves.csv` 和 `results.json`。`results.json` 包括区间、对数斜率、线性 R²、半区间斜率差、`early_window_warning`、自扩散系数（cm²/s）以及集体和 NE 电导率（mS/cm）。

**筛选规则**：默认滞后时间从 50 ps 起，最长为轨迹长度的 20%，且不超过 3000 ps；窗口长度至少 200 ps；双对数拟合斜率 0.8–1.2，线性拟合 R² ≥ 0.98，前后半段斜率为正且相对差异 ≤ 50%。符合条件时选最长窗口，再以更晚的开始时间打破平局。参数是为初步筛选设定的启发式，并非收敛或统计显著性的证明；起点早于 200 ps 会警告。长窗口和双对数接近 1 也可能是短暂表观行为或扫描选择偏差，需人工查看图、做窗口敏感性及轨迹长度检查。

集体结果用模拟轨迹中所有移动原子电荷和整个盒体积，包含中性 DMC 的部分电荷贡献；NE 用 Li 和 FSI 锚点的自扩散与模型中缩放电荷。FSI 锚点是长时间平移扩散的近似，不是分子质心；短时间区间会受分子转动影响。所有结果对应固定 MOF 框架和模型缩放电荷，不能直接等同实验电导率或形式电荷的结果。不同 lag 的时间起点高度重叠，拟合窗口不构成独立重复；只有不同初始构型/种子才能作为独立重复。`summary.csv` 给出通过筛选的重复数、均值和至少两个重复的样本标准差，但按通过筛选的重复求均值可能有选择偏差；缺失或失败的重复不可当零值，也不可默认为通过。15 ns、50 对离子、盒尺寸、平衡性、DCD 解包、有限尺寸和对电荷模型的敏感性仍需逐项检查。局部拟合误差或 R² 不代表整体物理不确定度。

`--cached-curves` 仅供用旧 `li_msd.csv`、`fsi_msd.csv` 和 `charge_moment_msd.csv` 进行开发对比：无法核验其来源和电荷约定，故不计算 NE，也不会计入汇总结果。正式分析请从原始轨迹运行。要求全部 15 个 DCD 和对应的完成 restart，且总帧数必须匹配元数据；程序不自动判断截断轨迹、周期解包正确性、平衡性或实际温度漂移。
