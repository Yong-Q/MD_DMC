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

## 推荐：整条轨迹的双对数区间筛选（v3）

对整条生产轨迹计算 Li、FSI 质心和集体电荷位移的 MSD，然后在**滞后时间**上扫描起点为 50 ps 步长、长度为 200、400、750、1000 ps 的窗口，滞后上限为 3000 ps 和轨迹长度 20% 中的较小者。滞后时间是两帧的时间差；整条轨迹中满足该时间差的时间起点都参与平均，并非截取模拟的前半段或后半段。对于每种量分别寻找双对数斜率 0.9–1.1、线性 R²≥0.98、正斜率和拟合区间前后两半斜率比为 0.5–2 的候选，从中选择斜率最接近 1 的窗口。这里的前后两半指拟合**滞后区间**，不是拆分模拟轨迹。找不到区间时标记 `no_candidate`。

```bash
python scripts/transport_v3.py --case runs/000 --output-dir analysis_v3/000
sbatch --array=0-3%4 slurm/analyze_array_v3.slurm
# 可选 Kinisi Li 交叉检查（先安装 requirements-kinisi.txt）
DMC_KINISI=1 sbatch --export=ALL --array=0-3%4 slurm/analyze_array_v3.slurm
python scripts/summarize_transport_v3.py --manifest manifest.csv --results-root analysis_v3 --output analysis_v3/summary.csv
```

输出双对数及拟合区间图 `diagnostics.png` / `.pdf`、曲线和 JSON；目录已存在时不覆盖。集体电导率以整个模拟盒体积和模型部分电荷计算；NE 用 Li 和 FSI 质心扩散估算。汇总只对 4 个不同种子重复都通过筛查的量报告均值和样本标准差；缺失或未通过时留空。Kinisi 是可选数值敏感性检查，不会替未进入扩散区的轨迹给出可信值。

扫描多个窗口并选择最接近 1 的结果会产生选择偏差，必须同时检查曲线、邻近区间和独立重复；这些候选窗口和阈值是探索性筛选，不保证真实扩散区、统计收敛或实验准确度。R² 和双对数斜率接近 1 不足以排除非平稳性、时间起点相关性或选择偏差；还须检查结构、能量、温度、周期解包与各重复的差异。需要完整 15 个 DCD/restart；数值检查和前后温度均值检查不能证明体系结构平衡。当前默认滞后上限为完整轨迹长度的 20%。

## 历史 v2（已被 v3 取代）：前后时间段对照

建议使用 `scripts/transport_v2.py`。默认 `first-half` 模式只用生产轨迹的前 7.5 ns 计算数值，以便在后半程可能发生变化时避免把两种状态混合。**这不是自动证明前半程已经平衡或处于长时间扩散区**，所得值一律标为 `first_half_provisional`。后半程不会进入该数值，但仍计算并显示为对照；前后差异可能来自结构/体系状态变化、采样噪声或处理误差，需结合能量、温度和结构诊断原因。

v2 从 Li 位置、FSI **质量加权质心**和全部移动原子的部分电荷分别计算 MSD 和集体电荷位移。它在前半程按预定顺序尝试滞后窗口 750–1500、500–1000、250–750、100–500 ps，选第一个满足双对数斜率 0.9–1.1、线性 R²≥0.98、斜率为正的窗口。没有候选窗口就标记 `no_candidate`，不输出该量的数值。然后在后半程同一区间做对照；它也满足上述条件且前后斜率比为 0.5–2 时记录 `candidate`，否则记录 `holdout_rejected`。不根据后半程重新挑窗口。时间段虽不重叠，仍可能受长相关时间影响，并不能保证统计独立。

这些阈值和时间窗口是**预设的工程筛查规则，而不是可靠度证明或经校准的统计检验**。默认要求完整 15 个 DCD 和 15 个 restart，全部帧数须与元数据一致；不会自动剔除生产初期。程序对生产阶段温度日志检查覆盖率，并比较前后半程平均温度与目标值（容差 10 K，仅作粗筛）。还必须额外检查前半程的能量、MOF/客体结构和离子分布是否平稳，以及周期边界解包是否正确。集体电导率使用整个模拟盒体积和模型中的部分电荷；NE 电导率用 Li、FSI 质心的自扩散和模型电荷，忽略离子相关性。固定 MOF、50 对离子、模拟时长、缩放电荷和有限尺寸都可能影响结果。

```bash
# 默认只对前半程估计；输出目录必须不存在。
python scripts/transport_v2.py --case runs/000 --output-dir analysis_v2/000
# 批量分析：在仓库根目录提交，先调整分区、内存、并发和时间。
sbatch --array=0-3%4 slurm/analyze_array_v2.slurm
# 安装可选依赖 `pip install -r requirements-kinisi.txt` 后，对 Li 前半程候选区间加做可选敏感性诊断：
DMC_KINISI=1 sbatch --export=ALL --array=0-3%4 slurm/analyze_array_v2.slurm
# 4 个不同重复均有前半程候选值时，计算其均值和样本标准差。
python scripts/summarize_transport_v2.py --manifest manifest.csv --results-root analysis_v2 --output analysis_v2/summary.csv
```

`diagnostics.png` / `.pdf` 展示前后半程双对数与线性曲线、选定窗口和前半程拟合；`results.json` 记录每个预设窗口的筛选指标、后半程状态以及结果；`curves.csv` 保存曲线。程序错误保留 `failure.json`。汇总要求对应 MOF 和比例的 4 个不同编号、不同速度种子重复全部有候选值，且目标温度、盒体积一致，所用时段的平均温度通过粗筛；否则均值留空。汇总保留各重复后半程对照状态，因此即使 4 个初步值均有，也要检查是否发生 `holdout_rejected`；样本标准差不能覆盖模型或力场偏差。

另一种 `--mode verified` 模式要求后半程也通过同窗口检查，才拟合整个 15 ns；它只适用于确认前后同一状态的情况，不应用于后半段异常的体系。例如先 `export DMC_ANALYSIS_MODE=verified DMC_ANALYSIS_ROOT="$PWD/analysis_verified"`，再通过 `sbatch --export=ALL ...` 提交；汇总时对该独立目录使用 `--mode verified`。

Kinisi 为可选 Li 交叉检查，使用选定窗口、16/31 个滞后点、两个随机种子和特征值裁剪；数值稳定不代表 MCMC 收敛、后验覆盖或扩散区已经获得证明。测试中其他矩阵处理方式和随机种子曾使结果大幅波动。混合正负部分电荷体系的集体电导率没有通过 Kinisi 独立验证。之前 v1 自动扫描的 10 条 Li 候选中 9 条斜率集中在 0.800–0.805，存在贴边选择风险，不应直接当作可靠结果。

## 旧版 v1：仅供历史对照

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
