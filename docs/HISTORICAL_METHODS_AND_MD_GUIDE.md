# MOF + LiFSI + DMC: 168 systems

## 统一后处理：唯一推荐入口

使用 [transport_framework.py](scripts/transport_framework.py) 和锁定配置 [config/transport.json](config/transport.json)。完整方法、近期论文/GitHub依据、质量标记和续跑说明见 [统一框架文档](docs/TRANSPORT_FRAMEWORK.md)。历史 v1–v5 保留作追溯，不能混合汇总。

```bash
python -m pip install -r requirements-transport.txt
sbatch slurm/transport_pipeline.slurm
```

它只处理已有轨迹：坐标缓存 → STACIE 电导率和频谱诊断 → 同一体系重复共享 MSD 平台区间 → 独立重复误差及 EH/NE 对照 → 图表。没有足够证据的结果会标为 review_required/no_plateau，描述均值与通过筛查的值分别输出。


Six MOFs × seven DMC/LiFSI ratios (2, 4, 6, 8, 10, 12, 14) × four independently prepared replicas. Each system has 50 LiFSI pairs; the number of DMC molecules is the ratio multiplied by 50. All 168 initial systems passed a 100-step LAMMPS check on the original cluster. `manifest.csv` maps array indices to compositions, replicas and SHA-256 of `system.data`.

These are **initial configurations**, including cases that had already started or finished elsewhere. The repository does not contain trajectories or restart files; derived local transport results are provided under `results/`. Running a case here starts it from its initial structure; it cannot continue the prior cluster's partial trajectory. `lists/cancelled_before_start_indices.txt` is the list of 115 jobs canceled while pending in the original array, not a live status list. Avoid duplicate runs by checking your own records.

The six materials were selected using historical Li self-diffusion classifications; these labels do not establish conductivity enhancement or suppression. All calculations use the original fixed framework and force-field inputs. Production is 15 chunks of 1 ns at 298 K, after the original relaxation/annealing protocol. Each replica has its own initial structure and velocity seed. Cell size is fixed across concentrations within each MOF.

## 本地后处理结果表

[168 个体系的已有扩散结果](results/existing_diffusion_20261115.csv) 和 [42 个已分析体系的已有电导率窗口表](results/existing_conductivity_windows_20261115.csv) 的字段、筛选条件与浓度定义见 [结果说明](results/README.md)。未计算或筛选不通过的数值留空，不表示零。

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

## 历史分析入口（用于复现旧结果）

Install the packages in `requirements.txt` into a suitable Python environment. After all 15 chunks finish:

```
python scripts/analyze_diffusion.py --case runs/000
PYTHONPATH="$PWD/src" python scripts/analyze_collective_conductivity.py --case runs/000
```

The first script computes Li and FSI MSD and log–log slopes, with linear MSD fitting rather than Kinisi to reduce memory use. Check the slope and fit window before interpreting diffusion. The second computes collective charge-displacement conductivity (Einstein–Helfand) and several time windows; it requires all chunks. Results belong to each run's `analysis/` directory. Use the four independent replicas to estimate variability; this package does not claim statistically converged transport coefficients. Analysis can still require substantial RAM and time for large trajectories. Run it as a separate appropriately sized cluster job.

The structures are ready to run; original source paths retained in metadata are provenance only and are not consulted for the MD run. Git does not include raw trajectories, restart files, existing outputs, or machine-specific LAMMPS binaries.

## 历史 v3：整条轨迹区间筛选

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

## 历史 v5：保守区间筛选

在检查 v4 和其结果后，v5 用于统一本地及另一台机器的体系：从完整轨迹计算 Li 自扩散、FSI 质心自扩散和 Einstein–Helfand 集体电导率；NE 电导率作为辅助。v4 在大量候选窗口中直接挑 α 最接近 1 的窗口，存在选区偏差；其允许的长滞后时间接近 15 ns 轨迹长度，采样更少。v4 历史结果仅供对照，不应与 v5 混合汇总。`v4_20261115_conductivity_by_state.csv` 还混用了旧固定窗口和 v4 结果。

v5 使用整条生产轨迹的全部可用时间起点计算 MSD，不按轨迹前后半段划分；拟合对象是滞后时间。预先按 2000、1000、500 ps 宽度和起点 100 ps 步长检查；最大滞后时间为完整轨迹长度的 50%。要求 α∈[0.9,1.1]、拟合 R²≥0.98、斜率为正、窗口两半的斜率比为 2/3–1.5，且窗口整体移动 ±100 ps 后斜率变化不超过 35%。按预定顺序选第一个通过的区间，避免以 α 离 1 最近作为优化目标。找不到时为 `no_candidate`，不输出该量的值。

这些均为启发式筛查参数，**不是严格的可靠性保证或实验验证**：时间起点高度相关，选择窗口仍有偏差，非扩散过程的交叉区也可能满足这些条件。还需要检查 MSD 图、拟合窗口敏感性、模拟平衡、结构、电荷约定、独立重复及尺寸效应。缩放电荷的集体电导率使用整个盒体积，输出 xyz 三方向的标量平均而不是电导率张量；FSI 是质量加权质心。Kinisi 可选，仅为 Li 的条件性数值交叉检查，不能替代采样与物理验证。

```bash
# 已有 MD 放在 runs/NNN 下；每个输出目录必须不存在。
git pull
python scripts/transport_v5.py --case runs/000 --output-dir analysis_v5/000
# 所有完整、尚未处理的体系：按目标集群选择分区和合适的并发。
ids=$(python scripts/list_completed_for_v5.py)
[ -n "$ids" ] && sbatch --array="${ids}%4" slurm/analyze_array_v5.slurm
# 可选安装 requirements-kinisi.txt，然后提交前设置 DMC_KINISI=1。
python scripts/export_transport_v5.py --manifest manifest.csv --results-root analysis_v5 --output analysis_v5/per_case.csv
python scripts/summarize_transport_v5.py --manifest manifest.csv --results-root analysis_v5 --output analysis_v5/summary.csv
```

程序输出每条轨迹的双对数与拟合区域图、曲线和 JSON。逐体系表额外标记 α 位于 0.9–0.92 或 1.08–1.1 的贴边结果，这只是根据已有数据添加的提醒、不是预校准的统计检验；汇总中含贴边重复的均值会标记 `complete_alpha_boundary_warning`。聚合表只有 4 个不同速度种子重复都通过筛选、温度粗筛、盒体积与目标温度一致时才给均值和样本标准差；缺失不能当零。`list_completed_for_v5.py` 跳过已有输出目录（包括失败目录），需检查并保留失败证据后再安排重跑。读取轨迹仍要求完整 15 个 DCD/restart、帧数正确；这些检查不保证周期解包或结构正确。

## 历史 v4（仅供对照，存在窗口选择与长滞后风险）

`scripts/transport_v4.py` 不再使用固定的 200/400/750/1000 ps 窗口或固定的 3 ns 滞后上限，而是对整条生产轨迹扫描：

- 滞后上限：`min(0.8 × 轨迹长度, 12 ns)`；
- 窗口宽度：500–6000 ps，步长 250 ps；
- 窗口起点：50 ps 起，步长 100 ps；
- 候选条件：双对数斜率 α∈[0.9,1.1]、斜率>0、线性 R²≥0.98、拟合区间前后两半斜率比 0.5–2；
- 选择 α 最接近 1 的窗口；找不到时输出 `no_candidate`，不强行给值。

每个 case 输出 Li 自扩散、9 原子质量加权 FSI 质心扩散、集体电荷位移电导率各一个值，以及 `curves.csv`、`results.json`、`diagnostics.png/pdf`。单作业批量运行示例：

```bash
sbatch slurm/analyze_all_v4.slurm
```

2026-11-15 的 80 个已完成体系结果见 `results/v4_20261115_*.csv`。与 `results/existing_*` 中的旧固定窗口结果相比，v4 的 Li/FSI 候选率约 96–98%，集体电荷约 81%；两种方法的电导率在比例 4 上相差约 10–25%，比较时须标注来源。
