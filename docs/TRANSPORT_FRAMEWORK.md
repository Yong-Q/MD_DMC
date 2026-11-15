# 统一输运后处理框架

唯一推荐入口为 `scripts/transport_framework.py`，统一配置为 `config/transport.json`。它处理已经完成的轨迹，不提交或增加 MD。历史 v1–v5 保留用于追溯，不应与本框架结果混表。

## 采用哪些近期方法

| 工作 | 核对的实现/说明 | 本框架采用的部分 |
|---|---|---|
| STACIE，JCIM 2025，1.3.0 软件/2026 文档 | https://github.com/molmod/stacie ; https://doi.org/10.1021/acs.jcim.5c01475 | 对电荷电流的低频谱估计电导率及条件误差，保留模型与有效频率诊断 |
| MDTransport，JCIM 2026 | https://github.com/ShahResearchGroup/MDTransport ; https://doi.org/10.1021/acs.jcim.6c00558 | 批量入口、分子质心、EH 与 NE 共同报告；没有复制或直接调用该项目的选窗实现 |
| GEMDAT，npj Computational Materials 2026 | https://github.com/GEMDAT-repos/GEMDAT ; https://doi.org/10.1038/s41524-026-02133-7 | 作为固体扩散分析工具对照；未把其 NE 估计当作混合电解质总电导率 |
| Kinisi，JCTC 2025 | https://github.com/kinisi-dev/kinisi ; https://doi.org/10.1021/acs.jctc.4c01249 | 比较了其条件协方差与贝叶斯误差；因本地测试存在数值敏感性，未强制作为默认估计器 |
| T-MSD，2025 预印本 | https://arxiv.org/abs/2503.04281 | 关注时间平均和重采样；未核实可直接采用的官方实现，因此未声称实现了该论文算法 |

STACIE 文档明确支持从已记录坐标构造**区间平均电流**：https://molmod.github.io/stacie/properties/electrical_conductivity.html 。坐标差除以帧间隔是该时段的平均电流，不是假装恢复 fs 尺度的瞬时速度。不得把每 1 ps 稀疏记录的瞬时速度等同此平均电流。它也不能恢复轨迹中丢失的跨周期信息。

## 统一计算与输出

1. **只读与缓存**：复用已核验的三斜盒 NoJump 处理；提取全部帧的 Li、九原子 FSI 质量加权质心、全部移动原子的电荷矩。核对 15 段、帧数、离子数、模型电荷与固定晶胞。压缩缓存用于后续重算，避免反复读取大轨迹；缓存以拓扑/元数据 SHA-256 和轨迹文件大小/mtime 做失效检查（没有对所有 DCD 内容逐字节校验）。能量、结构平衡与跨周期解包仍须额外审查。
2. **主要集体电导率**：对电荷矩相邻帧之差构造平均电流，转换为 SI；调用锁定版本 STACIE、固定 ExpPolyModel([0,1,2])、排除零频、自动边缘化频率截止。按 x/y/z 分别拟合，避免把 MOF 三方向强行当成同一分布。标量为三方向平均；同时保留每个方向的值。未知方向间相关性时，条件 SE 用三个方向 SE 之和/3 作为保守上界，不按平方和随意缩小。
3. **谱质量诊断**：保留 STACIE 的有效频率数和两个 Z 分数。有效频率数小于 20×模型参数数、|Z|>2，或相对条件 SE 超过 50% 时标记 needs_review。前两项来自软件的诊断建议；50% 是本框架的工程警示，不是普适验收标准。低频模型也可能不适合跳跃稀少、非平稳或受限运动。
4. **扩散与 EH 交叉检查**：不再逐条曲线最优化 α。对同一 MOF、比例的四次重复的平均 MSD，寻找局部双对数斜率持续在 [0.9,1.1] 的平台；要求跨度至少 500 ps、滞后时间比至少 2、最大滞后为轨迹长度的 1/4。只选择最长连续平台，全部重复使用同一个区间。每个重复都保留，逐条检查是否明显不一致。此平台检测是启发式，不是扩散证明。没有平台或重复不足时，用预先固定的 10%–50% 滞后区间给出 diagnostic 值，并明确 no_plateau / incomplete，不能当作已确认的长时间系数。
5. **独立重复误差**：四次重复是统计单位；不把 50 个 Li、XYZ 方向、重叠时间起点或多个拟合窗口当独立重复。均值、样本 SD、Student-t 95% 区间包含全部可用重复；不丢掉数值不好看的重复再求均值。缺少或失败的重复记录清楚，mean 为描述值，screened_mean 仅在完整、质量检查满足且区间下界为正时填写。即使 screened_mean 存在，也没有完成物理确认。
6. **一致性**：在同一重复上比较 EH 与 STACIE，用配对差值的重复间区间检查一致性；不把两者当作统计独立。NE 使用同一批 Li/FSI 扩散与模型电荷，作为忽略离子相关性的对照。

配置固定并有 SHA-256，所有结果携带配置。任意改参数、分析窗口或电荷约定，都应使用新输出目录，不能与其他机器的旧表混合。除前三方向的谱诊断外，有限轨迹、有限尺寸、力场和电荷缩放偏差没有被误差条覆盖。框架默认使用模型部分电荷（包括 DMC 的内部部分电荷）和整个模拟盒体积；不是用形式电荷替换模型电荷，不做未经验证的孔体积修正或 Yeh-Hummer 修正。

## 批量使用

在现有 Python 3.11+ 环境安装一次依赖：

```bash
python -m pip install -r requirements-transport.txt
# 在仓库根目录提交；按目标集群调整 partition、内存和时间。
sbatch slurm/transport_pipeline.slurm
```

默认读 `runs/NNN`，在 `transport_analysis` 中生成缓存、每个 case 结果、谱图和 MSD 图、体系级共同拟合区间图及 CSV。后续需要续跑时：

```bash
python scripts/transport_framework.py batch --runs-root runs --output transport_analysis --workers 2 --resume
```

也可直接 `DMC_RESUME=1 sbatch --export=ALL slurm/transport_pipeline.slurm` 续跑。续跑须在 Slurm 资源内执行；不能同时启动两个针对同一输出目录的批次。成功输出按配置和输入指纹校验后跳过；失败目录先移动到 failed_attempts 保留，再重试。新汇总用新目录保存，latest_summary.txt 指向最近一次汇总。配置变化则拒绝复用该输出根目录。工作进程数不得超过 Slurm 分配，两个进程默认申请 24 GiB；大体系应控制并发与内存。

已有缓存可以独立运行：

```bash
python scripts/transport_framework.py analyze --cache transport_analysis/cache/002.npz --output separate_analysis/002
python scripts/transport_framework.py summarize --manifest manifest.csv --runs-root transport_analysis --output separate_summary
```

`per_case.csv` 的 *_diagnostic_value 可能为有限采样的负值，不能解释为物理负扩散/负电导率。by_state.csv 的 mean 是全部重复的描述均值；**应同时看 status、complete、n、CI 和 screened_mean**。review_required 表示模型诊断、重复/方法一致性或采样存在问题。NE 字段带 diagnostic；不等于总集体电导率。内部 STACIE 误差与独立重复的置信区间分别报告，不能混为同一种误差。

## 验证范围

单元测试覆盖单位转换、全部重复参与、共同窗口、缺失值、配置一致性和异常曲线。合成电荷矩测试使用已知各方向增量方差的随机游走，核对 STACIE 电导率与解析值及重复间区间。真实本地轨迹测试会保留所有质量警示。通过这些测试验证的是实现及有限样例，不能证明每一种 MOF 的采样都足够；若没有足够扩散/电荷传输事件，正确输出应为不足以确认，而不是人为制造一个确定值。

配置变化后可保留同一原始输入缓存，例如在 Slurm 中运行：

```bash
python scripts/transport_framework.py batch --runs-root runs --cache-root transport_analysis/cache --output transport_analysis_newconfig --config config/transport.json --workers 2
```

这会重做统计估计，复用核验过的坐标缓存；不会重复读取所有 DCD。`cache` 子命令也可提前为每条轨迹准备缓存。请勿并行运行多个会写入同一个缓存或输出根目录的作业。

## 两台机器汇合

两台机器必须使用完全相同的配置文件。将各自 `cases/NNN/results.json` 和 `curves.csv` 复制至同一台机器上的两个目录后即可合并；不必转移 DCD 或坐标缓存：

```bash
python scripts/transport_framework.py merge --inputs local_analysis other_machine_analysis --output merged_analysis --manifest manifest.csv --config config/transport.json
```

框架拒绝重复 case 编号或不同配置的混合；合并后会重新确定同一体系的共同 MSD 平台和全部重复的误差，独立机器上生成的旧汇总均值不会被直接拼起来。若某机器仍用 v1–v5，必须先对已有轨迹运行本框架，不能从旧的几个系数重建频谱结果。
