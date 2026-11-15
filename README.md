# MOF + LiFSI + DMC

168 个已构建体系：6 个 MOF × 7 个 DMC/LiFSI 比例（2、4、6、8、10、12、14）× 4 次独立重复，每个体系 50 对 LiFSI。构建输入在 `cases/000`–`cases/167`，组成与重复映射见 [manifest.csv](manifest.csv)。MD 运行方法、原服务器状态快照和历史后处理记录保留在 [历史文档](docs/HISTORICAL_METHODS_AND_MD_GUIDE.md)。

## 统一后处理：唯一推荐入口

**脚本：[transport_framework.py](scripts/transport_framework.py)。固定配置：[config/transport.json](config/transport.json)。**

这套框架只处理已经完成的轨迹，不运行 MD。采用 STACIE 1.3.0 估计电导率和谱模型误差；用全部轨迹的 MSD、同一体系重复共享的扩散平台区间、独立重复误差和 EH/NE 对照诊断可靠程度。当前数据接口针对本项目完整 15 段 DCD/restart 和固定晶胞。

```bash
# 在已有 Python 3.11+ 环境中安装一次，不需要克隆其他仓库。
python -m pip install -r requirements-transport.txt
# runs/NNN 下已有完整轨迹后，在仓库根目录提交；先按集群设置分区。
sbatch slurm/transport_pipeline.slurm
```

默认输出到 `transport_analysis/`：坐标缓存、每条轨迹的频谱与 MSD 图和 JSON、按体系共享拟合区间的图、逐体系和汇总 CSV。默认两个进程、24 GiB 内存；可通过 Slurm 和环境变量指定输入/输出：`DMC_RUNS_ROOT`、`DMC_ANALYSIS_OUTPUT`。

```bash
# 中断后续跑；成功项校验后复用，失败尝试保留。
DMC_RESUME=1 sbatch --export=ALL slurm/transport_pipeline.slurm
# 不重读轨迹，重做统计汇总。
python scripts/transport_framework.py summarize --runs-root transport_analysis --output new_summary
# 两台机器的 cases/ 结果与曲线复制到同机后合并，不需要转移 DCD。
python scripts/transport_framework.py merge --inputs local_analysis other_analysis --output merged_analysis
```

两台机器必须使用同一代码、同一配置和一致的电荷/体积约定。框架拒绝不同配置或重复 case 编号的合并。修改配置时使用新输出目录，可用 `--cache-root` 复用原始坐标缓存。不要同时运行两个写同一输出目录的批次。

**解释输出时必须区分描述值与筛选值。** `by_state.csv` 的 `mean` 包含全部可用重复，`screened_mean` 只在完整与质量条件满足时填写。`review_required`、`no_plateau`、`incomplete` 不是零。双对数接近 1、谱拟合通过或误差条较小，都不能独立证明体系平衡或实验准确度。

完整公式、近期论文/GitHub依据、误差解释和测试边界见 [统一框架说明](docs/TRANSPORT_FRAMEWORK.md)。旧 v1–v5 只作追溯，不与本框架结果混表。

## 数据与旧结果

仓库不包含轨迹和重启二进制文件。全部168个体系的状态见 [统一结果索引](results/transport_framework_case_registry.csv)，本地42个结果已提供 [可合并包](results/transport_framework_local42_merge_bundle.tar.gz)。本地测试结果见 [统一框架逐体系表](results/transport_framework_local_20261115_per_case.csv)、[按体系汇总](results/transport_framework_local_20261115_by_state.csv) 和 [验证记录](results/transport_framework_validation_20261115.json)。解释方式和历史结果见 [results/README.md](results/README.md)；旧算法表不能与统一框架结果混用。初始结构在原集群通过了 100 步检查，但这不等于 15 ns 输运系数已收敛。
