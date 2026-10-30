# Toric Code 基态与任意子转移矩阵谱

本项目先用精确可解的 Toric Code 检查 PEPS 算法，再计算磁场中的基态，随后施加虚拟 Z₂ 对称约束，构造 transfer matrix，按对称性求解不同 sector 的 spectrum，最终研究任意子关联函数与关联长度。

**当前进展：零场 loop/cluster update 已从随机态得到基态；三个弱场点已有 loop 和 AD 结果。one-site/C₄/Z₂ AD 已得到弱场候选态及其四sector转移谱，部分谱支仍有χ依赖。** 后续对称基态计算从解析零场态出发；完整跨相变路径和任意子关联函数尚未完成。

## 模型与主要结果

```text
H = −Σₛ Aₛ − Σₚ Bₚ − hₓΣᵢXᵢ − h_zΣᵢZᵢ
Aₛ = ∏star X，Bₚ = ∏plaquette Z，Jₑ = Jₘ = 1
```

X、Z 为 Pauli matrices。物理自旋在无限平方晶格的 edge 上，d=2；现有部分计算使用含8个自旋的2×2顶点重复单元，one-site 约束实现则共享同一组张量，具体表示见相应算法说明。所有能量按 edge 自旋归一化。零场严格值 E/N=−1，所有 stabilizer 满足 ⟨Aₛ⟩=⟨Bₚ⟩=1。8-spin 周期小系统的 ED 仅用于模型检查，不代替无限平面计算。

| 计算 | 已得到的结果 | 适用范围 |
|---|---|---|
| 零场 loop/cluster update，D=2 | E/N≈−0.999999999913；最大 stabilizer 误差约 4.42×10⁻¹⁰ | 随机 seed=1，χ=8/16 一致 |
| 有限场 loop，(0.1,0)，D=2 | 能量误差约 1.76×10⁻⁶/spin，磁化误差约 0.143% | 达到既定弱场精度 |
| 有限场 loop，(0,0.1) 和 (0.1,0.1)，D=3 | 能量误差约 4.2–5.9×10⁻⁶/spin，Z 磁化误差约 0.77% | 来自 D=2 的增大 D 延续；双分量态另有 VUMPS 测量支持 |
| 有限场 split AD，三个弱场点，D=2 | 能量误差≤10⁻⁵/spin，非零磁化误差≤1% | 解析零场或近基态起点，采用 unrolled AD；随机全局寻优尚未通过 |

有限场 loop 的 D=2 Z 磁化误差约 1.16%，未达到 1% 标准；D=3 X 单轴尚有局部拟合和较大 χ 的问题；D=4 双分量计算未通过。详见 [零场总结](M2_REPORT.md)、[有限场总结](M3_REPORT.md) 和 [主要结果表](RESULTS.csv)。

严格零场的 [one-site/C₄/Z₂ SU](simple-update/one-site-c4-z2/README.md) 共享一组顶点张量 V、含边自旋的张量 P 与连接权重 λ，采用原生虚拟 Z₂ 分块。seed=5、6各完成700轮，均停在 $e=-0.5$，分别只有 star、plaquette 期望值接近1；seed=6的χ=8/16环境均收敛，仍未通过物理验收。这条SU路线没有提供合格的正式谱输入，当前不追加随机寻优；两种受困机制及完整逐轮数据见 [SU报告](simple-update/one-site-c4-z2/REPORT.md)。下述AD候选态的谱是另一条独立路线的结果。

## SU 与 AD 分别放在哪里

施加 one-site 平移、C₄ 和虚拟 Z₂ 的 AD 位于 [ad/one-site-c4-z2/](ad/one-site-c4-z2/README.md)，从两种初态得到的结果需分开理解：

| $(h_x,h_z)=(0,0.1)$、D=4 的初态来源 | 结果 | 后续安排 |
|---|---|---|
| 零场解析态加小扰动 | $e=-1.002515869800994$，与低场级数相差约 $8.16\times10^{-6}$/spin；梯度约 $4.14\times10^{-6}$。χ=8 测量重复一致，较大χ旧高精度环境检查未通过 | 作为后续对称 AD 的主线 |
| 独立随机态，多 seed | 最佳 seed=23 能量约 $-0.76707219$，梯度约0.214，最后线搜索失败，未达到 $10^{-4}$/spin 目标 | 已停止，不再追加随机寻优 |

各 seed 的汇总和数据位置见[AD报告](ad/one-site-c4-z2/REPORT.md)。固定解析来源D=4候选张量的四sector谱已完成：χ=8、64动量、每sector8模（2048点），另完成χ=12三个共同动量比较，总用时约30.4分钟。e、m、ε主导模在已比较动量上稳定，平凡sector与高支仍有χ依赖；图与结果见[谱报告](transfer-spectrum/REPORT.md)。

两种算法用文件夹区分，里面也记录了不同实现方式：

| 目录 | 内容 |
|---|---|
| [model/](model/README.md) | Hamiltonian、解析零场张量、ED、低场级数参考 |
| [simple-update/](simple-update/README.md) | SU 类方法总说明：旧 merged L 形 SU、split path/loop 对照、一般 D 的 loop 方法 |
| [simple-update/one-site-c4-z2/](simple-update/one-site-c4-z2/README.md) | one-site/C₄/原生虚拟 Z₂ 的零场 SU；两条失败轨迹、完整逐轮数据与结果报告 |
| [ad/](ad/README.md) | AD 总说明：旧 merged tied-tensor 零场优化，以及新 split V/P/λ 的弱场优化 |
| `tests/` | 保留算法的测试、历史张量读取检查及结果汇总脚本测试 |
| [benchmarks/filtered-toric-code/](benchmarks/filtered-toric-code/README.md) | 独立的解析过滤态 Fig.7(a) 复现、简短运行入口与代码阅读说明；与磁场 Hamiltonian 研究分开记录 |
| [benchmarks/filtered-toric-code-fig7/](benchmarks/filtered-toric-code-fig7/README.md) | 完整 Fig.7(a–h) 的三相谱计算与总图；保留原(a)实现，新增不同物理相的边界及 sector 构造 |

每套实现的 `results/` 保存本地计算数据，`figures/` 保存该实验的图，`frozen_states/` 保存可随代码携带的最终张量 JSON。项目根目录的 `RESULTS.csv` 只汇总主要结果，注明算法、初态、D/χ、状态和详细数据位置，不用它替代原始计算数据。

## 如何运行

命令均从仓库根目录运行。`julia-env/Project.toml` 和 `Manifest.toml` 固定 Julia 1.12.6 的依赖环境，使用 `--project=julia-env`，不需升级软件包。新计算请选择未使用的输出目录，避免覆盖已有结果。

```bash
# 零场：loop/cluster update，另带解析零场对照
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env simple-update/split/run.jl --init=both --plaquette=loop --seed=1 --schedule=0.05:200,0.01:200,0.002:200 --measure-every=40 --out=simple-update/split/results/new-zero-loop

# 有限场：从随机态求 (hₓ,h_z)=(0.1,0)，D=2
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env simple-update/loop/run.jl --D=2 --hx=0.1 --hz=0 --every=40 --out=simple-update/loop/results/new-x

# 全部保留的 Julia 实现测试；不做新的优化或磁场扫描
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env tests/runtests.jl

# 根据已有表格重新生成主要结果表，不执行物理计算
python3 scripts/summarize_results.py
```

Split AD 的弱场命令见 [ad/split/README.md](ad/split/README.md)，one-site/C₄/Z₂ 的解析初态入口见[对应说明](ad/one-site-c4-z2/README.md)。首次 Julia 编译显著慢于后续调用；运行前需估计实际时间和内存。

Python 数据汇总脚本需要 Python 3.11 或更新版本；`transfer-spectrum/plot.py` 使用标准库生成 SVG。其他绘图脚本按实际导入需要 NumPy、Matplotlib，运行前请确认所用 `python3` 环境具备对应依赖。

整体后续安排见 [PLAN.md](PLAN.md)。
