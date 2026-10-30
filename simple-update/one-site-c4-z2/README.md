# One-site/C₄/Z₂ 约束下的零场 loop/cluster update

这里尝试从随机对称初态求无限平方晶格 Toric Code 的零场基态，Hamiltonian 为

$$
H=-\sum_s A_s-\sum_p B_p,\qquad A_s=\prod_{i\in s}X_i,\quad B_p=\prod_{i\in p}Z_i.
$$

采用 Pauli convention，$J_e=J_m=1$，物理自旋位于 edge，$d=D=2$，能量按 edge 自旋归一化。全晶格共享一组 V/P/λ，保持 one-site 平移、C₄ 和原生虚拟 Z₂ 对称性。

**两条随机轨迹 seed=5、6 各完成700轮，均停在 $e=-0.5$，未得到零场基态。** 当前不追加随机寻优，也没有从这些失败态生成正式转移谱。具体结果与两种受困机制见 [REPORT.md](REPORT.md)，完整逐轮数据在 `results/seed5/` 和 `results/seed6/`。

## 如何运行

从仓库根目录执行，输出目录必须不存在：

```bash
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env simple-update/one-site-c4-z2/run.jl --seed=6 --steps=300 --dt=0.05 --als-sweeps=12 --budget=580 --out=simple-update/one-site-c4-z2/results/new-seed6-segment1
```

`--steps` 是本段新增轮数。一次完整 sweep 先做 star 更新，再做 plaquette 更新；`--als-sweeps` 是每次 plaquette 更新内的 ALS 轮数。`--budget` 的上限为600秒，独立 Julia 子进程接受总时间监督，导入、编译、更新和测量均占用预算；监督进程预留少量时间写入预算终止说明。预算较短时可能在首次编译期间结束。

每轮更新后保存 `latest.json` 并追加 `trajectory.jsonl`；每100轮做 χ=8 测量，段末用 χ=8/16 验收。正常结束写入 `final.json` 和 `summary.json`；预算终止保留最近完整状态。继续同一轨迹时，指定前段状态，其他数值设置不变，并使用新的输出目录：

```bash
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env simple-update/one-site-c4-z2/run.jl --seed=6 --steps=200 --dt=0.05 --als-sweeps=12 --budget=580 --resume=simple-update/one-site-c4-z2/results/new-seed6-segment1/latest.json --out=simple-update/one-site-c4-z2/results/new-seed6-segment2
```

正常结束的段优先从 `final.json` 继续，预算结束的段从 `latest.json` 继续。按保存状态的 `metadata.step` 决定后续轮数，不假设相同预算在不同机器上能完成相同轮数。验收以 `summary.json` 的 `passed` 和末态 metadata 的 `accepted` 为准，不能以进程正常结束代替物理验收。

## 张量表示与局部更新

顶点张量为 $V_{nesw}$，四腿依次指北、东、南、西；边张量为 $P^p_{ab}$，包含物理指标 $p=0,1$ 及 head、tail 两条虚拟腿。每条 V–P 连接带一个非负对角权重 λ。所有顶点共用 V，水平和竖直边共用 P，所有连接共用 λ。

虚拟空间为 $\mathbb C_0\oplus\mathbb C_1$，偶、奇 charge 各一维，群作用为 $G=\operatorname{diag}(1,-1)$。V 仅保留 $n+e+s+w$ 为偶数的分量；P 仅保留 $a=b$ 的块，物理空间的群作用取平凡表示。空间约束为

$$
V_{nesw}=V_{eswn}=V_{swne}=V_{wnes},\qquad P^p_{ab}=P^p_{ba}.
$$

V 的 C₄ 是循环旋转，并不要求任意指标交换对称。随机初始化仅在允许 charge 块内抽取实正态随机数，再施加上述空间投影。解析 parity V 与 X-basis copy P 仅用于代数检查。

**Star 更新**收缩中心 V 与四条 P，内部连接和外接腿各计入一次 λ，施加 $\exp(dt A_s)$。把每条臂的“物理腿、外侧虚拟腿”分组，对四组分别做 HOSVD，每组固定保留一偶、一奇，再重构中心 V 和四条 P。

**Plaquette 更新**保留四个 V、四个 P 组成的完整闭环，构造 $\exp(dt B_p)$ 作用后的目标。ALS 每次固定其余七个因子，在允许 charge 块内求当前因子的线性最小二乘解，循环更新八个因子；支撑上的伪逆处理秩亏问题。外环境仍由独立对角 λ 近似，因而这是 loop/cluster update，不是完整二维环境的 full update。

闭环的左右 cavity 密度矩阵由沿环传播并归一化得到；它们的正平方根与 SVD 确定成对变基和新的连接权重。对各位置张量作成对的 charge 块符号对齐后，程序检查加权区域张量是否保持不变，随后恢复共享参数并施加空间投影。纯 gauge 误差与共享恢复产生的变化分别记录。完整 sweep 在副本上执行，子更新失败时不替换原状态。

测量复用 `simple-update/split/`：将同一组 V/P/λ 复制到2×2测量布局，合并为物理维数4的 PEPS，再用 CTMRG 收缩。这不会增加变分晶胞。χ 是 CTMRG 环境维数；局部更新本身不使用这一完整二维环境。

## 代码与结果位置

| 文件 | 作用 |
|---|---|
| `run.jl` | 随机优化、分段续算、逐轮保存和时间预算监督入口 |
| `core.jl`、`tensors.jl` | 模块加载、Z₂ 张量表示、初态与数组转换 |
| `star.jl`、`plaquette.jl` | star HOSVD、闭环 ALS、cavity 变基与局部 stabilizer 残差 |
| `shared.jl` | 对称收缩、charge 块对齐、共享参数恢复与完整 sweep |
| `state_io.jl`、`measure.jl` | 状态读写、预算监督、物理验收与现有 CTMRG 测量适配 |
| `tests.jl` | 集中的轻量程序检查，不参与正式优化 |

`results/seed5/segment{1,2,3}/` 和 `results/seed6/segment{1,2,3}/` 保留全部原始逐轮 JSONL、分段设置、初态和段末状态。seed=5 的独立零更新测量在 `results/seed5/measurement/`，最终物理审查在 `results/seed5/audit.json`；seed=6 的最终测量直接见 `results/seed6/segment3/summary.json`。分段和保存点的具体对应见报告。

## 程序检查

```bash
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env simple-update/one-site-c4-z2/tests.jl
```

可用 `--group=tensors`、`--group=updates`、`--group=runner` 选择表示、局部更新或读写与运行控制检查。检查使用临时目录，不在 `results/` 中生成新材料。解析固定点和程序检查通过，不等于随机优化已经得到基态。
