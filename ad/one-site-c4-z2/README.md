# One-site/C₄/Z₂ 约束下的 Toric Code 能量 AD

这里计算无限平方晶格磁场 Toric Code 的 PEPS 基态，保持 one-site 平移不变性、C₄ 旋转对称性和虚拟 Z₂ 对称性。Hamiltonian 为

$$
H=-\sum_s A_s-\sum_p B_p-h_x\sum_i X_i-h_z\sum_i Z_i,
$$

采用 Pauli convention，$J_e=J_m=1$，物理自旋位于 edge，能量按 edge 自旋归一化。当前结果集中在 $(h_x,h_z)=(0,0.1)$、D=4。

**后续主线从零场解析态出发优化有限场基态。** 现有受小扰动解析初态得到 $e=-1.002515869800994$，与十阶低场级数相差约 $8.16\times10^{-6}$/spin；旧的较大χ高精度环境检查未通过，保留候选态身份。任意随机初态路线未取得基态，已停止研究。两条路线及全部 seed 的结果见 [REPORT.md](REPORT.md)，原始数据在 `results/`。

## 从解析初态运行

`run.jl` 是解析初态路线的入口。默认将零场 D=2 解析张量嵌入 D=4，再加入幅度0.03、seed=71 的对称参数扰动，使新增虚拟自由度参与优化。这仍是解析来源初态，不归入任意随机态寻优。

从仓库根目录执行，输出目录及同名 `.log` 应不存在：

```bash
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 python3 ad/one-site-c4-z2/supervise.py --init=exact --D=4 --hx=0 --hz=0.1 --noise=0.03 --seed=71 --budget=7200 --steps=1000 --out=ad/one-site-c4-z2/results/new-analytic
```

Python 只负责在指定总时间后停止 Julia；启动、编译、优化和最终测量均计入预算。需要直接运行 Julia 时，使用同样的参数调用 `julia --startup-file=no -O1 --project=julia-env ad/one-site-c4-z2/run.jl`。直接调用时程序在计算步骤之间检查时间，严格总时间限制由 Python 提供。`--noise=0` 表示不加扰动的解析嵌入态，与现有成功降能轨迹的起点不同。

每次接受更新后保存 `stepNNNN.json`，`trajectory.jsonl` 记录能量、梯度和拒绝原因；`final.json` 保存最后张量，`summary.json` 记录停止原因和接受状态。用户无需手动分段或操作续算。程序可选读入 `--probe=已有梯度记录.json` 核对设置，正常运行不依赖单独的梯度检查文件。

当前默认 CTMRG 环境维数为 χ=8、迭代容差 $10^{-8}$、最多300次；可用 `--chi`、`--tol`、`--ctm-maxiter` 指定。默认能量目标为相对低场级数的绝对误差 $10^{-4}$/spin，通过 `--energy-tolerance` 设置。复现已有结果时，应采用相应运行目录 `settings.json` 中的参数。

## 张量表示与优化方法

全晶格共享一个顶点张量 $V_{nesw}$、一个带物理自旋的 edge 张量 $P^p_{ab}$ 和对角连接权重 $\lambda$。$n,e,s,w$ 表示北、东、南、西方向的虚拟指标；$p=0,1$ 为物理指标，$a,b$ 为 P 的两条虚拟腿。每条 V–P 连接携带一个 $\lambda_a$。收缩时可把 V 与其北、东两个 P 合并为物理维数4的 PEPS 张量；2×2 测量布局中的张量均由同一组参数复制，并没有增加变分晶胞。

虚拟基矢的 charge 为 $q_a\in\{0,1\}$，定义 $G_{ab}=(-1)^{q_a}\delta_{ab}$。V 满足 $G^{\otimes4}V=V$，即只保留四腿 charge 和为偶数的分量，并满足

$$
V_{nesw}=V_{eswn}=V_{swne}=V_{wnes}.
$$

这施加的是 C₄ 循环旋转，不是任意指标交换。D=2 时，V 的允许分量分为0000、1111、相邻两个1、相对两个1四类循环轨道，各有一个系数。P 保留 $q_a=q_b$ 的块，并满足 $P^p_{ab}=P^p_{ba}$；D>2 时，同 charge 块内允许非对角分量。物理指标上的 Z₂ 表示取平凡作用，相同 P 用于水平和竖直 edge，连接权重与 G 对易。

允许分量用 Frobenius 内积下正交归一的张量基展开，V、P 分别归一化，$\lambda_a=\exp(u_a-\sum_bu_b/D)$。D=2 和偶奇 charge 各二维的 D=4 分别有10和54个原始实参数，含归一化及 gauge 冗余，不能理解为独立物理自由度。

能量通过双层 PEPS 与 CTMRG 环境收缩得到。AD 从已求出的参考环境出发，对20次固定展开的 CTMRG 迭代求导，再经链式法则把重复位置的贡献加回共享参数。L-BFGS 用梯度构造下降方向，Armijo 线搜索逐次减小步长，直到待更新 PEPS 的能量评价通过环境检查且满足下降条件。共同的收缩和求导实现复用 `ad/split/`，不执行 simple update。

`parameters.jl` 将共享参数转换为 `ad/split/` 所需的 V、北/东两个 P 和四个方向的连接参数，长度为 $D^4+4D^2+4D$，D=2、4 时分别为40、336。重复填入的是同一组物理参数，其导数须汇总回原来的10或54个系数。梯度方向差分使用同一参考环境和相同的 CTMRG 展开次数，以比较同一个有限计算图；20与40次展开的结果用于检查梯度对展开长度的敏感性。

AD 收缩把每条连接的 $\sqrt{\lambda}$ 分配到两端，谱模块则按其 blocking 约定吸收完整 λ。两者的边界 gauge 不同，单个 composite 张量的数组不必逐分量相同；应比较相同闭合网络的波函数收缩。保存的状态包含原始 V/P/λ，谱计算通过统一的 [V/P/λ 输入适配](../../transfer-spectrum/README.md#磁场数值态的输入与运行) 构造行算符。

当前每次评价使用两个随机环境初猜，要求主导范数网络收缩结果重复出现，并检查能量一致性与局部物理界。CTMRG 容差约束 corner/edge 奇异值在相邻迭代间的变化，不是能量误差界。不同初猜的网络范数相对差仍要求不超过 $10^{-8}$；这一判据对本项目较低精度目标是否必要，尚未校准。因此，未通过该检查只表示评价没有被程序接受，不能直接认定该 PEPS 不接近基态。

当能量达到目标或梯度足够小时，程序对冻结张量用新的两个环境初猜分别计算 χ=8 和12：要求能量均满足级数目标、两种χ能量差不超过 $10^{-5}$，并满足局部物理界；在 $h_x=0$ 时另检查所有 $|\langle B_p\rangle-1|\le10^{-3}$。现有解析来源候选态的能量满足目标，但未通过其较大χ的环境检查，具体结果见报告。

## 随机初态路线与代码结构

`random.jl` 是独立的随机路线入口，用于方法比较。固定一个 seed 就调用一次与解析路线相同的优化；给定多个 seed 时只增加种子遍历和总时间分配。随机路线已停止研究。

```bash
# 仅说明保留的用法；本轮不再执行随机寻优。
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 python3 ad/one-site-c4-z2/supervise.py --init=random --seeds=19,23,41 --D=4 --hx=0 --hz=0.1 --budget=21600 --steps=1000 --out=ad/one-site-c4-z2/results/new-random
```

只给一个 seed（如 `--seeds=19`）时，该初态使用整段预算。多个 seed 按剩余种子数分配剩余时间，失败的初态不占用尚未使用的预算。直接 Julia 入口为 `random.jl`；两种初态共享底层优化，不复制算法。

| 文件 | 作用 |
|---|---|
| `run.jl` | 解析初态运行入口，后续主要使用 |
| `random.jl` | 随机初态入口；单 seed 或多 seed 调度 |
| `parameters.jl` | one-site/C₄/Z₂ 参数构造、解析与随机初态、张量读写 |
| `objective.jl` | 参数到 PEPS 的能量与梯度适配 |
| `optimize.jl` | 两条路线共同的 L-BFGS、线搜索、保存和最终能量检查 |
| `core.jl` | 加载共同计算模块，供运行与谱程序使用 |
| `supervise.py` | 两条路线共用的外部时间预算管理 |
| `tests.jl`、`supervisor_tests.py` | 参数与程序控制的测试，不参与基态寻优 |
