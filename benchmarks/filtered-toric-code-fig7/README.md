# Filtered Toric Code：Fig. 7 的八个参数点

本目录计算 *Shadows of anyons* [Fig.7(a–h)](https://www.nature.com/articles/ncomms9284/figures/7) 的无限行转移谱。目标是定性复现低支、谱隙随过滤参数的变化，以及任意子凝聚后出现的禁闭特征。当前计算进度和最终图见 [REPORT.md](REPORT.md)。

**八图已完成。** [查看4行2列总图](results/figures/figure7.svg)，各子图及低支数值见报告。七组新增正式计算合计约62分钟；原(a)结果直接复用。

原来的 [Fig.7(a) 教学实现](../filtered-toric-code/README.md) 单独保留。它从顶点 V、边张量 P 和 bond tensor 推导四自旋张量 A，解释虚拟 Hadamard 变换、边界 MPS 与有效转移矩阵。本目录沿用完全相同的张量构造和块晶格，新增的是不同物理相中的边界表示和激发分类。

## 模型与八个参数点

无限方格的每条边有一个自旋，解析态为

$$
|\psi(\beta_z,\beta_x)\rangle\propto
\prod_i\exp[(\beta_x X_i+\beta_z Z_i)/4]|\mathrm{TC}\rangle.
$$

X、Z 为 Pauli 矩阵。每个选定棋盘 plaquette 的四个物理自旋合为一个 PEPS 张量 A，虚拟维数 D=2；每条边只属于一个块。块晶格基矢为 $(1,1)$ 和 $(-1,1)$，一个动量单位对应沿边界平移一个块。过滤参数 β 是波函数参数。

以下坐标按论文 [Fig.6](https://www.nature.com/articles/ncomms9284/figures/6) 的标记读取，用于定性复现：

| 子图 | βz | βx | 物理相 |
|---|---:|---:|---|
| a | 0.55 | 0.45 | topological |
| b | 0.80 | 0.45 | topological，接近 charge condensation |
| c | 1.05 | 0.45 | charge condensed |
| d | 1.30 | 0.45 | charge condensed |
| e | 0.45 | 0.55 | topological |
| f | 0.45 | 0.80 | topological，接近 flux condensation |
| g | 0.45 | 1.05 | flux condensed |
| h | 0.45 | 1.30 | flux condensed |

八图排列与论文一致：左列 a、b、c、d，右列 e、f、g、h。

## 同一个行转移算符，三种边界结构

沿物理指标收缩 A 与其共轭得到双层张量 Q，再沿水平方向连接成无限行 MPO，记作 T。每个边界位点的局部空间由 ket、bra 两个虚拟指标组成，维数为4。由于 A 的虚拟 Z₂ 对称性，T 与两个全边界操作分别对易：

$$
G_k=\bigotimes_j(Z\otimes I)_j,\qquad
G_b=\bigotimes_j(I\otimes Z)_j.
$$

VUMPS 用键维数χ的无限 MPS 近似 T 的一个主导右 fixed point r。同一个 T 的这些群作用可以产生其他 fixed points，但产生多少个不同的态由物理相决定。将边界限制在错误的对称性表示中，会遗漏破缺态或把本来相同的背景重复当作 domain wall。因此三个文件分别说明并实现其边界假设。

### `topological.jl`：a、b、e、f

选择的 r 保持 $C=G_kG_b$，但不保持 $G_k$ 或 $G_b$；另一个 fixed point 为 $r_1=G_kr$。边界 MPS 使用 C 的偶奇 charge 分块，和旧 Fig.7(a) 相同。

局部激发张量 B 的左侧使用 r，右侧可以使用 r 或 $r_1$，分别标记 g=0、1；后者为 domain wall。激发态在 C 下的本征值为 $(-1)^q$，q=0、1。四个 $(g,q)$ 对应 $1,e,m,\epsilon$，蓝/红表示 g，+ / ○ 表示 q。

### `charge_condensed.jl`：c、d

此时四个 fixed points 为 $r,G_kr,G_br,G_kG_br$，完整的 Z₂×Z₂ 对称性在所选的 r 中破缺。为允许这种结构，边界 MPS 使用普通复数空间，不强制任何 charge 分块。初始张量含 ket、bra 均沿虚拟 X 正方向的分量，用来选择一个破缺分支；随后由普通 T 进行 VUMPS 优化。

固定激发左侧背景为 r，计算三种右侧背景：

- r：普通局部激发，蓝色点。
- $G_kr$：单层 string 产生的 domain wall，红色点。bra 单层 string 与它由两层交换关联，图中采用一个代表。
- $G_kG_br$：双层 string 产生的 domain wall，绿色点。

这里 q 不再是所选破缺背景上可用的 charge 标签，因此使用点而不是 + / ○。双层 string sector 的最低谱量给出非零 string tension：在归一化重叠中，长度 L 的 string 带来随 L 指数衰减的因子，对应 flux 的禁闭。比较 c→d 时绿色低支的变化，是这一列的主要目标。

### `flux_condensed.jl`：g、h

r 同时保持 $G_k$ 和 $G_b$，两个单层群作用都回到同一个 fixed point。所有激发都使用相同的左右背景；相应地，边界 MPS 使用 Z₂×Z₂ 分块，每个激发标记两层独立 charge $(q_k,q_b)$。

$(0,0)$ 用蓝色 +，$(1,0)$ 和 $(0,1)$ 用蓝色 ○，$(1,1)$ 用蓝色 □。程序分别计算两种奇偶不同的 charge，保留其简并信息。双奇 charge sector 的谱隙对应 charge 的 string tension；它随 g→h 增大，体现 charge 禁闭。

程序将双层物理空间按 $(0,0),(0,1),(1,0),(1,1)$ 排序，并在 Q 的每条腿上施加同样的指标重排。这个操作仅改变数组存储顺序，T 所代表的物理收缩不变。

## 如何得到动量谱

在边界位置 j 插入 B，左右分别接上述选定的两个 MPS 背景，再取位置叠加

$$
|\Phi_k(B)\rangle=\sum_j e^{ikj}
|\cdots M_LM_L B_j M_RM_R\cdots\rangle.
$$

k 是沿边界的 Bloch 动量，图中横轴为 $k/\pi$。MPSKit 的左正交规范消去 B 的表示冗余。`sector_modes` 在剩余参数空间逐列建立 T 的有效矩阵并对角化，按本征值模长排序输出前 `num` 模。三个相分别传入 Z₂、普通空间的平凡 charge 和 Z₂×Z₂ charge，算法本身相同。

有效作用已按背景的主导每格尺度归一化，保存本征值比值 $\eta=\lambda/\lambda_0$，纵轴为 $\varepsilon=-\log|\eta|$。这是一维边界的激发变分谱，χ决定背景和局部激发参数空间的大小。

为与论文的横轴约定一致，单层 string 的右背景每格张量乘以−1；这等价于将该 domain-wall 态的动量平移π。双层 string 的右背景不变号，因此绿色低支的最低点出现在 k=0。无 domain wall 的 sectors 均不进行这种变换。

## 精度、保存内容和运行方法

默认χ=16、41个动量、每sector64模。边界环境容差 `--tol=1e-6`，VUMPS 内部停止目标为其十分之一；激发环境容差 `--spectrum-tol=1e-7`。容差的具体数学定义见旧 [README 第4节](../filtered-toric-code/README.md#4-精度参数和输出数值是什么意思)。结果中保存每个本征解的残差 $\|T_{\rm eff}b-\eta b\|$，它与有限χ的表示误差应分别理解。

`boundary.json` 保存边界张量、实际虚拟空间、局部残差和三个群作用的每格重叠，顺序为 ket、bra、两层同时作用。拓扑相预期前两者小于1、第三个为1；charge condensed 相三个都小于1；flux condensed 相三个都为1。这些量用于确认所使用的背景与相分类相符。

从仓库根运行，每次指定新的结果目录：

```bash
timeout --kill-after=10s 1820s env OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env benchmarks/filtered-toric-code-fig7/run.jl --panel=b --chi=16 --nk=41 --num=64 --budget=1800 --out=benchmarks/filtered-toric-code-fig7/results/b-chi16
```

每个动量和 sector 完成后立刻刷新 `spectrum.csv` 并打印进度。输出保存解析张量、边界、参数、sector 名称与绘图符号、运行源码、完整谱和实际运行时间。`source/` 是该次结果的复现记录，日常阅读只需本目录顶层程序。

`plot.py` 读取一个 JSON 文件，内容为子图字母到运行目录的映射，路径相对于该 JSON 文件。它可以直接读取旧(a)的结果，统一绘制八个独立 SVG 和 `figure7.svg`：

```bash
python3 benchmarks/filtered-toric-code-fig7/plot.py --runs=benchmarks/filtered-toric-code-fig7/results/runs.json --out=benchmarks/filtered-toric-code-fig7/results/figures
```

建议阅读顺序：先读旧(a)的张量推导，再读本说明的三相区别，随后查看对应相文件、`run.jl`、`plot.py`。
