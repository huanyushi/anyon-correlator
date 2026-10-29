# 从 V/P 网络理解 filtered Toric Code 的转移谱

本目录复现 Haegeman 等 [*Shadows of anyons*](https://www.nature.com/articles/ncomms9284) Fig.7(a)。输入是已知的解析过滤态；主项目中磁场 Hamiltonian 的数值基态研究另行记录。

最终图采用χ=24、61个动量、每个 sector 96个 modes。两个孤立低支的最低谱量约为1.13877、1.59674。图和χ=8、16的比较见 [REPORT.md](REPORT.md)。下面先从原始 V/P 网络说明张量构造，再解释谱计算中的符号。

## 1. 从顶点 V、边 P 和 bond tensor 出发

### 原始网络及解析张量

物理自旋在无限平方晶格的边上。每个 vertex 放置一个四虚拟腿张量 V；每条 edge 放置一个两虚拟腿、一物理腿张量 P。V 和 P 之间通过 bond tensor Λ 连接。这就是本项目使用的 V/P 图像。

对零场 Toric Code 的解析态，虚拟维数为 D=2。取所有二元指标为0或1，仓库中的张量可写成

$$
V_{i_1i_2i_3i_4}\propto
\delta_{i_1\oplus i_2\oplus i_3\oplus i_4,0},
\qquad
P^p_{ij}\propto\delta_{ij}(-1)^{pi},
\qquad \Lambda\propto I_2.
$$

$\oplus$ 表示模2加法，p 是边自旋在物理 Z 基底中的指标。这里省略了各张量的整体归一化常数。解析态的两个 bond 权重相等，因此 Λ 是单位矩阵乘一个标量；该标量可以吸收进张量的归一化。这个性质使下面的解析重组很简单。

### 将同一个网络改写成 copy V 和 XOR P

为了看清四自旋块的来源，先在每条 V–P 连接的两个张量端点同时变换虚拟基底，使用

$$
H=\frac{1}{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
$$

由于 $H^TH=I$ 且 $\Lambda\propto I$，同一条连接上的两个 H 在收缩时抵消，整个物理态保持不变。这个操作只改变虚拟指标的表示。

变换后的张量为

$$
\widetilde V_{abcd}\propto
\delta_{ab}\delta_{bc}\delta_{cd},
\qquad
\widetilde P^p_{ab}\propto\delta_{p,a\oplus b}.
$$

这时一个 vertex 的四个虚拟指标必须取相同的值，可以用一个变量 $a_v$ 表示。连接 vertices v、w 的 edge 上，P 要求物理指标满足 $p_{vw}=a_v\oplus a_w$。所以接下来出现的“角点变量”，就是**原来 vertex 上的共同虚拟变量**，没有增加新的物理自旋或自由度。

### 选择 plaquette，将网络精确重组为四自旋块 K

把原方格的 plaquettes 按黑白棋盘方式着色，只选黑色 plaquettes。每条物理边恰好属于一个黑色 plaquette，所以把每个黑色 plaquette 的四条边放进一个块，不会遗漏或重复计算物理自旋。

一个 vertex 同时被两个对角相邻的黑色 plaquettes 共享。它的 copy tensor 可以精确拆为两个三腿 copy tensors，通过一条新虚拟腿相连：

$$
\widetilde V_{i_1i_2i_3i_4}\propto
\sum_a C_{i_1i_2a}C_{i_3i_4a},
\qquad C_{ija}=\delta_{ia}\delta_{ja}.
$$

其中 $i_1,i_2$ 是朝向一个被选 plaquette 的两条腿，$i_3,i_4$ 朝向另一个被选 plaquette。每个三腿 C 分给其中一个块。对两个 C 收缩 a 后，仍然要求原来的四个指标相等。

现在，每个块包含四个 P 和四个 C。收缩它们内部的虚拟连接后，留下四个物理指标和四个朝向相邻块的虚拟指标，得到张量 K。用原方格的方向画一个块：

```text
a —— p₁ —— b
|           |
p₄    K     p₂
|           |
d —— p₃ —— c
```

图中的 a、b、c、d 是这个 plaquette 四个原始 vertices 的虚拟变量，p₁至p₄是四条边的物理指标，中心的 K 表示整个块。四个 P 分别给出

$$
K^{p_1p_2p_3p_4}_{abcd}=
\delta_{p_1,a\oplus b}\delta_{p_2,b\oplus c}
\delta_{p_3,c\oplus d}\delta_{p_4,d\oplus a}.
$$

相邻块在共享 vertex 处收缩同一个变量，正好恢复前面拆开的 copy V。**`tensor.jl` 直接写出这个收缩结果 K，而不是在运行时逐个建立 V、P、Λ、C 再收缩。** 上述推导说明了它与原始 V/P 网络的关系。

### 代码为什么还对 K 的虚拟腿作用 H？A 是什么？

同时将 $a,b,c,d$ 全部翻转，四个 XOR 都不变，所以 K 满足虚拟对称性 $X_{\rm v}^{\otimes4}K=K$。这里下标 v 表明算符作用于虚拟空间。

共享求解器采用 Z₂ charge 基底存储张量，要求群生成元在每条虚拟腿上表示为 $Z_{\rm v}=\operatorname{diag}(1,-1)$。因为 $HX_{\rm v}H=Z_{\rm v}$，代码对 K 的四条外部虚拟腿作用 H，得到

$$
A_0^{p_1p_2p_3p_4}{}_{nesw}
=\sum_{abcd}H_{na}H_{eb}H_{sc}H_{wd}
K^{p_1p_2p_3p_4}_{abcd}.
$$

于是 $A_0$ 满足 $Z_{\rm v}^{\otimes4}A_0=A_0$。新块之间每条连接的两个端点都使用 H，因此收缩后的物理态仍然不变。这一步是为了使用对称分块的数值表示，**H 没有作用在物理自旋上**。

前面把 V/P 改写为 copy/XOR 是解释等价构造用的中间步骤；代码从 K 的闭式表达式开始，只实际执行这里的四腿 H 变换。

接着，对每条物理边施加过滤

$$
F=\exp[(\beta_xX+\beta_zZ)/4],
\qquad (\beta_z,\beta_x)=(0.55,0.45),
$$

得到

$$
A^{p_1p_2p_3p_4}_{nesw}\propto
\sum_{t_1t_2t_3t_4}
\left(\prod_{\ell=1}^4 F_{p_\ell t_\ell}\right)
A_0^{t_1t_2t_3t_4}{}_{nesw}.
$$

**A 是最终用于构造 transfer matrix 的局部 PEPS 张量。** 它包含四条原始物理边；程序把四个物理指标合并成 $p=p_1+2p_2+4p_3+8p_4$，所以 A 的数组形状为 `(16,2,2,2,2)`，依次是物理、北、东、南、西指标。最后除以 A 的 Frobenius norm，只改变整体归一化。

过滤作用在物理空间，因此保留 A 的虚拟 Z₂ 对称性。这一波函数是 $|\psi\rangle\propto\prod_i F_i|\mathrm{TC}\rangle$，过滤参数 β 与磁场 Hamiltonian 的 h 是不同参数。

新块的中心构成对角方格。它的水平、竖直基矢可取 $\mathbf a_x=(1,1)$、$\mathbf a_y=(-1,1)$，单位为原 vertex 间距。新块的北、东、南、西虚拟腿分别指向上图的 a、b、c、d，也就是原方格的左上、右上、右下、左下角。这解释了代码的方向标签为什么不同于原始方格的轴向标签。

## 2. Transfer matrix、边界 MPS 和两个 fixed points

将 A 与其复共轭沿全部物理指标收缩，得到双层张量

$$
Q_{(n\bar n)(e\bar e)(s\bar s)(w\bar w)}
=\sum_p A^p_{nesw}\overline{A^p_{\bar n\bar e\bar s\bar w}}.
$$

无横线的指标来自 ket 层，带横线的来自 bra 层。沿新块晶格水平方向连接一整行 Q，得到行转移算符 T；代码用无限 MPO 表示它。T 沿 $-\mathbf a_y$ 传播一步。边界链上每个位点的局部空间是一个 ket/bra 虚拟指标对，维数为 $D^2=4$。

VUMPS 用一个无限 MPS 近似 T 的主导右本征态，记为 $r_0$。χ 是这个边界 MPS 的键维数，与 PEPS 的 D=2 不同。除去主导本征值对应的尺度后，可以形式地写成 $\widehat T r_0=r_0$，因此称 $r_0$ 为 fixed point。实际无限系统的尺度归一化通过每格本征值密度实现。

A 的虚拟对称性使 T 分别与以下两个全边界操作对易：

$$
G_{\rm ket}=\bigotimes_j(Z_{\rm v}\otimes I)_j,
\qquad
G_{\rm bra}=\bigotimes_j(I\otimes Z_{\rm v})_j.
$$

**“单层虚拟 Z₂ 作用”指在边界的每个位点，只变换 ket 虚拟指标，bra 虚拟指标保持不变。** 它不是对物理自旋施加一个新的操作，也不是只作用于边界的一个位点。

令 $r_1=G_{\rm ket}r_0$，则由对易关系可得

$$
\widehat T r_1=\widehat T G_{\rm ket}r_0
=G_{\rm ket}\widehat T r_0=G_{\rm ket}r_0=r_1.
$$

这代表**从一个 fixed point 得到另一个 fixed point。** 在本参数所属的拓扑相中，选取的 $r_0$ 在 $G_{\rm ket}$ 作用下不保持不变，因而 $r_1$ 与 $r_0$ 是两个不同的边界态。两者都保持 ket、bra 同时变换的对称性

$$
C=G_{\rm ket}G_{\rm bra},\qquad Cr_0=r_0,\quad Cr_1=r_1.
$$

文档以下直接称它们为“两个 fixed points”。程序只用 VUMPS 求 $r_0$，再由已知群作用构造 $r_1$，不需要第二次独立优化。$r_0,r_1$ 都是同一个 T 的右 fixed points。

## 3. k、g、q 分别标记什么？

在边界 MPS 的位置 j 用一个可变张量 B 替换原张量，B 左侧使用 $r_0$ 的 MPS 张量，右侧使用 $r_g$ 的 MPS 张量。再对 j 求 Bloch 叠加，构造

$$
|\Phi_{k,g}(B)\rangle
=\sum_j e^{ikj}|\ldots M_0M_0\,B_j\,M_gM_g\ldots\rangle.
$$

这里 $M_0,M_1$ 是分别表示边界态 $r_0,r_1$ 的局部 MPS 张量，和二维 PEPS 张量 A 是不同对象。三个标签的含义是：

- **k：沿边界链的 Bloch 动量。** j 每增加1，表示沿 $\mathbf a_x$ 移动一个四自旋块，叠加系数多一个 $e^{ik}$。k 的周期为 $2\pi$，图的横轴是 $k/\pi$。
- **g：B 两侧是否使用不同的 fixed points。** g=0 时两侧都用 $r_0$，得到同一个边界态上的局部扰动；g=1 时左侧用 $r_0$、右侧用 $r_1$，B 位于两种边界态的界面，因此是一个 **domain wall**。
- **q：上述边界激发在 C 对称性下的 Z₂ charge。** 将 B 限制在满足 $C|\Phi\rangle=(-1)^q|\Phi\rangle$ 的子空间；q=0 为偶，q=1 为奇。C 同时作用于 ket、bra，因而这个 charge 表示两层 charge 相同还是不同。

在这里的拓扑相，四个 $(g,q)$ 子空间对应：

| g | q | 边界激发 | 任意子标签 | 图中符号 |
|---|---|---|---|---|
| 0 | 0 | 同一 fixed point 上的偶 charge 扰动 | 1 | 蓝色加号 |
| 0 | 1 | 同一 fixed point 上的奇 charge 扰动 | e | 蓝色空心圆 |
| 1 | 0 | 偶 charge 的 domain wall | m | 红色加号 |
| 1 | 1 | 奇 charge 的 domain wall | ε | 红色空心圆 |

这些分类没有改变行 MPO 本身；改变的是用于表示其本征态的边界条件和 charge 子空间。

固定 k、g、q 后，程序用 MPSKit 计算 T 在 B 的变分参数空间中的有效作用 $T_{\rm eff}(k)$，逐列建立矩阵并对角化。按本征值模长由大到小取前 `num` 个解，每个解就是这里所说的一个 mode。输出比值 $\eta_j(k)=\lambda_j(k)/\lambda_0$ 以及谱量 $\varepsilon_j(k)=-\log|\eta_j(k)|$。这里 λ₀ 表示主导尺度，实际有效作用已按 fixed point 的每格尺度归一化。

对于 g=1，两个 fixed points 的局部 MPS 张量还有相对相位约定。将右侧每个 $M_1$ 乘以−1，并同时把 k 平移π，表示同一个 domain-wall 态：在有限长度 L 的表示中，位置 j 右侧产生 $(-1)^{L-j}=(-1)^L e^{-i\pi j}$，恰好补偿 Bloch 相位变化。当前入口采用这个变号约定，与论文中 m 最低点出现在 $k=\pm\pi$ 的画法一致。

## 4. 精度参数和输出数值是什么意思？

### `--tol`：求边界 fixed point 时的数值精度设置

VUMPS 反复更新边界 MPS。更新一个张量时，需要把它左右的无限网络收缩成有限张量，这些收缩结果称为 **environments**。`--tol` 控制这一步的环境求解精度，并决定 VUMPS 的停止精度。

保存的 `boundary.json` 中 `residual` 是 VUMPS 的局部 Galerkin residual。具体地，把左右 environments 和行 MPO 作用到 MPS 中心张量，所得结果归一化为 Y；将左正交 MPS 张量写成满足 $M_L^\dagger M_L=I$ 的矩阵，则程序测量

$$
R_{\rm boundary}=\|(I-M_LM_L^\dagger)Y\|_F.
$$

它表示当前局部固定点方程在该 MPS 子空间之外还剩多少分量。在给定χ下，这个量越小，说明局部迭代越接近自洽。当前 `--tol=1e-6`；代码把 VUMPS 内部停止目标设为它的十分之一，随后重建 environments 并报告上述 residual。

### `--spectrum-tol`：构造有效转移作用时的环境求解精度

计算 $T_{\rm eff}(k)B$ 时，同样需要收缩 B 两侧的无限网络，并求和 B 位于不同位置时的贡献。无限尾部可转化为固定点问题或线性方程，例如形式上的 $G=b+e^{ik}\mathcal E G$；这里 $\mathcal E$ 是对应的混合转移作用，b 是含有一次 B 插入的贡献，G 代表其无限求和。

`--spectrum-tol` 控制这些环境和线性方程的迭代精度，当前取 `1e-7`。它不是“最终谱量 ε 的误差保证为10⁻⁷”，也不是完整矩阵对角化的独立精度参数。

### `projected_residual`：每个谱解的本征方程残差

对求出的本征值 η 和归一化的 B 参数向量 b，程序重新计算一次有效作用，记录

$$
R=\|T_{\rm eff}(k)b-\eta b\|,\qquad\|b\|=1.
$$

这是 CSV 中的 `projected_residual`。`summary.json` 的 `max_projected_residual` 是所有 k、g、q、mode 中这个 R 的最大值，也就是报告所说的**本征方程残差的最大值**。它衡量有效变分问题的求解情况；增大χ改变变分空间所产生的谱差异，需要另行比较。

### 实际运行时间和内存

`wall_seconds` 是从程序开始到结束的实际经过时间，包含首次编译、数值计算和输出，报告统一称“实际运行时间”。`peak_rss_bytes` 是进程使用的物理内存峰值，报告换算成 GiB；1 GiB为 $2^{30}$ bytes。

## 5. 阅读代码、运行和保留的数据

现在可以按以下顺序阅读，符号与上面的推导对应：

1. **`tensor.jl`**：直接建立 K，变换虚拟基底，再对四个物理自旋过滤，返回 A。
2. **`run.jl`**：由 A 构造行 MPO；VUMPS 求 $r_0$，群作用构造 $r_1$；遍历 k、g、q，调用 `solve_modes` 并保存结果。
3. **`plot.py`**：读取 CSV，用颜色和符号画四个 sector。

共享数值实现在 [transfer-spectrum](../../transfer-spectrum/README.md)。`network.jl` 构造 Q 和行 MPO，`boundary.jl` 的 `boundary_solve` 求 fixed point，`spectrum.jl` 的 `effective_action` 和 `_solve_modes` 构造并对角化有效矩阵。

从仓库根运行，选择不存在的新输出目录：

```bash
# 较快的χ=8设置
timeout --kill-after=10s 1800s env OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --heap-size-hint=4G --project=julia-env benchmarks/filtered-toric-code/run.jl --out=benchmarks/filtered-toric-code/results/new-chi8 --chi=8 --nk=41 --num=24 --tol=1e-6 --spectrum-tol=1e-7 --budget=1780
python3 benchmarks/filtered-toric-code/plot.py --run=benchmarks/filtered-toric-code/results/new-chi8 --out=benchmarks/filtered-toric-code/results/new-chi8/figure7a.svg

# 最终主图的χ=24设置
timeout --kill-after=10s 7200s env OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --heap-size-hint=6G --project=julia-env benchmarks/filtered-toric-code/run.jl --out=benchmarks/filtered-toric-code/results/new-chi24 --chi=24 --nk=61 --num=96 --tol=1e-6 --spectrum-tol=1e-7 --budget=7180
python3 benchmarks/filtered-toric-code/plot.py --run=benchmarks/filtered-toric-code/results/new-chi24 --out=benchmarks/filtered-toric-code/results/new-chi24/figure7a.svg
```

`--nk` 是动量采样数，`--num` 是每个 k、g、q 下输出的 modes 数。`--budget` 以秒设置数值计算时间上限，外层 `timeout` 对整个进程限制时间。每完成一个 k、g、q，程序立即将新结果写入 CSV 并刷新到文件。

`results/` 只保留当前有效的三组计算：

- [`fig7a-chi8/`](results/fig7a-chi8/)：快速参考，41动量、每sector24模。
- [`fig7a-chi16/`](results/fig7a-chi16/)：中间χ的比较，41动量、每sector64模。
- [`fig7a-chi24/`](results/fig7a-chi24/)：最终主图，61动量、每sector96模。

每组保存输入张量、边界 MPS、参数、谱、运行日志和当次代码。图显示 $0\le\varepsilon\le6$，CSV保留该范围以外的 modes；图中额外的零谱量点是归一化 fixed point，不占用 `num` 个边界激发解。
