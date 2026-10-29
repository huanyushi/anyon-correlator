# 无限行转移算符的动量谱

这里保存磁场 Toric Code 数值态的转移谱入口，以及独立 benchmark 共用的数值求解器。[Filtered Toric Code benchmark 的说明](../benchmarks/filtered-toric-code/README.md) 从 V—bond tensor—P 网络推导四自旋 PEPS 张量，并详细定义 fixed points、动量和 charge，可作为算法的入门说明。

## 算法主干与阅读顺序

计算流程是 **PEPS 张量 → 双层行 MPO → VUMPS 边界 → 四 sector 准粒子谱**。`numerics.jl` 只加载以下数值部分，不加载 SU、AD 或 Hamiltonian 优化程序。

1. **`network.jl`：构造转移算符。** 输入 $A^p_{nesw}$，物理指标 p 可以包含一个块内的多个自旋。将 A 与其共轭沿 p 收缩得到 $Q=\sum_p A^p\otimes\overline{A^p}$，每条双层腿的维数为 $D^2$。ket/bra 腿融合到按对角 Z₂ charge 分块的空间，Q 作为一位点无限 MPO 的局部张量。程序保留独立 ket、bra 群作用。一个转移步对应一个输入 PEPS 块；晶格方向由输入张量的 blocking 决定。

2. **`boundary.jl`：求主导边界 MPS。** `boundary_solve` 用 VUMPS 求行 MPO 的主导右本征态 $r_0$，χ 是这个无限 MPS 的键维数。拓扑相中，单层虚拟群作用给出另一个固定点 $r_1$。这里求解的是给定 PEPS 的边界，PEPS 张量本身保持固定。`solve_boundaries` 是磁场数值态使用的完整边界检查入口；简单 benchmark 直接调用 `boundary_solve`。

3. **`sectors.jl`：选择背景与 charge。** g=0 表示扰动左右均为 $r_0$，g=1 表示左右分别为 $r_0,r_1$，即一个 domain wall。q 是 ket、bra 同时进行虚拟 Z₂ 变换时的 charge：q=0/1 对应激发态的变换本征值 +1/−1。拓扑相中的 $(g,q)=(0,0),(0,1),(1,0),(1,1)$ 分别对应 $1,e,m,\epsilon$。

4. **`spectrum.jl`：在边界扰动空间中求谱。** 用局部张量 B 替换边界 MPS 的一个张量，并对位置 j 叠加相位 $e^{ikj}$。左右背景及 q 确定变分空间；左正交规范消去 B 的表示冗余。MPSKit 提供归一化有效转移作用，`solve_modes` 逐列构造其在 B 参数基中的矩阵并对角化，再按本征值模长取所需 modes。

输出谱量为 $\varepsilon(k)=-\log|\lambda(k)/\lambda_0|$。这里的有限矩阵是无限 MPS 的准粒子变分矩阵，不是有限圆柱的转移矩阵。给定 χ 下的投影残差反映求解精度；增加 χ 或扩大扰动范围会改变变分空间，因此应按所需精度判断谱的稳定性。当前完整矩阵算法在同一χ下增加输出 modes 的主要开销较小，增加χ则会明显扩大矩阵。

这里 k 是沿边界块晶格平移的 Bloch 动量，周期为 $2\pi$。domain-wall 背景具有相对逐格相位约定：$r_1$ 的每格张量变号，同时令 $k\mapsto k+\pi$，表示同一个 domain-wall 态。与文献比较时需要写清这个约定及转移晶格的步长。

## 磁场数值态的输入与运行

`core.jl` 在数值主干上加载项目适配层：`input.jl` 将 V/P/λ 张量精确合并为东、北两个边自旋的轴向块；`output.jl` 读取已保存 SU/AD 态并处理结果输出。它们服务于磁场 Hamiltonian 的数值态研究。

在 V/P 表示中，合并张量为

$$
A^{p_Ep_N}_{nesw}=\sum_{a,b}V_{basw}P^{p_E}_{ae}P^{p_N}_{bn}
\lambda_a\lambda_e\lambda_b\lambda_n,
$$

物理融合次序为 $p=p_E+2p_N$，每条 V–P 连接的 λ 恰好计入一次。该表示沿原顶点晶格轴向转移，不能将谱值直接当作对角四自旋块的谱。

从仓库根运行，输出目录须不存在：

```bash
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env transfer-spectrum/run.jl --state=已验收final.json --out=新谱目录 --chis=8 --nk=32 --num=8 --budget=3600
python3 transfer-spectrum/plot.py --input=新谱目录/spectrum.csv --out=新图目录
```

边界 MPS 的χ与基态优化中 CTMRG 的环境维数属于不同计算步骤。解析过滤态的运行命令和χ比较见 [benchmark README](../benchmarks/filtered-toric-code/README.md) 和 [报告](../benchmarks/filtered-toric-code/REPORT.md)。

已有测试入口为 `tests.jl`，可用 `--group=network`、`--group=boundary` 或 `--group=spectrum` 选择相应部分。
