# 无限行转移算符的动量谱

这里保存磁场 Toric Code 数值态的转移谱入口，以及独立 benchmark 共用的数值求解器。[Filtered Toric Code benchmark 的说明](../benchmarks/filtered-toric-code/README.md) 从 V—bond tensor—P 网络推导四自旋 PEPS 张量，并详细定义 fixed points、动量和 charge，可作为算法的入门说明。

## 算法主干与阅读顺序

计算流程是 **PEPS 张量 → 双层行 MPO → VUMPS 边界 → 四 sector 准粒子谱**。`numerics.jl` 只加载以下数值部分，不加载 SU、AD 或 Hamiltonian 优化程序。

1. **`network.jl`：构造转移算符。** 输入 $A^p_{nesw}$，物理指标 p 可以包含一个块内的多个自旋。将 A 与其共轭沿 p 收缩得到 $Q=\sum_p A^p\otimes\overline{A^p}$，每条双层腿的维数为 $D^2$。ket/bra 腿融合到按对角 Z₂ charge 分块的空间，Q 作为 one-site 无限 MPO 的局部张量。程序保留独立 ket、bra 群作用。一个转移步对应一个输入 PEPS 块；晶格方向由输入张量的 blocking 决定。

2. **`boundary.jl`：求主导边界 MPS。** `boundary_solve` 用 VUMPS 求行 MPO 的主导右本征态 $r_0$，χ 是这个无限 MPS 的键维数。拓扑相中，单层虚拟群作用给出另一个固定点 $r_1$。这里求解的是给定 PEPS 的边界，PEPS 张量本身保持固定。`solve_boundaries` 是磁场数值态使用的完整边界检查入口；简单 benchmark 直接调用 `boundary_solve`。

3. **`spectrum.jl`：选择背景与 charge。** g=0 表示扰动左右均为 $r_0$，g=1 表示左右分别为 $r_0,r_1$，即一个 domain wall。q 是 ket、bra 同时进行虚拟 Z₂ 变换时的 charge：q=0/1 对应激发态的变换本征值 +1/−1。拓扑相中的 $(g,q)=(0,0),(0,1),(1,0),(1,1)$ 分别对应 $1,e,m,\epsilon$。

4. **`spectrum.jl`：在边界扰动空间中求谱。** 用局部张量 B 替换边界 MPS 的一个张量，并对位置 j 叠加相位 $e^{ikj}$。左右背景及 q 确定变分空间；左正交规范消去 B 的表示冗余。MPSKit 提供归一化有效转移作用，`solve_modes` 逐列构造其在 B 参数基中的矩阵并对角化，再按本征值模长取所需 modes。

输出谱量为 $\varepsilon(k)=-\log|\lambda(k)/\lambda_0|$。这里的有限矩阵是无限 MPS 的准粒子变分矩阵，不是有限圆柱的转移矩阵。给定 χ 下的投影残差反映求解精度；增加 χ 或扩大扰动范围会改变变分空间，因此应按所需精度判断谱的稳定性。当前完整矩阵算法在同一χ下增加输出 modes 的主要开销较小，增加χ则会明显扩大矩阵。

这里 k 是沿边界块晶格平移的 Bloch 动量，周期为 $2\pi$。domain-wall 背景具有相对逐格相位约定：$r_1$ 的每格张量变号，同时令 $k\mapsto k+\pi$，表示同一个 domain-wall 态。与文献比较时需要写清这个约定及转移晶格的步长。

## 磁场数值态的输入与运行

`core.jl` 在数值主干上加载项目适配层：`input.jl` 将 V/P/λ 张量精确合并为东、北两个边自旋的轴向块；`output.jl` 读取已保存 SU/AD 态并处理结果输出。`run.jl` 是唯一计算入口，先求完整动量网格，再按剩余预算比较另一χ，不另设验证驱动。

在 V/P 表示中，合并张量为

$$
A^{p_Ep_N}_{nesw}=\sum_{a,b}V_{basw}P^{p_E}_{ae}P^{p_N}_{bn}
\lambda_a\lambda_e\lambda_b\lambda_n,
$$

物理融合次序为 $p=p_E+2p_N$，每条 V–P 连接的 λ 恰好计入一次。该表示沿原顶点晶格轴向转移，不能将谱值直接当作对角四自旋块的谱。

从仓库根运行，输出目录须不存在：

```bash
timeout --signal=TERM --kill-after=10s 7200s env OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --heap-size-hint=2G --project=julia-env transfer-spectrum/run.jl --state=ad/one-site-c4-z2/results/d4-exact-noise-cont1-20261029/final.json --candidate=true --out=新谱目录 --chi=8 --compare-chi=12 --nk=64 --num=8 --budget=7080
python3 transfer-spectrum/plot.py --input=新谱目录/spectrum.csv --out=新图目录
```

`--chi` 指定主网格的边界维数；`--compare-chi` 指定可选比较维数，设为0即禁用。主网格包含 $k=-\pi+2\pi j/N_k$、$j=0,\ldots,N_k-1$，不重复周期端点。剩余时间超过300秒时尝试较大χ的 $k=0,\pi/2,-\pi$；比较失败或超时不废弃已保存主谱。外部 `timeout` 限制整个进程，内部预算用于在收缩步骤间主动停止。

`--candidate=true` 明确允许读取历史未验收张量，不改变原态的 `accepted=false`。本轮固定解析来源张量，不做基态优化。边界 MPS 的χ与基态 CTMRG 的χ属于不同计算步骤。

默认边界投影残差容差为 $10^{-6}$，左右固定点一致性容差为 $10^{-5}$；激发环境容差为 $10^{-10}$，有效本征方程残差容差为 $10^{-7}$。这是面向定性谱的设置，不是有限χ误差界。尤其小模长本征值的对数谱量会放大环境误差，不能仅凭残差小就认定高支可靠。

每完成一个sector即追加JSONL，整个动量完成后更新CSV与 `progress.json`，避免后续sector超时时丢弃已算完的sector。`settings.json` 记录参数、输入和源码哈希，`input.json` 保存冻结张量，`boundaries.json` 保存边界张量，`source/` 保存本次谱程序。源码副本用于核对，独立复现仍需本仓库及根 `julia-env/Manifest.toml`。`summary.json` 区分主网格完成情况和可选比较状态。

## 当前磁场结果与已有检查

D=4冻结张量的χ=8、64动量、2048模主谱及χ=12三个共同动量比较均已完成，实际约30.4分钟，结果和图见 [REPORT.md](REPORT.md)。e、m、ε主导模在已比较动量上稳定，但平凡sector及高支仍有χ依赖。输入来自[one-site/C₄/Z₂ AD](../ad/one-site-c4-z2/README.md)的受扰动解析初态优化，不来自已停止的随机路线。

实现检查集中在 `tests.jl`（张量合并、虚拟对称性、sector、输入记录）与 `plot_tests.py`（已保存数据的绘图）；不作为每次正式计算的前置长任务。运行方式分别为 `julia --project=julia-env transfer-spectrum/tests.jl` 和 `python3 -m unittest discover -s transfer-spectrum -p plot_tests.py`。物理可靠性主要通过所关心谱支的χ比较判断，而不是增加测试文件数量。
