# One-site/C₄/Z₂ SU：零场随机优化结果

## 研究对象与结论

这里研究无限平方晶格 Toric Code 的零场基态，$H=-\sum_sA_s-\sum_pB_p$，$J_e=J_m=1$，X、Z 采用 Pauli convention，物理自旋位于 edge。全网络共享一组顶点张量 V、含边自旋的 P 和对角连接权重 λ，保持 one-site 平移、V 的 C₄ 旋转、P 的端点交换及虚拟 Z₂。物理维数与 bond dimension 均为2，虚拟偶、奇 charge 各一维。

算法交替执行 star HOSVD 和完整 plaquette 环的逐张量 ALS，区域外使用独立对角 λ 近似环境；具体收缩与截断见 [README](README.md#张量表示与局部更新)。正式初态只从允许 charge 块内随机生成，没有使用解析态或对失败态添加扰动。

**seed=5、6 两条轨迹各完成700轮，均未得到零场基态。** 两者每边自旋能量均为 $e=-0.5$，但分别满足不同的一类 stabilizer；严格目标为 $e=-1$ 且全部 $\langle A_s\rangle=\langle B_p\rangle=1$。末态均保持 `accepted=false`，没有用作正式 SU 转移谱输入。当前不追加随机 SU 寻优，整体后续安排见 [PLAN.md](../../PLAN.md)。

## 参数与最终物理量

两条轨迹均取 $dt=0.05$，每次 plaquette 更新做12轮 ALS。CTMRG 只用于测量，不参与局部张量优化；χ 是其环境维数。

| 轨迹 | 累计 sweep | 每边自旋能量 | 平均 ⟨A⟩ | 平均 ⟨B⟩ | 开放局部残差 $(r_A,r_B)$ |
|---|---:|---:|---:|---:|---|
| seed=5 | 700 | −0.5 | ≈1 | ≈0 | $(2.90\times10^{-15},\sqrt2)$ |
| seed=6 | 700 | −0.5 | $5.31\times10^{-85}$ | ≈1 | $(\sqrt2,9.72\times10^{-10})$ |

能量按 $e=-(\overline{\langle A\rangle}+\overline{\langle B\rangle})/2$ 计算。开放局部残差使用区域内部 λ、去掉外接环境 λ 后的 star 张量 $\Theta_A$ 与 plaquette 张量 $\Theta_B$：

$$
r_A=\frac{\|A_s\Theta_A-\Theta_A\|_F}{\|\Theta_A\|_F},\qquad
r_B=\frac{\|B_p\Theta_B-\Theta_B\|_F}{\|\Theta_B\|_F}.
$$

这两个量检查保留开放虚拟指标的区域张量是否满足 stabilizer 恒等式，不是 CTMRG 收缩的期望值，也不是能量误差界。

| 轨迹 | χ=8 的 CTMRG 残差 | χ=16 的 CTMRG 残差 | 环境情况 |
|---|---:|---:|---|
| seed=5 | $7.13\times10^{-11}$ | 0.04633 | χ=8 收敛，χ=16 未收敛 |
| seed=6 | $3.33\times10^{-16}$ | $4.44\times10^{-16}$ | 两者均收敛，能量一致 |

seed=5 的两种χ都给出约−0.5，但χ=16的读数不能当作已收敛结果。seed=6 两种环境均收敛，仍有明显的 star 违例，因此它的失败不能归因于 CTMRG 未收敛。

当前零场验收要求两种χ的环境均收敛且残差不超过 $10^{-9}$，各自 $|e+1|\le10^{-6}$、所有 stabilizer 误差不超过 $10^{-6}$，两种χ能量差不超过 $10^{-6}$；另要求开放局部残差不超过 $10^{-6}$，张量空间/虚拟对称残差不超过 $10^{-12}$，并且程序正常完成。两条轨迹都明确失败，不只是未达到较严格的环境标准。

## 为什么停在两个错误态

### seed=5：X 极化与逐张量 ALS 驻点

P 的两个虚拟 charge 对应的物理列都趋向 $|+\rangle$，第700轮归一化列重叠为1，$|-\rangle$ 振幅约为 $10^{-15}$ 或零。物理态因而趋近全体边自旋的 X 极化乘积态，解释了 $\langle A\rangle=1$、$\langle B\rangle=0$。

对该极限，plaquette 门的目标为

$$
e^{dt B_p}|++++\rangle
=\cosh(dt)|++++\rangle+\sinh(dt)|----\rangle.
$$

固定另外三个物理因子、只改变一个 P 的 ALS 方向，与第二项正交，无法产生四个物理因子同时翻转的分量。逐张量最小二乘因此存在驻点。归一化目标中遗漏分量的平方范数为

$$
w=\frac{\sinh^2(dt)}{\cosh^2(dt)+\sinh^2(dt)}
=0.002489625523386754\quad(dt=0.05).
$$

实际 plaquette 拟合平方残差与此吻合。轨迹约从第40轮起 $r_B$ 已接近 $\sqrt2$，而 $r_A$ 继续下降；第700轮仍停在相同物理状态。

### seed=6：Z-basis copy 极限与 star 截断

P 的两个 charge 列趋近物理 Z 的两个本征态，其极限形式可写为 $P^p_{ab}\propto\delta_{ab}\delta_{pa}$。此时 V 的偶 parity 条件变成额外的物理约束

$$
C_s=\prod_{i\in\operatorname{star}(s)}Z_i=+1.
$$

相邻顶点的 $A_{s'}$ 与 $C_s$ 共享一条边并反对易。因此在严格满足这个附加约束的态中，$\langle A_{s'}\rangle=0$。这来自本条轨迹中 P 的极限形式，并不是一般允许的 Z₂ 张量都被固定为这种形式。

对严格 Z-basis copy P，star 门生成原物理构型和四自旋翻转构型。每个 charge 块仅保留一个奇异向量时，较小的翻转分量被 HOSVD 截去；相对丢弃权重同样为上式的 $w$，与保存的 star 更新记录一致。seed=6 最终满足 plaquette，却无法消除 star 违例。

两种机制解释了已观察到的有限精度轨迹，不是关于任意初态或无限迭代次数的失败定理，也不能据此证明受约束的表示容纳不了 Toric Code 基态。解析参考态通过完整 sweep 的固定点检查；实际问题在于这两条随机轨迹没有找到它。

## 完整逐轮数据与复现

每条轨迹的数据分为三段：

| 数据位置 | 累计轮数区间 | 该段末态 | 停止情况 |
|---|---|---|---|
| `results/seed5/segment1/` | 1–300 | `final.json` | 更新完成，物理验收失败 |
| `results/seed5/segment2/` | 301–600 | `final.json` | 更新完成，物理验收失败 |
| `results/seed5/segment3/` | 601–700 | `latest.json` | 更新完成，随后测量超时；保留 `recovery.json` |
| `results/seed6/segment1/` | 1–292 | `latest.json` | 预算结束 |
| `results/seed6/segment2/` | 293–500 | `final.json` | 正常结束，物理验收失败 |
| `results/seed6/segment3/` | 501–700 | `final.json` | 正常结束，物理验收失败 |

每段的 `trajectory.jsonl` 保留全部 update、residual 和 measurement 记录；update 包含 star/plaquette 的拟合或丢弃误差、纯 gauge 残差、共享恢复变化及空间/虚拟对称残差。不是每轮都做 CTMRG 测量，逐轮更新记录不能理解为逐轮全环境能量。

`settings.json` 保存参数、起始累计轮数、父状态 SHA256、Julia/包版本及当时源码哈希；`initial.json` 保存每段起始张量。相邻段的张量逐分量一致，父状态哈希相符。700轮共有1400个子更新，两条轨迹记录的最大纯 gauge 残差分别约 $4.21\times10^{-16}$、$7.88\times10^{-16}$，最大共享恢复变化分别约0.001856、0.006684；共享恢复不是纯 gauge 变换，不能因为 gauge 误差小就忽略它。

seed=5 的第700轮末态为 `results/seed5/segment3/latest.json`。该态的 χ=8/16 测量见 `results/seed5/measurement/trajectory.jsonl`，完整物理验收结果见 [results/seed5/audit.json](results/seed5/audit.json)，其中 `accepted=false`。测量进程在汇总阶段触及预算，`measurement/summary.json` 因而标为 `budget`；物理结论依据已完成的测量和局部残差。

seed=6 的最终验收直接见 [results/seed6/segment3/summary.json](results/seed6/segment3/summary.json)，末态为同目录 `final.json`。两条失败态均不满足正式谱输入的接受条件。
