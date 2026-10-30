# Automatic differentiation（AD）算法

AD 通过能量对张量参数的导数优化基态，与施加虚时门并截断的 simple update 分开。七月和十月的实现都归在本目录，按张量表示区分，不再以月份建立平行项目。

| 目录 | 参数化与方法 | 完成情况 |
|---|---|---|
| [merged/](merged/README.md) | 一个 dense complex D=2 merged 张量共享在四个位置，固定点求导、投影梯度、Armijo descent | 七月零场成功；作为已完成对照 |
| [split/](split/README.md) | 实数 V/P/λ，位置独立或平移共享；fixed-point / unrolled AD | 十月三个弱场点从解析零场或近基态起点成功；随机全局寻优未通过 |
| [one-site-c4-z2/](one-site-c4-z2/README.md) | 一组实数 V/P/λ，保持 one-site、C₄、虚拟 Z₂，一般 D | 受扰动解析初态得到 e≈−1.00251587 的候选态，作为后续主线；任意随机初态寻优未成功，已停止 |

这些实现的张量表示与对称约束不同：一般 merged 张量允许的自由度比固定内部 D 的 split 网络更宽，相同 D 下的可变参数也不同。项目总体仍保留 loop/cluster update 与 AD 比较；对 one-site/C₄/Z₂ 约束 AD，后续重点是从零场解析态出发，而不是要求任意随机态寻优成功。解析和随机路线保留独立入口，共享参数化、收缩和优化实现。

有限场的具体结果见 [M3](../M3_REPORT.md)。
