# One-site Z₂ AD Implementation Plan

> 历史实施计划，不再作为待办执行。当前目录已改为 `ad/one-site-c4-z2/`，解析初态为后续主线，随机寻优已停止。当前运行方式与汇总见该目录 README 和 REPORT；下文旧路径和步骤仅保留实施背景。

> **For agentic workers:** REQUIRED SUB-SKILL: Use `executing-plans` to implement this plan task-by-task.

**Goal:** 独立实现 one-site/C₄/虚拟 Z₂ 的有限场能量 AD，先验证参数化和梯度，再运行双起点对照。

**Architecture:** `ad/one-site-z2/` 保存模块、参数化、能量适配、优化运行器和独立结果。只读复用 `ad/split/`，通过 10→40 参数映射及其 pullback 连接现有能量接口；不调用 SU 更新。

**Tech Stack:** 原 `julia-env/` 的 Julia、TensorKit、PEPSKit、Zygote、JSON。

**Spec:** `docs/2026-10-29-one-site-z2-ad-design.md`。

## 修订与执行状态（用户已批准增大 D 和延长计算）

Ruling: 参数化改为一般 D，偶、奇表示均允许重数；D=2 作校准，D=4（各二维）作首个较大表示。下文原10→40接口相应推广为 `ParameterBasis` 显式参数：D=4 为54→336。若较大表示仍不足，需要进一步做 D 收敛研究。

- [x] 已建立 `ad/one-site-z2` 分支，保留全部原有改动。
- [x] D=2/4 参数、对称性及映射差分测试26项通过。
- [x] Task 1：状态读写、D=2/4闭合网络规范等价性及全部零场稳定子验证。
- [x] Task 2：D=2/4共享梯度拉回、真实 CTMRG 方向导数及20/40次稳定性通过。
- [ ] Task 3：预算运行器、证明设置绑定和继续计算来源检查已实现；双起点计算进行中。
- [ ] Task 4：集成、文档与独立审查已完成，待双起点计算最终结果回填。

2026-10-29 验证记录：新模块71项轻量检查通过，旧 split AD 两组测试75项通过。600秒全仓回归已通过新增测试所在组，在旧 `FailureChecks` 编译期间超时，强制终止退出137；不能标记全仓通过。D=2真实CTMRG六个方向差分最大误差 $5.2732\times10^{-9}$，20/40次梯度差 $6.9410\times10^{-11}$，通过。首次梯度约509.6秒，后续40次约4.0秒。D=4前向收缩137.19秒、残差 $2.05\times10^{-11}$。原始日志保存在 `ad/one-site-z2/results/software-and-d2-probes-20261029/`。

Ruling: 用户已允许延长计算，真实梯度诊断单段预算从600秒增至1800秒；保持独立超时监督与逐阶段输出。成本是单次失败诊断最长占用30分钟，不改变物理验收标准。

Final review: 一次独立只读审查完成。三项 Important 均已修复：Python 独立监督在评价/验收期间终止真实子进程并保存预算结论；最终验收对状态与汇总统一回写；预期非有限导数和无效试探权重转换为可保存的数值失败。各项均 RED→GREEN，新模块90项 Julia 检查与 Python 超时测试通过。

Final: minor (deferred): D=2相邻/相对C₄轨道的逐轨道隔离系数测试；当前已有严格C₄/charge检查、参数计数及闭合网络验证。

Final: minor (deferred): `latest.json` 和独立 settings 文件；当前保留编号 checkpoint、初态来源、父态哈希，以及监督进程的完整命令和预算。

Final: Ruling: D4数值有效性由真实梯度诊断确认；两种起点的收敛性由正在运行的独立轨迹判定；转移谱兼容性仍待合格态出现后单独验证。提前将这些未完成项解释为成功会产生错误物理结论，因此报告明确保留待定。

D=4诊断：六个方向差分最大误差 $3.2133\times10^{-10}$，梯度差 $1.6527\times10^{-12}$，首次求导532.94秒、40次后续求导21.50秒。剩余旧 `FailureChecks` 和状态加载检查单独运行116项通过，保留600秒全仓单次超时记录。

## Global Constraints

- 一般 D，默认偶奇尽量均分、物理两个平凡表示；C₄ 仅循环置换。
- 初始有限场 (0,0.1)，一组 V/P/λ，全程满足约束。
- 保留现有未提交改动、旧源码、历史结果与 Julia 环境；不提交或推送。
- 生产计算单段含编译≤600秒，RSS预算<16GiB，逐步保存；先完成软件验证和成本测量。
- 使用独立 AD provenance，不接受旧 SU 成功标志作为 AD 验收。

## Review Focus

1. C₄ 相邻/相对轨道不可混合；任务1逐分量验证。
2. λ 平衡规范与谱输出规范不同；任务1独立波函数检查。
3. 共享参数的所有复制位置均进入梯度；任务2有限差分检查。
4. 梯度依赖环境和有限 unroll 图；任务2固定参考环境、任务3拒绝失效梯度。
5. 超时、小梯度、解析起点与随机起点不能混淆；任务3输出/验收回归。

所有 Julia 命令前缀：`OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --startup-file=no -O1 --project=julia-env`。

### Task 1: 共享对称参数化与独立状态格式

**Files:** `ad/one-site-z2/{core.jl,parameters.jl,state_io.jl,tests.jl,README.md}`。

**Interfaces:** `initial_parameters(seed;exact=false)::Vector{Float64}` 长度10；`tensors(theta)` 返回 V/P/λ；`split_parameters(theta)::Vector` 长度40；`save_state(path,theta;metadata)`、`load_state(path)`。

- [ ] 先写 `--group=parameters`：四个 V 轨道逐分量、P forbidden 分量严格零、非均匀 λ、随机可复现、零参数拒绝、JSON往返与拒绝覆盖；新接口缺失时确认失败。
- [ ] 实现纯函数参数构造和10→40映射；λ的整体尺度只影响未归一化波函数。
- [ ] 独立小闭合网络核对平衡/谱规范的归一化波函数，误差≤1e-12；解析参考检查全部稳定子。
- [ ] 运行 `ad/one-site-z2/tests.jl --group=parameters`，预期全部通过。

### Task 2: 受约束能量和梯度适配

**Files:** `ad/one-site-z2/{objective.jl,gradient_check.jl,tests.jl}`。

**Interfaces:** `contract(theta,H;kwargs...)` 包装原AD的D=2收缩；`energy_gradient(theta,H,reference;unroll_steps=20)` 返回能量与10维梯度；`check_gradient(theta,hx,hz;options)` 返回设置、方向导数、误差和通过标志。

- [ ] 添加参数映射标量目标的AD/差分测试，观察缺失接口失败；方向差分 atol=1e-9、rtol=1e-6。
- [ ] 实现链式pullback，复用原AD的40维能量导数，不重复实现CTMRG求导。
- [ ] 运行轻量导数组；再以有界独立命令测量真实梯度成本，并在相同参考环境下对非零方向导数比较步长1e-4、1e-5、1e-6，要求存在稳定区间且误差≤max(1e-7,1e-3|方向导数|)。
- [ ] 比较 N=20/40 梯度差≤2e-6；不满足时记录失败原因，不授权正式优化。

### Task 3: 独立有预算的AD运行与验收

**Files:** `ad/one-site-z2/{optimize.jl,run.jl,tests.jl}`。

**Interfaces:** runner接受 `--init=exact|random --seed --hx --hz --steps --budget --resume --out`；输出settings、initial/latest/final JSON、trajectory JSONL、summary。

- [ ] 先测试预算耗尽、过期截止时间、错误梯度设置、环境失败、不同场点resume、已存在路径拒绝；确认失败后实现。
- [ ] 复用L-BFGS方向，线搜索每次真正评价共享受约束态；独立硬预算监督，保存最后完整步骤。不要照搬旧optimizer将budget改为stationary的逻辑。
- [ ] 验收按设计中的物理/梯度/环境标准；解析来源与随机来源分别输出，AD格式不得混作SU格式。
- [ ] 运行 `ad/one-site-z2/tests.jl --group=runner`，预期全部通过；确认成本后在新目录进行正式双起点计算，失败原样保存。

### Task 4: 集成、说明与结果边界

**Files:** `tests/runtests.jl`、`ad/README.md`、根`README.md`、`PLAN.md`、`M3_REPORT.md`、新目录`README.md`与`REPORT.md`。

- [ ] 总入口隔离命名空间，仅接入轻量参数/导数/运行器检查。
- [ ] 新模块全组及受影响旧AD测试通过，再运行有预算总回归；超时如实记录，不混同生产失败。
- [ ] 同步算法说明、实际可运行命令和验证记录；新有限场结果取得前不写成已成功。
- [ ] 一次独立最终审查并核对`git diff --check`。转移谱生产接入须在AD态通过验收后单独验证，保留尚未完成项。

## 执行记录

- [x] 已读取当前AD网络、目标函数、优化器、梯度测试及AGENTS；旧优化器不能直接将10维变量传入，必须通过映射和提供者适配。
- [x] 设计与任务接口自查：10维参数、40维复用接口、平衡规范、独立AD来源与预算失败均已指定。
- [ ] 计划审阅与工作分支确认后开始Task 1；当前main存在上一轮未提交研究数据，不能静默迁移或覆盖。
