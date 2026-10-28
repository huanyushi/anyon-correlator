# 独立仓库迁移记录

## 已确认方案

本项目从 `/home/huanyushi/code/quantum.harness/tracks/peps/solutions/anyon-correlator/` 完整复制到 `/home/huanyushi/code/anyon-correlator`，以当前状态建立首次 Git 提交。包含此前尚未提交的程序、文档，以及全部原始 `results/` 数据和最终张量。原仓库保留作历史档案。

Julia 的 Project.toml 和 Manifest.toml 从原仓库 `julia-env/` 复制，保留原版本；复用当前 WSL 用户的软件包缓存。独立仓库不需要安装 quantum.harness 的技能集合或 Ion。

## 执行与检查

- [x] 完整复制课题目录和 Julia 环境配置。
- [x] 更新各级 README 的运行路径和项目约定。
- [x] 逐文件核对原始数据、程序与环境配置；检查 Julia 环境可用性。
- [x] 全部研究文件加入首次提交的暂存区，并逐文件检查；共 1354 个文件，仅排除两个 Python 字节码缓存。

历史运行目录中的源代码副本、参数、日志与哈希保持原样。迁移只改变当前文档与仓库配置，不重新计算基态、运行磁场扫描或更新依赖。

完整性检查覆盖原目录的 1352 个文件，其中 1149 个 `results/` 文件内容完全一致；12 个 README/AGENTS 文档更新了路径或协作说明，5 个子目录 `.gitignore` 取消了对研究数据、日志和图片的忽略，其余原文件逐字节一致。两个 Julia 环境文件也与原仓库逐字节一致。

在新仓库中使用 `julia --startup-file=no --project=julia-env` 成功加载 TensorKit、PEPSKit、MPSKit、KrylovKit、JSON、JLD2 和 Zygote，运行版本为 Julia 1.12.6。本次检查验证环境加载和迁移完整性，没有重跑数值测试或研究计算。

历史报告和旧计划保留原文。复用其中的命令时，删除 `tracks/peps/solutions/anyon-correlator/` 路径前缀；旧文档中的 `.venv/bin/python` 指原仓库环境，需根据具体绘图脚本的依赖另行选择解释器。新仓库的全部原始研究数据纳入版本管理，Python 缓存等可再生文件由 `.gitignore` 忽略。
