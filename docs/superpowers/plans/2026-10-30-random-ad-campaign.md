# Random AD ground-state campaign implementation plan

> 本随机计算已结束且路线已停止，不再执行后续诊断或重算。代码现位于 `ad/one-site-c4-z2/`，各 seed 的合并结果见该目录 REPORT。以下是历史计划和执行记录，不是当前运行说明。

> **For agentic workers:** Use executing-plans to implement the following tasks inline.

**Goal:** Test random-origin D=4 one-site/C₄/virtual-Z₂ AD at hx=0, hz=0.1 within six hours, using a 10⁻⁴ per-spin energy target against the existing low-field series.

**Architecture:** Retain shared parameterization and unrolled CTMRG AD. Make contraction precision configurable; replace mandatory high-precision stationarity auditing with a short energy-based check. Run multiple attempts in one Julia process to amortize compilation, with external wall-time supervision.

**Tech Stack:** Julia repository environment, PEPSKit/Zygote, Python process supervisor.

**Spec:** User-approved conversation: continue seed 19 first, switch seeds according to progress, total budget extended to six hours; no spectrum and no extensive validation.

## Global constraints

- Preserve all historical results and unrelated working-tree changes; no package changes or cloud/cluster access.
- Keep one-site, C₄ and virtual Z₂ constraints exact; retain random provenance and parent hashes.
- Energy reference is a truncated series, not an exact ground-state certificate.
- Single Julia process, single BLAS/Julia thread, 2G GC heap hint; save every accepted update and flush progress.
- No automatic commits. No full Julia regression or repeated gradient probes.

## Review focus

- Budget termination must leave a complete checkpoint, never label incomplete work accepted.
- Reject unconverged/nonphysical contractions and inconsistent independent environments.
- Do not choose a contraction branch by lowest energy.
- Use identical accuracy settings across random seeds; keep continuation distinct from fresh starts.
- Passing the energy target must still require fresh-environment and larger-χ consistency.

### Task 1: Accuracy settings and energy-based completion

Files: `objective.jl`, `optimize.jl`, `runner.jl`, `tests.jl` in `ad/one-site-z2/`.

- [x] Expose CTMRG tolerance (default 10⁻⁸), maxiter (300), and energy target (10⁻⁴); retain 20-step AD and existing gradient evidence without rerunning it.
- [x] Use exactly two independent initializations per evaluation, retaining dominant-norm selection, not minimum-energy selection.
- [x] Add an optimizer energy-target stop and keep stationarity as a diagnostic, not a success prerequisite.
- [x] Final check: fresh χ=8 pair plus χ=12 pair, tolerance 10⁻⁸, energies within 10⁻⁵ of each other and 10⁻⁴ of the series, physical bounds and all local B within 10⁻³ of +1. Save all check rows. No gradient audit.
- [x] Run only targeted lightweight optimizer/control tests and syntax checks.

### Task 2: Six-hour adaptive campaign

Files: new `ad/one-site-z2/campaign.jl`; existing `supervise.py`.

- [x] Run a single persistent Julia process, continuing seed19's last complete checkpoint for up to one hour first.
- [x] Allocate later attempts 45 minutes each, interleaving fresh seeds with continuation of the best progressing trajectory. Stop on energy-based acceptance or total deadline. Limit repeated continuation without progress.
- [x] Persist settings, per-attempt summaries, latest complete-state paths and progress after every attempt; outer Python supervisor enforces six-hour wall budget including compilation.
- [x] Launch locally in background after targeted checks; report actual startup and progress, not speculative success.

### Task 3: Clear framework and current documentation

- [x] Remove current one-off high-precision `ctm_trace.jl`, `audit_state.jl` and its dedicated tests, plus the repeated standalone gradient probe driver; retain their historical source snapshots.
- [x] Keep focused parameterization/optimizer and supervisor tests; retain shared final-state synchronization tests because the energy check still uses that control path.
- [x] Update directory README/REPORT and root README/PLAN to explain the energy-based campaign and its scope; preserve previous numerical conclusions.
- [x] Record campaign running status accurately; inspect changes without full-regression reruns. Numerical outcome is pending.

## Execution ledger

- User explicitly requested minimal testing and immediate execution: implement inline in the existing working directory, whose one-site AD code and many unrelated changes are not yet committed. No worktree migration, automatic commits, full-suite runs, or repeated gradient checks.
- Target-stop test first failed on unsupported keywords, then the focused directory tests passed (94 checks). The supervisor subprocess timeout test passed. Campaign load check returned the series reference and `CAMPAIGN_LOAD_OK`.
- Launch: external supervisor background task `sh_1895805f9001KEWHE7ioGPop8d`, six-hour hard budget, output `ad/one-site-z2/results/random-campaign-20261030/`.
- Original diagnostic scripts were verified present in the campaign source snapshot before deletion from current code. No historical results were edited.
- Independent focused read-only review found one important reporting issue: outer supervisor looked for checkpoints only at campaign root. Added a real nested-child timeout regression (observed failure then pass) and recovery of per-attempt checkpoint paths, completed attempts and best trajectory in future supervisor invocations. The already-running Python supervisor loaded the previous version; at a hard timeout, its root summary may still have a null last_checkpoint, while child states and progress.json remain authoritative and intact. Do not restart a compiling numerical job merely to improve timeout summary formatting.
- Confirmed live progress.json names attempt01, seed19, continuation=true; log records the 3600-second first allocation. Process is alive with 2G GC hint. No first new gradient/energy result at that observation; do not infer convergence.
- First campaign completed normally after 16676.73 supervised seconds, because six fresh seeds and the limited continuations were exhausted, not because of the six-hour cap. No accepted state; best seed23 e=-0.7669864963106053, 28 accepted updates across three segments. Ruling: spend remaining authorized compute on best seed23 with unchanged settings rather than stop at the arbitrary list cap. A separate supervised continuation has budget4900 seconds (including startup/compilation), total supervised allowance 21576.73 seconds. Background task `sh_18a56de04001ZTklMME2pzvHVZ`, output `ad/one-site-z2/results/random-campaign-20261030-best23-tail/`. No automatic further compute beyond this allowance.
- Tail completed normally after4223.62 supervised seconds, status=line_search_failed, two additional updates, final e=-0.7670721872726077, gradient norm0.21376308586145004. Last16 trials rejected by contraction checks, not by budget or a final high-precision audit. Total actual supervised compute20900.35s (5h48m20s); accepted=false. Final states and all logs preserved. Documentation synchronized; no more computations authorized automatically by this plan.
