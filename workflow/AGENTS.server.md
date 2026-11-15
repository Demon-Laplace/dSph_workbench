# Server data workspace

All server mutations must remain inside this data directory. Existing solver
installations outside data may be read and executed, but never edited.

Keep server production models in model/ and isolation models separately in
isoFornax/model/. Do not consolidate server isolation into production folders.
New automated isolation runs use isoFornax/model/runs/run/. Only the LOCAL
workspace consolidates ICs and isolation under each production model directory.
Never move active Slurm work or overwrite snapshots, particle ICs, model images,
videos, or frozen run inputs. Existing frozen validation cases keep their paths.
Equal model numbers alone do not prove matching IC provenance.

GitHub main is the source of truth for code. Root workflow*.py files correspond
to workflow/ in the Git checkout; dSph_workbench/ corresponds to repo/. After a
validated code change, synchronize that exact change to GitHub main. Use local
authenticated synchronization when server credentials are unavailable; never
copy credentials or force-push. Stop on divergent history and report it.
Machine-specific workflow.json, particle data, logs, and private configuration
must not be committed. IC_MODELS.csv and MODEL_SOURCES.csv are small summaries.

Respect the more specific dSph_workbench/AGENTS.md. Immutable workspace/releases/
are never edited in place; publishing analysis requires a new validated release.
Use float32 for inspection/analysis where verified, and preserve simulation
precision. Run focused tests and distinguish short integration checks from
scientific equilibrium or production-scale validation.

Existing user .sh files are protected: do not rewrite, rename, move, replace,
centralize, or delete them unless the user explicitly requests that operation.
Preserve their invocation paths when organizing Python dependencies. Obsolete
standalone temporary-check shell scripts may be archived only within an
explicitly authorized cleanup; model simulation launchers are not such scripts.
The legacy isoFornax/python/ files are byte-identical shared copies tracked in
workflow/isolation_python/. Model .py compatibility links preserve existing .sh
commands. Preserve isoFornax/python/ at this depth: the density script derives
its workspace paths from __file__. Do not mix these with a different analysis
version or edit scientific behavior as part of file organization.

The user explicitly authorized adapting the unnumbered isolation template on
2026-10-13. Its five versioned sub*.sh launchers now call shared Python and
prepare-legacy. This exception does not authorize changing existing numbered
models' shells or frozen run.sh files. Future template changes still require
explicit user authorization. Template Python duplicates are archived; do not
reintroduce per-model copies. TREECOOL is materialized from the configured
template, and new manual inputs must be named after their model directory.

Local model results use a single <model>/isolation/ directory. Do not recreate
isolation/legacy/<model>/ or isolation/<run>/client/ locally. Local transport
receipts belong in workspace/receipts/isolation/<model>/<run>/. Server isolation
paths and frozen cases remain separate and are not flattened by this convention.

The user requested lightweight command-line cleanup on 2026-11-04, superseding
background replacement cleanup. Use isolation-clean-failed MODEL --name FAILED_RUN
from the login node; --name identifies the failed directory itself. Query once,
require terminal FAILED/OUT_OF_MEMORY within 140 seconds, memory/segfault evidence,
verified inputs and no active job using the directory, then delete only that case.
Keep a small deletion receipt. Never apply this policy to successful cases or
long-running interrupted simulations. No background scans, repeated Slurm queries,
or extra scheduler client configurations. Accepted upload staging copies are
removed immediately after receipt verification; keep the simulation's IC.
Old frozen monitors receive automatic_cleanup=disabled from the central cleanup
entry point inside Slurm, ending their retries without changing the solver.
Density checks use the shared isoFornax/python/check_density_profile.py and
central workflow_isolation_ops; preserve existing user shell launchers.
Checks for multiple simulations must share one Slurm job on one node, including
different isolation run names (use repeated --run MODEL=RUN arguments). Process
models sequentially with the configured worker count, normally 16, including
video-only checks. Do not submit one job per model or run concurrent worker
groups. Preserve already submitted tasks rather than cancelling or repeating them.

Production orbital simulations use astro_fat via existing subastrofat.sh,
one separate job per model with the launcher's current resource settings.
Fornax7004/7005 must use snap2csv exports from isolation snapshot 50 and retain
timescale_fd=16.6. Reuse the already verified source exports for orbit-only reruns;
do not regenerate or independently revalidate the same source IC. Keep the
existing one-pass merge workflow. Preserve source isolation snapshots and user
shell scripts.
The user superseded the old Fornax5028 orbit/orientation on 2026-10-20:
current Fornax7004/7005 use joint_GC_radius_target6D_PA43 from the GitHub memo,
not the PA-only candidate. See workflow/REORBIT_7004_7005_20261020.md in the checkout.
The user explicitly authorized permanent deletion of their retired production
outputs, restarts, old analyses and videos. Only local snapshot 282/289 and
elinfo remain in each model's reference/orbit5028-20261017/ archive, with metadata.
Do not reuse old prepared analyses, CSV, PNG or cohort identity with the new
snapshot numbering. This cleanup authorization does not extend to other models.
Production startup monitoring is automatic: exit after 140 running seconds and
nonzero integration progress, leaving the simulation running. Avoid repeated
manual checks or confirmation requests.
After detection succeeds, generate simulation-details.json with IC/model and
Milky Way/background inputs, particle counts and mass resolution, orbit and
orientation, physics parameters and job resources. Return the report and stop
monitoring while leaving the simulation running.
The monitor runs as a sleeping child in the existing batch allocation, checking
every 15 seconds, with no extra Slurm job, MPI task or reserved CPU. Queue time
uses no monitor node resources. Preserve user launchers; freeze a generated
wrapper and monitor sources instead. Local report reception is independent.

New agents should first read workspace/code/github-sync/workflow/AGENT_HANDOFF.md
for the complete local/server/GitHub workflow and the user-approved boundaries.

## 模型变更备注与 GitHub 汇总（2026-11-15）

每次实际变更都要在对应模型根目录的 `MODEL_CHANGES.txt` 追加简洁的 UTF-8 备注，
并在同一轮工作中同步到 GitHub `Demon-Laplace/dSph_workbench` 的
`model_changes/<模型名>.txt` 和 `SIMULATION_HISTORY.md`。
这包括 IC、轨道、姿态、物理或数值参数、程序、分析结果、图件及运行的启动、
停止、重试、清理；同一任务的相关改动合为一条，纯读取或状态查询不另记。
Fornax5028 等正式模型的备注放在 `/travail/xhuang/<模型名>/MODEL_CHANGES.txt`；
孤立模型放在 `/travail/xhuang/isoFornax/<模型名>/MODEL_CHANGES.txt`，
并使用 `model_changes/isolation/<模型名>.txt` 镜像，避免同编号混淆。
共用程序或模板变更在所属目录留同名 TXT；影响具体模型时在模型备注中引用。

每条只需：唯一记录号、日期、模型/轨道/运行名、前序来源、改动原因与实际内容、
验证到哪一步，以及必要的回执路径。使用清楚的名称，不凭模型编号推断来源。
区分候选、已写入 IC、已提交、启动通过与完整科学验证；图形预测不能写成演化结果。
删除前先保存必要的小型结构化结果、来源、参数与校验回执，核实无活动作业使用目标，
只清理用户明确授权的范围。备注记录保留位置和实际删除范围。

GitHub 汇总保留既有条目，并追加同一记录号、简短内容及 TXT 镜像链接；
镜像保留模型 TXT 的完整备注。修正历史时追加更正，不抹去旧记录。
不上传粒子数据、完整日志、缓存或凭据。遵守既有 Git 状态检查和并发更新保护，
不得覆盖他人改动或强推；同步失败时保留本地备注并明确标记“待同步”，
恢复后补齐，不得声称已经同步。只在两端回读一致后报告同步完成。
