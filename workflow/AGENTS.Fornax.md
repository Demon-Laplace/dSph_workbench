# Fornax 服务器目录约定

本文件位于 `/travail/xhuang/Fornax/AGENTS.md`，适用于此轨道模板目录。
同时遵守上级 `/travail/xhuang/AGENTS.md`；带编号的同级模型也由上级规则覆盖。
先读 `/travail/xhuang/workspace/code/github-sync/workflow/AGENT_HANDOFF.md`，
并查看 GitHub 的 `SIMULATION_HISTORY.md` 与对应模型的 TXT 备注，核实来源和当前状态。

仅在 `/travail/xhuang` 内写入。现有用户启动脚本、冻结输入、原始快照及其他模型
保持保护；本文件不自动授权修改模板、提交作业或删除数据。
轨道模拟沿用既有一次合并流程和 `subastrofat.sh` 的资源配置。
新分析使用已选定的 `astropy_native_direct_v1` 观察者合同，显式固定太阳常数。

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

简洁备注示例：

```text
记录: F5028-PA43-20261115
日期: 2026-11-15
对象: Fornax5028 / orbit_PA43
来源: 独立 dwarf snapshot50 IC；既有 MW/背景；注明轨道候选及前序运行
变更: 为什么改、具体改变了什么
验证: IC 校验/提交/启动/完整演化的实际状态
依据: 参数、校验和、启动或清理回执路径
```
