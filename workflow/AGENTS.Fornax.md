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


## 粒子性质与空间成员优先（2026-11-17 用户明确更新）

- 模拟分析按粒子性质、当前空间位置及必要的运动学信息选择矮星系成员。粒子 ID 或初始来源不能作为主成员选择的门槛；不同来源的气体进入矮星系孔径后也应计入。
- ID 重复本身不视为物理错误，不据此终止计算、去重或排除粒子。保留各物理粒子行；只有位置、速度、质量、类型及数值/物理一致性证据明确说明数据有问题时，才针对实际问题报告。不得重新引入 ID 唯一性/来源复核作为空间统计的前置条件。
- 矮星系气体同时输出两个三维口径：严格 `r<5 kpc`，宽松预筛选 `r<15 kpc`；另记 `5<=r<15 kpc` 外层。中心与每个孔径、温度/密度条件、速度参考系均显式保存。宽松空间预筛选不自动解释为引力束缚。
- 恒星年龄直接读取粒子记录的年龄或形成时间字段，并按该字段的定义和单位换算，不依赖 ID 或来源标签；不得将预存恒星统一视为年龄未知。有明确形成时间时使用同一快照时间减形成时间。
- 旧流程 `30 kpc` 只能作为预筛选，不直接当作最终气体成员或 HI fraction 的积分范围。`f_HI=MHI/(Mstar+MHI)` 的分子、分母必须取同一最终三维孔径，分别输出 5 kpc 和 15 kpc 版本；不把外围预筛选气体混入分子而高估比例。
- 本节替代前文中将固定 uint64 ID/cohort 或来源认证作为新主分析必需条件的部分；canonical 物理观察者、太阳、逐星射线和孔径/估计器定义仍保留。源 ID 可仅作附加元数据，不得决定主样本。
- 此前仅因 ID 重复或来源/parent 编码失败的计算须按空间方法重新计算。旧结果保留原定义与历史身份，不能改列名后冒充新空间成员统计。
- 本次 `sigma_spatial_parameters_v1` 结构化补算统一一个 Slurm 任务、单节点顺序处理多模型，16 workers；不依赖视频编码。每模型 CSV/JSON/NPZ 回传到 `F:/Fornax_data/Fornax处理/<实际模型编号目录>/` 直接平铺，轨道/分支用文件名前缀区分。暂不绘图或制作 PDF，不监控作业；自动回传只读取交付清单。

## Ordinary scientific transfers (owner update 2026-11-17)

Use normally readable data directly. Do not routinely compare file sizes,
SHA256/content, reread outputs, or add acceptance/verification gates before
ordinary scientific processing and transfers. Diagnose actual errors only.
This supersedes conflicting older verification requirements below/above.
The current receiver takes an explicit project data root and writes each
model's machine data into its single flat _support directory. Keep only
necessary resume state while transfer is active, not per-file audit reports.
