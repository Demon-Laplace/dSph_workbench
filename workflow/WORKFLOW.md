# Fornax 本地、服务器与 GitHub 工作流

日常入口是本目录的 `workflow.py` 和 `workspace/`。服务器保存完整模拟数据并执行 Slurm 任务，本地维护参数、检查版本、运行交互式分析和查看回传成果。服务器写入范围固定为 `workflow.json` 中的远程数据根，即 `~/data`。

## 目录与数据归属

两端均使用以下分类。历史模型采用目录链接保留原路径，避免破坏 `mergegalaxy.par`、已有脚本和笔记本的相对路径；新分析任务和新生产案例实际写入新结构。

| 目录 | 用途 |
| --- | --- |
| `workspace/models/production/` | 生产模拟，历史模型索引及新案例 |
| `workspace/models/isolation/` | 服务器隔离模拟索引，原始存储为 `isoFornax` |
| `workspace/models/milky-way/` | 银河系模型索引 |
| `workspace/code/production/` | 当前生产分析代码；本地为服务器源码镜像 |
| `workspace/code/github/` | 独立 GitHub 工作副本 |
| `workspace/inputs/parameters/` | GitHub 参数、旧参数库和本地参数草稿 |
| `workspace/inputs/parameter_batches/` | 服务器已校验的上传参数批次 |
| `workspace/inputs/observations/` | 观测表和参考数据索引 |
| `workspace/inputs/shared/` | 大型共享输入索引 |
| `workspace/templates/` | 原有生产与隔离模板索引 |
| `workspace/analysis/<任务名>/` | 固定代码、输入记录、统计、图片、视频及运行状态 |
| `workspace/results/remote/<任务名>/` | 本地回传成果，不含原始快照 |
| `workspace/results/models/<模型名>/` | 从历史模型直接取回的成果 |
| `workspace/notebooks/` | 旧笔记本索引；运行旧笔记本仍应使用原工作目录 |
| `workspace/catalog/` | 目录分类与数据放置清单 |
| `workspace/archive/` | 旧版脚本或服务器待删归档入口 |

快照和粒子 IC 的完整副本以服务器为主。参数 INI 是小型模型定义，不能与数十 MB 到数 GB 的同扩展名粒子 IC 混淆。默认成果同步排除所有 `.ini` 和 `.hdf5`，不会自动上传或下载大型原始输入，也不使用 `rsync --delete`。

本地 42 个快照文件约 40.63 GiB，其中 14 个找到同尺寸远程候选。`workspace/catalog/local_snapshot_placement.json` 记录匹配结果。同尺寸尚不能证明内容相同，因此本次没有删除本地快照；释放这些副本前应做逐文件 SHA256 校验。

## 连接与配置

配置集中在 `workflow.json`，不保存密码。当前 SSH 连接复用 `/tmp/fornax-mesopsl-v2.sock`。若连接过期，在本机终端执行并在终端内输入密码：

```bash
ssh -M -S /tmp/fornax-mesopsl-v2.sock -o ControlPersist=1h -o ProxyJump=xhuang@rubicon.obspm.fr xhuang@mesopsl
```

保持该终端连接，另一个终端执行工作流命令。服务器侧 Python 和 ffmpeg 使用已有 module 环境；没有修改服务器 home、系统安装或调度器配置。

## GitHub、源码同步与发布

GitHub 保存代码、参数定义、测试和工作流文档；本地和服务器保存机器配置及数据。GitHub 远端只保留 main 分支。正确工作副本为 `tools/dSph_workbench/`，同步分支为 `main`。`repo/` 是统一分析源码目录，包含已核对的分支功能与服务器兼容修复；`workflow/` 管理运行、同步与发布。

服务器 `dSph_workbench/AGENTS.md` 要求智能体在修改源码、完成相关验证后自动同步。当前推荐在本地工作目录执行（复用本机 GitHub 凭据）：

```bash
python3 workflow.py sync-server-github --message "说明服务器修改与验证结果"
```

该入口先导出允许的服务器源码并验证 SHA256，再提交和推送。服务器当前没有独立 GitHub 推送凭据；直接在服务器操作的智能体必须报告认证失败，不能伪报同步完成。若服务器已有认证，可在数据根直接执行：

```bash
python3 workflow.py sync-github --message "说明修改内容与验证结果"
```

在服务器数据根运行这条命令，源码来自 `github.source`，Git 工作副本由 `github.checkout` 指定。命令检查 Python 语法、仅复制允许的源码/文档/测试，提交并推送main 分支。它不会推送数据，也不会强制覆盖远程历史；发现 GitHub 与源码都改动的文件会停止并要求先合并；认证、网络或分支冲突失败时保留本地提交并返回失败。AGENTS.md 是智能体操作约定，不是监控所有编辑器的后台服务。

发布前确保本地工作副本已获取并快进到所需提交。发布只读取干净 Git 提交中的文件：

```bash
git -C tools/dSph_workbench fetch origin
git -C tools/dSph_workbench merge --ff-only origin/main
python3 workflow.py release --name analysis-v1
python3 workflow.py deploy --name analysis-v1
```

`release` 将提交号和逐文件 SHA256 写入 `workspace/releases/<版本>/release.json`；`deploy` 上传到相同分类并验证完整性，验证通过后才切换服务器 `analysis_code`。服务器原始可编辑源码保留，已有分析任务继续使用准备时的代码副本。回退只需激活仍保留的旧版本：

```bash
python3 workflow.py remote activate-release --name analysis-v1
```

发布包包含工作流程序与参数定义用于版本追溯；当前 `deploy` 自动切换的是**分析模块**，工作流入口程序的升级需从该提交的 `workflow/` 安装到数据根后验证。机器专属 `workflow.json` 不进入 GitHub，也不被代码发布覆盖。

## 参数维护和上传

```bash
python3 workflow.py doctor
python3 workflow.py check-ini 7002-7003
python3 workflow.py new-ini 9001 --base 7003 --set Density_stellar.axisRatioZ=0.7
python3 workflow.py push-params 9001 --batch experiment-9001-v1
```

`new-ini` 的新编号必须尚未使用；草稿写入 `workspace/inputs/parameters/drafts/`，不修改 GitHub 原件。它复用现有模型工具的轴比约束同步。上传前先校验参数，按不可覆盖的新批次上传，服务器逐文件验证 SHA256。重复批次名会报错，上传失败的批次应检查后用新名称重试。

当前 GitHub 的 `IC_Fornax4001.ini` 有轴比不一致：`Component_stellar=0.5`、`Density_stellar=0.4`。本次没有擅自更改该物理设定。服务器独有的统计、输入路径和绘图修复已与其他分支逐项比对后整合；保留历史源码备份，避免用某个分支整体覆盖。

## 自动分析与成果回传

从本地执行以下两个命令即可提交任务，并等待成功后自动回传结果：

```bash
python3 workflow.py remote pipeline Fornax7001 --name run7001-example --range 100,102 --processes 2 --submit
python3 workflow.py collect --name run7001-example --wait
```

快照区间为闭区间，`100,102` 包含 100、101、102。省略 `--submit` 时只准备任务并执行 `sbatch --test-only`，不会消耗计算节点运行分析；随后可用 `remote submit --name ... --range ...` 正式提交。任务名必须唯一。

流水线依次执行统计、验证 CSV、画图、验证完整 PNG，并在服务器合成视频。失败时停止后续步骤。每次任务复制当前生产模块并记录 SHA256，快照和粒子 IC 使用链接并记录大小及修改时间。大输入并非冻结副本；源数据改变时会拒绝继续，需重新准备任务。

```bash
python3 workflow.py remote status
python3 workflow.py remote status --name run7001-example --json
python3 workflow.py pull --name run7001-example --dry-run
python3 workflow.py pull --name run7001-example
python3 workflow.py fetch-model Fornax5026 --dry-run
```

`collect --wait` 默认最多等待 7200 秒，可用 `--timeout` 调整。它只在任务完成后下载；Slurm 超时、取消或失败会返回错误。`pull` 可手工下载未完成任务的诊断记录。`fetch-model` 从已有模型获取图片、视频、统计表和小型运行参数，不获取快照。

同一任务目录有进程锁，提交回执防止误重复提交。失败后或希望补做视频时，可以显式续跑：

```bash
python3 workflow.py remote submit --name run7001-example --range 100,102 --stage plot --resume
```

续跑前检查旧任务已离开队列。原始模型的 `elinfo` 和 `output_img` 不受影响。每个分析目录保存 `manifest.json`、`status.json`、`submission.json` 和分阶段日志；旧统计表不直接混入新字段定义。

## 新模拟案例

已有隔离阶段粒子 IC 后，可从明确模板创建新生产案例：

```bash
python3 workflow.py remote new-case Fornax9001 --template Fornax7001 --dwarf-ic isoFornax/Fornax9001/IC_snapshot50_Fornax9001.ini
```

这会复制模板的物理参数、轨道和方向，规范输入路径，在 `workspace/models/production/` 建立实际目录，并保留数据根下的兼容路径。它不会复制模板快照，也不会运行模型。提交前应确认 `case.json` 及物理参数符合本次实验。

```bash
python3 workflow.py remote simulate Fornax9001 --name sim9001-check --launcher subastrothin.sh
python3 workflow.py remote simulate Fornax9001 --name sim9001-run --launcher subastrothin.sh --submit
```

模拟提交复用该模型已有的 Slurm 资源配置和 GIZMO 启动命令。默认仅预检；只有 `--submit` 才提交。存在快照或 restart 数据时拒绝重跑，避免覆盖已有模拟。这里没有自动替代原有重启程序。

AGAMA 初始模型构建尚未纳入可执行的自动流水线：当前检查的两端 Python 环境均未找到 `agama`。粒子 IC 合并工具位于服务器已有 `zimpack` 安装中，但新案例的合并步骤和长时间 GIZMO 模拟本次没有执行验证。参数上传不会伪装成已生成粒子初始条件；`simulate` 会检查合并后的 `InitCondFile` 是否已存在。

## 整理与清理记录

```bash
python3 workflow.py organize
python3 workflow.py remote organize
python3 workflow.py list
python3 workflow.py audit
```

`organize` 更新分类入口，不搬动历史模拟数据。旧笔记本和旧脚本仍可通过原路径使用。五个本地旧脚本变体已移入 `workspace/archive/local_legacy_20261012/`，原位置保留链接。

服务器归档位置：

```text
~/data/observer_review_delete_20261006/review_before_delete/workflow_reorg_20261012/
```

本次迁移包含 32 个与模型并列的检查目录、11 个明确旧源码备份及 2 个完全空的模型目录。模型内部图片、视频、快照和统计结果保留；论文主目录保留。归档是可逆移动，没有删除这些检查成果。

已按用户要求直接删除 741 个已结束任务的 Slurm 输出，文件逻辑大小合计 107940127534 字节，约 100.53 GiB。扫描覆盖整个 `data/`，包括 `isoFornax`、旧代码实验目录及待删归档。识别依据为标准 Slurm 文件名或提交脚本明确指定的输出，未泛删所有 `.out` 文件。删除时重新查询所有集群，校验文件 inode、大小和修改时间；活动任务及变化文件会跳过。删除的日志未做内容备份。

迁移计划、逐项记录和日志删除清单位于：

```text
~/data/observer_review_delete_20261006/_receipts/workflow_reorg_20261012/
```

需要整体恢复本次归档移动时，在服务器已加载 Python 3.9 模块的终端执行：

```bash
python3 ~/data/workflow_archive.py --root ~/data --plan ~/data/observer_review_delete_20261006/_receipts/workflow_reorg_20261012/archive_plan.json --restore
```

恢复不会恢复已经明确删除的 Slurm 日志。归档中的历史清理脚本仅作记录，不能作为日常命令运行。`workflow_cleanup.py` 的删除操作需要实时 Slurm 查询；计算节点只能使用队列快照做盘点，不能凭旧队列快照删除。

## 已完成的验证

已实测参数 7002、7003 的批次上传和 SHA256 校验，以及 Fornax7001 快照 100–102 的 Slurm 提交、统计、三张 PNG、视频合成和成果回传。测试成果位于 `workspace/results/remote/workflow-smoke-20261012/`。这是工作流冒烟测试，三帧历史不足以用于完整轨道或恒星形成史推断。

`tests/test_workflow.py` 覆盖输入与输出隔离、原数据变化检测、禁止覆盖、路径逃逸、快照位置歧义、统计完整性及日志识别。完整新模型生成、IC 合并和长时间模拟不在已验证范围内。

## 分支整理与精度约定

原有 5 个功能分支已逐文件比对；演化和论文图分支被恒星质量分支包含，圆孔径速度弥散改动已等价包含。统一源码合并这些功能、服务器输入路径/绘图修复及分块合并工具，并保留新旧恒星质量分量，避免把新生恒星质量重复计入总质量。各原分支用 `archive/20261012/` 标签及本地 Git bundle 保留历史；日常远端分支仅使用 `main`。

检查分析使用 float32 粒子和质量数组；不改模拟求解器精度，不统一启用 float16。真实 Fornax7001 的 100–102 帧回归已验证原有 36 列一致，新增 6 个诊断字段。分支科学测试及工作流保护测试各 12 项通过。论文图整套生成和完整历史演化任务未重新运行。

工作目录根部原先指向其他项目的空 `.git` 已可逆归档到 `workspace/archive/local_legacy_20261012/empty_git_metadata`。日常 Git 命令应在 `tools/dSph_workbench/` 执行。
