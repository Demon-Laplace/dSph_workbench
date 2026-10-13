## 本地归并、服务器分离（2026-10-13）

本地新 IC、隔离运行回执和正式模拟按同一模型 ID 归属。以下本地布局取代旧 generated/ 和
workspace/models/isolation/ 的新建路径；正式模拟原有文件位置保持可用：

```text
Fornax/Fornax7003/
  ic/<运行名>/            # AGAMA 结果、冻结参数、回执
  isolation/<运行名>/     # 隔离模拟；本地 client/ 保存上传回执
  isolation/legacy/Fornax7003/ # 历史文件；保留末级模型名供旧脚本识别
  isolation/source       # 正式模拟实际引用的隔离模型（有依据时创建）
  isolation/source.json  # 来源与编号复用记录
  GZWJL.PARAM             # 正式模拟参数
  output/                 # 正式模拟输出
```

服务器保持分离：正式模拟在 `data/Fornax7003/`，历史隔离模拟在
`data/isoFornax/Fornax7003/`，新隔离运行在 `data/isoFornax/Fornax7003/runs/<运行名>/`。
本地目录与服务器目录由工作流映射，上传时不改变服务器已有模型组织。
分析工作区配置 `models_root: "."`，IC 工具工作区配置 `models_root: "../../Fornax"`。
参数库、总表和 Git 代码仍在 IC 工具工作区。旧 generated/ 结果仍可复用，已迁移目录的旧位置保留兼容链接。

```bash
python3 workflow.py build-ic --model Fornax7003 --name ic-v2
python3 workflow.py isolation Fornax7003 --ic-run ic-v2 --name isolation-v2 --submit
python3 workflow.py isolation-status --model Fornax7003 --name isolation-v2
python3 workflow.py ic-catalog --model Fornax7003 --run ic-v2
```

`new-case` 可在已含 ic/、isolation/ 的模型目录准备正式输入，已有正式文件时拒绝覆盖。
同一运行名可在不同模型中使用，查询时用 --model 消除歧义。

2026-10-13 已按服务器 `mergegalaxy.par` 和 IC 链接盘点 101 个正式模拟、89 个隔离模型。
没有一个可靠的统一编号起点：2028 已使用同号 IC，但 2037 引用 2035、2042 引用 158，
2052/2053/2055 引用 2047。[MODEL_SOURCES.csv](MODEL_SOURCES.csv) 记录声明的来源，
未能从参数和链接确定的 10 个条目明确保留 unknown，不推测对应关系。

本地 Fornax2047+ 的 91 个目录（含 7002、7003）已归入同号模型的
`isolation/legacy/<模型名>/`，旧位置保留兼容链接。本地 `source` 链接另指向正式模拟
实际使用的隔离模型；这是当前参数管理范围，不表示 2047 是编号统一起点。
较早的 207 个本地目录保持原位，没有按相同编号推断对应关系。

服务器此前移动的 89 个目录已全部恢复到 `data/isoFornax/`，并复核原目录 inode 和
顶层文件大小；新增的归并来源链接、共享输入别名及空目录已清理，来源记录留在回滚归档。
恢复回执在 `workspace/review/unified-models-20261013/server-restoration-summary.json`。
已完成且冻结绝对路径的测试运行保留原位置，重试会识别旧运行，避免重复提交。

工作流支持 `remote.proxy_jump` 显式配置跳板。公钥授权并解锁到 ssh-agent 后，可不依赖临时共享
连接。远端认证配置由用户处理，数据迁移继续限于 data/。

# Fornax 本地、服务器与 GitHub 工作流

日常入口是本目录的 `workflow.py` 和 `workspace/`。服务器保存完整模拟数据并执行 Slurm 任务，本地维护参数、检查版本、运行交互式分析和查看回传成果。服务器写入范围固定为 `workflow.json` 中的远程数据根，即 `~/data`。

## 目录与数据归属

两端代码、任务和结果使用以下分类；模型存储按上节分别处理。兼容索引避免破坏 `mergegalaxy.par`、已有脚本和笔记本的相对路径；新分析任务和新生产案例实际写入新结构。

| 目录 | 用途 |
| --- | --- |
| `<模型名>/` | 本地归并模型目录；服务器仅正式模拟 |
| `workspace/models/production/` | 正式模拟兼容索引 |
| `workspace/models/isolation/` | 服务器隔离模型索引，实际数据仍在 isoFornax/ |
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

配置集中在 `workflow.json`，不保存密码。`remote.proxy_jump` 指定跳板，公钥授权并
在 ssh-agent 中解锁后，工作流可直接连接；已验证无需保留 Tabby 或临时共享连接。
现有 `control_socket` 仅在对应 socket 可用时复用，失效时使用正常公钥连接。
服务器使用既有 module 环境。服务器修改范围限定为配置的 data 根目录。

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

本机 `simulation` Conda 环境已有 AGAMA，当前默认 Python 未安装该模块。`workflow.json` 的 `ic_python` 指定建模解释器，也可用 `--python` 显式覆盖。统一建模入口会拒绝已有输出目录，冻结输入参数、记录随机种子和求解器 SHA256，生成独立粒子 IC；不会把粒子数据推入 GitHub：

参数与生成工具在 `model/IsolationModel/Fornax/`；粒子 IC 统一生成到
`model/Fornax/<模型名>/ic/<运行名>/`，与同模型的 isolation 和正式模拟归在一起。
在分析目录的 `workflow.json` 配置 `ic.root` 和 `ic.repository`，二者均为
`../IsolationModel/Fornax`；在 IC 工作区配置为 `.`，`github.checkout` 为 `.`。
`ic_python` 指向安装了 AGAMA 的 Python。路径相对于各自工作区，服务器需要单独配置到 data/ 内。

```text
IsolationModel/Fornax/
  IniFile/              # 已确认的参数；原有文件保留
  drafts/               # new-ini 产生的参数草稿
  generated/<旧运行名>/ # 兼容旧 IC 运行，新任务写入模型目录
  Fornax2047+/          # 指向统一模型内历史目录的兼容链接
  IC_MODELS.csv         # 自动维护并同步 GitHub 的参数总表
  local_tools/ic_generation/  # 版本化生成程序
```

在分析工作区运行：

```bash
python3 workflow.py build-ic --model Fornax7003 --name Fornax7003-v1 --threads 2 --seed 42
python3 workflow.py new-ini 7004 --base 7003 --set Density_stellar.axisRatioZ=0.6
python3 workflow.py build-ic --ini ../IsolationModel/Fornax/drafts/IC_Fornax7004.ini --name Fornax7004-v1
```

在 IC 工作区也可使用 `python3 workflow.py ...`，或者批量生成：

```bash
python3 model_tool.py build "7002-7003" --processes 2 --threads 2 --seed 42
```

旧 `makeini.py` 兼容入口转交给 `model_tool.py`。批量命令生成带时间和随机标识的独立目录；
`--overwrite` 已停用。旧流程的关联扁率修正和 Fornax2047+ 参数同步保留，在批量生成前执行。
新 `build-ic` 使用指定参数的原样快照，不自动修正科学参数。`--ini` 路径相对于工作区；
`--model` 从 IC 工作区的 `IniFile/` 读取。草稿没有自动成为正式参数文件。

每次成功生成并检查粒子数量、有限值和正质量后，会追加 `IC_MODELS.csv`，自动提交并推送
GitHub `main`。重复登记同一运行不会增加重复行，不同运行保留各自记录。锁串行化同一
Git 工作区的并行登记；先快进获取远端变化，再提交精确的总表文件。代码工作区存在未提交
修改或分支分歧时停止同步，不强制覆盖。生成失败不会登记为成功。

网络或 Git 推送失败时，退出码为 3，已完成的粒子文件和 `catalog-entry.json` 保留；
`build.json` 中生成状态仍为 `complete`，`catalog_sync.state` 为 `pending`。恢复连接或
解决 Git 工作区问题后，仅重试总表，不重新计算：

```bash
python3 workflow.py ic-catalog --run Fornax7003-v1
python3 workflow.py ic-catalog --refresh  # 扫描历史 IniFile/ 和原模型目录，补充历史记录
```

总表包含名义气体/恒星质量（M_sun）、粒子数、尺度半径及厚度（kpc）、Sersic 指数、恒星
扁率、beta/kappa、随机种子、UTC 生成时间及文件 SHA256。质量单位来自求解器的 1e10 M_sun；
它们是参数的名义质量，并非对旧粒子文件重新积分。气体 scaleHeight 的负号保留原 AGAMA
参数约定。历史多气体分量的质量和粒子数求和，尺度按分量列出，未知分量保存在 other_components。

`existing_ic` 表示保存的参数和 IC 文件存在，且文件头气体/恒星数量匹配；没有逐粒子复验。
`existing_ic_check_needed` 表示文件头数量需要检查；`parameters_only` 表示仅有该版本参数，
不推定已经成功生成；`generated` 是新流程验证完成的运行。历史冻结参数与当前参数不同
时分别记行，不把新版参数错误地归给旧 IC。测试模型以 smoke 命名，不用于物理结论。
历史生成时间未知时留空。Git 只保存总表和源码，粒子、日志、图片和机器配置不上传。

该命令保留原 Schwarzschild 求解器及计算精度。每次输出包含 `Fornax.ini`、`parameters.ini`、
`build.json`、`build.log` 和分量文件。失败目录保留供排查；下一次使用新名称。

两级粒子 IC 合并工具仍使用服务器已有 `zimpack` 安装。2026-10-12 已在独立测试目录实测合并与短时模拟；`new-case` 仍只准备输入，未自动执行合并。`simulate` 会检查合并后的 `InitCondFile` 是否存在。完整规模及长时间稳定性不在本次验证范围内。

## 整理与清理记录

```bash
python3 workflow.py organize
python3 workflow.py remote organize
python3 workflow.py list
python3 workflow.py audit
```

`organize` 更新分类入口，不搬动历史模拟数据。本地根目录原有 23 个旧 Python 脚本或链接已清除：12 个仍被历史笔记本或专用检查使用的工具迁到 `tools/dSph_workbench/local_tools/analysis/`，通过 `workspace/code/local-analysis/` 访问；建模核心统一放到 `local_tools/ic_generation/`。生产分析仍使用 `repo/` 中经过验证的当前版本，两套源码用途明确区分。

本地按用户要求使用系统回收站：31 个旧脚本/旧链接、12 个迁移前原件、15 个笔记本修改前原件，以及 4 个缓存/旧分类链接目录，均可在回收站恢复。15 个笔记本已加入显式路径初始化，保留原工作目录及已存输出；已打开的旧内核需重启后从首单元运行。两个旧 `makeini.ipynb` 已迁至独立 IC 工作区的 `notebooks/`，接入统一生成命令；分析工作区只在 `workspace/notebooks/` 保留访问链接。原件在系统回收站，原诊断单元与已存输出保留，旧粒子重导出单元停用。`test_pyscript` 的专用研究程序保留，其空脚本已清理。

逐文件清单及恢复位置：`workspace/catalog/local-cleanup-receipt-20261012.json`。不要恢复旧入口覆盖统一工作流；如需历史代码，可从回收站恢复到单独目录。模型快照、IC、图片和视频未参与这次本地脚本清理。

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

`tests/test_workflow.py` 覆盖输入与输出隔离、原数据变化检测、禁止覆盖、路径逃逸、快照位置歧义、统计完整性及日志识别。另已实测统一 `build-ic` 入口：1536 个粒子与旧流程输出逐值一致，并验证了拒绝覆盖已有输出。长时间模拟尚未验证。

### 2026-10-12 建模与短时演化验证

完整记录位于本地及服务器的 `workspace/validation/ic-smoke-20261012/`，摘要为 `workspace/catalog/ic-validation-20261012.json`。这是缩小轨道库和粒子数的流程验证模型，轨道约束存在残差，不用于物理结论。

- AGAMA 生成 512 个气体和 1024 个恒星粒子；质量与 AGAMA 实际密度积分匹配。Sersic 的实际积分质量与 INI 名义质量约差 0.0102%，没有人为修改粒子质量。
- 两级 IC 合并得到 3136 个粒子，独立读取二进制块验证了数量、唯一 ID、位置、速度及质量。
- 原隔离求解器成功推进到 `TimeMax=1e-5`，写出三帧有效快照。
- 原周期边界求解器在 `PartAllocFactor=5.5`、4 或 2 个 MPI 任务时崩溃。GDB 显示 `drift_particle` 中的 SIMD 指令访问仅按 8 字节对齐的地址；旧分配器的四处申请/重申请都只向 8 字节取整。
- 在独立 `data/workspace/validation/` 副本中将四处取整改为 16 字节，复用既有其他目标文件，仅重编译分配器并重新链接；原参数下 4 和 2 个 MPI 任务均通过，Slurm 作业 236563 为 `COMPLETED`。候选没有替换原安装或正式模型启动命令。
- 原程序仅改变测试内存预留系数到 6.0 也能完成，验证了布局相关性；该数值不是通用修复。修复候选的 4 任务结果与这个对照的三帧坐标、速度、质量、密度及内能逐值一致。2 任务结果有微小并行差异。

补丁保存在 `workflow/solver_patches/gizmo_allocator_alignment16.patch`，只适用于经过匹配检查的旧分配器源码。不要直接修改服务器 `data/` 外的程序。补丁没有改变物理方程或数值类型；生产规模、长时间演化及重启读写仍未验证，不能据短时测试直接替换所有正式求解器。

`workflow/tests/verify_simulation.py` 验证快照确实推进到目标时间、粒子数和 ID、质量及有限值；它针对固定粒子数的短时测试，不适用于发生星形成/粒子合并的通用生产验证。Slurm `COMPLETED` 本身不足以判定成功，调试器可在被调试程序崩溃后返回零退出码。

## 分支整理与精度约定

原有 5 个功能分支已逐文件比对；演化和论文图分支被恒星质量分支包含，圆孔径速度弥散改动已等价包含。统一源码合并这些功能、服务器输入路径/绘图修复及分块合并工具，并保留新旧恒星质量分量，避免把新生恒星质量重复计入总质量。各原分支用 `archive/20261012/` 标签及本地 Git bundle 保留历史；日常远端分支仅使用 `main`。

检查分析使用 float32 粒子和质量数组；不改模拟求解器精度，不统一启用 float16。真实 Fornax7001 的 100–102 帧回归已验证原有 36 列一致，新增 6 个诊断字段。分支科学测试及原工作流保护测试各 12 项通过，新增 6 项 IC 总表及目录隔离测试也已通过。论文图整套生成和完整历史演化任务未重新运行。

工作目录根部原先指向其他项目的空 `.git` 已可逆归档到 `workspace/archive/local_legacy_20261012/empty_git_metadata`。日常 Git 命令应在 `tools/dSph_workbench/` 执行。

### 独立 IC 工作区验证

2026-10-12 在 `IsolationModel/Fornax/generated/ic-catalog-smoke-20261012/` 实际生成了
1536 个粒子，`Fornax.ini` SHA256 与先前验证输出完全相同；自动登记与 GitHub main 推送成功。
历史总表包含 299 个模型的 308 条参数记录（288 条 existing_ic、20 条 parameters_only），
另有两次明确标记为 smoke 的流程测试记录。旧测试 IC 已迁至同一 generated/，原路径保留兼容链接。
18 项定向测试通过，覆盖重复登记、远端快进、推送失败重试、文件篡改拒绝、质量单位和目录逃逸。
两个笔记本所有代码单元通过语法检查，目录初始化实测通过；整套诊断绘图未重跑。


## 一条命令运行 isolation

新入口把本地 AGAMA 生成、总表同步、粒子与参数校验上传、服务器 case 准备、IC 格式转换、
Slurm 提交和结束时间检查串起来。服务器首次需要配置 `isolation_profiles.default`；
配置必须从实际可用的隔离模板和启动脚本核对，不能猜测转换器参数或复制周期边界求解器。
服务器已配置 default 和 smoke 两套 profile。default 沿用隔离模板与 astro_thin 的
24 进程、99 小时配置。smoke 使用 mesopsl1 的 4 进程、8 分钟上限及短积分参数。
真实端到端测试作业 236569 已完成：1536 粒子、3 个快照、最终时间 1e-5，重复命令仍返回
原 JobID。总表自动提交为 ae63a71；这验证传输、转换和短时积分，不替代长期稳定性检查。

在分析工作区或 IC 工作区使用：

```bash
# 完全只读预览，不连接服务器、不生成 IC
python3 workflow.py isolation Fornax7003 --name Fornax7003-iso-v1 --dry-run
# 一条命令：生成 -> 登记总表 -> 校验上传 -> 提交隔离模拟
python3 workflow.py isolation Fornax7003 --name Fornax7003-iso-v1 --submit
# 复用模型 ic/ 中已完成的 IC（兼容旧 generated/），不重复计算
python3 workflow.py isolation Fornax7003 --ic-run Fornax7003-v1 --name Fornax7003-iso-v2 --submit
python3 workflow.py isolation-status --name Fornax7003-iso-v1
```

不传 `--submit` 会完成生成/上传/准备，然后仅运行 `sbatch --test-only`。随后加 `--submit`
重复同一命令即可提交。相同运行名绑定固定 IC、参数与 profile：已生成的 IC 校验后复用，
已接收的 case 不重复上传，已提交的作业返回原 JobID，不再提交。参数或种子改变需要新运行名。
上传中断时重试使用新的 staging 目录，完整接收的 case 原子发布；未完成的 staging 保留供检查。
提交响应不确定时停止盲目重试，需要先核对 Slurm，防止产生两个同时写入同一目录的作业。
失败/取消任务的重启不是这个命令的隐式行为，需要检查后使用新的运行名或既有重启流程。

本地回执在 `<模型名>/isolation/<name>/client/`，服务器输出在
`data/isoFornax/<模型名>/runs/<name>/`。`isolation-status` 读取服务器保存的状态；
如果 Slurm 在启动前取消作业或强制杀死进程，应根据回执中的 JobID 查 `sacct`，因为进程可能
没有机会更新状态。生成、转换、求解日志分别保留。最后快照必须达到 `TimeMax` 才标记完成；
这不代替平衡性、收敛性或科学结果检查。服务器执行验证的 Python 需要 h5py。

服务器 `isolation_root` 配置为 `isoFornax`（默认值），与本地 `models_root` 分开。
服务器 profile 保存在 data 根目录的 `workflow.json` 的 `isolation_profiles` 字段，键为 profile 名：

- `template`：data 内已有隔离模板目录，包含 `GZWJL.PARAM` 和可选 `TREECOOL`。
- `converter` / `solver`：根据既有可运行脚本填写的命令参数数组。占位符有 `{ascii}`、
  `{binary}`（预期 `Fornax.IC_000`）、`{prefix}`（Fornax 路径前缀）、`{parameters}`。
- `setup`：既有模块加载命令，写入作业脚本；`python` 指定支持 h5py 的 Python。
- `sbatch`：资源配置，可用 clusters、qos、partition、account、nodes、ntasks、cpus-per-task、time、mem、constraint。
- `parameter_overrides`：明确指定的参数覆盖；默认保留模板物理参数。输入和输出路径自动重写，
  `ResubmitOn` 置零，由本流程负责记录提交，避免旧脚本自动重提绕开防重复保护。

模板、脚本、配置和输入 SHA256 一同冻结；修改已准备的参数会被拒绝。求解器及转换器可以
只读调用现有安装，所有新文件限制在指定 data 根目录内。GitHub 仍只同步源码和 IC 总表。

## IC 工作区旧文件清理

IC 根目录的 24 个历史检查/分析工具迁入 `local_tools/isolation_analysis/`，其中 6 个逐字节相同
的模块用相对链接复用公共工具，18 个有不同逻辑的文件独立保留；不以日期判断科学程序是否无用。
从 IC 根目录运行这些工具，例如 `python3 local_tools/isolation_analysis/CheckDensityProfile.py --help`。
旧工具原件、退役 benchmark 启动脚本、重复 Schwarzschild 求解器、根缓存和 5 张独立比较图放入
系统回收站。模型内部图片、粒子、轨道及模拟数据保留。

8 个其余笔记本迁到 `notebooks/`，添加工作目录和旧工具导入初始化，原件在回收站。打开后重启
内核并先执行首单元。`Test.ipynb` 原有一个未完成的 `df[]` 草稿单元保持原样。
历史统计表格和演示文件集中在 `reports/historical/`；包含粒子文件的旧性能测试目录保存在
`archive/benchmark_results/`。观测输入仍保留原路径，避免旧检查程序失去输入。
逐项清理记录在分析工作区 `workspace/catalog/ic-cleanup-20261012.json`。

## 服务器历史 isolation 的公共 Python 脚本

有效历史模型的 14 个 Python 文件集中于 `data/isoFornax/python/`，源码对应
GitHub 的 `workflow/isolation_python/`。各模型保留同名 `.py` 符号链接，所以
原有 `.sh` 的内容、位置和运行命令保持不变。仍应在模型目录内提交，例如：

```bash
cd data/isoFornax/Fornax7001
sbatch subchecdensity.sh
```

这次只合并逐字节相同的历史版本，不改变分析精度或科学算法。公共目录与模型
目录处于同一层，保留密度检查程序推导 `isoFornax` 和 `data` 路径的行为。
原始副本在 `observer_review_delete_20261006/isolation-python-20261013/`；
核对后删除这些备份不会破坏公共脚本。审计清单和验证结果存于
`workspace/review/isolation-python-20261013/`。

现有用户 `.sh` 未经明确要求不得改写、移动或替换。新 isolation 继续使用前述
`workflow.py isolation ... --submit` 自动流程，不需要手动向每个模型复制检查脚本。
历史公共代码不覆盖已冻结运行内的代码或启动文件。
