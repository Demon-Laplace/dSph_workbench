# dSph 工作流：智能体交接与操作说明

本文件供没有对话上下文的智能体使用。先读本文件，再读相关代码；不要重新设计已稳定的目录结构。
最后更新：2026-10-16。示例中的运行状态是历史事实，执行前必须重新查询。

## 1. 三端职责与边界

- 本地分析工作区：`/home/x-huang/simulation/model/Fornax`，统一入口 `python3 workflow.py`。
- IC 工具与独立 Git checkout：`../IsolationModel/Fornax`。IC 参数在其 `IniFile/`；新模型生成结果落入分析工作区 `<model>/ic/<build>/`。
- 本地源码 checkout：`tools/dSph_workbench`，GitHub `Demon-Laplace/dSph_workbench`，只使用远端 `main`。先检查 status/fetch，禁止覆盖用户修改或 force-push。
- 服务器：`xhuang@mesopsl`，经 `xhuang@rubicon.obspm.fr` 跳转。使用现有 SSH 配置和凭据，绝不复制/打印密钥。连接示例：`ssh -J xhuang@rubicon.obspm.fr xhuang@mesopsl`。
- **服务器写入仅限 `/travail/xhuang`，它是用户的 data 目录。** `/obs/...` 下转换器、求解器只能读取/执行，不得修改。
- 服务器源码 checkout：`data/workspace/code/github-sync`。根目录 `workflow*.py` 是已安装副本；更新代码必须同步 GitHub main 和安装副本。
- `data/workspace/releases/` 是不可变分析版本，禁止直接编辑。`data/dSph_workbench/` 是另一套生产分析源码，不能用 isolation 工具整体覆盖。
- 机器配置 `workflow.json`、大数据、IC 粒子文件、日志不进入 Git；代码和小型 IC 总表 `IC_MODELS.csv` 进入 Git。

服务器与本地结构不同，这是用户的明确选择：

```text
本地 <model>/
  ic/<build>/                 新生成且冻结的 IC
  isolation/                  旧 IC、孤立模型结果；保留一级
    density_output/           下载的检查结果
本地 workspace/receipts/isolation/<model>/<run>/  上传与运行回执
服务器 data/isoFornax/<model>/runs/<run>/
  output/                     HDF5 快照
  density_output/             检查图、统计表、视频
服务器 data/<model>/          正式轨道模拟，与 isolation 分开
```

不要重新创建本地 `isolation/legacy/<model>/` 或 `isolation/<run>/client/`。
不要按编号相同就认定旧 IC 与正式模拟同源。

## 2. 接手先检查

1. 读当前目录及目标服务器目录的 AGENTS.md；检查 Git 状态和 `workflow.json`。
2. `python3 workflow.py --help` 查看实际安装命令；必要时与源码版本比对。
3. 查 `squeue -M all --me`、目标运行 `status.json` 和 `sacct`，避免重复提交或修改活动目录。
4. 任何需要删除的运行必须先验证其身份、Slurm 终态和保留副本。查询失败时停止清理。
5. 不读/展示凭据文件。服务器 Python 使用 `. /etc/profile.d/z00_lmod.sh && module load python/3.9.1`；Git 模块为 `git/2.23.0`。

## 3. 生成或接入 IC 并运行

在本地分析工作区执行：

```bash
# 从已配置的 IniFile 参数生成新 IC、登记总表、上传并提交
python3 workflow.py isolation Fornax7004 --name fd25-v1 --profile fd25-paf40 --submit
# 复用新流程已有的冻结 IC
python3 workflow.py isolation Fornax7004 --ic-run BUILD --name fd25-v2 --profile fd25-paf40 --submit
# 接入用户以前生成、缺少 build.json 的 IC；不能重新生成替代用户文件
python3 workflow.py isolation Fornax7002 Fornax7003 --existing --name fd25-paf40-20261013 --profile fd25-paf40 --submit
python3 workflow.py isolation-status --model Fornax7002 --name fd25-paf40-20261013
```

`--existing` 读取 `<model>/isolation/Fornax.ini` 和冻结的 `IC_<model>.ini`，校验全部粒子、计数、质量、唯一 ID、有限值和 SHA256。生成历史/随机种子未知时不得编造。
新生成的粒子文件带模型编号，兼容名 `Fornax.ini` 是链接。上传按模型身份、校验和分发，不靠手动重命名。
不传 `--submit` 会上传、准备并执行 Slurm 预检。转换二进制 IC 与运行在计算节点执行。
TREECOOL 自动从配置模板复制为普通文件，禁止依赖复制后失效的相对链接。

当前已配置 `fd25-paf35` 和 `fd25-paf40`：feedback 对应 `timescale_fd=24.9`，**不是 25.0**；BufferSize=120；PartAllocFactor 分别为 3.5/4.0。仅 7002、7003 已明确授权这个 feedback，其他模型按用户参数选择配置，不擅自套用。
普通 `default` 模板的 feedback 是 16.6，不能误用于这两个模型。
模拟精度保持原样；检查工具可沿用已经验证的 float32，不能顺手把模拟数据降为 float16。

## 4. 启动失败、重试和自动替换清理

用户经验标准：正常运行超过 **140 秒（2 分 20 秒）且已进入积分**，判断启动成功；完整模拟仍要求快照达到 TimeMax。
7002、7003 在 PartAllocFactor=3.5 时均于初始快照后退出码 11；提高到 4.0 后成功。退出码 11 是段错误，不能断言系统物理内存耗尽。

处理新失败：

1. 确认作业已停止，检查 solver.log/slurm 和 sacct。内存/分配相关启动失败按用户授权以 0.5 步长调整 PartAllocFactor；保留 feedback、BufferSize 和其他物理参数。
2. 基于原配置新增明确的 profile（不要修改其他运行的冻结配置），用相同 IC 和新运行名重试。不要向有快照的目录盲目重启，也不要手改冻结 PARAM 的校验和。
3. 本地和服务器 `workflow.json` 都启用 `cleanup_failed_isolation: true`。新运行的作业内监视器在运行 140 秒后检查积分进展，并自动删除同模型、同 IC 的旧启动失败目录，包括旧 slurm 输出、初始快照、重复 IC 和自动生成的 run.sh。本地会话退出不影响监视器。
4. 这是对失败尝试的自动替换：最终只保留成功运行数据；不会永久堆积每次尝试。旧版本已冻结的作业不含监视器时，用以下命令补执行：

```bash
python3 workflow.py isolation-clean-failed Fornax7002 Fornax7003 --name fd25-paf40-20261013 --dry-run
python3 workflow.py isolation-clean-failed Fornax7002 Fornax7003 --name fd25-paf40-20261013
```

`--name` 是保留的替代运行。新运行尚未正常启动时延后清理；查询 Slurm 失败时不删。
旧运行须为 FAILED/OUT_OF_MEMORY、运行不超过 140 秒，并有内存/退出码 11 证据；正常完成或长时间运行后中断的结果不能按此规则删除。
仅保留 `workspace/receipts/isolation-cleaned/` 中的小型删除记录；不归档失败的大文件。已验证接收的上传暂存副本也清除。
已删除的运行名有删除回执保护，不能误用旧命令重新提交；按保留运行名查询当前状态。

## 5. 完成后处理 isolated 模拟

```bash
# 默认处理全部快照并生成视频；与历史密度检查流程一致
python3 workflow.py isolation-check Fornax7002 Fornax7003 --name fd25-paf40-20261013
# 相同命令不会重复提交，返回原检查 JobID 和当前调度状态
# 检查结束后取回图和统计表，不搬运全部快照
python3 workflow.py isolation-fetch Fornax7002 Fornax7003 --name fd25-paf40-20261013
# 可选：显式 quick 只抽查首、中、末帧，不生成视频
python3 workflow.py isolation-check Fornax7002 Fornax7003 --name fd25-paf40-20261013 --range quick
```

默认不传 --range 会处理全部快照并生成视频；服务器按原检查脚本使用 16 个进程。
`--range quick` 只抽查首、中、末三帧，不生成视频。
`--range 0,300,601` 只检查列出的快照；`--range 0-601` 检查区间。
共享程序：服务器 `data/isoFornax/python/check_density_profile.py`。它依赖 `isoFornax/` 下的共享 Python 模块，保持该安装深度。
历史 shell 实际为 `data/isoFornax/Fornax/subchecdensity.sh`；没有要求每个新运行复制一份。在成功运行目录 `cd` 后，可以直接 `sbatch` 该脚本的绝对路径进行全量检查。
统一入口使用服务器 `isolation_check` 配置提交独立分析任务；既有 .sh 不修改。
分析结果不能证明科学平衡性已通过：应看初/中/末密度、结构、速度弥散等趋势，再决定是否用于正式轨道模拟。本流程不会自动挑选某快照并启动正式轨道模拟。

## 6. 修改、测试、发布

主要文件：
- `workflow/workflow.py`：命令入口。
- `workflow/workflow_isolation.py`：传输、冻结、转换、提交、启动监视。
- `workflow/workflow_isolation_ops.py`：失败清理、密度检查和结果下载。
- `workflow/isolation_python/check_density_profile.py`：共享密度分析；模型名从 request.json 读取，旧模型保留兼容。
- `local_tools/ic_generation/catalog.py`：IC 校验和 GitHub 总表同步。

在源码 `workflow/` 目录运行 `python3 -m unittest discover -s tests`，避免从安装根目录误导入旧 workflow 模块。先跑相关测试，再做适当的真实服务器验证；禁止为了验证而重复跑完整模拟。

通过后提交并推送 GitHub main，再更新独立 IC checkout 和服务器 checkout。服务器无 GitHub 凭据时，从本地传输增量 git bundle 到 data/workspace，再在服务器 `git fetch BUNDLE main`、`git merge --ff-only FETCH_HEAD`。不能把 bundle 路径直接当 merge 参数，不能复制本机凭据。
安装变更文件前，与服务器 Git 旧版本逐字节比对，发现额外改动先核对。根目录 workflow 模块只安装变更部分；共享密度脚本部署到 isoFornax/python，AGENTS.server.md 对应 data/AGENTS.md。不要编辑已提交运行里的冻结程序。

## 7. 用户明确偏好

- 不擅改已有用户 `.sh`。2026-10-13 仅授权过无编号 Fornax 模板五个 sub*.sh 的适配；本次仅额外允许删除失败运行目录内自动生成的 run.sh。
- 模型自身正常的 img/video 保留；不要把检查图和模拟数据混同。
- 服务器历史归档 `observer_review_delete_20261006` 留给用户检查，不顺手永久清空。
- 服务器历史完成作业的无用 Slurm 输出曾获授权清理；本次自动策略仅针对被新尝试替代的短时失败，不扩展成删除所有日志。
- 编号、物理参数、路径不明确时先查现有记录，常规可验证操作直接执行，不反复索要已经给出的授权。

### 多模型密度检查共用一次 Slurm 提交

一次 `isolation-check MODEL1 MODEL2 ... --name RUN` 只调用一次 sbatch，所有模型
共用一个 JobID。在该作业内依次检查各模型，每个模型仍使用配置的 16 个进程；
不把多组进程同时挤在同一份 CPU/内存配额里。申请时长按待检查模型数乘以单模型时长。
每个模型的 density_output/ 与 density-<JobID>.log 分开保存，图和视频不会混合。
某模型失败时停止本批后续模型；已生成结果保留。相同命令重复执行返回原作业号。
已有独立任务不撤回、不重提；只把尚未提交的模型合并为一个新任务。
同模型不同检查范围也会检查已有批任务回执，避免共用输出目录的任务并发写入。
本次 7002/7003 的全量检查已分别开始，保留原作业 236608、236609；批量合并用于后续提交。
