# dSph 工作流：智能体交接与操作说明

本文件供没有对话上下文的智能体使用。先读本文件，再读相关代码；不要重新设计已稳定的目录结构。
最后更新：2026-10-17。示例中的运行状态是历史事实，执行前必须重新查询。

### 正式模拟资源与本次 IC 选择（2026-10-17）

用户指定正式轨道模拟使用 `astro_fat`：通过原有 `subastrofat.sh`，每个模型
分别提交一个作业，保持启动脚本现有资源配置（当前模板为单节点、96 MPI 进程）。
此规则与多模型检查共用单节点、顺序处理的规则不同；不要改写既有用户 `.sh`。
7004、7005 从各自成功孤立运行的 snapshot 50 用服务器 `isoFornax/snap2csv.py`
转换 IC，保留其 `timescale_fd=16.6`；用户更正：正式轨道与朝向参照 Fornax5028，
不得沿用 Fornax7001 的轨道和角度。其他物理设置保持现有正式模拟约定。
新孤立流程的快照在 `runs/<run>/output/`，可调用 snap2csv 的原有读取、居中、
导出函数，显式传入运行目录；输出放 `isoFornax/<model>/IC_snapshot50_<model>.ini`。
导出后核对全部粒子 ID、类型映射、质量、速度、气体内能和恒定居中位移，再用
既有 `merge-ic` 完成合并及二进制逐粒子验证，成功后提交正式模拟。

`simulate --submit` 默认自动检测启动：作业正常运行至少 140 秒且 info.txt 有非零
时刻积分记录后，写入作业目录的 `startup.json`，退出检测，模拟继续运行。
初始化或只有初始快照不算通过；失败作业直接报错。不要手动反复轮询或询问确认。
仅需提交后立即返回可加 `--no-wait-startup`；`--startup-timeout` 默认 7200 秒，超时
不会取消模拟。已提交作业可用 `remote simulation-startup MODEL --name JOB_NAME`
自动检查，支持 `--timeout`，不会重复提交。
检测通过后生成 `simulation-details.json`，包含模型、snapshot 来源、银河系和背景气体
输入、分类型粒子数、实际合并 IC 的粒子质量范围/均值（太阳质量单位）、轨道角度、
物理参数和申请资源。远程命令自动回传到本地
`workspace/receipts/production/<model>/<job-name>/`，校验 SHA256 后结束，不继续监测。

新提交作业使用独立生成的 `run-monitored.sh`，完整保留原 `subastrofat.sh` 等脚本及
其资源指令，冻结监测 Python 源码。监测是同一 batch 作业内的后台轻量进程，不额外
sbatch/srun、不多申请节点或 CPU；每 15 秒读一次 info.txt，其余时间休眠。
从作业实际执行开始计时，排队期间不启动监测。生成报告后后台进程自行退出，模拟
继续；模拟提前退出时 wrapper 停止监测，并保留原模拟退出码。本地只等待报告文件并
下载，断线不影响服务器检测；恢复后用 `simulation-startup` 取回。旧作业保留原样，
独立检测入口兼容旧提交回执，不为安装监测而取消、重提已有模拟。

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
    IC_density_<model>_*.mp4 检查视频，直接放在 isolation/
    density_evolution.csv     密度统计表
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

当前已配置 `fd25-paf35` 和 `fd25-paf40`：feedback 对应 `timescale_fd=24.9`，**不是 25.0**；BufferSize=120；PartAllocFactor 分别为 3.5/4.0。7002、7003 使用此 feedback；用户已另行指定 7004、7005 使用 16.6，不得混用。
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
# 检查结束后取回视频和统计表，不搬运全部快照
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

多个模拟的检查必须合并为一个 Slurm 作业、一个节点，包括运行名不同的模拟，
不能因 `--name` 不同就拆成多个独立作业。相同运行名仍可使用
`isolation-check MODEL1 MODEL2 ... --name RUN`；不同运行名用重复的 `--run MODEL=RUN`
明确指定，也可用 `--name` 提供默认值、`--run` 覆盖个别模型。例如：

```bash
python3 workflow.py isolation-check Fornax7004 Fornax7005 \
  --run Fornax7004=fd16p6-paf60-20261016 \
  --run Fornax7005=fd16p6-paf40-20261016 --fetch
```

`isolation-fetch` 支持相同的映射与 `--wait`；`--video-only` 也遵守此规则。
一次批量检查只调用一次 sbatch，所有待提交模型共用一个 JobID，强制 nodes=1、ntasks=1。
在该作业内依次检查各模型，每个模型仍使用配置的 16 个进程；
不把多组进程同时挤在同一份 CPU/内存配额里。申请时长按待检查模型数乘以单模型时长。
每个模型的 density_output/ 与 density-<JobID>.log 分开保存，图和视频不会混合。
某模型失败时停止本批后续模型；已生成结果保留。相同命令重复执行返回原作业号。
已有独立任务不撤回、不重提；只把尚未提交的模型合并为一个新任务。
同模型不同检查范围也会检查已有批任务回执，避免共用输出目录的任务并发写入。
本次 7002/7003 的全量检查已分别开始，保留原作业 236608、236609；批量合并用于后续提交。

### 视频完整性与仅重建视频

服务器原 FFmpeg 对未指定时间戳的 PNG concat 输入会丢帧：602 张图曾只编码出 4 帧。
现改用临时的连续编号图片链接序列，25 fps，每张检查图保留一帧，不复制图片数据。
完成后用 ffprobe 验证视频帧数等于所选快照数，不能只凭 MP4 存在或作业退出码判断成功。
已有完整检查图时可以只重建视频，多个模型仍只提交一个 Slurm 作业：

```bash
python3 workflow.py isolation-check Fornax7002 Fornax7003 --name fd25-paf40-20261013 --video-only
```

这不重新计算密度。检查完成后仍使用 isolation-fetch 取回视频和统计结果。

### 正式轨道 IC 一次合并（2026-10-16）

使用 `workflow.py remote merge-ic MODEL --name NAME --submit`，完成后再 `remote simulate`。
读 [MERGE_IC.md](MERGE_IC.md)：它记录原程序的实际旋转/舍入/ID 规则与验证方法。
不要恢复两次手动合并，也不要合并矩阵后跳过旧中间舍入。保留原 .par 与 .sh。
源码新增 workflow_merge.py 和 ic_merge.c，两者必须一起安装到本地/服务器工作区根目录。
新作业冻结两者及配置；仅在逐粒子二进制校验、原子发布成功后删除自己生成的 merged.dat。
现存 IC 与有快照的模型不重建；合并作业不是新的模拟作业。
历史合并 .dat 已移入两端 observer_review_delete_20261006/merge-dat-20261016/，
恢复清单在 workspace/receipts/merge-dat-20261016/moved.jsonl。其他 .dat 保留。

### 检查完成后自动下载到对应本地模型（2026-10-16）

在本地分析工作区运行以下命令，批量检查仍只提交一个 Slurm 作业；完成后自动下载
每个模型的视频和统计表到 `<model>/isolation/`：

```bash
python3 workflow.py isolation-check Fornax7002 Fornax7003 --name fd25-paf40-20261013 --fetch
# 已经提交了检查或视频任务：只等待并取回，不重新提交
python3 workflow.py isolation-fetch Fornax7002 Fornax7003 --name fd25-paf40-20261013 --wait
```

`--video-only --fetch` 可用于仅重建视频后自动下载。`--timeout` 设置每个模型最多等待
秒数，默认 7200；本地命令需保持运行。终端断开或超时不会取消服务器任务，之后用
`isolation-fetch --wait` 继续取回。失败作业直接报错，不把旧视频当作本次成功结果下载。
全量检查在下载前再次验证视频帧数等于快照数，下载后比对视频大小和 SHA256；校验
通过才写入本地 `density-result.json` 完成回执。CSV 随同同步；逐帧 PNG 留在服务器用于生成视频，不回传，也不搬运快照。
7002、7003 的视频作业 236610 已完成；各 602 帧、25 fps、24.08 秒，已同步至本地。

本地回传仅允许当前已验证的视频文件名及 `density_evolution.csv`，直接写入 isolation/。
IC 文件及原建模诊断图保留；回执仍集中存放在 workspace/receipts/。
2026-10-16 已将 7002/7003 的视频和 CSV 上移，旧 density_output/ 中的逐帧图片移入本地回收站。
服务器 density_output/ 结构与检查程序保持原样。

### 正式模拟回传与本机运行（2026-10-16）

正式模拟也采用精简回传：fetch-model MODEL、pull --name RUN、collect --wait 都只取
顶层视频和 elinfo_<model>.csv，直接存入本地 <model>/，逐帧 PNG 留在服务器；回执在
workspace/receipts/production/。下载检查 SHA256 与视频帧数。不会再创建回传专用目录。
Canonical 观测同样遵守此本地路径规则，不把服务器的 output/analysis/<observer>/
层级复制到本地。observer_snapshot_*.json 存入回执的 <model>/<run>/observer/，
视频和 CSV 直接放在 Fornax 对应编号目录顶层。
旧 workspace/results/remote/ 是历史验证产物，不能随意将不同检查范围的同名视频互相覆盖。

本机原生 isolation 已配置，读 [LOCAL_ISOLATION.md](LOCAL_ISOLATION.md)，使用
isolation-local / isolation-local-status。本机源码和环境在 workspace/software/ 与
workspace/environments/，计算数据在 workspace/local-isolation/，和服务器结果分开。
原生构建使用同一套物理宏，应用已有 16 字节分配器对齐补丁；不替换服务器原二进制。
运行时使用短相对参数路径，旧程序 ParameterFile/InitCondFile 缓冲区只有 100 字符。
本机时间估计读取 progress.json 的实际积分速度；初始化阶段不得给出确定耗时。
7004、7005 本次用户指定 timescale_fd=16.6；服务器 profile 命名为 fd16p6-pafXX。

### 正式模拟的 HI 定义与旧结果迁移（2026-10-16）

用户明确选择 HI fraction = M_HI / M_star，可大于 1；分子、分母必须使用相同孔径。
求解器 NeutralHydrogenAbundance 是 n_HI/n_H，所以单粒子 HI 质量为
0.76 * m_gas * NeutralHydrogenAbundance；冷、热两种温度成分都要统计。
新 elinfo 包含 hi_mass、hi_half_mass、hi_fraction、hi_half_fraction、gas_mass、
gas_half_mass。旧 coldgas_* 质量列保留为修正后 HI 质量的兼容别名；hotgas_*
仍是历史定义，不能与 HI 相加充当总气体质量。重子理论弥散使用恒星加总气体。

GetInfo、PlotFig、CheckFig 和 IC 报告读取相同定义。旧 CSV 缺少新字段时应在
新的分析目录重算，不能更改列名后继续，也不能混用旧统计与新帧。既有分析和
release 保持冻结。真实验证为 7001 快照 100–102 的统计和绘图，以及快照 100
的独立质量/孔径计算：39 项无关统计与旧版完全一致。机器上的验证回执在
workspace/validation/hi-fix-20261016/，不会上传粒子数据或图片到 GitHub。

7004/7005 此次服务器反馈 16.6；7004 在 PAF=4.0/4.5/5.0 的短时失败已验证
并清理，PAF=6.0 的替代作业号 305984 保留；7005 为 PAF=4.0，作业号 305981。
这仅为历史记录，接手时必须重新查询状态，不能重复提交。
本机耗时估计的实测依据和局限见 LOCAL_ISOLATION.md。

### 日心观测与结果目录规范（2026-10-18，覆盖旧正式模拟回传规则）

新正式观测分析固定 `astropy_native_direct_v1`：原生银心轴直接对应显式
Astropy Galactocentric，太阳在负X，参数为R0=8.122 kpc、zSun=0.0208 kpc、
VSun=(12.9,245.6,7.78) km/s、GC ICRS=(266.4051,-28.936175)度、roll=0。
同快照MW H/B与单位各应用一次，恒星和冷热气体共用完整六维转换。
相机只换表示；逐星LOS、RA/Dec、距离与PM必须来自同一观察者。
这套方法本身不证明模型绝对六维已匹配真实Fornax。

当前显示合同为`heliocentric_tan_xi_eta_2p1_v1`，真实WCS-TAN ξ/η两轴各
±2.1度、等比例方形、北上东左、40×40分箱，有星即着色、色标±12 km/s；
日心自行箭头从(0,0)出发，球面方向长度先取tan后显示，不再除cosDec。
图内不标梯度PA或角动量轴。Walker61仍用真实球面角距与投影圆R<5 kpc，
不能按显示方窗裁切拟合。测量定义及字段见repo/SERVER_WORKFLOW.md。

固定cohort使用snapshot目录旁的
`analysis/astropy_native_direct_v1/cohort_ids.npz`（唯一uint64 particle_id），
或显式DSPH_PREEXISTING_ID_PATH。无文件直接报错，不猜成员。标记
selected_then_frozen，不冒称全部初始IC恒星。

物理CSV与observer JSON保存在对应快照目录旁的
`analysis/astropy_native_direct_v1/`；TAN图与视频使用其独立子目录
`display/heliocentric_tan_xi_eta_2p1_v1/`（帧再入`output_img/`），视频文件名
包含`tan_xi_eta_2p1_v1`，保留旧raw产物。prepare冻结cohort、observer模块SHA，
并另记display合同、目录及模块SHA；每帧sidecar绑定PNG SHA，显示manifest列出
本次帧/视频及物理身份。run/collect/fetch验证身份、SHA和明确视频列表，不从
磁盘glob取未声明TAN视频，也不回退旧顶层elinfo。回传视频和CSV直接写本地
`<model>/`；observer JSON进入`workspace/receipts/production/<model>/<run>/observer/`，
显示manifest进入相邻`display/`。逐帧PNG不自动回传。此处本地回传规则与
2026-10-20 WORKFLOW.md一致，替代旧的本地嵌套observer目录描述。
旧prepared分析保留原代码、旧路径并明确legacy身份，不重写冻结运行。

更改正式源码后仍须同步GitHub main、创建并激活新不可变release；只改editable
source不会改变workflow.json中已经选择的analysis_code。历史release、用户.sh、
模拟/隔离快照、IC与论文不修改。多模型验证仍合并一个Slurm、顺序执行。
