# Fornax 生产天空投影迁移交接说明

记录日期：2026-10-20。本文供继续维护 Linux—GitHub—服务器工作流的人或智能体使用，记录本次 TAN 绘图迁移的已完成内容、验收依据和客户端待办。它是交接记录；观测方法以 [生产分析约定](../repo/AGENTS.md) 为准，操作入口见 [工作流说明](WORKFLOW.md)。

服务器已发布 `tan-sky-20261020-v1`，5026/300 和 5028/285 的真实快照验证通过。用户实际 Linux 客户端尚未连接或安装，本机隔离环境中的安装器测试不能代替该步骤。论文及正式论文图未修改。

## 1. 版本与三端状态

| 对象 | 本次交付状态 |
| --- | --- |
| 已部署科学代码 | [`4ff4884a1a84fc98e834ecce6c5b9446ad951e4f`](https://github.com/Demon-Laplace/dSph_workbench/commit/4ff4884a1a84fc98e834ecce6c5b9446ad951e4f)；对比基线为 `e8f3bb76d01ab02611ab310e2c6053ad93871132` |
| 服务器活动 release | `tan-sky-20261020-v1`；发布完成时间为 2026-10-20 11:43:48 UTC |
| 服务器 `workflow.json.analysis_code` | `workspace/releases/tan-sky-20261020-v1/repo` |
| 服务器生产导航 | `/travail/xhuang/workspace/code/production` 指向上述 release 的 `repo` |
| 服务器可编辑分析源码 | `/travail/xhuang/dSph_workbench`；发布时其 59 个清单文件与新版本一致 |
| 服务器 Git checkout | `/travail/xhuang/workspace/code/github-sync`；发布时为干净的 `4ff4884` |
| 服务器根目录工作流 | 新版 `workflow_runtime.py`、`workflow_management.py` 已安装并回读 SHA |
| 用户 Linux 客户端 | **待安装、待实机核验**，见第 6 节 |

这些是发布时的状态。本文及后续文档提交会使 GitHub main 超过 `4ff4884`，不代表服务器科学 release 已随之更新。开始下一次修改前，应重新查询活动配置、源码、已安装副本和任务 manifest，不能只看 GitHub HEAD。

发布回执和备份在服务器：

```text
/travail/xhuang/workspace/validation/tan-production-20261020-v1/publication/
  publication_receipt.json
  backups/
```

该回执中的 `jobs_submitted=false`、`scientific_analysis_executed=false` 只描述**发布阶段**；此前真实快照分析确已在作业 236660 中完成。旧 release、冻结任务和历史产物保留。

## 2. 固定物理观察者，独立版本化显示

物理观察者保持 `astropy_native_direct_v1`。原生 MW 银心轴直接作为显式 Astropy Galactocentric 轴，`native_to_GC=I`。参数为 `galcen_distance=8.122 kpc`、`z_sun=0.0208 kpc`、`galcen_v_sun=(12.9,245.6,7.78) km/s`、GC ICRS `(266.4051,-28.936175) deg`、`roll=0`；太阳在原生负 X 侧：`(-sqrt(8.122²-0.0208²),0,0.0208) kpc`。

先以同模型的 MW 中心 H、整体速度 B 和单位换算形成银心六维，再做一次完整 Galactocentric→ICRS 六维变换，得到一致的日心 RA/Dec、距离、LOS 和自行；平移与整体速度各只扣一次。原生轴的无 R 点积另作 LOS 不变量核验，不能把已日心化的数据再次声明为 Galactocentric 输入。逐星 LOS 按唯一 uint64 ID 绑定。相机只改变表示，不能改变太阳、逐星 LOS、距离、成员或物理天东/天北。本次显示迁移没有重新选择太阳、翻轴、拟合滚转、主动注册模型或调整轨道。实现见 [canonical_observer.py](../repo/canonical_observer.py)。

新显示身份为 `heliocentric_tan_xi_eta_2p1_v1`，与物理 observer/config 身份分别保存。物理配置 SHA 为 `e7fd3ff559ba72ba8fed4bc4b149a403cbdd137a4ed4a7b986e5ee357ae9fde2`。该历史物理配置仍含旧显示字符串；不能凭这个字符串认定新图采用 raw RA/Dec，应以新的显示合同、源码 SHA 和逐帧记录判断。

| 项目 | 当前生产定义 |
| --- | --- |
| 天空图坐标 | 同一实际图心的 WCS-TAN ξ、η；各 ±2.1°，等比例方形、北上东左 |
| TAN 公式 | `ξ=(180/π)(n·E)/(n·n0)`，`η=(180/π)(n·N)/(n·n0)`；不是 raw ΔRA/ΔDec，也不是每轴单独 atan |
| 速度场 | 40×40 切平面分箱，有星即着色；星数等权平均冻结 LOS，减同目录 Walker 拟合截距 v0；色标 ±12 km/s |
| 自行箭头 | 日心，标 `Heliocentric proper motion`；从 `(0,0)` 出发，沿同一 E/N 方向；端点为 `(180/π)tan(0.65°)×unit_EN`，不再除 cosDec |
| 图内标注 | 不标速度梯度 PA 或角动量轴 |
| 面亮度 | 在 TAN 上重新分箱，使用逐像元精确立体角；保留原 M/L、零点和距离约定 |
| H I 与轮廓 | 仍按物理 kpc 像元面积算柱密度，保留冷/热气氢质量定义；±2.1 TAN 平面度对应物理半宽 `D0*radians(2.1)`，网格用 `rad2deg(x_kpc/D0)` 映到 TAN；椭圆和气体中心按同一映射核验 |
| 显示有效域 | 各成分均检查前向半球；显示窗不代替测量孔径 |

当前生产动力学成员是矮星系**当前 old+new 恒星**；光度参考仍为冻结 old cohort。Walker 理想模拟测量使用前向的 `Rproj=D0*tan(rho)<5 kpc`，真实球面 rho 以度计，北向东扫描 `0,3,...,180°` 共 61 个方向，各向拟合并选最大 `|k|`，单位 km/s/deg。这里等权、Pi=1、测量误差为零。2.1° 方窗、R5 圆孔径和有限抽样足迹是三个不同定义。

主要绘图实现为 [PlotFig.py](../repo/PlotFig.py)、[canonical_display.py](../repo/canonical_display.py) 和 [display_outputs.py](../repo/display_outputs.py)。此外，[snapshot_io.py](../repo/snapshot_io.py) 的临时标记改用分析 `TMPDIR`，无此变量时回退到工作目录 `.dsph_tmp`，以适配只读 release；物理数据加载定义保持。

## 3. 实际入口与输出位置

```text
Linux workflow.py remote pipeline MODEL
  → SSH 调用服务器根目录 workflow.py
  → 显式 --code，或服务器 workflow.json.analysis_code
  → prepare 冻结源码到 run/code、cohort 到 run/inputs，SHA 记录于 run/manifest.json
  → run 使用冻结的 GetInfo.py、PlotFig.py
  → collect / pull 校验 manifest、来源相对路径和 SHA，再回传成果
```

普通 `remote pipeline` 不会上传 Linux 本机的 PlotFig；本机 `analysis_code` 也不会替换服务器的配置。GitHub 源码同步、安装工作流模块、构建并激活 release、准备冻结任务，是各自独立的步骤。新 release 不会改变已准备任务的 `run/code`；从旧 raw 任务获得 TAN 图需要用新名称重新 prepare。

历史 `.sh` 直接引用 `../dSph_workbench/PlotFig.py`，走可编辑源码而非 release。本次该源码也已同步，现有用户 `.sh` 保留。

正常生产输出布局如下；独立验证使用第 5 节的新分析目录：

```text
<model>/output/analysis/astropy_native_direct_v1/
  CSV、cohort、逐快照 observer JSON
  display/heliocentric_tan_xi_eta_2p1_v1/
    output_img/check_nature_NNN.png
    output_img/check_nature_NNN.png.display.json
    display_manifest.json
    FornaxNNNN_tan_xi_eta_2p1_v1_v5.mp4
```

每帧 sidecar 和显示 manifest 记录物理观察者、cohort、显示合同、源码 SHA 及 PNG/视频 SHA、大小。完整 PNG 如果身份记录缺失或不匹配，不能作为续跑成功依据；新视频名称避免覆盖旧 raw 视频。视频合成已去除固定 3416×1886 拉伸；本次 3150×2250 PNG 输出为 3416×2440、SAR 1:1。

新版 [workflow_runtime.py](workflow_runtime.py) 和 [workflow_management.py](workflow_management.py) 识别这些路径：视频、CSV 直接回传本地 `<model>/`；observer JSON 和显示 manifest 分别放入 `workspace/receipts/production/<model>/<run>/observer/`、相邻 `display/`。逐帧图片留在服务器。历史冻结 manifest 没有显示身份时仍按历史 raw 规则处理；完整、显式命名的 info-only 任务有单独的 CSV 回传规则，不允许退回未版本化的默认 CSV。

## 4. 已完成验证及结果口径

实际分析作业为 **236660，mesopsl1，COMPLETED，退出码 0:0**，一个节点、16 核，两个模型顺序处理。验证比较上述基线和候选代码，并作独立数值重建。

| 模型 / 快照 | 当前成员 N | 冻结 old N | R5 N | Walker \|k\| (km/s/deg) | Walker v0 (km/s) | σ raw MLE (km/s) | σ after Walker MLE (km/s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5026 / 300 | 635802 | 607492 | 635689 | 3.977061 | 131.262677 | 8.905824 | 8.879907 |
| 5028 / 285 | 618981 | 605980 | 618033 | 4.829982 | 127.267294 | 8.299276 | 8.237061 |

表中 v0 是所选 Walker 方向的拟合截距；σ 是全部 R5 当前成员的 MLE，残差值为 `sqrt(SSE/N)`。它们不是论文 Fig.1(a) 的 age>2 Gyr、96 径向分箱结果，也不是历史 CSV 中 Re 孔径的平面扣除 `sigma`。

验收结果：

- 基线与候选按 ID 对比，恒星及冷/热气六维、LOS、R5 mask、61 向扫描、Walker v0 和日心自行的最大数值差均为 0。
- 独立 WCS 坐标误差小于 `5e-14 deg`；像元立体角相对误差约 `2.27e-12`，面亮度差小于 `1.7e-12 mag/arcsec²`。速度分箱的数量、支持域和均值一致；H I 物理面密度、同孔径质量比、轮廓、气体中心和自行端点均通过核对。
- 52 项针对性 Linux 检查通过，覆盖观测/显示、临时目录及工作流路径。Windows 完整 96 项测试未全部通过，不能宣称全套测试通过。
- 两模型多面板 PNG 已目视检查；各验证了一个单帧视频，`nb_frames=1`、SAR 1:1。未验证完整时间序列或 7004/7005 的 TAN 重绘。

本次没有重跑有限抽样、旋转搜索或 GIZMO 模拟。旧观察者下的概率不自动成为当前合同的结果。论文、正式图、旧数据产物和冻结任务保留。

## 5. 可追溯证据

完整独立 QA 位于各输入快照旁的新分析目录：

```text
/travail/xhuang/Fornax5026/output/analysis/snapshot_300/tan_migration_validation_20261020_v3/
/travail/xhuang/Fornax5028/output/analysis/snapshot_285/tan_migration_validation_20261020_v3/
```

每个目录包含 `independent_real_snapshot_QA.json`（输入身份/SHA 与科学对比）、`_validation/candidate/renderer_independent_QA.json`、`_validation/video_QA.json` 及版本化图像。v1/v2 的失败或部分结果保留为历史记录；**本次完整验收指 v3**。

| 文件 | SHA256 |
| --- | --- |
| 5026 `independent_real_snapshot_QA.json` | `7a612b73995fcbef9c8e72d866fa66077bec84f573435a23391935cc216cb663` |
| 5028 `independent_real_snapshot_QA.json` | `f4d0252e21651a662af222b444eaac5277cc96c9addee789f5cea098cac89f33` |

Windows 汇总与发布回执在工作区 `analysis/server_sky_projection_repair_20261020_v1/implementation/` 下的 `validation_summary.json`、`DELIVERY_REPORT.txt`、`publication/publication_receipt.json`。这些生成文件、原始数据、图片、日志和私有 `workflow.json` 未随本交接文稿上传 GitHub。

观测面亮度表由 PlotFig 源码目录的父目录解析，不是由当前 shell 目录解析。正常 prepare 会把源码放在 `run/code/`、输入表放在 `run/`；单独搭建验证目录时需保持该关系，不能靠更改观测数值修复路径。此次沿用表的 SHA 为 `71ee8205036a5f0ef6dfcd2e2e9467497d14a211bb680f286a7a067103824105`。

## 6. Linux 客户端接手步骤

服务器发布已完成。Linux 仍需把**实际执行 workflow.py 的目录**中的以下两个模块成对更新；只更新 checkout 中 PlotFig 或只切服务器 release，不能更新正在运行的本机 collect/pull。

| 安装文件 | `4ff4884` 对应 SHA256 |
| --- | --- |
| `workflow_runtime.py` | `dcbcc0d9e996c6dbf78737560c721371c77c4d430f635ebdeccf54b4f65c7b78` |
| `workflow_management.py` | `9e88005040028314a854bac3881a97250655bf41b35e3ef7482ecbdb1080dc0d` |

表中 SHA 对应 Git 中的 LF 文件及服务器安装字节；Windows checkout 的自动 CRLF 转换可能产生不同文件 SHA，安装时以核对过的 Git 版本或更新包为来源。

1. 在 Linux 先确认源码 checkout 与实际工作流根目录，检查 Git 状态和已安装模块差异，备份现有两个模块。正常 fetch 并快进 GitHub main，保留个人修改。
2. 从已核对版本的 `workflow/` 安装这两个文件到实际工作流根目录，回读 SHA；保留机器自己的 `workflow.json`、输入路径、模型数据和旧任务。如遇其他修订或个人改动，先比较再合并，不能直接覆盖。
3. 重新启动 collect/pull；旧进程不会热加载。取回已有 TAN 任务不需要重新提交 Slurm。若任务本身冻结的是旧 raw 源码，则另建新名称准备 TAN 分析。
4. 用指定 run 的 `pull --dry-run` 核对来源，再回传并检查视频、CSV 与两类回执。Linux 若直接运行本地分析，还要单独确认本机 `analysis_code`；更新客户端模块不会替它激活本地科学 release。

已另备更新包 `analysis/server_sky_projection_repair_20261020_v1/implementation/linux_client_update_4ff4884.zip`，包含只替换上述两模块的安装器。它先核对旧/新 SHA、备份、成对安装，未知版本保留并停止；已在本机隔离 fixture 验证，**尚未在用户 Linux 实机执行**。压缩包未上传 GitHub，GitHub 中的同版模块也可作为安装来源。

发布后还以默认入口准备了 `tan-entry5026-20261020-v1`、`tan-entry5028-20261020-v1`，核对其冻结代码和显示路径；两者**仅 prepare，未提交、未执行**，不能引用为新的科学结果。

## 7. 后续维护约束

继续维护时分别核对物理 observer、显示合同、cohort 和来源 SHA。新快照分析结果放在输入快照旁的新 analysis 子目录，保护快照与历史产物；冻结 release 和旧 run 不原位修改。此次 GitHub 交接文稿是文档提交，不触发新 release 或模拟。

后续若改变物理观察者、成员、测量孔径、抽样或噪声，应另起明确命名的分析并重新测量；纯显示修改也必须从同一逐星目录重新分箱并保留数值不变量。针对性回归入口见 [TAN 测试](../repo/tests/test_tan_display.py)、[显示身份测试](../repo/tests/test_display_outputs.py)、[临时目录测试](../repo/tests/test_snapshot_temp_path.py) 和 [工作流路径测试](tests/test_observer_paths.py)。
