# 模拟变更记录

每次实际变更须在对应模型的 `MODEL_CHANGES.txt` 留简洁备注，并在同一轮工作同步到
`model_changes/` 的对应 TXT 镜像和本汇总；同一任务的相关改动合记，纯读取或查询不另记。
范围包括 IC、轨道、内部姿态、物理/数值配置、程序、分析、图件及运行启动、停止、重试、清理。
镜像保存完整备注，本表记录原因、改动、来源、验证状态及同一记录号和镜像链接。
IC 来源、轨道参考和运行重试分别注明，不能仅凭模型编号推断继承关系。
候选、已应用、启动通过、完整演化验证分开表述；状态只代表记录当时，不代表实时作业状态。
后续修正追加新行，不覆盖旧结论；未知来源写“待核实”。详细参数、日志和校验和留在原回执。

以下先补录已核实的近期变更，不作为全部历史模拟的完整谱系。

| 日期 | 模拟 / 版本 | 前序来源 | 原因与改动 | 验证状态 / 依据 |
| --- | --- | --- | --- | --- |
| 2026-10-20 | Fornax7004、7005 / `joint-gc-pa43-20261020` | 各自 isolation snapshot 50；替换各自 `orbit5028-20261017` 正式轨道 | 按用户新目标采用 `joint_GC_radius_target6D_PA43`，改变轨道及内部朝向；复用已验证源 IC，保留 feedback=16.6 | 已应用并启动；末态匹配尚待验证。[换轨记录](workflow/REORBIT_7004_7005_20261020.md) |
| 2026-11-07 | Fornax7006 / `anchor139p4-pose-density` 初始试验 | IC：自身 isolation snapshot 50；轨道：Fornax5028 `anchor139p4_orbitonly_20261102_v1`；姿态：snapshot 266 的 `pose_01` 静态候选 | 继承已指定质心六维，回推初始内部姿态以试验投影形态；feedback=16.6 | 条件姿态试验；后续启动版本见下一行，不能视为已演化形态匹配。回执索引① |
| 2026-11-11 | Fornax7006 / `paf80-20261111-v1` | 同一7006试验的失败启动版本（最近作业19062） | SSE 粒子数组对齐错误；PartAllocFactor 提高到8.0，保持同一二进制 IC、求解器和物理配置 | 作业19064已实际演化；未重新导出或合并 IC。回执索引① |
| 2026-11-14 | Fornax7006 / `fixed_program` 轨道修正候选 | 基线：上述作业19064的 snapshots 000–315；方法：既有 Fornax5028 轨道修正程序及 MW53 缓存 | 原轨道已过近心点（R_GC≈156 kpc），未达目标 R_GC≈141.44 kpc；固定比较时刻2.540276 Gyr，以 D_sun=139.4 kpc、v_LOS=55.3 km/s及指定天空位置/自行修正初始质心六维；内部姿态和太阳约定不改 | 原模拟按用户要求停止；候选闭合及步长检查通过，尚未写新 IC 或重跑。11-15已生成条件轨道预期图。回执索引② |
| 2026-11-14 | Fornax5029 / 双区 β 参数候选 | 原始 [Fornax5028 IC](IniFile/IC_Fornax5028.ini)，并非其正式末态 | 试验中心切向、外围径向偏置：r≤0.5 kpc 的 β=−0.2，1–3 kpc 的 β=+0.8；保留原恒星/气体结构、质量和粒子数 | 参数候选；[模型说明](IniFile/IC_Fornax5029.md) |
| 2026-11-15 | Fornax5029 / `Fornax5029-20261114T231856Z-38156f9b` 构建 | 上述双区 β 参数；同模型失败构建 | 旧入口不能读取 `two_zone`，且外区约束未归一化导致 β 偏离；接入专用求解器并归一化两区约束，目标和验收阈值不变 | 静态构建通过，共614,500粒子；不代表孤立平衡。[构建与验证](IniFile/IC_Fornax5029.md#使用与当前状态) |
| 2026-11-15 | Fornax5029 / `fd16p6-paf45-mesopsl1-20261115-fix1` isolation | 同一已验证5029 IC；前次 `fd16p6-paf40-mesopsl1-20261115-retry1` | GDB 确认 SSE 访问8字节错位；PartAllocFactor 4.0→4.5改变分配布局，保持原求解器、IC及物理参数（feedback=16.6） | 作业236841运行超过140秒且有非零积分进度，启动通过；完整演化与平衡未验收。回执索引③ |
| 2026-11-15 | Fornax 模板及全模型备注规则 / `FORNAX-NOTES-20261115` | 既有服务器与 GitHub 工作流约定 | 新建 Fornax/AGENTS.md；要求模型 TXT、GitHub TXT 镜像与本汇总同步记录每次实际变更 | 两份服务器 AGENTS 已与 main 逐字节核对；不改变模拟。[完整备注](model_changes/Fornax.txt) |

| 2026-11-15 | Fornax5028 / `orbit_PA43`，`F5028-PA43-20261115` | 自身孤立第50帧既有独立IC；原anchor139p4轨道的PA/面亮度候选 | 应用minor_axis_transport初始姿态和新质心六维，目标PA43并近似保持面亮度；保留feedback=24.9、3.6 Gyr终点；旧anchor139大文件按授权保留结构化结果后清理 | IC及既有合并校验通过，已提交但受账户CPU时数额度限制排队，末态未验证；额外成分检查程序按用户要求移除。[完整备注](model_changes/Fornax5028.txt) |
| 2026-11-15 | Fornax5029 / `orbit_PA43`，`F5029-PA43-20261115` | 用户选定自身孤立模拟snapshot_030；轨道/初始朝向参考本次5028 PA43 | 按自身源质心和bulk调整放置以共用目标质心六维；保留feedback=16.6、3.6 Gyr终点，原孤立模拟继续 | 既有合并流程完成并提交；受同一CPU时数额度限制排队，尚未启动积分，演化形态未验证。[完整备注](model_changes/Fornax5029.txt) |
| 2026-11-15 | Fornax7006 / 状态核查，`F7006-STATUS-20261115` | 上述fixed_program候选及旧模拟的退役回执 | 补记旧模型已退役，纠正“原始输出仍保留”的过期导航说明 | 当前无活动作业；轨道修改仍为候选和条件预期图，未写新IC或重跑。本次只核查并补记状态。[完整备注](model_changes/Fornax7006.txt) |

| 2026-11-15 | Fornax7006 / `orbit_corrected`，`F7006-CORRECTED-20261115` | 自身孤立第50帧既有独立IC；此前fixed_program目标；用户授权本地复算 | 用既有MW53固定场程序重新求初始质心六维，只改轨道；原内部姿态和完整物理/数值配置保持，feedback=16.6、PartAllocFactor=8.0；不声称逐字恢复已转走的旧候选 | 缓存、共享观察者、闭合/半步及既有一次合并校验通过；正式模拟已提交，受账户CPU时数额度限制排队，尚未积分，末态目标未验证。[完整备注](model_changes/Fornax7006.txt) |

回执索引（相对本地模拟工作区，不将大数据上传 GitHub）：

- ① `workspace/receipts/production/Fornax7006/fornax7006-anchor139p4-pose-density-paf80-20261111-v1/`：`inputs/case.json`、启动/取消回执。
- ② `Fornax7006/analysis/`：`s272_orbit_correction__fixed_program_candidate.json`、`s272_orbit_correction__fixed_program_verification.json`、`Fornax7006_expected_orbit_images.json`。
- ③ `workspace/receipts/isolation/Fornax5029/fd16p6-paf45-mesopsl1-20261115-fix1/startup-repair.json`。

追加格式：`日期 | 模型/唯一版本 | IC来源、轨道参考或前次运行 | 为什么改、改了什么 | 已验证到哪一步、依据`。

| 2026-11-20 | 共享 isolation 工作流 / `WORKFLOW-RUN-VERSIONS-20261120` | main `02da0be`；服务器 11-17 非例行核验接收路径 | 按用户要求，以声明的 IC/参数版本自动绑定运行并生成 JSON/TXT 台账；同版本复用、不同版本新建，不扫描粒子或做哈希认证 | 合成版本、并发和客户端/接收器测试；未改真实 IC 或提交模拟。[共享程序备注](workflow/MODEL_CHANGES.txt)、[使用说明](workflow/RUN_VERSIONS.md) |

| 2026-11-20 | Fornax8001 / `F8001-ROTATION-SEED-20261120` | 原Fornax5028 AGAMA参数；旧模型保持原样 | 用户指定8001进行恒星旋转修正测试，仅新增stellar.icbeta=0.0，保留密度、质量、ickappa及最终beta-second；INI静态检查发现当前96份旧恒星配方同类遗漏 | 参数候选已准备，尚未生成IC或演化；[完整备注](model_changes/isolation/Fornax8001.txt)、[参数排查](workflow/reports/IC_ROTATION_INI_AUDIT_20261120.md) |
| 2026-11-20 | 共享IC生成器 / `IC-ROTATION-SAMPLER-20261120` | main `8a8829f`；声明版本绑定规则 | 记录采样分支、提示缺失icbeta；新构建和正常目录登记不再计算粒子哈希，本机兼容新build_id及旧运行身份 | 43项定向单测通过；旧production-fetch断言不在本任务修正范围；[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-20 | Fornax8001 / `F8001-ROTATION-SEED-20261120` 进展 | 同任务首次构建v1；无旧模型变化 | v1因本地命令120秒执行上限中断，保留部分输出并改用新构建名v2；加入本征旋转科学测量及同作业内后处理入口 | 测量程序3项合成样本测试通过；完整IC仍待构建，不能将命令超时解释为物理失败；[追加备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8001 / `F8001-ROTATION-SEED-20261120` IC完成 | 5028解析配方+stellar.icbeta=0；build 4258fa52-a693-4112-99f7-08ca903f94ac | v2完整生成600000旧星+14500气体；最终IC旧星Vphi=5.556 km/s、V/sigma_phi=0.5951、顺行74.18%，r50相对解析目标差−0.114% | 仅IC实测确认，中心100pc及beta定义局限保留；到代码时间0.5的isolation已发起上传/提交，未提前宣称演化通过；[完整备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8001 / `F8001-ROTATION-SEED-20261120` 已提交 | 上述完整旋转IC；IC001/P001 | rotation-test提交JobID313990，单节点24MPI、TimeMax=0.5、feedback=16.6；同作业内自动输出分量旋转和结构演化统计 | 提交成功；启动、达到目标时间及平衡性未提前认定；[完整备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8001 / `F8001-ROTATION-SEED-20261120` 启动完成 | 同一Job313990及IC001/P001；未重新提交 | node0251运行约155秒，非零积分到代码时间0.00268555；结束启动检查、模拟继续至所设0.5 | 仅启动标准满足；演化保持性及平衡结论仍待自动科学后处理；[完整备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8002 / `F8002-COLD-DISK-20261120` 候选 | 用户授权的新冷盘测试；8001/5028质量及径向配方，旧模型不动 | 恒星q=.12、seed kappa=1、beta_z=.65；移除球坐标beta-second约束，细化薄盘网格；新增Q、圆度与径向冷度科学测量 | 60000点初步采样用于选参，不冒充完整IC；完整构建与隔离演化尚待执行；[完整备注](model_changes/isolation/Fornax8002.txt)、[物理方案](workflow/reports/FORNAX8002_COLD_DISK_DESIGN_20261120.md) |
| 2026-11-20 | IC默认图与冷盘诊断 / `IC-COLD-DISK-PLOTS-20261120` | model_tool现有build→workflow→build_ic链；旧funcini图未接入新构建器 | 恢复默认初态图，优化质量/密度及旋转弥散内容；新增Q、圆度、径向冷度/厚度及同作业内轻量启动检测 | 53项定向单测通过；IC图片与物理测试仍需实际执行；[程序备注](workflow/MODEL_CHANGES.txt) |
