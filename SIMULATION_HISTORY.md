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
| 2026-11-20 | Fornax8001 / `F8001-IC-PLOTS-20261120` | 既有完整旋转IC build4258fa52；Job313990不动 | 按用户反馈补齐并优化5张初态图，保存于本地isolation；共用model_tool下游构建器已恢复默认出图 | 真实8001图片已生成；没有改IC、重读演化快照或重提模拟；[完整备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8002 / `F8002-COLD-DISK-20261120` IC完成及提交 | build27ace52b；上述q=.12、kappa=1、beta_z=.65配方 | 614500粒子与5图已生成；主体84%旧星V/sigma_R=2.30、V/sigma_1d=2.55、顺行96.98%，但Qmin=1.14未达预设1.2；提交Job313997（IC001/P001）到代码时间.5 | 初态4/5指标满足，保留暖中心与边缘Q限制；不声称稳定性或演化通过；[完整备注](model_changes/isolation/Fornax8002.txt)、[设计与实测](workflow/reports/FORNAX8002_COLD_DISK_DESIGN_20261120.md) |
| 2026-11-20 | Fornax8002 / `F8002-COLD-DISK-20261120` 启动通过 | 同一Job313997、IC001/P001，未重提 | 同分配内自动检测在150秒记录非零积分Time=.000732422，自行退出检测、保留模拟继续 | 仅启动通过，未达到.5或确认冷盘稳定；结束本次监测；[完整备注](model_changes/isolation/Fornax8002.txt) |
| 2026-11-20 | IC面亮度绘图 / `IC-FORNAX-STELLAR-SB-20261120` | 用户明确不加WLM；既有Fornax V波段观测表与分析M/L=2.6 | 轮廓图分质量/光度面板，加入Fornax观测误差及既有Muñoz标定点；空间图明确整体与中心尺度；未调IC/参数、未拟合归一化 | 58项定向测试及两模型真实出图通过；[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-20 | Fornax8001 / `F8001-STELLAR-SB-20261120` | 既有build4258fa52、Job313990不动 | 仅刷新3张初态图，加入恒星mu_V观测而非WLM，澄清整体/放大与成分/投影含义 | 图已更新、旧图备份；不是演化拟合结果；[完整备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8002 / `F8002-STELLAR-SB-20261120` | 既有build27ace52b、Job313997不动 | 同样更新初态观测面亮度参照，保持固定M/L与原始观测标定，未改冷盘IC | 图已更新、旧图备份；未重新运行或查询演化；[完整备注](model_changes/isolation/Fornax8002.txt) |
| 2026-11-20 | Muñoz发表轮廓 / `IC-MUNOZ-PUBLISHED-FIT-20261120` | 原刊ApJ860,66表1/5/6与图6；既有Yang参考保留 | 用n=.71、theta_e=18.4arcmin、mu0V=23.59重建完整发表Sersic拟合；明确非新增bins、外推与归一化误差范围，统一D139.4显示角半径 | 63项定向测试通过；[共用程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-20 | Fornax8001 / `F8001-MUNOZ-PROFILE-20261120` | 同一旋转IC、Job313990不动 | 只刷新面亮度轮廓，加入Muñoz发表拟合和来源，不重标定观测或拟合M/L | 真实图已更新、旧图保存；[完整备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8002 / `F8002-MUNOZ-PROFILE-20261120` | 同一冷盘IC、Job313997不动 | 同样补充完整Muñoz拟合，保留非独立标定/外推说明，物理参数不变 | 真实图已更新、旧图保存；[完整备注](model_changes/isolation/Fornax8002.txt) |
| 2026-11-20 | Fornax8003 / `F8003-LOW-MSTAR-20261120` | 8002初态配方、build27ace52b；不是演化快照 | 用户确认只将恒星质量6e7降至2e7 Msun，气体1.45e8及q=.12等其余配方、粒子数保持；重新构建势与两种分量速度，不改旧模型 | 参数比较仅恒星质量一项变化；完整IC已开始构建，冷度/Q/稳定性待实测；未提交演化。[完整备注](model_changes/isolation/Fornax8003.txt) |
| 2026-11-20 | Fornax8003 / `F8003-LOW-MSTAR-20261120` IC完成 | 同一配方，build20d56906、seed42；非演化快照 | 已完成614500粒子与5张默认初态图；按用户最终约定合并到本地Fornax8003/isolation平铺 | IC构建及出图完成，总表4349bc9已发布；Q/冷度及演化稳定性未验收，未提交模拟。[完整备注](model_changes/isolation/Fornax8003.txt) |
| 2026-11-20 | 本地目录 / `LOCAL-FLAT-ISOLATION-20261120` | 原model/ic/build与model/isolation并存的布局 | 按用户明确要求合并8个模型的136个文件到model/isolation平铺；重名版本用文件前缀保留，工具链接集中到workspace；后续构建同样采用此布局 | 38项定向测试、真实8003路径dry-run通过；不重算IC、不读粒子做认证、不移动服务器运行。[说明](workflow/LOCAL_LAYOUT.md)、[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-20 | Fornax8001 / `F8001-AVAILABLE-DENSITY-20261120` | rotation-test及其P002续跑的000–304 | 沿用原密度/视频脚本，与8002共用Job236887；不中断模拟 | 视频与CSV已回传，范围约到2.97 Gyr；不代表3 Gyr完成或稳定性通过。[备注](model_changes/isolation/Fornax8001.txt) |
| 2026-11-20 | Fornax8002 / `F8002-AVAILABLE-DENSITY-20261120` | cold-disk-test及其P002续跑的000–225 | 同一Job236887顺序处理，共226帧 | 视频与CSV已回传，范围约到2.20 Gyr；不代表演化完成或冷盘稳定。[备注](model_changes/isolation/Fornax8002.txt) |
| 2026-11-20 | 共享工作流 / `ISOLATION-INPLACE-PARTIAL-20261120` | 原版本绑定、二进制restart及密度分析流程 | 新增原目录续跑，保留版本历史；现有分析入口支持未完成模拟的固定快照范围及模拟运行中回传 | 65项定向测试通过，无额外真实演化；既有分段作业不动。[说明](workflow/ISOLATION_CONTINUATION.md)、[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-21 | Fornax8001正式模拟 / `F8001-S001-7006-THIN-20261121` | 自身rotation-test的snapshot001；7006当前orbit_corrected轨道及朝向、同一MW/背景 | 按源当前恒星质心换算轨道偏移；astro_thin 4节点96MPI、3800M/核，MaxMemSize3800，保留反馈16.6，终点3.6 Gyr | IC准备236889完成，正式Job314007已提交，启动待检测；热能导出不等于继承独立反馈库。[完整备注](model_changes/Fornax8001.txt) |
| 2026-11-21 | Fornax8002正式模拟 / `F8002-S001-7006-THIN-20261121` | 自身cold-disk-test的snapshot001；同一7006修正轨道及朝向、MW/背景 | 同一资源与物理约定；保留原孤立运行 | 首次314008节点启动失败且挂起，无快照；原样重发为314009，启动待检测，不宣称冷盘稳定。[完整备注](model_changes/Fornax8002.txt) |
| 2026-11-21 | Fornax8001正式启动 / `F8001-S001-7006-THIN-20261121` | 同一Job314007，未重提 | 运行275秒、Sync-Point16且代码时间5.61792e-05，自动监测退出 | 启动通过，实际总粒子14144922；详情已回传，不代表长期科学验证。[备注](model_changes/Fornax8001.txt) |
| 2026-11-21 | Fornax8002启动失败 / `F8002-S001-7006-THIN-20261121` | 同一Job314009，原IC与资源不动 | 24分13秒后MPI通信失败；仅t=0快照、无restart。定位并回传实际Slurm日志，增加模型根目录入口 | TimeBinCount同步涉及rank58；无明确OOM或core证据，根因未定。按用户要求先排查、不再重提。[备注](model_changes/Fornax8002.txt) |
| 2026-11-21 | Fornax8002继续排查 / `F8002-S001-7006-THIN-20261121` | 同一失败日志及实际GIZMO二进制，只读分析 | 确认8字节分配器与16字节SSE要求冲突；从分配表重建8001的P对齐、8002偏8字节，可解释两者启动差异 | 首要嫌疑，尚缺rank58栈最终验证；PartAlloc8.5仅为未执行规避建议，3800MB上限/现有IC与作业均未改。[备注](model_changes/Fornax8002.txt) |
| 2026-11-21 | Fornax8003 / `F8003-RES100-ISOLATION-20261121` | 8002解析配方；8003保留原build20d56906，采用新版build08392e71（seed42） | 2e7 Msun恒星以200000粒子保持100 Msun/粒子；恒星轨道库20000、气体14500不变。按8002原始cold-disk-test配置提交314028（IC001/P001，TimeMax=0.5） | 新版IC及图已完成，作业已提交并开始运行；非零积分及演化稳定性尚未验收。[完整备注](model_changes/isolation/Fornax8003.txt) |
| 2026-11-21 | Fornax8003 / `F8003-ISOLATION-3GYR-20261121` | 同一build08392e71；前序314028退出码11，仅t=0，无restart，保留原目录 | 从零运行至3 Gyr（代码时间3.0681360790075956）；PAF=5.5作对齐规避，反馈16.6不变；本地isolation增加res100初态链接 | 314029（IC001/P002）运行至代码时间0.00299623；未达到启动观察窗或目标时间。[完整备注](model_changes/isolation/Fornax8003.txt) |
| 2026-11-21 | 共享初态生成 / `IC-PUBLISHED-PATH-20261121` | 查阅7004/8002/8003发布记录，8003新版build08392e71已自动平铺 | 生成器打印最终发布路径而非工作视图；后续常规约3 Gyr采用代码时间TimeMax=3.0 | 31项定向测试通过；未改IC与作业314029。[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-21 | 本地工作区 / `LOCAL-ISOLATION-CONSOLIDATION-PILOT-20261121` | 旧IsolationModel/Fornax/Fornax100隔离数据，非主工作区正式模拟输入 | 15份隔离文件原位移动到主Fornax100/isolation，旧路径改链接，无数据复制 | 旧新路径指向同一文件，正式输入和快照保留，主工作流可发现；未批量迁移或改服务器。[布局](workflow/LOCAL_LAYOUT.md)、[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-21 | 本地工作区 / `LOCAL-ISOLATION-CONSOLIDATION-20261121` | 原Fornax100试点、旧IsolationModel/Fornax工作区 | 迁入205个旧隔离目录、208份历史参数及旧资料；Fornax1004空目录归档，保留主目录视频；统一源码和命令入口，旧根仅作兼容 | 319个链接及迁移元数据、旧构建定位、主入口dry-run与32项定向测试通过；无重建、覆盖或服务器变更。[布局](workflow/LOCAL_LAYOUT.md)、[程序备注](workflow/MODEL_CHANGES.txt) |
| 2026-11-21 | Fornax8004 / `F8004-MSTAR4E7-IC-20261121` | 8003新版解析配方及build08392e71，非演化快照 | 恒星质量增至4e7 Msun、400000粒子，保持100 Msun/粒子；气体、几何和轨道采样不变，重新构建两分量速度 | build e913d247 的414500粒子及5张图完成；未提交演化，冷度/Q/稳定性未验收。[完整备注](model_changes/isolation/Fornax8004.txt) |
| 2026-11-21 | Fornax8002 / `F8002-ISOLATION-RESUME-T6-20261121` | build27ace52b，Job314006在cold-disk-test-IC001-P002-01的二进制检查点，旧日志末时刻4.0 | 使用现有流程原地续跑，`TimeMax=6.0`；单独保留上一段末帧，不重建IC或改物理 | Job314033已提交（IC001/P003）；配置终点非实际终止实测。[完整备注](model_changes/isolation/Fornax8002.txt) |
| 2026-11-21 | Fornax8004 / `F8004-ISOLATION-T6-20261121` | 自身完整build e913d247，初态t=0 | 用现有IC运行至配置`TimeMax=6.0`；反馈16.6，初始PAF6.5作布局规避 | Job314032已提交并进入RUNNING（IC001/P001）；未认定完成或稳定。[完整备注](model_changes/isolation/Fornax8004.txt) |
| 2026-11-21 | Fornax8002 / `F8002-REMOVE-PRE-RESUME-SNAPSHOT-20261121` | Job314006在续跑前另存的末帧 | 按用户要求仅删除`prior-output-job314006/snapshot_306.hdf5` | 确认文件已移除且不在314033输出目录内；其余快照及作业不动。[完整备注](model_changes/isolation/Fornax8002.txt) |
| 2026-11-21 | Fornax8001 与 Fornax5028 / `F8001-VS5028-IC-20261121` | 归档5028原始IC与8001完整build4258fa52，非演化快照 | mesopsl1 单作业236890做配方、轨道、质量结构、速度矩及投影对照；5 kpc恒星Vφ 0.013→5.556 km/s、顺行49.98→74.18%，恒星半质量半径相差约0.20% | Slurm退出0:0，报告和表图已回传；未修改IC、启动演化或验证长期稳定性。[完整备注](model_changes/isolation/Fornax8001.txt) |
