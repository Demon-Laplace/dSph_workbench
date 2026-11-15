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

回执索引（相对本地模拟工作区，不将大数据上传 GitHub）：

- ① `workspace/receipts/production/Fornax7006/fornax7006-anchor139p4-pose-density-paf80-20261111-v1/`：`inputs/case.json`、启动/取消回执。
- ② `Fornax7006/analysis/`：`s272_orbit_correction__fixed_program_candidate.json`、`s272_orbit_correction__fixed_program_verification.json`、`Fornax7006_expected_orbit_images.json`。
- ③ `workspace/receipts/isolation/Fornax5029/fd16p6-paf45-mesopsl1-20261115-fix1/startup-repair.json`。

追加格式：`日期 | 模型/唯一版本 | IC来源、轨道参考或前次运行 | 为什么改、改了什么 | 已验证到哪一步、依据`。
