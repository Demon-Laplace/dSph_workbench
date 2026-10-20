# 从真实演化残差修正初态（2026-10-21）

本入口解决此前“银心总速接近，但日心 LOS 和绝对天区明显不符”的问题。所有观测量来自
固定 `astropy_native_direct_v1`。程序只读取演化结果、报告残差或输出新初态候选；不会提交
模拟、覆盖生产参数、改变太阳、旋转最终快照或给最终速度加常数。

## 本版的实质变化

旧候选生成仍可作为初始猜测。新版静态预览使用 `R @ r_measured` 和 `R @ v_measured`，
保留变换后的真实总速度；不再把预览的质心速度直接设为目标速度。
`rigid_preview_measured_COM` 保存实际预览六维，`static_target_velocity_injected=false`，
`evolved_6D_match_validated=false`。原候选 ID/文件名保留用于兼容历史来源，名称中的
`target6D` 不能作为演化匹配的证据。目标六维往返通过、候选构造完成、真实演化达标分别报告。
静态候选汇总的 `state` 改为 `UNEVOLVED_CANDIDATES_READY`，不再使用含糊的 `PASS`；
历史冻结报告保留原值。本版唯一的六维验收标志来自真实轨迹 `evaluate` 的逐项容差检查。

新增三个阶段：

| 阶段 | 输入 | 输出 |
| --- | --- | --- |
| `evaluate` | 冻结原生轨迹索引、目标、共同时间或搜索时间窗 | 各实际输出帧的六维残差、共同时间的验收结果 |
| `plan-feedback` | 同模型基准轨迹、明确步长与固定比较时间 | 12 个初态中心差分候选；等待真实演化 |
| `update-feedback` | 基准及 6 对真实扰动演化的冻结索引 | 有限差分响应、受步长限制的新初态；等待再次演化 |

更新采用真实演化测得的局部响应，不再假定旧轨道的经验能量/角动量系数在新轨道上不变。
它仍是局部迭代，不能保证全局最优。缺少扰动运行、时间支持或独立响应时直接拒绝更新。

## 残差、权重与时间

六维残差顺序固定为：

```text
sky_east_deg, sky_north_deg, D_kpc, vlos_kms,
mu_alpha_cosdec_masyr, mu_delta_masyr
```

天空误差采用目标方向处的球面 logarithmic map，其两分量模长等于大圆角距离。
它不是 raw ΔRA，也不是用于图片显示的 TAN。距离和各速度量由同一个原生质心六维生成。
轨道质心统一使用冻结 old_full 成员的位置/速度质量均值，不以 Walker 零点或 R5 自行均值代替。
因此验收不代表 Walker、形态 PA、气体分布或观测选择函数同时达标。

`residual_scales` 是用户明确给出的工程归一化尺度，`acceptance_tolerances` 是逐项容差。
示例数值只演示接口，不能当作测量误差、完整协方差、χ²、似然或统计置信水平。
没有容差时只输出测量，绝不认证匹配。

可指定一个 `comparison_time_Gyr`，在该同一时刻同时比较全部六项。相邻快照之间仅插值
同一成员定义的原生质心 r/v，再做完整观测转换；不逐项插值观测量，不跨缺失帧插值，不外推。
该线性插值没有给出时间分辨率误差界；严格验收应优先使用实际输出时刻或加密时间采样。

也可在 `feedback` 中用 `search_time_window_Gyr: [t_min,t_max]` 替换固定时间，按全部六项
归一化残差平方和遍历每个相邻区间、寻找共同时间。窗口必须已被真实演化完整覆盖。
部分覆盖不产生“最佳今日帧”；边界解有专门标记。选定时间后，先显式冻结该时间再生成响应
试验，不能让各扰动或各观测分量分别选择一个有利时刻。下一轮实际演化后可重新评估时间窗。

## 准备可复核的真实轨迹

复制 `config.feedback.example.json`，核对目标与尺度。提取阶段为每个模型写明确的
`initial/latest/extract_snapshots`，包括 snapshot 0、比较时刻的两个相邻端点，以及搜索窗
内所有所需输出。`extract` 可用于早期演化，不再要求尚不存在的目标半径交点。

```bash
python -m orbit_calibration_v1 --config my_feedback.json --stage extract \
  --workspace-root /data --source /verified/repo \
  --run-directory /data/workspace/orbit-feedback-extract --tag orbit-feedback-v1
```

启用 `feedback.freeze_physics_inputs=true` 时，额外冻结 MW、dwarf、cube 三份实际输入的
SHA256；共同文件在一次提取命令中只哈希一次。此项会读取大输入文件，应在计算节点执行。
跨多个模型的提取按原工作流共用一个检查作业、顺序处理，不能各自新增检查作业。

新的索引为每份 `state.json` 保存 SHA；状态还绑定原快照 SHA、原生缓存 SHA、输入参数、
源代码、cohort 及完整 uint64 old ID 集。读取时验证缓存中的质心与状态一致，成员和提取方法
跨帧一致。旧索引缺少绑定信息时必须重新核查提取；不能补标签或改 SHA 绕过拒绝。
反馈阶段验证冻结状态与缓存，原 HDF5 SHA 继承自提取；不声称每次反馈都重新读过原 HDF5。

正常运行的身份来自提取时捕获的 `case.json`，不能仅凭模型名和帧号对应。对于没有
`case.json` 的历史模型（如原 Fornax5026），可在新分析目录建立明确的身份 sidecar，
在该模型配置项中指定 `run_identity_file`；它不写回生产模型目录。sidecar 必须声明：

```json
{
  "identity_kind": "content_addressed_legacy_run",
  "run_id": "legacy-Fornax5026-actual-snapshot0-hash-prefix",
  "initial_raw_snapshot_sha256": "完整的实际snapshot000 SHA256",
  "mw_ic": "/actual/mw-input.ini",
  "dwarf_ic": "/actual/dwarf-input.ini",
  "cube": "/actual/cube-input.ini"
}
```

初始快照身份必须与提取到的 snapshot 0 一致，路径从实际原合并输入核实。
这只命名可追溯的历史运行，不补造过去的作业、生成史或验证结论。

## 验收与响应更新

```bash
python -m orbit_calibration_v1 --config my_feedback.json --stage evaluate \
  --index /actual/extraction/frame_state_index.json --model Fornax5026 \
  --run-directory /data/workspace/feedback-evaluate-v1

python -m orbit_calibration_v1 --config my_feedback.json --stage plan-feedback \
  --index /actual/extraction/frame_state_index.json --model Fornax5026 \
  --run-directory /data/workspace/feedback-plan-v1 --physical-candidates
```

六个初态控制依次为原生 MW 相对 x/y/z（kpc）、vx/vy/vz（km/s）。每次只改变一个控制的
正/负步长；保持 dwarf 内部朝向、cube 变换、三个固定 IC 输入、物理参数和 old ID 集一致。
输出 `.CANDIDATE.par` 只改初始 placement/boost，经过当前合并规则的线性质心核查。
实际粒子 IC 验证及正向演化仍是后续必要步骤；程序不启动这些运行。

基准加 12 个扰动演化到共同时间后，在独立 JSON 中列出六对索引：

```json
{
  "schema_id": "canonical_orbit_feedback_trials_v1",
  "pairs": [
    {"axis": 0,
     "minus": {"index": "axis0_minus/frame_state_index.json", "sha256": "完整索引SHA256"},
     "plus": {"index": "axis0_plus/frame_state_index.json", "sha256": "完整索引SHA256"}}
  ]
}
```

实际文件必须包含 axis=0..5 全部六对；不得把 7004 与 7005 等不同内部模型当成轨道微扰对。
每个运行须有独立 run_id，并验证实际 snapshot 0 的扰动与计划相符。

```bash
python -m orbit_calibration_v1 --config my_feedback.json --stage update-feedback \
  --index /actual/baseline/frame_state_index.json --model Fornax5026 \
  --feedback-spec actual_trials.json \
  --run-directory /data/workspace/feedback-update-v1 --physical-candidates
```

程序由实际初态差值及演化残差求局部 Jacobian，检测秩与条件数，再做有阻尼、受信赖步长
限制的更新。归一化初态步长的欧氏长度不超过 `trust_radius<=1`，避免向未测过的区域大步
外推。输出始终为 `AWAITING_FORWARD_EVOLUTION`；预测残差不是实际测量。

新候选经过独立演化后用 `evaluate --proposal /previous/FEEDBACK_RESULT.json` 核验。
程序检查实际初态、物理上下文、目标、共同时间与尺度，比较实际改善和预测改善。
实际变差会拒绝本次更新，并建议减小步长；不能凭预测通过就把候选升级为演化达标。
重复迭代时重新测量当前轨道附近的响应；不把旧 Jacobian 当作全局动力学定律。

## 验证边界

`tests/test_orbit_feedback.py` 覆盖真实输入格式、SHA/质心绑定、球面残差、时间覆盖、
共同时间搜索、缺少响应、秩亏、跨模型拒绝、受限更新及独立正向结果拒绝。
其中正向响应采用明确标为 synthetic 的自由飞行模型，只验证数值反馈流程。
真实缓存回归另检查 Fornax5026/300 的 607492 颗星；不冒称已完成新的 N-body 轨道标定。
