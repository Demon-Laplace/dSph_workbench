# Canonical Fornax orbit and morphology candidates

程序与文档现纳入 GitHub main。入口为 `python -m orbit_calibration_v1`；
具体修改事项和验证层级见 [CHANGELOG.md](CHANGELOG.md)。此前“轨道程序仅保留本地”的约定已撤销。

独立真实演化可用 `python -m orbit_calibration_v1.rephase_ic --config CONFIG.json`
复用原始合并 IC，仅对全部矮星系粒子施加一致的初始平移和速度增量。配置显式绑定原 IC、
原 snapshot 000、冻结 cohort、候选绝对物理六维、原 cohort 初始 COM、单位及各文件 SHA。
目标目录与原目录必须不同且 IC/报告均不存在；按完整 ID 核对原初态并验证所有保护字段。
新初态的实际 COM、量化误差和新 IC SHA 写入报告。启动仍使用现有工作流，程序不提交作业。
原IC可能已有隔离阶段的恒星形成generation位；必须提供真实源ASCII粒子文件及历史ID偏移，
保留完整ID，不按连续编号或粒子数上界选择成员。原snapshot000的初始气体同步输出及周期卷回
可使用显式、独立审计过的微小数值容差（见修改记录）；恒星速度和全部质量仍精确核对，
原/新IC中未修改字段仍逐字节核验。此容差不能用于后续轨道终点的匹配。

2026-10-21：新增真实演化六维残差、共同时间搜索、有限差分初态响应与独立重演化验收。
操作与科学边界见 [真实演化反馈说明](FEEDBACK_GUIDE_ZH.md)。旧静态预览已移除目标质心速度置入；
候选构造完成不能作为演化匹配通过的证据。三个新阶段为 `evaluate/plan-feedback/update-feedback`。

本工具从相同冻结 old 恒星 ID 的原生银心六维，生成轨道速度、轨道平面及内部朝向的修正候选。输出是**尚未重新演化的候选**；不会启动 GIZMO、替换生产参数、修改太阳、修改旧快照或更新论文。

## 固定的物理输入

2026-10-20 参考系修正：目标必须显式声明 `target.input_frame=icrs_heliocentric`。
Table 1 的日心六维通过 `target_frame.py` 做一次完整反变换；Li (2021)
Table 2 的银心球坐标不能混入这一输入。派生银心数值使用当前太阳常数，
不能标成 Li Table 2 的原值。`--stage frame-check` 从实际共享观测入口往返
核对六个观测量，不读取模拟、不旋转模型、不生成候选。`validate` 也执行此核查。

候选求解须另加 `--physical-candidates`；它会改变物理初始条件，不是参考系修正。
历史配置、冻结源码及旧结果保留。旧 CSV 缺少 H/B 和观察者身份时不能作为新标定输入。

配置 schema 为 `canonical_orbit_recalibration_v1`，观察者固定为 `astropy_native_direct_v1`。原生 MW 银心轴直接对应显式 Astropy Galactocentric，太阳位于负 X：`(-sqrt(8.122²−0.0208²),0,0.0208) kpc`，速度为 `(12.9,245.6,7.78) km/s`。相机滚转不能代替物理轨道修正。

从原快照读取同一模型 MW 中心 H、bulk B 后，只做一次 `(p−H)L`、`(v−B)V`。抽取的 native `r_grf_kpc/v_grf_kms` 已完成这些步骤，不能再扣 H/B 或太阳。目标日心六维由同一固定映射反解为原生银心位置和速度，并通过实际 `canonical_observer.production_observables` 往返核验。

轨道标定使用 `old_full`：冻结 old cohort 的质量加权位置与质量加权速度，二者来自相同成员。`current_full/old_core5/current_core5` 另名保存供诊断。不能把旧 CSV 的 `distance_gal`、区域平均 `vr/vtheta/vphi`、R5 日心自行混成同一 COM 六维。

形态仍用冻结 old 恒星的 TAN 天东/天北坐标，质量乘 30 邻居局部密度加权协方差；`PA=(90−degrees(pa_math)) mod180`，标准北→东。示例的 **43° 是形态长轴 PA**，不是 Walker 梯度方向。完整日心自行、Walker 测量及显示继续使用共同观测规范；本工具不改这些规范。

## 交点与经验速度修正

主锚点在目标**银心半径**处。示例目标 `Dhel=139.4 kpc` 对应 `RGC≈141.440346 kpc`，二者不能替换。另保存日心距离交点作为敏感性对照。相邻帧按唯一 uint64 ID 及 native 三维位置/速度插值，再对矢量长度求根；不对 CSV 标量半径直接线性插值来冒充真实三维交点。

同模型初态与锚点给出：

```text
C_E = |v_anchor|² − |v_initial|²
k_L = |r_anchor × v_anchor| / |r_initial × v_initial|
```

再沿种子初态的径向/切向基底反解满足目标能量和角动量尺度的初始速度。这里 `C_E` 是速度平方差的经验标定，`k_L` 是实测角动量比；它们不是保证在新轨道上不变的动力学定律。

轨道平面的 proper SO(3) 注册是对模型的**物理候选变换**，不是更改观察者。银河系盘、热气及非球形引力场保持，因此原轨道的旋转不保证新轨道有相同动力学。内部 PA 修正也需传回初态，并以重新演化检验；静态图上达到 43° 不等于已证明未来仍达到 43°。

## 一次性合并与近心点边界

当前一次性 merger 的线性部分以 `Q` 表示外层 cube 旋转、`G` 表示 dwarf 旋转：

```text
r_initial = Q (G c_dwarf + placement_position − c_MW)
v_initial = Q (G b_dwarf + placement_velocity − b_MW)
```

四个输入中心/bulk 偏置均须保留；共同 cube 平移/boost 在 MW 相对量中消去。不能直接套用旧二次并合角度顺序，或假定所有输入质心和 bulk 严格为零。参数角的符号、顺序、`.017453` 常数、float32 运算及两次 `%.7g` 输出需按当前 `ic_merge.c` 核验；数值 helper 的往返测试不能替代实际 C 探针。

近心点未知不妨碍已覆盖的入落锚点标定，但本工具不把最后一帧的最小半径称为已测近心点，也不强制某个未知 pericenter。默认另检查当前最新完整快照及前一帧；如仍向内运动，结果必须保留“尚未观测到转折”的范围。近心点、轨道时相和 PA 的最终结果需要更长的原演化或新的完整模拟，而不是当前静态候选证明。

## 调用

Python 3.10 或更新版本。`frame-check/validate` 需要 NumPy、Astropy；
候选数学还需要 SciPy。真实快照提取复用当前生产分析环境（含 h5py、pandas、
Matplotlib、SciPy，以及生产源码实际导入的依赖），不能仅安装前两项便当作完整模拟环境。

最小的独立参考系核查不需要任何模拟文件：

```bash
cd repo
python -m orbit_calibration_v1 --config orbit_calibration_v1/config.frame_check.example.json --stage frame-check
```

成功时标准输出为包含 `state=PASS`、目标六维、派生银心六维和往返误差的 JSON。
可以把标准输出另存为报告；程序在此阶段不创建分析目录。

从含本包的 `repo/` 目录运行。先复制并核对 [config.example.json](config.example.json)，尤其模型、帧号、两个相邻交点区间和目标六维；示例帧号只对应本次 7004/7005，不是任意模拟的自动搜索设置。

```bash
python -m orbit_calibration_v1 --config my_orbit_config.json --stage validate
```

`validate/frame-check` 核验配置、完整 ICRS/GC 往返及实际共享观察者，不读取模拟。`extract/solve/all` 均应在批准的计算节点中执行；两模型放同一 Slurm 作业、按配置顺序处理，不为每个模型单独提交任务。

```bash
python -m orbit_calibration_v1 --config my_orbit_config.json \
  --workspace-root /path/to/data --source /path/to/verified/repo \
  --run-directory /path/to/new/run --tag orbit-check-v1 --stage all --physical-candidates
```

`extract` 读取快照并生成 `frame_state_index.json`；`solve` 使用同一 run 的经核验状态和缓存生成候选。为避免混用已有产物，物理输入、配置或代码改变时采用新的 tag/run。不要通过改 SHA 绕过拒绝，也不要给旧 raw 缓存补一个 canonical 标签。

| 阶段 | 输入与执行内容 | 是否求物理候选 |
|---|---|---|
| `frame-check` / `validate` | 配置、固定观察者、目标完整六维往返；允许 `cases=[]` | 否 |
| `extract` | 明确的模型、帧号、原快照、冻结 cohort、生产源码及合并参数；提取 native 状态 | 否 |
| `solve` | 同一 run 的冻结索引、缓存、原输入与源码 SHA；验证后求解 | 需要 `--physical-candidates` |
| `all` | 顺序执行提取和求解；默认阶段，但仍要求显式候选开关 | 需要 `--physical-candidates` |

`--workspace-root` 是含 `FornaxNNNN/` 的数据根；每个模型需有
`output/snapshot_NNN.hdf5`、`output/analysis/astropy_native_direct_v1/cohort_ids.npz`、
`mergegalaxy.par`、`mergecube.par` 和 `GZWJL.PARAM`。`--source` 必须指向一套完整且
通过核对的生产 `repo/`，不能只复制本包后使用缺少依赖的目录。

`config.example.json` 的 7004/7005 帧号是旧轨道历史示例。后续重跑状态另见
[重跑交接](../../workflow/REORBIT_7004_7005_20261020.md)；使用新模拟前必须重新确认
输入快照与交点区间，不可将示例帧号套到同名的新运行。

目标 `ra_deg/dec_deg` 为 ICRS 度，`D_kpc` 为日心 kpc，`vlos_kms` 为日心 km/s，
`mu_alpha_cosdec_masyr/mu_delta_masyr` 为日心 mas/yr（第一项已含 cos Dec），
`morphology_PA_E_of_N_deg` 为形态长轴角，单位度。示例是项目采用的参考值，
不代表把文献的全部参数和协方差完整移植。`physical_to_icrs` 是固定物理坐标映射，
不是可调相机；更改后实际共享观察者核查必须通过，否则直接拒绝。

每帧状态、恒星缓存及候选放在输入快照旁的新 `output/analysis/snapshot_NNN/<tag>/`。保存原始快照、cohort、配置、源码及缓存 SHA，raw H/B、单位、同成员 COM 六维和插值来源。ID 保持 uint64；缓存先尝试 float32 offset，并对 LOS、R5 ID、PA、q、半光半径与中心作回归，达不到门限时保留较高精度。原产物不覆盖。

| 主要输出 | 含义 |
|---|---|
| `state.json` / `stars.npz` | 每帧 native 状态、成员定义和缓存及来源 SHA |
| `frame_state_index.json` | 本次配置、源码身份、各帧输出位置与输入参数的冻结索引 |
| `PERICENTRE_STATUS.json` | 已采样径向转折区间；无转折不外推近心点 |
| `ORBIT_CANDIDATES_SUMMARY.json/.csv` | 跨模型候选摘要，保存在 run 目录 |
| `candidates/mergegalaxy.*.CANDIDATE.par` | 独立候选参数，保存在相关快照分析目录，不能当成已演化结果 |

常见拒绝信息分别对应缺少 `input_frame`、观察者/SHA 不符、GC 球量混入、
缺少显式候选开关、cohort 缺失，或冻结输入/源码变化。应核实来源并采用新的配置和 tag，
不要修改冻结报告的 SHA 来强行继续。

## 旧入口迁移

仓库内 `calibrate_orbit.py`、`analyze_orbit_orientation.py` 是新版 CLI 的薄适配器，
就地导入同级包。`python run_calibration_from_config.py my_orbit_config.json --stage frame-check`
提供位置参数形式的配置入口。仓库外的旧安装适配器可用 `DSPH_CANONICAL_ORBIT_CODE`
显式指向含本包的 `repo/`。旧正 X 太阳、旧相机、旧 CSV 系数和配置不静默改译，也不提供自动执行旧备份的后门。

`orientation_Fornax5028_current/reproduce_density_PA.py` 只保留显式 `--legacy-observer-reproduction` 的历史重现；新目录与结果必须标 legacy，不能覆盖或认证为 canonical。该目录的 `original_core` 科学函数是历史 helper 的依赖，不能用新的薄适配器替换其函数定义。

## 验证范围

```bash
python -m unittest discover -s tests -p test_orbit_calibration_v1.py
python -m unittest discover -s tests -p test_orbit_frame_corrections.py
```

第一组使用小型 native 矢量与合并数学，需要 NumPy/SciPy，不需要 Astropy、真实快照或 Slurm；
覆盖配置拒绝、目标矢量、形态角符号、经验反解、合并偏置和矢量交点。
第二组明确需要 Astropy，覆盖实际观察者与参考系语义及兼容入口。
实际 C merger 探针、真实快照提取/缓存回归及候选图检查属于不同层级的验收，不能由单元测试替代。
此次发布的具体测试结果见 [验证记录](VALIDATION_20261020.md)。
GitHub 源码更新不会替换服务器或 Linux 客户端的已安装副本、不可变 release 或已准备作业；
后续安装需连同共享生产源码核对，并重新选择版本化输出，保留历史结果。
