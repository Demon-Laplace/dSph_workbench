# Canonical Fornax orbit and morphology candidates

本工具从相同冻结 old 恒星 ID 的原生银心六维，生成轨道速度、轨道平面及内部朝向的修正候选。输出是**尚未重新演化的候选**；不会启动 GIZMO、替换生产参数、修改太阳、修改旧快照或更新论文。

## 固定的物理输入

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

从含本包的 `repo/` 目录运行。先复制并核对 [config.example.json](config.example.json)，尤其模型、帧号、两个相邻交点区间和目标六维；示例帧号只对应本次 7004/7005，不是任意模拟的自动搜索设置。

```bash
python -m orbit_calibration_v1 --config my_orbit_config.json --stage validate
```

`validate` 只检查配置结构及声明的合同；正式执行还会核验实际共享观察者。`extract/solve/all` 均应在批准的计算节点中执行；两模型放同一 Slurm 作业、按配置顺序处理，不为每个模型单独提交任务。

```bash
python -m orbit_calibration_v1 --config my_orbit_config.json \
  --workspace-root /path/to/data --source /path/to/verified/repo \
  --run-directory /path/to/new/run --tag orbit-check-v1 --stage all
```

`extract` 读取快照并生成 `frame_state_index.json`；`solve` 使用同一 run 的经核验状态和缓存生成候选。为避免混用已有产物，物理输入、配置或代码改变时采用新的 tag/run。不要通过改 SHA 绕过拒绝，也不要给旧 raw 缓存补一个 canonical 标签。

每帧状态、恒星缓存及候选放在输入快照旁的新 `output/analysis/snapshot_NNN/<tag>/`。保存原始快照、cohort、配置、源码及缓存 SHA，raw H/B、单位、同成员 COM 六维和插值来源。ID 保持 uint64；缓存先尝试 float32 offset，并对 LOS、R5 ID、PA、q、半光半径与中心作回归，达不到门限时保留较高精度。原产物不覆盖。

## 旧入口迁移

历史 `calibrate_orbit.py`、公共 `analyze_orbit_orientation.py` 应只作新版 CLI 的薄适配器：使用 `DSPH_CANONICAL_ORBIT_CODE` 指向含本包的 repo，再传入新 `--config`。旧正 X 太阳、旧相机、旧 CSV 系数和配置不静默改译，也不提供自动执行旧备份的后门。

`orientation_Fornax5028_current/reproduce_density_PA.py` 只保留显式 `--legacy-observer-reproduction` 的历史重现；新目录与结果必须标 legacy，不能覆盖或认证为 canonical。该目录的 `original_core` 科学函数是历史 helper 的依赖，不能用新的薄适配器替换其函数定义。

## 验证范围

```bash
python -m unittest discover -s tests -p test_orbit_calibration_v1.py
```

本测试使用小型 native 矢量与合并数学，不需要 Astropy、真实快照或 Slurm。它覆盖配置拒绝、目标矢量、形态角符号、经验反解、合并偏置和矢量交点。正式交付还须完成实际 Astropy 观察者往返、当前 C merger 探针、真实快照提取/缓存回归及候选图检查；这些不同层级的验收不得互相冒充。
