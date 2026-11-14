# Fornax5028：中心 −0.2、外围 +0.8 的双区 β 候选

[完整 INI](IC_Fornax5028_betaInner_m02_betaOuter_p08.ini) 基于原始
[IC_Fornax5028.ini](IC_Fornax5028.ini)，保留恒星与气体结构、质量及粒子数，
将恒星运动学改为两区约束。它是**尚未实际运行构模的参数候选**。

## 密度与数值参数

AGAMA 单位为 **kpc、km/s、10¹⁰ M☉**；恒星和气体各为一个物理成分，无暗物质成分。

| 参数 | 恒星 | 气体 |
|---|---:|---:|
| 密度模型 | 单 Sérsic | 盘 |
| 总质量 | 6×10⁷ M☉ | 1.45×10⁸ M☉ |
| `sersicIndex` | 0.6 | 1 |
| `scaleRadius` | 0.54 kpc | 0.9 kpc |
| 形状参数 | `axisRatioZ=0.55` | `scaleHeight=-0.45` kpc，沿用原参数符号 |
| 输出粒子数 | 600,000 | 14,500 |
| 名义单粒子质量（等质量采样） | 100 M☉ | 10,000 M☉ |
| 轨道库大小 `numOrbits` | 30,000 | 10,000 |
| 积分时间 `intTime` | 100×各轨道的 `Tcirc` | 100×各轨道的 `Tcirc` |
| 轨道初值抽样 `ickappa` | 0.70 | 0.85 |

`ickappa` 仅用于轨道库初值抽样，不是最终旋转速度或旋转能比例；气体另保留
`icbeta=0`。两种势均为轴对称 `CylSpline`（`mmax=0`、30×25 网格）：恒星
R=0.02–20、z=0.02–6 kpc，气体 R=0.02–40、z=0.02–12 kpc。
恒星密度目标为 `DensityClassicLinear`，r=0.01–10 kpc、24节点、`stripsPerPane=2`；
气体密度目标为 `DensityCylindricalLinear`，完整网格表达式保留在 INI 中。

原 INI 的 `trajsize=400/1000` 字段保留，但当前原求解器及本兼容读取器未将它传入
`agama.orbit`，不能把它当作已生效的轨迹采样分辨率。

## 双区各向异性

本体静止系下的球坐标二阶矩定义为

\[
\beta=1-\frac{\langle v_\theta^2+v_\phi^2\rangle}{2\langle v_r^2\rangle}.
\]

| 真实三维半径 r | 约束 | 含义 |
|---|---:|---|
| 0–0.5 kpc | β=−0.2 | 轻微切向偏置 |
| 0.5–1 kpc | 自由过渡 | 不直接指定 β，不预先保证平滑 |
| 1–3 kpc | β=+0.8 | 强径向偏置 |
| r>3 kpc | 不直接约束 | 不能外推为 β=+0.8 |

```ini
beta=two_zone
kinemDegree=2
kinemgrid_inner=0, .1, .2, .35, .5
beta_inner=-0.2
kinemgrid_outer=1, 1.4, 2, 3
beta_outer=0.8
```

两区使用**同一个恒星成分、同一组非负轨道权重**，分别建立 `KinemShell` 目标。
约束作用于二阶 B 样条加权的径向、切向二阶矩，不是每个半径点精确恒定的 β，
也不是扣除每壳平均流后的方差 β。求解后须逐基函数检查实际 β；中心约束按动能矩
尺度逐行归一化，避免小尺度约束因绝对量级过小而被忽略。

## 使用与当前状态

原 `schwarzschild.py` 只支持单个标量 β。本配置需使用配套的
[schwarzschild_two_beta.py](../local_tools/ic_generation/schwarzschild_two_beta.py)，
以及其同目录的原求解器。在装有 NumPy、支持轨道权重优化的 AGAMA 环境中，
从仓库根目录运行：

```bash
python local_tools/ic_generation/schwarzschild_two_beta.py \
  IniFile/IC_Fornax5028_betaInner_m02_betaOuter_p08.ini \
  --output /path/to/simulation-workspace/Fornax5028/ic/beta_m02_p08 \
  --seed 42
```

输出目录必须尚不存在；应位于模拟工作区。这里直接使用双区读取器，
`build_ic.py` 尚未集成本配置的双区入口。`beta=two_zone` 会使旧单 β 读取器明确报错，
避免静默漏掉外围约束。兼容读取器检查非负权重、总质量、密度残差及两区逐基函数 β；
未达标则退出，部分输出不能视为完整成功的 IC。

本次发布只包含配置、读取器及软件检查。**没有生成实际粒子IC，没有确认该密度下
β外围=+0.8的合格平衡解，也未验证孤立或轨道演化稳定性。** 气体热状态、feedback
与外部轨道由后续模拟配置指定。

这份配置保留的是原始5028的单 Sérsic与初始气体盘；它与此前6.69×10⁷ M☉、
s282气体近似的静态模型不同，不能沿用那组静态 σ 数值或承诺复现其末态轮廓。
初始输入 β、孤立交接帧的实测 β(r)、正式末态 β(r) 应分别记录。
