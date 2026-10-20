# 轨道修正入口迁移（2026-10-20）

7004/7005 的轨道与恒星长轴方向现在使用
[`repo/orbit_calibration_v1`](../repo/orbit_calibration_v1/README.md)。
`repo/calibrate_orbit.py` 和 `repo/analyze_orbit_orientation.py` 为兼容入口，
共同转交同一 CLI。它不自动启动 GIZMO，也不替换现有 `.par`、IC 或快照。

## 原因与边界

旧轨道工具采用正 X 太阳、旧速度符号、固定相机映射及旧标量目标。
这些约定与已经确定的 `astropy_native_direct_v1` 不同。
不能仅修改 PA 的符号或把旧配置改名后继续使用。

新方法固定太阳、MW 中心/bulk、显式单位和原生六维转换；相机只负责显示。
PA 是恒星长轴从天北向天东的角度，不是速度梯度方向。
当前项目目标为 RA=39.9971°、Dec=-34.4492°、D=139.4 kpc、
vlos=55.3 km/s、日心自行=(0.384,-0.364) mas/yr，用户指定形态 PA=43°。
完整目标必须同时通过实际 canonical 六维入口验收。

CE/kL 是当前模拟的经验响应，不是精确守恒量。主要标定点使用目标六维对应的
相同银心半径；相同日心距离仅为敏感性对照。
两端必须来自同一模型同一套质量加权质心六维，不把 CSV 的逐星球坐标均值当质心速度。
近心点尚未出现时仍可生成经验候选；程序不外推近心点，不把末帧最小半径标为近心点。

当前服务器一次性合并程序仍逐步模拟历史两次变换及中间舍入。
新工具使用实际 C 的角度顺序、负号、外层变换和中心/bulk 偏置反解新放置量。
改变内部朝向时，要同步补偿输入质心偏置，才能保留要求的初始 COM 轨道。
具体参数是静态及经验候选，必须通过新 IC 全粒子验证和重新演化才能认定科学结果。

特别地，改变初始径向／切向速度也会改变到达同一半径时的累计轨道方位；
即使球形静态势且 CE/kL 精确，沿用旧轨道锚点求出的 P 也不是新的边值问题解。
下一轮须重新积分到目标半径、更新轨道内相位与 P。
精确目标六维对应 old_full COM 单点，不等于 R5 目录平均自行、平均 LOS 或 Walker v0 全部精确命中。
合并偏置是由实际 snapshot-0 反推的有效 old-tracer 偏置，不能称为独立测量的整份原始 IC 质心。

## 代码冻结与旧入口

`workflow_runtime.prepare` 和 GitHub 源码导出共同允许并核验
`orbit_calibration_v1` 的递归源码；每个冻结文件进入 SHA 清单。
不把其数据、缓存、生成图片或运行输出放入 Git。
参见 [源码包冻结规则](SOURCE_PACKAGE_FREEZING.md)。

Windows 研究归档中的两个旧公共入口转交新包；原始程序完整备份。
`orientation_Fornax5028_current` 的旧几何重现保留科学公式，默认拒绝；只有
`--legacy-observer-reproduction` 可在独立新目录输出显式 legacy 结果。
旧缓存和旧候选不认证为当前观察者的结果。

分析结果保存在输入快照旁的新 `output/analysis/snapshot_NNN/<tag>/`。
本次验证与回执根目录为服务器
`workspace/validation/orbit_correction_7004_7005_20261020_v1/`。
实际部署和验证状态以该目录的提交、作业与发布回执为准；本文不代替作业成功证据。

Linux 客户端需先更新 GitHub main 并安装更新的 runtime/github 模块和新分析 release。
既有已准备作业继续使用其冻结代码；不得为迁移而修改或重启活动模拟。
