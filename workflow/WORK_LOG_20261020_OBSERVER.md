# 工作日志：旧太阳观测代码修正（2026-10-20）

本次更新的范围是旧太阳观测处理的入口和坐标约定。
最初交付时轨道分析程序仅保留本地；用户随后明确要求将更新程序同步 GitHub。
当前程序、配置示例、测试及说明文档纳入 main；对应快照的数据、报告和图片继续保存在分析目录。
新轨道候选只以简明交接数据记录在 [下一步工作](NEXT_WORK_7004_7005_ORBIT.md)。

## 修改与核查

- 旧轨道工具默认入口改走 `astropy_native_direct_v1`，不再默认执行正 X 太阳、旧速度符号与旧相机映射的混合处理。
- 统一太阳位置为 `(-sqrt(8.122²-0.0208²),0,0.0208) kpc`，太阳速度为 `(12.9,245.6,7.78) km/s`。
- MW 中心和整体速度各扣一次；原生银心六维经一次完整转换产生日心 RA/Dec、距离、LOS 和自行。相机 R 只负责置中，不改变固定物理观察者的 LOS。
- 旧配置不自动改译成新配置。历史几何重现保留原始代码与备份，必须显式标为 legacy；不能再充当当前分析的默认观测方法。

本仓库相应的共享实现为 [canonical_observer.py](../repo/canonical_observer.py)，
调用入口为 [GetInfo.py](../repo/GetInfo.py) 和 [PlotFig.py](../repo/PlotFig.py)。
正确的日心六维实现已在提交 `3921bb0` 引入；本轮核对后保留其公式，没有再次换太阳。
显示采用 [canonical_display.py](../repo/canonical_display.py) 中的真实 TAN 切平面，
详见 [既有 TAN 迁移说明](TAN_SKY_MIGRATION_20261020.md)。

## 验证与交付

- 7004/291、7005/306 的真实快照复核及补充最新完整帧检查已完成。
- 38 项相关测试通过；实际一次性合并 C 程序的 4608 个 float 字段按位比较通过。
- 分析作业 `236671`（mesopsl1）完成，退出码 `0:0`。未启动新的 GIZMO 模拟。
- 新候选没有应用到生产 `.par`、IC 或活动模拟，论文未修改。
- 结果目录为 `<model>/output/analysis/snapshot_NNN/orbit_correction_7004_7005_20261020_v2/`；小型结果已回传对应 F 盘模型目录。

此前提交 `0ccde06` 曾额外纳入轨道分析实现。本次按用户澄清，用后续提交从主分支撤出这些实现及详细分析说明；本地副本与结果保留，既有 Git 历史保留。

## 后续范围变更：轨道程序纳入 GitHub（2026-10-20）

上述撤出是当时的历史决定，已被用户最新指令取代。此次重新纳入的是完成完整六维
参考系修正后的程序，不是直接恢复 `0ccde06` 的旧版本。入口为
`repo/orbit_calibration_v1/`，兼容命令为 `repo/calibrate_orbit.py`、
`repo/analyze_orbit_orientation.py` 和 `repo/run_calibration_from_config.py`。
使用方法见 [程序说明](../repo/orbit_calibration_v1/README.md)，修改事项与验证范围见
[修改记录](../repo/orbit_calibration_v1/CHANGELOG.md)。本次 GitHub 发布不启动模拟或更改论文。
