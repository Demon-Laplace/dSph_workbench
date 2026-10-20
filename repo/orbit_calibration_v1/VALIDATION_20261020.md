# GitHub 同步前验证记录（2026-10-20）

测试在已快进到 `7c350f61b17f291f5372bd08bee1b7b8a71d4a7e` 的仓库基础上运行，
使用本机 Python 3.12 与已安装科学依赖。共享观察者合同 SHA 为
`e7fd3ff559ba72ba8fed4bc4b149a403cbdd137a4ed4a7b986e5ee357ae9fde2`。

| 检查 | 本次结果 |
|---|---|
| `tests/test_orbit_calibration_v1.py` | 16 项通过，包括实际启动两个 worker 子进程、矢量插值、经验反解、合并偏置与形态角符号 |
| `tests/test_orbit_frame_corrections.py` | 14 项通过，包括完整六维往返、24 组随机天空探针、太阳共速零点、混合参考系拒绝、三个兼容入口和显式物理候选开关 |
| workflow 轨道包同步白名单 | 1 项通过；源码/说明/config 示例被包含，快照、缓存、图片、候选参数及隐藏目录被排除 |
| 5028/285 真实气体保留相空间输入 | 14,011 粒子按同一 uint64 ID 对照；RA、Dec、距离、日心 LOS 最大绝对差均为 0 |
| 5028/285 真实旧恒星保留相空间输入 | 605,980 颗按同一 uint64 ID 对照；RA、Dec、距离、日心 LOS、两个日心自行分量最大绝对差均为 0 |
| 本机路径依赖 | 三个兼容入口在仓库外临时目录运行 `frame-check` 成功；临时目录未生成模拟或分析文件 |

真实输入回归调用本次 GitHub checkout 的 `canonical_observer.production_observables`，
先用同一配置 H/B 和单位把保存的原始相空间输入变为 MW 相对位置/速度，再做完整观测变换。
比较的是上一次 5028 原快照重算交付保存的目录，原始粒子和结果没有修改。

这是保留真实输入的回归，本次没有重新读取整份 HDF5、重新求 Walker、生成新的 H I 图，
也没有重演化轨道。用于审计的原始输入与参考目录 SHA、逐项日志和机器可读回执保存在
本机 `analysis/orbit_github_publication_20261020_v1/VALIDATION.json`；这些快照数据不上传 Git。

完整 workflow 测试依赖 POSIX `fcntl`，不能直接在 Windows 运行。本次只执行受改动影响的
同步白名单测试，导入时提供空的 `fcntl` 占位模块；未执行或模拟文件锁、Slurm 或服务器部署。
初次测试遇到 Windows 临时目录权限限制，改用工作区内临时目录后完成回归。
因此本记录不代表 Linux 全工作流、服务器已安装版本或新模拟的科学验收。

没有启动服务器作业、修改生产 `.par`、粒子 IC、论文或不可变 release。
