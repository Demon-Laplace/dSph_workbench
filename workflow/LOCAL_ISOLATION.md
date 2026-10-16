# 本机 isolation 运行

本机可运行与服务器相同物理配置的 isolation。当前独立安装复用系统 GCC、OpenMPI、GSL、
HDF5，Python 环境在 `workspace/environments/isolation/`，不更改 IC 生成环境。
只读取并复制服务器非周期 isolation 求解器源码；服务器原安装与既有用户 .sh 均不修改。

## 当前配置与命令

```bash
python3 workflow.py isolation-local Fornax7004 \
  --ic-run Fornax7004-20261016T185019Z-fc66b3c6 \
  --name fd16p6-native8-20261016 --processes 8
python3 workflow.py isolation-local-status Fornax7004 --name fd16p6-native8-20261016
```

启动器后台运行，终端关闭后继续。相同名字与相同配置不会重复启动；失败或完成后不会
盲目覆盖。修改配置或 IC 需要新名称。`--dry-run` 只预览，`--time-max` 可显式设置短时
验证范围；未传时沿用模板的 TimeMax=6。正常完整运行仍以快照达到 TimeMax 为准。

本机配置在不进入 Git 的 `workflow.json` 的 `isolation_local` 中：模板、二进制路径、
MPI 参数列表、Python 路径、进程数及参数覆盖均显式配置。当前 feedback=16.6、
PartAllocFactor=4.0、BufferSize=120，MaxMemSize=2048 MiB/进程。
默认 4 进程；本次实际为 8 进程，以低调度优先级运行，避免占满全部桌面资源。
本机 30 GiB 内存，8 个进程的内存上限合计约 16 GiB；实际 RSS 以运行时采样为准。

计算文件与快照放在 `workspace/local-isolation/<model>/<name>/`，与服务器结果分开。
其中 request.json 记录冻结文件 SHA256；status.json 记录运行/失败/完成；progress.json
每约 15 秒记录积分进度，并在至少两次有效采样后给出按近期积分速度外推的剩余小时数。
初始化不适合估计整个模拟时间，早期时间步也可能与后期不同，不能用单核频率推断总耗时。
本机输出不能用服务器 `isolation-fetch` 覆盖；该命令仍专用于服务器成果。

## 本机构建记录

来源是服务器 `gizmowjl_MFM_gas_Nstar100` 的只读源码副本，保存于
`workspace/software/gizmo-isolation/source/`；编译配置宏与原隔离求解器相同。
本机原生编译命令记录在 `workspace/software/gizmo-isolation/build.json`：

```bash
make -j6 SYSTYPE=Local CC=mpicc CXX=mpicxx \
  OPTIMIZE='-O2 -g -fcommon -ffp-contract=off' \
  MPICHLIB= GRACKLELIBS= GSL_INCL= GSL_LIBS= \
  HDF5INCL='-I/usr/include/hdf5/serial -DH5_USE_16_API' \
  HDF5LIB='-lhdf5_serial -lz'
```

该旧源码的分配器只按 8 字节取整，本机 GCC/SIMD 在读取粒子时会因此段错误。
本机副本应用了仓库已有的 `solver_patches/gizmo_allocator_alignment16.patch`，
四处分配/重分配统一取整为 16 字节；只改变内存布局，不修改物理方程或数值类型。
原始源码归档 SHA256、编译器/MPI 版本、配置宏和最终二进制哈希保存在 build.json。
配置文件和程序文件名使用短相对路径，绕开原程序 100 字符路径缓冲区限制。

跨编译器、MPI 实现和进程数的演化结果不承诺逐位相同。保留原 IC 与物理参数，分别
检查快照和密度演化；不能把“运行正常”直接等同于科学平衡性验证。

## 首次实测与时间估计（2026-10-16）

7004 本机 8 MPI 进程运行到模拟时刻 0.00357056，累计计算约 1400 秒；同模型
服务器 24 MPI 进程约 169 秒。本机日志中密度计算约占 47%，密度负载不均等待
约占 46%。当前配置没有表现出比服务器更快；单核频率不能代表整套模拟吞吐量。

扣除初始化后，用已完成 7002/7003 的同阶段速度和完整耗时校准，估算本机
7004 全程约 33–36 小时，暂用 24–60 小时范围。模型、feedback、后期时间步及
桌面负载都可能改变该估计；这不是运行时限或完成保证。证据在本地
workspace/receipts/local-isolation/Fornax7004/runtime-estimate-20261016.json。
progress.json 的近期线性外推不包含后期加速，不应把初始化附近的数百小时
直接当作完整运行预测。
