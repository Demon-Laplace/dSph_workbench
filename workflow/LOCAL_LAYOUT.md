# 本地模型目录（2026-11-20）

用户最终约定：初态和isolation结果放在一起，不分成ic和isolation两个目录。

```text
Fornax8003/
  isolation/
    Fornax8003.ini
    Fornax.ini
    parameters.ini
    build.json
    model_component_stellar.nbody
    model_component_gas.nbody
    model_potential_stellar
    model_potential_gas
    Fornax8003_components.png
    Fornax8003_central_components.png
    Fornax8003_column_density_profile.png
    Fornax8003_Mstar_Mgas_relation.png
    Fornax8003_kinematics.png
    ...
```

- 不再创建 `<model>/ic/<build>/`；isolation内不创建构建子目录。
- 原文件名尽量不变。不同版本重名时使用 `<build>__<原文件名>`，不覆盖旧版本。
- `workflow.py build-ic`、`model_tool.py build` 自动创建并发布到对应模型的
  `isolation/`，不需要另跑同步命令。成功时打印发布后的粒子 INI 和参数快照路径；
  已有无前缀文件仍属于旧构建，不作为最新版入口。
- 工具需要固定文件名和独立工作目录，内部链接视图放在
  `workspace/ic-builds/<model>/<build>/`；实际数据只在isolation存一份。
- `--ic-run BUILD` 仍用原构建名，build_id、随机种子和科学参数不变。
- 旧服务器运行、快照、restart、启动脚本及正式轨道模拟目录不搬动。

## 旧本地工作区合并试点（2026-11-21）

目标是日常只使用 `model/Fornax/`，每个模型的本地隔离产物只保存一份。
`model/IsolationModel/Fornax/Fornax100/` 已原位重命名至
`model/Fornax/Fornax100/isolation/`，旧路径保留为相对符号链接。
主目录下原有的正式模型粒子 INI 和快照保留，不与隔离文件混用。

这是单模型试点，不代表旧工作区已经退役。尚余 206 个实体模型目录（其中 17 个
在主工作区有同名模型）、208 份未跟踪的历史参数，以及旧笔记本和构建记录。
后续批量迁移先按模型检查目标路径、保留冲突版本，再原位移动并保留旧链接。
源码入口应最终归一到主工作区的 `tools/dSph_workbench`；切换前须保留历史参数和构建的访问。
不要直接把整个旧工作区链接到主根目录：旧入口下的 `Fornax.ini` 应继续解析到隔离输入，
而不是主模型目录同名的正式模拟输入。服务器和正在进行的运行不变。
