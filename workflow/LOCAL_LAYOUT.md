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
