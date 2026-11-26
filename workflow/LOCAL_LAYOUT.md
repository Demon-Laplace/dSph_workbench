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
    Fornax8003_initial_diagnostics.pdf
    initial_figures.json        仅 plot-only 原地更新时记录新版定义和统计
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
- 初态诊断默认合为三页 PDF（空间/轮廓/运动学）。原地更新只替换派生报告，
  不重建粒子；成功后清理该前缀的旧诊断 PNG。旧 IC 版本和冻结构建记录保留。
- 旧服务器运行、快照、restart、启动脚本及正式轨道模拟目录不搬动。

## 旧本地工作区合并试点（2026-11-21）

目标是日常只使用 `model/Fornax/`，每个模型的本地隔离产物只保存一份。
`model/IsolationModel/Fornax/Fornax100/` 已原位重命名至
`model/Fornax/Fornax100/isolation/`，旧路径保留为相对符号链接。
主目录下原有的正式模型粒子 INI 和快照保留，不与隔离文件混用。

试点后完整迁移了剩余 206 个历史模型目录中的 205 个；另一个 `Fornax1004` 原目录为空，
保留在旧工作树归档中，主模型 `isolation/` 中现有视频不动。208 份仅本地的历史参数和
笔记本、旧构建等资料也已迁入。参数按原 Git 规则保留，不将旧编号参数新增到上游。

```text
model/Fornax/                              # 唯一本地工作目录
  workflow.py, workflow.json
  IniFile/ -> tools/dSph_workbench/IniFile/
  model_tool.py, makeini.py -> tools/dSph_workbench/ 中同名入口
  FornaxNNNN/
    isolation/                             # 初态、隔离结果，平铺
    analysis/                              # 正式分析（若有）
    …                                      # 正式结果、快照保持原位
  tools/dSph_workbench/                     # 唯一活动源码库和 IC 参数总表
  workspace/
    ic-builds/<model>/<build>/              # 指向平铺文件的工具链接
    legacy-isolation/                      # 旧路径视图；含历史笔记本及验证构建
    archive/legacy-ic-checkout-20261121/    # 原 Git 工作树、旧配置和回滚资料
    catalog/legacy-isolation-merge-20261121/migration.json
model/IsolationModel/Fornax -> ../Fornax/workspace/legacy-isolation
```

在 `model/Fornax/` 下，先用 `python3 workflow.py check-ini 8003` 核对配方；
要新生成一次，可用 `python3 workflow.py build-ic --model Fornax8003 --name build-v2`，
或改用 `python3 model_tool.py build 8003`，不要两条都执行。提交已有初态时
用 `workflow.py isolation MODEL --ic-run BUILD --name RUN --profile PROFILE --submit`；
不提交时可用 `--dry-run`。后续常规约 3 Gyr 的测试选用 `TimeMax=3.0` 的服务器配置。

不要直接把整个旧工作区链接到主根目录：旧入口下的 `Fornax.ini` 应继续解析到隔离输入，
而不是主模型目录同名的正式模拟输入。服务器和正在进行的运行不变。
