# 自动 IC／参数版本与孤立运行绑定

2026-11-20 用户确认：版本判断交给程序，不依赖人或智能体查看 TXT，
不以粒子 ID、SHA、文件大小比对或全量回读作为正常操作的门槛。
本次接入 `isolation` 的本机客户端、服务器接收器及其新建运行；
不更改正式轨道、成员选择、绘图行为或既有冻结运行。

## 台账与版本的含义

服务器 `isoFornax/<model>/MODEL_VERSIONS.json` 是唯一权威台账，
`MODEL_VERSIONS.txt` 是从同一记录生成的易读视图。台账由程序加锁、
原子更新，不需要手工维护；TXT 写入失败不改变已经登记的运行绑定。
同一模型的 IC、参数、运行编号及作业状态放在一起，不新建逐文件核验报告。

- `IC001` 等是该模型内的 IC 来源版本别名。AGAMA 已有的 `build.json.build_id`
  作为声明的构建版本；本程序不重读粒子文件来认证它。
- `P001` 等绑定实际生效的 profile、模板参数、参数覆盖及执行/资源配置。
  只读取小型配置；生成的 InitCondFile、OutputDir、ResubmitOn 不造成参数版本变化。
  选择版本时保存模板和配置快照，上传期间修改模板不改变已经选择的版本。
- 首次使用的 `--name` 保持原名；同一名字下更换 IC 或配置后，自动使用
  `<name>-IC002-P003-01` 这样的新名字。重发相同请求会找到同一运行，
  不会每次再追加一个编号。
- 版本比较不保证绕过入口的手工文件改动可被发现。不要将新内容写入旧构建目录后
  继续宣称它是同一构建；生成新 IC 使用新的构建目录，保留旧文件。

## 使用

```bash
# 首次准备，未加 --submit 不调用 sbatch 或 sbatch --test-only。
python3 workflow.py isolation Fornax9001 --ic-run ic-v1 --name iso

# 明确提交这个版本。重复相同请求只返回原作业。
python3 workflow.py isolation Fornax9001 --ic-run ic-v1 --name iso --submit

# 换 IC；自动分配新的运行名，旧 IC、快照、脚本和作业保持不动。
python3 workflow.py isolation Fornax9001 --ic-run ic-v2 --name iso --submit
```

输出明确列出实际运行名、IC/参数版本、是否复用以及真实的作业号。
同一版本的失败运行仍返回失败状态，不擅自启动重试；明确重试时指定新的 `--name`。
结果取回及密度检查使用输出中的实际运行名。本机
`isolation-status --model MODEL --name 原请求名` 可查询本机最近一次成功请求所选的运行；
传实际运行名则直接查询对应运行。纯状态查询不额外写版本日志。

`--existing` 导入的目录没有独立构建记录时，以声明的来源作为版本。
如果在同一导入位置更换了 IC，应显式给出新的来源版本：

```bash
python3 workflow.py isolation Fornax9001 --existing --ic-version imported-v2 --name iso --submit
```

这不是粒子 ID，也不是要求人工核验数据；它是导入这份新 IC 时声明的版本名称。
程序不会为推测导入目录是否被修改而重新读取粒子。

## 历史兼容与并发

- 旧运行的来源及小型配置已明确且与请求相同时，可直接复用，仍标
  `legacy-unversioned`，不改旧 request.json、不补造版本、不重提交。
- 明确换了 IC，或旧运行来源不明而本次明确指定了输入时，自动创建新运行；
  不把旧结果认作新版本，也不因缺少历史字段阻挡新请求。
- 保留现有模型目录、写入边界、已清理运行名以及活动/已有输出保护。
- 同时请求同一版本时，由模型版本锁分配一个运行，再由运行锁保护提交；
  已处于 submitting 且提交响应丢失时，仍暂停自动重提，不能凭版本号猜测作业未创建。
- `reserved` 仅表示已选版本，`prepared` 表示已准备输入；
  `submitted` 必须有实际 sbatch 返回的作业号。`complete` 仅表示求解器进程正常退出，
  不声称达到目标终点、平衡、收敛或科学验收。

安装时将 `workflow_versions.py` 与 `workflow_isolation.py` 一起放到工作区根；
CLI 选项来自同版本的 `workflow.py`。新建的本机/服务器运行会一并复制所需模块，
不得向已经提交或冻结的运行补装。

测试（合成输入，调度器替身，不加载真实快照或提交 Slurm）：

```bash
cd workflow
PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest -v \
  tests.test_run_versions tests.test_isolation
```
