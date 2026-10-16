# 正式轨道模拟：一次合并与二进制 IC

在本地分析工作区执行（也可在服务器省略 `remote`）：

```bash
# new-case 准备输入后，先预检，再提交同一个准备名称
python3 workflow.py remote merge-ic Fornax7004 --name orbit-v1
python3 workflow.py remote merge-ic Fornax7004 --name orbit-v1 --submit
# 合并作业成功后，按原 launcher 提交模拟
python3 workflow.py remote simulate Fornax7004 --name orbit-v1 --submit
```

`merge-ic` 是正式轨道模拟的准备步骤，不用于重新生成 AGAMA 模型，也不重跑 isolation。
原来的 `mergegalaxy.par`（或 `mergegalaxye.par`）与 `mergecube.par` 继续保存输入路径、
轨道、朝向和盒子平移；不再手动运行两次合并。`new-case` 仍负责复制模板参数，物理参数
应在提交前核对。已有用户 `.sh` 原样保留。

一个计算节点作业内：

1. 编译冻结的 `ic_merge.c`，一次读取 MW、dwarf、cube 三份原始粒子文件。
2. 直接写最终 `merged.dat`，不生成 `mergegalaxy.dat`。
3. 调用配置的原 `IC_ascii2snap`，生成 `GZWJL.PARAM` 中 `InitCondFile` 指定的二进制文件。
4. 校验头信息、全部数据块、所有粒子位置/速度/ID/质量及气体能量与密度块；记录 SHA256。
5. 在模型准备锁保护下原子重命名发布 IC，随后删除本次生成的 `merged.dat`；不保留二进制暂存副本。

当前参数习惯使用 `mergecube.IC_000`；不强制重命名为 `IC_000`，因此原模拟 launcher
可以继续使用。已存在的 IC、快照和 restart 文件受到保护，不覆盖。转换或验证失败时保留
本次中间文件，且不发布 IC。重新运行需检查失败原因并使用新 `--name`。
相同已提交名称只返回原回执，不重复提交；相同模型的活动准备任务会阻止另一次提交。

回执、冻结代码、日志位于 `data/workspace/merge-ic/<model>/<name>/`。
成功标志为 `completed.json`，包含 `all_particle_fields_verified`、二进制 SHA256 和
`intermediate_deleted: true`；不能仅凭转换器退出码判定成功。

## 为什么与原两步结果一致

以实际服务器 `ICmerger_matrixMW` 可执行文件为参考核对其行为，并通过原程序运行作对照。
未找到原 C 源码，新文件是按已核对行为独立实现。必须保留以下兼容规则：

- 原合并器用 float 运算和固定角度系数 **0.017453**，并非精确的 pi/180。
- 只对第二个输入按 X、Y、Z 顺序旋转位置与速度，再添加偏移。第一输入的三个角度在
  实际旧程序中没有应用。新流程保留这一行为；纠正它属于另一个会改变物理结果的任务。
- 两步流程中银河系和矮星系的数据会先经过 `%.7g` 输出再读回。新程序在内存中显式
  复现这次舍入，然后应用盒子变换；不能直接合成矩阵跳过舍入。
- 最终粒子按类型稳定排列，每种类型内为 cube、MW、dwarf；ID 分别偏移 0、Ncube、
  Ncube+Nmw，质量和额外气体字段沿用原规则。输入 ID 冲突或不合法类型会直接报错。
- 编译固定使用 `-ffp-contract=off -fno-fast-math`，不允许为提速开启 FMA/fast-math。

这些精度是原 IC 合并器本来使用的精度；没有把模拟求解器或原始粒子输入降精度。

参考程序 SHA256：

```text
ICmerger_matrixMW 8657a9c232e53809a81862b5651b62f2dca79d47b40ae95d5bf0ebe7105be2de
IC_ascii2snap     217b2c6928dbef1d2f0a0f59be66b192733b1c089d68a348a5920e181d817e5b
```

本地 3,136 粒子真实 IC 样本的最终 ASCII 和二进制均逐字节一致；20 组独立测试覆盖
Fornax7001/1071 轨道参数、任意三轴旋转/平移、混合粒子类型、无气体输入、非零气体能量
及乱序输入。回归命令（从源码 workflow/ 运行，二进制路径由本机提供）：

```bash
IC_LEGACY_MERGER=/path/to/ICmerger_matrixMW \
IC_LEGACY_CONVERTER=/path/to/IC_ascii2snap \
python3 -m unittest discover -s tests -p test_merge.py
```

生产规模对照入口为 `tests/validate_merge_reference.py`。必须在独立的计算节点 scratch
目录运行，传入原参数、原两个可执行文件和新 `--work`；它比较完整 ASCII 与二进制
SHA256，逐粒子验证二进制，成功后自动清除自己产生的所有临时粒子数据，只保留报告。
不要以模型目录作为 scratch，也不要为验证重提完整模拟。

## 配置与部署

服务器 `workflow.json` 的 `ic_merge` 配置包括 `converter`（参数列表，使用 `{ascii}`、
`{binary}`）、`setup`（模块环境）、`python`、`compiler`（默认 gcc）和 `sbatch`。
若未显式指定 converter，复用 isolation 默认 profile 的转换器。Slurm 资源必须显式配置；
此次服务器采用 mesopsl1、1 个 CPU、8 GiB、1 小时。参数和配置在提交时冻结，执行前后
校验输入文件身份；执行时记录原始输入哈希。禁止修改排队任务的输入。

部署时同步 `workflow.py`、`workflow_models.py`、`workflow_merge.py` 和 `ic_merge.c`；
C 源码必须位于 workflow_merge.py 同目录。不要只安装 Python 文件。

## 旧合并文件归档

2026-10-16 清单覆盖本地工作区和服务器 data/，跳过已有待删归档和符号链接。
服务器 5 个、本地 3 个旧合并 `.dat` 已按原相对路径移入各自的：

```text
observer_review_delete_20261006/merge-dat-20261016/
```

逐文件恢复位置可从 `workspace/receipts/merge-dat-20261016/moved.jsonl` 查询。
移动前检查文件 inode/大小/修改时间与活动 Slurm 作业。归档不会立即释放磁盘空间，
用户检查后删除归档才释放。Ewald 表和其他数据 `.dat` 不属于本次归档。

## 生产规模实测结果（2026-10-16）

服务器作业 236611 已完成：使用 Fornax7001 的全部 14,150,204 个粒子输入，新旧最终
ASCII 与二进制分别 SHA256 完全一致，所有二进制粒子字段复验通过。旧两次合并合计
165.26 秒，新一次合并 147.57 秒（单次实测，不作为通用性能承诺）。验证产生的临时
粒子文件已删除，原模拟目录和原始输入保持原样。

作业 236613 另验证了统一提交入口、冻结输入、原转换器、原子发布和自动删除全过程：
3,136 粒子结果与历史二进制逐字节一致。服务器不允许跨目录硬链接，因此发布使用
准备锁保护下的原子重命名；不要改回硬链接发布，也不要手工并发修改准备中的模型。
本地全部 65 项测试通过。机器可读证据见 [验证报告](reports/merge-onepass-20261016.json)。
