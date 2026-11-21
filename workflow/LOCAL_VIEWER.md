# 本地 Fornax 模拟查看器

## 启动

在 `/home/x-huang/simulation/model/Fornax/` 用同时装有 `numpy`、`h5py`、
`plotly` 的桌面 Python 运行（打开选取窗口另需 Tkinter）：

```bash
python3 visualize.py
```

本机也可打开工作区的「Fornax 模拟查看器」桌面启动器。
选择或输入模型编号、选择一个或多个本地快照/粒子初态，点击“生成并打开”。
列表只查模型根、`isolation/`、`output/` 及本地轨道分支；没有快照也可选粒子
INI，或手动指定文件。没有粒子文件也可查看登记的配方摘要。头部会显示来源路径、
时间、各类型计数、已登记的初态构建等。
多选文件生成一个含帧滑块的网页（默认最多 8 帧）。网页内置 Plotly，可离线打开；
不启动 HTTP 服务、不上传数据、不提交模拟。默认输出在 `analysis/YYYY-MM-DD_3D_<模型>/`。

无桌面时可直接导出：

```bash
python3 visualize.py Fornax8004/isolation/Fornax8004.ini \
  --units fornax --axis short \
  --output analysis/2026-11-21_3D_Fornax8004/viewer.html

python3 visualize.py Fornax81/isolation --skip-latest --frames 4 \
  --units fornax --output analysis/2026-11-21_3D_Fornax81/viewer.html
```

支持单文件和常规分片 GIZMO HDF5、工作流的十列 ASCII 粒子 INI；不支持
求解器 restart、Gadget 二进制或参数 INI。结果不会覆盖既有 HTML。

## 看什么

- 鼠标旋转和缩放；快捷切正视/侧视；开关气体、旧星、新星等。
- 按速度分量、速度大小、带符号切向速度和高度着色；查看质量投影、旋转与弥散曲线，
  并读取球孔径内的质量、半质量半径、厚度、顺行比例及单位质量角动量。
- 拖动 3D 图或切正/侧视时，右侧按当前相机方向正交投影，重新计算恒星二维面亮度
  和圆环径向轮廓。使用恒定 `--mass-to-light`（默认 V 波段 2.6）、选中的恒星类型
  和可调的独立恒星样本 `--profile-points`；样本总质量归一到该组孔径总质量。
  使用 `--units native` 则显示恒星面质量，不把原生单位误标成物理面亮度。

固定轴的质量投影与动力学统计使用整个选定球孔径，不用显示点数计算。
`--units fornax` 明确按 kpc、km/s、\(10^{10} M_\odot\) 初态/快照单位解释，
不根据 `HubbleParam` 或红移自动乘除；如果文件使用其他换算，可显式设置乘数。
默认 TP2 作定心和定轴参考，正式轨道或包含银河系时请核对并按需指定 `--center`
和 `--bulk`。此工具的空间孔径不代表束缚成员；视角变化不是日心坐标变换，也不用于
观测拟合、H I 估计或内部姿态搜索。
