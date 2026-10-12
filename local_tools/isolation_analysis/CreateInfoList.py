#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_model_args(args):
    """解析命令行参数，例如 45, 47, 52-56"""
    models = []
    for arg in args:
        if '-' in arg:
            start, end = map(int, arg.split('-'))
            models.extend(range(start, end+1))
        else:
            models.append(int(arg))
    return models

def process_model(dw_name, modelnum):
    filename = f'./{dw_name}{modelnum}/IC_{dw_name}{modelnum}_evolution.csv'
    if not os.path.exists(filename):
        print(f"❌ 文件未找到: {filename}")
        return None

    # --- 读取头部信息 ---
    with open(filename, 'r') as f:
        header1 = f.readline().strip()  # e.g., "#Timescale_fd: 16.6"
        header2 = f.readline().strip()  # e.g., "#sfr_effiency: 0.015"

    try:
        timescale_fd = float(header1.split(':')[1].strip())
    except Exception:
        timescale_fd = np.nan
    try:
        sfr_eff = float(header2.split(':')[1].strip())
    except Exception:
        sfr_eff = np.nan

    feedback_value = timescale_fd / 8.3 if pd.notna(timescale_fd) else np.nan
    if pd.notna(feedback_value):
        feedback_str = (
            f"{int(feedback_value)} median"
            if feedback_value.is_integer()
            else f"{feedback_value:.2f} median"
        )
    else:
        feedback_str = "N/A"

    # --- 读取数据 ---
    df = pd.read_csv(filename, skiprows=2)

    # 删除 snapshot 和 hotgas 列
    drop_cols = [c for c in df.columns if "hot_gas" in c or c == "snapshot"]
    df = df.drop(columns=drop_cols, errors="ignore")

    # 重命名 cold_gas_mass_rhalf -> M_Gas_rhalf
    if "cold_gas_mass_rhalf" in df.columns:
        df = df.rename(columns={"cold_gas_mass_rhalf": "M_Gas_rhalf"})

    # 新建 M_Gas = cold_gas_mass + total_gas_mass，并删除原列
    if "cold_gas_mass" in df.columns and "total_gas_mass" in df.columns:
        df["M_Gas"] = df["cold_gas_mass"] + df["total_gas_mass"]
        df = df.drop(columns=["cold_gas_mass", "total_gas_mass"])

    # 定义时间点并插值到固定网格
    time_points = np.arange(0, 8.5, 0.5)  # 0, 0.5, ..., 8.0
    df = df.sort_values("time_gyr")
    df = df.groupby("time_gyr", as_index=False).mean()  # 去掉重复时间点
    df_interp = df.set_index("time_gyr").reindex(time_points).interpolate(method="index")


    # 计算 SFR (/1e9, 保留4位小数)
    sfr = []
    if "stellar_mass_5kpc" in df_interp.columns:
        stellar_mass = df_interp["stellar_mass_5kpc"].values
        for i in range(len(time_points)):
            if i == 0:
                sfr.append("-")
            else:
                delta_mass = stellar_mass[i] - stellar_mass[i-1]
                delta_time = time_points[i] - time_points[i-1]
                sfr.append(round((delta_mass / delta_time) / 1e9, 5))
    else:
        sfr = ["-"] * len(time_points)
    df_interp["SFR"] = sfr

    # 调整列顺序（仅保留存在的列）
    desired_order = [
        "sigma_x", "sigma_y", "sigma_z",
        "stellar_mass_5kpc", "M_Gas", "gas_fraction", "half_light_radius", "stellar_mass_rhalf",
        "M_Gas_rhalf", "gas_fraction_rhalf",
        "SFR",
    ]
    ordered = [c for c in desired_order if c in df_interp.columns]
    df_interp = df_interp[ordered]

    # 添加元信息
    df_interp.attrs["Model"] = f"{dw_name}{modelnum}"
    df_interp.attrs["Feedback"] = feedback_str
    df_interp.attrs["sfr_eff"] = sfr_eff

    return df_interp

# 简短表头名映射
SHORT_NAME_MAP = {
    "sigma_x": "σx",
    "sigma_y": "σy",
    "sigma_z": "σz",
    "stellar_mass_5kpc": "M_star",
    "half_light_radius": "R_half",
    "stellar_mass_rhalf": "M_star_rhalf",
    "M_Gas_rhalf": "M_Gas_rhalf",
    "M_Gas": "M_Gas",
    "gas_fraction": "f_gas",
    "gas_fraction_rhalf": "f_gas_rhalf",
    "SFR": "SFR",
}

def apply_short_names(df):
    """将列名映射为简短版本（仅映射存在的列）"""
    mapping = {k: v for k, v in SHORT_NAME_MAP.items() if k in df.columns}
    return df.rename(columns=mapping)

def format_numbers(df, threshold=6, precision=3):
    """根据位数决定是否使用科学计数法"""
    def smart_fmt(x):
        if isinstance(x, (int, float, np.floating)) and x != "-":
            int_part = str(int(abs(x))) if x != 0 else "0"
            if len(int_part) >= threshold:
                return f"{x:.{precision}e}"
            else:
                return f"{x:.{precision}f}"
        return str(x)
    return df.applymap(smart_fmt)

def save_to_pdf(df, dw_name="Fornax", modelnum=0):
    """保存单个模型结果为 PDF 文件"""
    model = df.attrs["Model"]
    feedback = df.attrs["Feedback"]
    sfr_eff = df.attrs["sfr_eff"]

    # 应用简短列名，并格式化数值
    df_short = apply_short_names(df)
    df_display = format_numbers(df_short, threshold=6, precision=3)

    # 输出路径
    folder = f"./{dw_name}{modelnum}"
    os.makedirs(folder, exist_ok=True)
    output = os.path.join(folder, f"{model}_Info.pdf")

    with PdfPages(output) as pdf:
        fig, ax = plt.subplots(figsize=(11, 6))
        ax.axis("off")

        # 标题（包含 Feedback 与 sfr_eff）
        title = f"Model: {model}\nFeedback = {feedback}, sfr_eff = {sfr_eff}"
        ax.set_title(title, fontsize=14, loc="left")

        # 表格：左侧第一列为 Time
        row_labels = [f"T={t:.1f}" for t in df.index]
        col_labels = ["Time"] + list(df_display.columns)
        cell_text = [[row_labels[i]] + list(df_display.iloc[i].values) for i in range(len(df_display))]

        table = ax.table(cellText=cell_text, colLabels=col_labels, loc="center")

        # 删除内部网格线
        for _, cell in table.get_celld().items():
            cell.set_edgecolor("white")

        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.2, 1.2)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"✅ 已保存: {output}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python3 CreateInfoList.py 45 47 52-56")
        sys.exit(1)

    dw_name = "Fornax"
    models = parse_model_args(sys.argv[1:])

    for m in models:
        df = process_model(dw_name, m)
        if df is not None:
            save_to_pdf(df, dw_name=dw_name, modelnum=m)



