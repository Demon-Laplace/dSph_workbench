#!/usr/bin/env python3
import os
import sys
import configparser
import pandas as pd
import numpy as np
from odf.opendocument import OpenDocumentSpreadsheet
from odf.style import Style, TableCellProperties, TextProperties
from odf.table import Table, TableRow, TableCell, TableColumn
from odf.text import P

ini_model_file = './modelfile'

# --------------------------
# Reference RPF orbit dataset
# --------------------------
RPF_REFERENCE = {
    0.0: [-9.2, 177.8, -71.8, -26.3, -32.7, -85.8],
    0.2: [-9.9, 186.1, -76.7, -25.7, -39.8, -82.0],
    0.4: [-10.6, 194.4, -81.3, -25.2, -46.6, -78.5],
    0.6: [-11.3, 202.8, -85.8, -24.7, -53.4, -75.1],
    0.8: [-12.0, 211.3, -90.3, -24.2, -60.3, -71.7],
    1.0: [-12.7, 220.1, -94.9, -23.7, -67.2, -68.4]
}

def match_rpf(orbit_data):
    """Find closest RPF value for the given Orbit Data."""
    arr = np.array(orbit_data)
    best_rpf, best_diff = None, np.inf
    for rpf, ref in RPF_REFERENCE.items():
        diff = np.linalg.norm(arr - np.array(ref))
        if diff < best_diff:
            best_rpf, best_diff = rpf, diff
    return best_rpf

def read_feedback_rpf(galaxy_name, num):
    """Read feedback and Orbit Data from elinfo csv."""
    folder = f"./{galaxy_name}{num}"
    file_path = os.path.join(folder, f"elinfo_{galaxy_name}{num}.csv")
    feedback, rpf = None, None
    if not os.path.isfile(file_path):
        print(f"⚠️ Missing elinfo file: {file_path}")
        return feedback, rpf

    with open(file_path, 'r') as f:
        lines = [l.strip() for l in f.readlines()[:2]]

    orbit_data = None
    for l in lines:
        if "Orbit Data" in l:
            nums = l.split(":")[1].split(",")
            orbit_data = [float(x) for x in nums]
        elif "Timescale_fd" in l:
            feedback = float(l.split(":")[1].strip())

    if orbit_data:
        rpf = match_rpf(orbit_data)

    return feedback, rpf

def parse_ini_file(file_path):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(file_path)

    model_data = {}
    gas_count, star_count = 0, 0

    for section in config.sections():
        if not section.startswith("Density_"):
            continue

        component_type = section.split("_")[1].lower()
        if component_type.startswith("gas"):
            gas_count += 1
            comp_name = "Gas" + (str(gas_count) if gas_count > 1 else "")
        elif component_type.startswith("stellar"):
            star_count += 1
            comp_name = "Stars" + (str(star_count) if star_count > 1 else "")
        else:
            continue

        data = {}
        for key, val in config[section].items():
            if key.lower() == "type":
                continue
            k = key.lower()
            try:
                val = float(val)
            except ValueError:
                pass

            if k == "sersicindex":
                data["sersic"] = val
            elif k == "scaleradius":
                data["radius"] = val
            elif k == "scaleheight":
                data["height"] = val
            elif k == "axisratioz":
                data["ratioZ"] = val
            elif k == "mass":
                # Convert to 1e10 scale and scientific notation
                val = val * 1e10
                data["mass"] = f"{val:.1e}"

        model_data[comp_name] = data

    return model_data

def collect_models(galaxy_name, start, end):
    base_dir = ini_model_file
    summary = {}

    for num in range(start, end + 1):
        file_name = f"IC_{galaxy_name}{num}.ini"
        file_path = os.path.join(base_dir, file_name)
        if not os.path.isfile(file_path):
            print(f"⚠️ Missing INI file: {file_path}")
            continue

        feedback, rpf = read_feedback_rpf(galaxy_name, num)
        model_info = parse_ini_file(file_path)
        summary[num] = {"feedback": feedback, "rpf": rpf, "data": model_info}

    return summary

def build_structured_df(galaxy_name, summary_dict):
    """Reorganize DataFrame with grouped component blocks."""
    components = set()
    for info in summary_dict.values():
        components.update(info["data"].keys())
    components = sorted(components, key=lambda x: (x.startswith("Gas"), x))

    rows, index_labels = [], []

    index_labels.append("Feedback")
    rows.append([summary_dict[n]["feedback"] for n in summary_dict])
    index_labels.append("RPF")
    rows.append([summary_dict[n]["rpf"] for n in summary_dict])

    params = ["mass", "sersic", "radius", "height", "ratioZ"]

    for comp in components:
        index_labels.append(comp)
        rows.append(["" for _ in summary_dict])
        for p in params:
            vals = [summary_dict[n]["data"].get(comp, {}).get(p, "") for n in summary_dict]
            # Skip rows that are entirely blank
            if all(v == "" for v in vals):
                continue
            index_labels.append(f"   {p}")
            rows.append(vals)

    df = pd.DataFrame(rows, index=index_labels,
                      columns=[f"{n}" for n in summary_dict])
    return df

def save_as_ods(df, output_file):
    doc = OpenDocumentSpreadsheet()
    table = Table(name="Galaxy Models")

    for _ in range(len(df.columns) + 1):
        table.addElement(TableColumn())

    header_style = Style(name="HeaderCell", family="table-cell")
    header_style.addElement(TableCellProperties(backgroundcolor="#DCE6F1"))
    header_style.addElement(TextProperties(fontweight="bold", fontsize="12pt"))
    doc.styles.addElement(header_style)

    title_style = Style(name="TitleCell", family="table-cell")
    title_style.addElement(TextProperties(fontweight="bold", fontsize="11pt", color="#000000"))
    doc.styles.addElement(title_style)

    # Header row
    header_row = TableRow()
    first = TableCell(stylename=header_style, valuetype="string")
    first.addElement(P(text=f"{galaxy_name}"))
    header_row.addElement(first)
    for col in df.columns:
        cell = TableCell(stylename=header_style, valuetype="string")
        cell.addElement(P(text=col))
        header_row.addElement(cell)
    table.addElement(header_row)

    # Data rows
    for idx, row in df.iterrows():
        tr = TableRow()
        style = title_style if not idx.startswith(" ") else None
        cell = TableCell(valuetype="string", stylename=style)
        cell.addElement(P(text=str(idx)))
        tr.addElement(cell)
        for val in row:
            val_str = "" if pd.isna(val) else str(val)
            cell = TableCell(valuetype="string")
            cell.addElement(P(text=val_str))
            tr.addElement(cell)
        table.addElement(tr)

    doc.spreadsheet.addElement(table)
    doc.save(output_file)
    print(f"✅ Saved as {output_file}")

def parse_model_argument(arg):
    model_numbers = set()
    parts = arg.split(",")
    for part in parts:
        part = part.strip()
        if "-" in part:
            try:
                start, end = map(int, part.split("-"))
                model_numbers.update(range(start, end + 1))
            except ValueError:
                print(f"❌ Invalid range format: '{part}'. Use something like 199-202")
                sys.exit(1)
        else:
            try:
                model_numbers.add(int(part))
            except ValueError:
                print(f"❌ Invalid number format: '{part}'. Use integers like 187,189")
                sys.exit(1)
    return sorted(model_numbers)

def flatten_model_info(model_num, info):
    """Flatten model info into one-line dict for CSV output."""
    row = {
        "model": model_num,
        "feedback": info["feedback"],
        "rpf": info["rpf"],
    }

    params = ["mass", "sersic", "radius", "height", "ratioZ"]

    for comp, comp_dict in info["data"].items():
        for p in params:
            key = f"{comp}_{p}"
            row[key] = comp_dict.get(p, "")

    return row


def save_as_csv(summary_dict, galaxy_name, model_numbers):
    """Save summary_dict as one-line-per-model CSV."""
    flattened = []

    for num in model_numbers:
        if num in summary_dict:
            flattened.append(flatten_model_info(num, summary_dict[num]))

    df = pd.DataFrame(flattened)
    output_file = f"{galaxy_name}_{model_numbers[0]}-{model_numbers[-1]}_summary.csv"
    df.to_csv(output_file, index=False)
    print(f"✅ Saved CSV: {output_file}")


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python3 GetModelInfo.py 187,189,199-202")
        sys.exit(1)

    arg = sys.argv[1]
    model_numbers = parse_model_argument(arg)

    # auto-detect galaxy name
    galaxy_name = os.path.basename(os.getcwd())

    summary_dict = {}
    for num in model_numbers:
        summary_dict.update(collect_models(galaxy_name, num, num))

    save_as_csv(summary_dict, galaxy_name, model_numbers)

