#!/usr/bin/env python3
"""Summarize isolation-model physical parameters from INI files.

This script is designed to run inside a galaxy directory such as ``Fornax``.
It scans ``IniFile/IC_<Galaxy><model>.ini`` files, extracts the gas and stellar
component parameters, computes mass resolution from ``nbody``, and exports a
summary table as CSV.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from odf.opendocument import OpenDocumentSpreadsheet
from odf.style import Style, TableCellProperties, TextProperties
from odf.table import Table, TableCell, TableColumn, TableRow
from odf.text import P


INI_DIR = "IniFile"
FILENAME_PATTERN = re.compile(r"^IC_(?P<galaxy>[A-Za-z]+)(?P<model>\d+)\.ini$")
COMPONENTS = ("gas", "stellar")
COLUMNS = [
    "model",
    "gas_mass",
    "gas_sersic_index",
    "gas_radius",
    "gas_height",
    "gas_ratio_z",
    "gas_resolution",
    "stellar_mass",
    "stellar_sersic_index",
    "stellar_radius",
    "stellar_height",
    "stellar_ratio_z",
    "stellar_resolution",
]
COLUMN_LABELS = {
    "model": "Model",
    "gas_mass": "Gas Mass",
    "gas_sersic_index": "Gas SersicIndex",
    "gas_radius": "Gas Radius",
    "gas_height": "Gas Height",
    "gas_ratio_z": "Gas RatioZ",
    "gas_resolution": "Gas Resolution",
    "stellar_mass": "Stellar Mass",
    "stellar_sersic_index": "Stellar SersicIndex",
    "stellar_radius": "Stellar Radius",
    "stellar_height": "Stellar Height",
    "stellar_ratio_z": "Stellar RatioZ",
    "stellar_resolution": "Stellar Resolution",
}
GROUP_LABELS = {
    "full": "Full Parameter Groups",
    "gas": "Gas Parameter Groups",
    "stellar": "Stellar Parameter Groups",
}


def matching_sections(config: configparser.ConfigParser, prefix: str) -> List[str]:
    prefix_lower = prefix.lower()
    matches = []
    for section in config.sections():
        section_lower = section.lower()
        if section_lower == prefix_lower or section_lower.startswith(prefix_lower + "_"):
            matches.append(section)
    return matches


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract gas and stellar physical parameters from isolation-model "
            "INI files and save them as a CSV summary table."
        )
    )
    parser.add_argument(
        "models",
        nargs="?",
        help=(
            "Model numbers to include, for example '1,5,8-12'. "
            "If omitted, all matching INI files in IniFile/ are used."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output CSV path. Default: <Galaxy>_isolation_summary.csv for all "
            "models, or <Galaxy>_<selection>_summary.csv for a subset."
        ),
    )
    parser.add_argument(
        "--no-preview",
        action="store_true",
        help="Do not print a terminal preview table.",
    )
    return parser


def parse_model_selection(selection: str) -> List[int]:
    model_numbers = set()
    for chunk in selection.split(","):
        part = chunk.strip()
        if not part:
            continue
        if "-" in part:
            try:
                start, end = (int(x.strip()) for x in part.split("-", 1))
            except ValueError as exc:
                raise ValueError(f"Invalid range: {part}") from exc
            if start > end:
                start, end = end, start
            model_numbers.update(range(start, end + 1))
        else:
            try:
                model_numbers.add(int(part))
            except ValueError as exc:
                raise ValueError(f"Invalid model number: {part}") from exc
    return sorted(model_numbers)


def get_galaxy_name(base_dir: Path) -> str:
    return base_dir.resolve().name


def find_ini_files(ini_dir: Path, galaxy_name: str) -> Dict[int, Path]:
    model_files: Dict[int, Path] = {}
    for path in ini_dir.glob(f"IC_{galaxy_name}*.ini"):
        match = FILENAME_PATTERN.match(path.name)
        if not match:
            continue
        if match.group("galaxy") != galaxy_name:
            continue
        model_files[int(match.group("model"))] = path
    return dict(sorted(model_files.items()))


def load_ini(path: Path) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.optionxform = str
    with path.open("r", encoding="utf-8") as handle:
        config.read_file(handle)
    return config


def to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def sci_mass_string(mass_1e10: Optional[float]) -> str:
    if mass_1e10 is None:
        return ""
    return f"{mass_1e10 * 1e10:.3e}"


def compact_number(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:g}"


def resolution_string(mass_1e10: Optional[float], nbody: Optional[float]) -> str:
    if mass_1e10 is None or nbody in (None, 0):
        return ""
    return f"{(mass_1e10 * 1e10) / nbody:.3e}"


def get_section(config: configparser.ConfigParser, section_name: str) -> Dict[str, str]:
    if not config.has_section(section_name):
        return {}
    return {key.lower(): value.strip() for key, value in config.items(section_name)}


def has_multiple_components(config: configparser.ConfigParser) -> bool:
    for component in COMPONENTS:
        density_sections = matching_sections(config, f"Density_{component}")
        component_sections = matching_sections(config, f"Component_{component}")
        if len(density_sections) > 1 or len(component_sections) > 1:
            return True
    return False


def extract_component_row(model: int, config: configparser.ConfigParser) -> Dict[str, str]:
    row: Dict[str, str] = {"model": str(model)}

    for component in COMPONENTS:
        density = get_section(config, f"Density_{component}")
        comp = get_section(config, f"Component_{component}")

        mass = to_float(density.get("mass"))
        sersic = to_float(density.get("sersicindex"))
        if sersic is None and density:
            sersic = 1.0
        radius = to_float(density.get("scaleradius"))
        height = to_float(density.get("scaleheight"))
        ratio_z = to_float(density.get("axisratioz"))
        if ratio_z is None:
            ratio_z = to_float(comp.get("axisratioz"))
        nbody = to_float(comp.get("nbody"))

        prefix = "gas" if component == "gas" else "stellar"
        row[f"{prefix}_mass"] = sci_mass_string(mass)
        row[f"{prefix}_sersic_index"] = compact_number(sersic)
        row[f"{prefix}_radius"] = compact_number(radius)
        row[f"{prefix}_height"] = compact_number(height)
        row[f"{prefix}_ratio_z"] = compact_number(ratio_z)
        row[f"{prefix}_resolution"] = resolution_string(mass, nbody)

    return row


def build_rows(model_files: Dict[int, Path], selected_models: Iterable[int]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    missing: List[int] = []
    skipped_multi_component: List[int] = []

    for model in selected_models:
        path = model_files.get(model)
        if path is None:
            missing.append(model)
            continue
        config = load_ini(path)
        if has_multiple_components(config):
            skipped_multi_component.append(model)
            continue
        rows.append(extract_component_row(model, config))

    if missing:
        print(
            "Warning: missing INI files for model(s): "
            + ", ".join(str(num) for num in missing),
            file=sys.stderr,
        )
    if skipped_multi_component:
        print(
            "Skipped multi-component model(s): "
            + ", ".join(str(num) for num in skipped_multi_component)
        )
    return rows


def active_columns(rows: List[Dict[str, str]]) -> List[str]:
    columns = ["model"]
    for column in COLUMNS:
        if column == "model":
            continue
        if any(row.get(column, "") != "" for row in rows):
            columns.append(column)
    return columns


def sanitize_selection_label(selection: str) -> str:
    return re.sub(r"[^0-9A-Za-z,_-]+", "_", selection)


def default_output_path(base_dir: Path, galaxy_name: str, selection: Optional[str]) -> Path:
    if selection:
        return base_dir / f"{galaxy_name}_{sanitize_selection_label(selection)}_summary.csv"
    return base_dir / f"{galaxy_name}_isolation_summary.csv"


def write_csv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_ods(path: Path, rows: List[Dict[str, str]], columns: List[str], sheet_name: str) -> None:
    doc = OpenDocumentSpreadsheet()
    table = Table(name=sheet_name[:31] or "Summary")

    for _ in columns:
        table.addElement(TableColumn())

    header_style = Style(name="HeaderCell", family="table-cell")
    header_style.addElement(TableCellProperties(backgroundcolor="#DCE6F1"))
    header_style.addElement(TextProperties(fontweight="bold"))
    doc.styles.addElement(header_style)

    header_row = TableRow()
    for column in columns:
        cell = TableCell(valuetype="string", stylename=header_style)
        cell.addElement(P(text=COLUMN_LABELS[column]))
        header_row.addElement(cell)
    table.addElement(header_row)

    for row in rows:
        table_row = TableRow()
        for column in columns:
            cell = TableCell(valuetype="string")
            cell.addElement(P(text=row.get(column, "")))
            table_row.addElement(cell)
        table.addElement(table_row)

    doc.spreadsheet.addElement(table)
    doc.save(str(path))


def group_rows(rows: List[Dict[str, str]], group_columns: List[str]) -> List[Dict[str, str]]:
    grouped: Dict[tuple[str, ...], List[Dict[str, str]]] = {}
    for row in rows:
        key = tuple(row.get(column, "") for column in group_columns)
        grouped.setdefault(key, []).append(row)

    sorted_groups = sorted(
        grouped.items(),
        key=lambda item: (-len(item[1]), [int(entry["model"]) for entry in item[1]]),
    )

    result: List[Dict[str, str]] = []
    for index, (key, members) in enumerate(sorted_groups, start=1):
        record: Dict[str, str] = {
            "group_id": str(index),
            "model_count": str(len(members)),
            "models": ",".join(member["model"] for member in members),
        }
        for column, value in zip(group_columns, key):
            record[column] = value
        result.append(record)
    return result


def parameter_value_groups(rows: List[Dict[str, str]], columns: List[str]) -> List[Dict[str, str]]:
    result: List[Dict[str, str]] = []
    for column in columns:
        if column == "model":
            continue
        buckets: Dict[str, List[str]] = {}
        for row in rows:
            value = row.get(column, "")
            if value == "":
                continue
            buckets.setdefault(value, []).append(row["model"])
        for value, models in sorted(
            buckets.items(),
            key=lambda item: (-len(item[1]), item[0]),
        ):
            result.append(
                {
                    "parameter": COLUMN_LABELS[column],
                    "value": value,
                    "model_count": str(len(models)),
                    "models": ",".join(models),
                }
            )
    return result


def sort_models(models: List[str]) -> List[str]:
    return sorted(models, key=lambda value: int(value))


def sort_values(values: List[str]) -> List[str]:
    def key(value: str) -> tuple[int, float | str]:
        try:
            return (0, float(value))
        except ValueError:
            return (1, value)

    return sorted(values, key=key)


def single_parameter_sweeps(rows: List[Dict[str, str]], columns: List[str]) -> List[Dict[str, str]]:
    physical_columns = [column for column in columns if column != "model"]
    sweeps: List[Dict[str, str]] = []
    seen_signatures = set()

    for varying_column in physical_columns:
        fixed_columns = [column for column in physical_columns if column != varying_column]
        grouped: Dict[tuple[str, ...], List[Dict[str, str]]] = {}

        for row in rows:
            varying_value = row.get(varying_column, "")
            if varying_value == "":
                continue
            key = tuple(row.get(column, "") for column in fixed_columns)
            grouped.setdefault(key, []).append(row)

        for fixed_key, members in grouped.items():
            distinct_values = sort_values(
                list({member.get(varying_column, "") for member in members if member.get(varying_column, "") != ""})
            )
            if len(distinct_values) < 2:
                continue

            models = sort_models([member["model"] for member in members])
            signature = (varying_column, tuple(models))
            if signature in seen_signatures:
                continue
            seen_signatures.add(signature)

            fixed_parts = []
            for column, value in zip(fixed_columns, fixed_key):
                if value != "":
                    fixed_parts.append(f"{COLUMN_LABELS[column]}={value}")

            sweeps.append(
                {
                    "varying_parameter": COLUMN_LABELS[varying_column],
                    "model_count": str(len(models)),
                    "models": ",".join(models),
                    "sweep_values": ",".join(distinct_values),
                    "fixed_parameters": "; ".join(fixed_parts),
                }
            )

    sweeps.sort(
        key=lambda row: (
            row["varying_parameter"],
            -int(row["model_count"]),
            [int(model) for model in row["models"].split(",")],
        )
    )
    for index, row in enumerate(sweeps, start=1):
        row["family_id"] = str(index)
    return sweeps


def add_sheet(
    doc: OpenDocumentSpreadsheet,
    table_name: str,
    rows: List[Dict[str, str]],
    columns: List[str],
    header_style: Style,
) -> None:
    table = Table(name=table_name[:31] or "Sheet")
    for _ in columns:
        table.addElement(TableColumn())

    header_row = TableRow()
    for column in columns:
        cell = TableCell(valuetype="string", stylename=header_style)
        cell.addElement(P(text=column))
        header_row.addElement(cell)
    table.addElement(header_row)

    for row in rows:
        table_row = TableRow()
        for column in columns:
            cell = TableCell(valuetype="string")
            cell.addElement(P(text=row.get(column, "")))
            table_row.addElement(cell)
        table.addElement(table_row)

    doc.spreadsheet.addElement(table)


def write_comparison_text(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    full_columns = [column for column in columns if column != "model"]
    gas_columns = [column for column in columns if column.startswith("gas_")]
    stellar_columns = [column for column in columns if column.startswith("stellar_")]

    sections: List[str] = []
    intro_lines = [
        "Model Classification Report",
        f"Total models included in this report: {len(rows)}",
        "This file groups models by shared physical structure and highlights likely single-parameter scan families.",
        "Please compare these groups with the summary table for the detailed parameter values.",
    ]
    sections.append("\n".join(intro_lines))

    sweep_rows = single_parameter_sweeps(rows, columns)
    if sweep_rows:
        lines = [
            "Single-Parameter Scan Families",
            "The following families are likely designed to vary one physical parameter while keeping the others fixed.",
        ]
        for sweep in sweep_rows:
            lines.append(
                f"Family {sweep['family_id']}: models {sweep['models']} "
                f"({sweep['model_count']} models) appear to scan {sweep['varying_parameter']}."
            )
        sections.append("\n".join(lines))
    else:
        lines = [
            "Single-Parameter Scan Families",
            "No clear single-parameter scan family was identified in the selected models.",
        ]
        sections.append("\n".join(lines))

    for title, group_columns, description in [
        (
            "Full Physical Configuration Groups",
            full_columns,
            "Models in the same group share the same full set of extracted physical parameters.",
        ),
        (
            "Gas Configuration Groups",
            gas_columns,
            "Models in the same group share the same gas-component properties.",
        ),
        (
            "Stellar Configuration Groups",
            stellar_columns,
            "Models in the same group share the same stellar-component properties.",
        ),
    ]:
        if not group_columns:
            continue
        lines = [title, description]
        groups = group_rows(rows, group_columns)
        groups = [group for group in groups if int(group["model_count"]) > 1]
        if not groups:
            lines.append("No multi-model group was identified in this category.")
            sections.append("\n".join(lines))
            continue
        for group in groups:
            lines.append(
                f"Group {group['group_id']}: models {group['models']} "
                f"({group['model_count']} models)"
            )
        sections.append("\n".join(lines))

    path.write_text("\n\n".join(sections) + "\n", encoding="utf-8")


def preview_table(rows: List[Dict[str, str]], limit: int = 20) -> str:
    preview_columns = [
        "model",
        "gas_mass",
        "gas_radius",
        "gas_height",
        "gas_resolution",
        "stellar_mass",
        "stellar_radius",
        "stellar_ratio_z",
        "stellar_resolution",
    ]
    headers = {
        "model": "Model",
        "gas_mass": "GasMass",
        "gas_radius": "GasR",
        "gas_height": "GasH",
        "gas_resolution": "GasRes",
        "stellar_mass": "StarMass",
        "stellar_radius": "StarR",
        "stellar_ratio_z": "StarZ",
        "stellar_resolution": "StarRes",
    }
    preview_columns = [
        column
        for column in preview_columns
        if any(row.get(column, "") != "" for row in rows)
    ]

    preview_rows = rows[:limit]
    widths: Dict[str, int] = {}
    for column in preview_columns:
        widths[column] = max(
            len(headers[column]),
            *(len(row.get(column, "")) for row in preview_rows),
        )

    def render_line(values: Dict[str, str]) -> str:
        return "  ".join(values[col].ljust(widths[col]) for col in preview_columns)

    header_line = render_line(headers)
    divider = "  ".join("-" * widths[col] for col in preview_columns)
    body = [
        render_line({col: row.get(col, "") for col in preview_columns})
        for row in preview_rows
    ]
    if len(rows) > limit:
        body.append(f"... ({len(rows) - limit} more rows)")
    return "\n".join([header_line, divider, *body])


def main() -> int:
    args = make_parser().parse_args()
    base_dir = Path.cwd()
    galaxy_name = get_galaxy_name(base_dir)
    ini_dir = base_dir / INI_DIR

    if not ini_dir.is_dir():
        print(f"Error: missing directory: {ini_dir}", file=sys.stderr)
        return 1

    model_files = find_ini_files(ini_dir, galaxy_name)
    if not model_files:
        print(
            f"Error: no files matching IC_{galaxy_name}<model>.ini were found in {ini_dir}",
            file=sys.stderr,
        )
        return 1

    if args.models:
        try:
            selected_models = parse_model_selection(args.models)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
    else:
        selected_models = sorted(model_files)

    rows = build_rows(model_files, selected_models)
    if not rows:
        print("Error: no model data was collected.", file=sys.stderr)
        return 1

    columns = active_columns(rows)
    output_path = Path(args.output) if args.output else default_output_path(
        base_dir, galaxy_name, args.models
    )
    ods_path = output_path.with_suffix(".ods")
    compare_txt_path = output_path.with_name(output_path.stem + "_compare.txt")
    write_csv(output_path, rows, columns)
    write_ods(ods_path, rows, columns, f"{galaxy_name}_summary")
    write_comparison_text(compare_txt_path, rows, columns)

    if not args.no_preview:
        print(preview_table(rows))
        print()
    print(f"Saved summary CSV: {output_path}")
    print(f"Saved summary ODS: {ods_path}")
    print(f"Saved comparison TXT: {compare_txt_path}")
    print("Resolution is computed as component_mass / nbody in solar-mass units.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
