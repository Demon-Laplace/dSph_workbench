import argparse
import csv
import importlib.util
import math
import re
import sys
from functools import lru_cache
from multiprocessing import get_context
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_CASE_ROOT = WORKSPACE_ROOT / "fixtures_realistic" / "mock_server"
DEFAULT_ORBIT_ALIAS_FILE = Path(__file__).with_name("orbit_aliases.py")

MODEL_SPEC_TOKEN_PATTERN = re.compile(r"^\d+(?:-\d+)?$")
IC_INI_PATTERN_TEMPLATE = r"^IC_snapshot(?P<snapshot>\d+)_{prefix}(?P<ic_model>\d+)\.ini$"
MW_MODEL_PATTERN = re.compile(r"^(MW\d+)")
ELINFO_REQUIRED_COLUMNS = (
    "distance",
    "sigma",
    "star_mass",
    "rhalf",
    "hotgas_mass",
    "coldgas_mass",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Inspect Fornax production cases and report which initial condition "
            "each model uses."
        )
    )
    parser.add_argument(
        "model_specs",
        nargs="+",
        help=(
            'Model selectors such as "1035", "1035-1087", or "1035,1087,2001-2027". '
            "Ranges are inclusive."
        ),
    )
    parser.add_argument(
        "--case-root",
        type=Path,
        help=(
            "Root directory containing Fornax model folders. "
            "Default: auto-detect current working directory first, then local fixture path."
        ),
    )
    parser.add_argument(
        "--prefix",
        default="Fornax",
        help='Model directory and IC prefix. Default: "Fornax".',
    )
    parser.add_argument(
        "--format",
        choices=("pretty", "csv", "tsv"),
        default="pretty",
        help="Output format. Default: pretty.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Optional output file. If provided, parent directories are created. "
            "Write outside fixtures, e.g. under sandbox_runs/."
        ),
    )
    parser.add_argument(
        "--orbit-alias-file",
        type=Path,
        default=DEFAULT_ORBIT_ALIAS_FILE,
        help=(
            "Python file defining orbit aliases as top-level names mapped to six-value "
            f"strings. Default: {DEFAULT_ORBIT_ALIAS_FILE}."
        ),
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of worker processes to use. Default: 1.",
    )
    return parser.parse_args()


def expand_model_specs(specs):
    model_ids = set()
    for spec in specs:
        for token in spec.split(","):
            token = token.strip()
            if not token:
                continue
            if not MODEL_SPEC_TOKEN_PATTERN.fullmatch(token):
                raise ValueError(
                    f'Invalid model selector "{token}". Use forms like "1035" or "2001-2027".'
                )
            if "-" in token:
                start_text, end_text = token.split("-", maxsplit=1)
                start = int(start_text)
                end = int(end_text)
                if start > end:
                    raise ValueError(
                        f'Invalid model selector "{token}": range start must be <= end.'
                    )
                model_ids.update(range(start, end + 1))
            else:
                model_ids.add(int(token))
    if not model_ids:
        raise ValueError("No model ids were provided.")
    return sorted(model_ids)


def parse_ic_token(token, prefix):
    ic_pattern = re.compile(IC_INI_PATTERN_TEMPLATE.format(prefix=re.escape(prefix)))
    match = ic_pattern.fullmatch(Path(token).name)
    if match is None:
        return None
    snapshot = match.group("snapshot")
    raw_ic_model = f"{prefix}{match.group('ic_model')}"
    return {
        "snapshot": snapshot,
        "raw_ic_model": raw_ic_model,
        "ic_model": f"IC_{raw_ic_model}",
    }


def parse_orbit_text(value):
    parts = str(value).split()
    if len(parts) != 6:
        raise ValueError(f'Expected six orbit values, got {len(parts)} in "{value}".')
    return tuple(float(part) for part in parts)


def parse_orbit_alias_value(value):
    if isinstance(value, str):
        return parse_orbit_text(value)
    if isinstance(value, (tuple, list)) and len(value) == 6:
        return tuple(float(part) for part in value)
    raise ValueError(
        "Orbit alias values must be either a six-value string or a six-value tuple/list."
    )


@lru_cache(maxsize=None)
def load_alias_module(config_path):
    config_path = config_path.resolve()
    if not config_path.exists():
        return None

    spec = importlib.util.spec_from_file_location("build_production_ic_orbit_aliases", config_path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Could not load orbit alias file: {config_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_orbit_aliases(config_path):
    module = load_alias_module(config_path)
    if module is None:
        return []

    aliases = []
    for name, value in vars(module).items():
        if name.startswith("_"):
            continue
        if isinstance(value, (str, tuple, list)):
            try:
                orbit_values = parse_orbit_alias_value(value)
            except ValueError:
                continue
            aliases.append((name, orbit_values))
    return aliases


def load_dwarf_targets(config_path):
    module = load_alias_module(config_path)
    if module is None:
        return {}

    targets = {}
    for name, value in vars(module).items():
        if not name.startswith("d_"):
            continue
        dwarf_name = name[2:]
        try:
            targets[dwarf_name] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f'Invalid target distance for "{name}": {value!r}') from exc
    return targets


def load_scalar_config(config_path, variable_name):
    module = load_alias_module(config_path)
    if module is None:
        return None
    value = getattr(module, variable_name, None)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'Invalid scalar config for "{variable_name}": {value!r}') from exc


def load_mw_model_aliases(config_path):
    module = load_alias_module(config_path)
    if module is None:
        return {}

    aliases = {}
    for name, value in vars(module).items():
        if name.startswith("_") or not isinstance(value, str):
            continue
        if MW_MODEL_PATTERN.fullmatch(name) is None:
            continue
        try:
            parse_orbit_text(value)
        except ValueError:
            aliases[name] = value.strip()
    return aliases


def resolve_orbit_alias(orbit_values, orbit_aliases):
    if orbit_values is None:
        return ""
    for name, alias_values in orbit_aliases:
        if all(
            math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-6)
            for observed, expected in zip(orbit_values, alias_values)
        ):
            return name
    return ""


def parse_mergegalaxy_metadata(case_dir, prefix, orbit_aliases, mw_model_aliases):
    merge_path = case_dir / "mergegalaxy.par"
    if not merge_path.exists() or not merge_path.is_file():
        return {
            "mw_model": "",
            "merge_ic": None,
            "orbit_values": None,
            "orbit_label": "",
        }

    with merge_path.open("r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    metadata = {
        "mw_model": "",
        "merge_ic": None,
        "orbit_values": None,
        "orbit_label": "",
    }

    if lines:
        tokens = lines[0].split()
        if len(tokens) >= 1:
            mw_source = Path(tokens[0]).stem
            mw_match = MW_MODEL_PATTERN.match(mw_source)
            mw_model_name = mw_match.group(1) if mw_match else mw_source
            metadata["mw_model"] = mw_model_aliases.get(mw_model_name, mw_model_name)
        if len(tokens) >= 2:
            metadata["merge_ic"] = parse_ic_token(tokens[1], prefix)

    if len(lines) >= 2:
        orbit_parts = lines[1].split()
        if len(orbit_parts) >= 6:
            orbit_values = tuple(float(part) for part in orbit_parts[:6])
            metadata["orbit_values"] = orbit_values
            metadata["orbit_label"] = resolve_orbit_alias(orbit_values, orbit_aliases)

    return metadata


def parse_float_field(row, column_name):
    value = row.get(column_name, "")
    if value is None:
        return math.nan
    text = str(value).strip()
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def parse_feedback_value(param_path):
    if not param_path.exists() or not param_path.is_file():
        return None

    with param_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("%"):
                continue
            if not line.startswith("timescale_fd"):
                continue
            value_text = line.split(";", maxsplit=1)[0].split()[1]
            try:
                return float(value_text)
            except (IndexError, ValueError):
                return None
    return None


def parse_softening_values(param_path):
    if not param_path.exists() or not param_path.is_file():
        return {}

    values = {}
    target_keys = {"SofteningGas", "SofteningStars"}
    with param_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split(";", maxsplit=1)[0].split()
            if len(parts) < 2 or parts[0] not in target_keys:
                continue
            try:
                values[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return values


def format_feedback_strength(fd_value, median_fd):
    if fd_value is None or median_fd in (None, 0):
        return ""

    ratio = fd_value / median_fd
    quantized = round(ratio * 2.0) / 2.0
    if math.isclose(quantized, round(quantized), rel_tol=0.0, abs_tol=1e-9):
        number_text = str(int(round(quantized)))
    else:
        number_text = f"{quantized:g}"
    return f"{number_text} m"


def summarize_softening(param_path):
    values = parse_softening_values(param_path)
    gas_value = values.get("SofteningGas")
    stars_value = values.get("SofteningStars")

    if gas_value is None and stars_value is None:
        return "", "softening_missing"
    if gas_value is None or stars_value is None:
        return "", "softening_incomplete"
    if math.isclose(gas_value, stars_value, rel_tol=0.0, abs_tol=1e-12):
        return f"{gas_value:g}", ""
    return f"Gas={gas_value:g},Stars={stars_value:g}", "softening_mismatch"


def find_elinfo_file(case_dir, prefix, model_id):
    expected = case_dir / f"elinfo_{prefix}{model_id}.csv"
    if expected.exists() and expected.is_file():
        return expected

    candidates = sorted(case_dir.glob("elinfo_*.csv"))
    if len(candidates) == 1:
        return candidates[0]
    return None


def summarize_elinfo(case_dir, prefix, model_id, target_distance):
    empty = {
        "elinfo_path": "",
        "target_distance": "",
        "matched_distance": "",
        "distance_delta": "",
        "matched_sigma": "",
        "matched_star_mass": "",
        "matched_rhalf": "",
        "matched_gas_fraction": "",
        "elinfo_status": "",
        "elinfo_message": "",
    }

    if target_distance is None:
        empty["elinfo_status"] = "target_distance_missing"
        empty["elinfo_message"] = f"no d_{prefix} configured in orbit alias file"
        return empty

    elinfo_path = find_elinfo_file(case_dir, prefix, model_id)
    if elinfo_path is None:
        empty["target_distance"] = target_distance
        empty["elinfo_status"] = "elinfo_missing"
        empty["elinfo_message"] = "no unique elinfo csv file found"
        return empty

    best_summary = None
    row_count = 0
    with elinfo_path.open("r", encoding="utf-8", newline="") as handle:
        filtered_lines = (line for line in handle if not line.lstrip().startswith("#"))
        reader = csv.DictReader(filtered_lines)
        if reader.fieldnames is None:
            empty["elinfo_path"] = str(elinfo_path)
            empty["target_distance"] = target_distance
            empty["elinfo_status"] = "elinfo_invalid"
            empty["elinfo_message"] = "elinfo file has no header row"
            return empty

        missing_columns = [name for name in ELINFO_REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing_columns:
            empty["elinfo_path"] = str(elinfo_path)
            empty["target_distance"] = target_distance
            empty["elinfo_status"] = "elinfo_invalid"
            empty["elinfo_message"] = (
                "elinfo file is missing required columns: " + ",".join(missing_columns)
            )
            return empty

        for csv_row in reader:
            row_count += 1
            distance = parse_float_field(csv_row, "distance")
            if math.isnan(distance):
                continue

            delta = abs(distance - target_distance)
            if best_summary is not None and delta >= best_summary["distance_delta"]:
                continue

            star_mass = parse_float_field(csv_row, "star_mass")
            hotgas_mass = parse_float_field(csv_row, "hotgas_mass")
            coldgas_mass = parse_float_field(csv_row, "coldgas_mass")
            gas_mass_total = hotgas_mass + coldgas_mass + star_mass
            gas_fraction = (
                coldgas_mass / gas_mass_total
                if not math.isnan(gas_mass_total) and gas_mass_total > 0.0
                else math.nan
            )

            best_summary = {
                "elinfo_path": str(elinfo_path),
                "target_distance": target_distance,
                "matched_distance": distance,
                "distance_delta": delta,
                "matched_sigma": parse_float_field(csv_row, "sigma"),
                "matched_star_mass": star_mass,
                "matched_rhalf": parse_float_field(csv_row, "rhalf"),
                "matched_gas_fraction": gas_fraction,
                "elinfo_status": "ok",
                "elinfo_message": "",
            }

    if best_summary is None:
        empty["elinfo_path"] = str(elinfo_path)
        empty["target_distance"] = target_distance
        empty["elinfo_status"] = "elinfo_no_valid_rows"
        empty["elinfo_message"] = f"no valid distance rows found in {row_count} elinfo entries"
        return empty

    return best_summary


def inspect_model(case_root, prefix, model_id, orbit_aliases, dwarf_targets, median_fd, mw_model_aliases):
    case_dir = case_root / f"{prefix}{model_id}"
    row = {
        "model": f"{prefix}{model_id}",
        "model_id": str(model_id),
        "case_dir": str(case_dir),
        "mw_model": "",
        "orbit_label": "",
        "ic_snapshot": "",
        "ic_model": "",
        "target_distance": "",
        "matched_distance": "",
        "distance_delta": "",
        "matched_sigma": "",
        "matched_star_mass": "",
        "matched_rhalf": "",
        "matched_gas_fraction": "",
        "fd_value": "",
        "fd_strength": "",
        "softening": "",
        "softening_status": "",
        "elinfo_path": "",
        "elinfo_status": "",
        "elinfo_message": "",
        "status": "",
        "message": "",
    }

    if not case_dir.exists():
        return None

    if not case_dir.is_dir():
        row["status"] = "case_invalid"
        row["message"] = "model path exists but is not a directory"
        return row

    merge_metadata = parse_mergegalaxy_metadata(case_dir, prefix, orbit_aliases, mw_model_aliases)
    row["mw_model"] = merge_metadata["mw_model"]
    row["orbit_label"] = merge_metadata["orbit_label"]

    fd_value = parse_feedback_value(case_dir / "GZWJL.PARAM")
    row["fd_value"] = "" if fd_value is None else fd_value
    row["fd_strength"] = format_feedback_strength(fd_value, median_fd)
    row["softening"], row["softening_status"] = summarize_softening(case_dir / "GZWJL.PARAM")

    elinfo_summary = summarize_elinfo(
        case_dir=case_dir,
        prefix=prefix,
        model_id=model_id,
        target_distance=dwarf_targets.get(prefix),
    )
    row.update(elinfo_summary)

    matched_ics = []
    for path in sorted(case_dir.iterdir()):
        if not path.is_file():
            continue
        parsed_ic = parse_ic_token(path.name, prefix)
        if parsed_ic is not None:
            matched_ics.append(parsed_ic)

    if not matched_ics:
        row["status"] = "ic_ini_missing"
        row["message"] = "no IC_snapshot ini file found"
        return row

    selected_ic = None
    if len(matched_ics) == 1:
        selected_ic = matched_ics[0]
    else:
        unique_ic_models = {ic_info["raw_ic_model"] for ic_info in matched_ics}
        merge_ic = merge_metadata["merge_ic"]
        if (
            len(unique_ic_models) == 1
            and merge_ic is not None
            and merge_ic["raw_ic_model"] in unique_ic_models
        ):
            selected_ic = merge_ic
        else:
            row["status"] = "ic_ini_ambiguous"
            row["message"] = (
                f"found {len(matched_ics)} IC_snapshot ini files with "
                f"{len(unique_ic_models)} IC models"
            )
            row["ic_model"] = ";".join(sorted(ic_info["ic_model"] for ic_info in matched_ics))
            return row

    row["ic_snapshot"] = selected_ic["snapshot"]
    row["ic_model"] = selected_ic["ic_model"]
    row["status"] = "ok"
    row["message"] = f"{row['model']} uses {row['ic_model']} from snapshot {row['ic_snapshot']}"
    return row


def has_model_dirs(case_root, prefix):
    if not case_root.exists() or not case_root.is_dir():
        return False
    return any(
        path.is_dir() and path.name.startswith(prefix)
        for path in case_root.iterdir()
    )


def resolve_case_root(case_root_arg, prefix):
    if case_root_arg is not None:
        return case_root_arg.resolve()

    cwd_root = Path.cwd().resolve()
    if has_model_dirs(cwd_root, prefix):
        return cwd_root

    if has_model_dirs(FIXTURE_CASE_ROOT, prefix):
        return FIXTURE_CASE_ROOT.resolve()

    return cwd_root


def inspect_model_task(task):
    return inspect_model(*task)


def collect_rows(case_root, prefix, model_ids, orbit_aliases, dwarf_targets, median_fd, mw_model_aliases, processes):
    tasks = [
        (case_root, prefix, model_id, orbit_aliases, dwarf_targets, median_fd, mw_model_aliases)
        for model_id in model_ids
    ]
    if processes <= 1:
        rows = [inspect_model_task(task) for task in tasks]
    else:
        with get_context("fork").Pool(processes=processes) as pool:
            rows = pool.map(inspect_model_task, tasks)
    return [row for row in rows if row is not None]


def render_pretty(rows):
    headers = [
        "model",
        "mw_model",
        "orbit_label",
        "ic_model",
        "ic_snapshot",
        "matched_distance",
        "matched_sigma",
        "matched_star_mass",
        "matched_rhalf",
        "matched_gas_fraction",
        "fd_strength",
        "softening",
        "status",
    ]
    widths = {header: len(header) for header in headers}
    for row in rows:
        for header in headers:
            widths[header] = max(widths[header], len(str(row[header])))

    lines = [
        "  ".join(header.ljust(widths[header]) for header in headers),
        "  ".join("-" * widths[header] for header in headers),
    ]
    for row in rows:
        lines.append(
            "  ".join(str(row[header]).ljust(widths[header]) for header in headers)
        )

    notes = [f"{row['model']}: {row['message']}" for row in rows if row["status"] != "ok"]
    notes.extend(
        f"{row['model']}: {row['elinfo_message']}"
        for row in rows
        if row["elinfo_status"] not in ("", "ok") and row["elinfo_message"]
    )
    notes.extend(
        f"{row['model']}: {row['softening_status']} ({row['softening']})"
        for row in rows
        if row["softening_status"]
    )
    if notes:
        lines.append("")
        lines.append("Notes:")
        lines.extend(notes)
    return "\n".join(lines)


def write_delimited(rows, delimiter, output_stream):
    fieldnames = [
        "model",
        "model_id",
        "case_dir",
        "mw_model",
        "orbit_label",
        "ic_snapshot",
        "ic_model",
        "target_distance",
        "matched_distance",
        "distance_delta",
        "matched_sigma",
        "matched_star_mass",
        "matched_rhalf",
        "matched_gas_fraction",
        "fd_value",
        "fd_strength",
        "softening",
        "softening_status",
        "elinfo_path",
        "elinfo_status",
        "elinfo_message",
        "status",
        "message",
    ]
    writer = csv.DictWriter(output_stream, fieldnames=fieldnames, delimiter=delimiter)
    writer.writeheader()
    writer.writerows(rows)


def ensure_output_parent(output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)


def main():
    args = parse_args()
    case_root = resolve_case_root(args.case_root, args.prefix)
    try:
        orbit_aliases = load_orbit_aliases(args.orbit_alias_file)
        dwarf_targets = load_dwarf_targets(args.orbit_alias_file)
        median_fd = load_scalar_config(args.orbit_alias_file, "median_fd")
        mw_model_aliases = load_mw_model_aliases(args.orbit_alias_file)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    try:
        model_ids = expand_model_specs(args.model_specs)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    rows = collect_rows(
        case_root=case_root,
        prefix=args.prefix,
        model_ids=model_ids,
        orbit_aliases=orbit_aliases,
        dwarf_targets=dwarf_targets,
        median_fd=median_fd,
        mw_model_aliases=mw_model_aliases,
        processes=args.processes,
    )

    if args.output is not None:
        output_path = args.output.resolve()
        ensure_output_parent(output_path)
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            if args.format == "pretty":
                handle.write(render_pretty(rows))
                handle.write("\n")
            elif args.format == "csv":
                write_delimited(rows, delimiter=",", output_stream=handle)
            else:
                write_delimited(rows, delimiter="\t", output_stream=handle)
    else:
        if args.format == "pretty":
            print(render_pretty(rows))
        elif args.format == "csv":
            write_delimited(rows, delimiter=",", output_stream=sys.stdout)
        else:
            write_delimited(rows, delimiter="\t", output_stream=sys.stdout)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
