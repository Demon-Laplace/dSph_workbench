import argparse
import importlib.util
import math
import os
import re
import sys
import warnings
from pathlib import Path
from typing import Optional

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")
warnings.filterwarnings(
    "ignore",
    message=r"mergecube\.par not found under .*; falling back to legacy cube_path=.*",
    category=RuntimeWarning,
)

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent if (SCRIPT_DIR.parent / "data_processing.py").exists() else SCRIPT_DIR
WORKSPACE_ROOT = PROJECT_ROOT.parent


def load_check_density_module():
    candidate_paths = [
        SCRIPT_DIR / "check_density_profile.py",
        SCRIPT_DIR / "check_density_profile" / "check_density_profile.py",
        PROJECT_ROOT / "check_density_profile" / "check_density_profile.py",
        Path.cwd().resolve() / "check_density_profile.py",
        Path.cwd().resolve() / "check_density_profile" / "check_density_profile.py",
    ]
    seen = set()
    for candidate in candidate_paths:
        if candidate in seen or not candidate.exists():
            continue
        seen.add(candidate)
        spec = importlib.util.spec_from_file_location("_compare_ic_check_density_profile", candidate)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    tried = "\n".join(f"  - {path}" for path in candidate_paths)
    raise ModuleNotFoundError(
        "Could not locate check_density_profile.py. Tried:\n"
        f"{tried}"
    )


CHECK_DENSITY_MODULE = load_check_density_module()
compute_summary = CHECK_DENSITY_MODULE.compute_summary
load_centered_dwarf_dataframe = CHECK_DENSITY_MODULE.load_centered_dwarf_dataframe
load_density_history = CHECK_DENSITY_MODULE.load_density_history
load_wlm_reference = CHECK_DENSITY_MODULE.load_wlm_reference


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot IC density-profile panels for one or more models. "
            "Reproduces the original fig[0,1] and fig[1,1] content in a shared figure."
        )
    )
    parser.add_argument(
        "models",
        nargs="+",
        help="Model specs, e.g. 119-500 131-300. Format: <model>-<snapshot> or plain model id/path.",
    )
    parser.add_argument(
        "--dwarf-name",
        default="Fornax",
        help="Dwarf prefix used when positional models are numeric. Default: Fornax",
    )
    parser.add_argument(
        "--model-root",
        help=(
            "Optional base directory for model lookup. "
            "Useful when models are not under the default fixture/current-workdir search roots."
        ),
    )
    parser.add_argument(
        "--snapshot",
        type=int,
        help="Snapshot number used for the radial surface-density comparison. Defaults to the latest available snapshot per model.",
    )
    parser.add_argument("--radius-max", type=float, default=10.0, help="3D radial cut in kpc for analysis samples.")
    parser.add_argument("--nbin", type=int, default=100, help="Number of radial bins.")
    parser.add_argument("--sigma-radius", type=float, default=1.0, help="Projected radius in kpc for velocity-dispersion metrics.")
    parser.add_argument(
        "--output-dir",
        help="Directory for the comparison figure. Defaults to ../sandbox_runs/compare_ic",
    )
    parser.add_argument(
        "--output-name",
        help="Optional custom filename for the generated figure.",
    )
    return parser.parse_args()


def build_model_label(model_spec: str, dwarf_name: str):
    if Path(model_spec).exists():
        return Path(model_spec).resolve().name
    if model_spec.isdigit():
        return f"{dwarf_name}{int(model_spec)}"
    return model_spec


def parse_model_spec(raw_spec: str):
    model_part = raw_spec
    snapshot_part = None
    if not Path(raw_spec).exists():
        match = re.match(r"^(?P<model>.+)-(?P<snapshot>\d+)$", raw_spec)
        if match:
            model_part = match.group("model")
            snapshot_part = int(match.group("snapshot"))
    return {
        "raw_spec": raw_spec,
        "model_spec": model_part,
        "snapshot": snapshot_part,
    }


def candidate_search_roots(dwarf_name: str, model_root: Optional[str]):
    roots = []
    if model_root:
        root = Path(model_root).resolve()
        roots.extend([root, root / dwarf_name])
    else:
        cwd = Path.cwd().resolve()
        roots.extend(
            [
                cwd,
                cwd / dwarf_name,
                cwd / "IsolationModel" / dwarf_name,
                WORKSPACE_ROOT / "fixtures_realistic" / "local_pc" / "IsolationModel" / dwarf_name,
                WORKSPACE_ROOT / "fixtures_realistic" / "local_pc" / dwarf_name,
                WORKSPACE_ROOT / "fixtures_realistic" / "mock_server",
            ]
        )
    return [root for root in roots if root.exists()]


def resolve_model_dir(model_spec: str, dwarf_name: str, model_root: Optional[str]):
    model_path = Path(model_spec)
    if model_path.exists():
        if not model_path.is_dir():
            raise NotADirectoryError(f"Model path is not a directory: {model_path}")
        return model_path.resolve()

    model_label = build_model_label(model_spec, dwarf_name)
    tried = []
    for root in candidate_search_roots(dwarf_name, model_root):
        candidate = root / model_label
        tried.append(candidate)
        if candidate.exists():
            if not candidate.is_dir():
                raise NotADirectoryError(f"Resolved model path is not a directory: {candidate}")
            return candidate.resolve()

    tried_text = "\n".join(f"  - {path}" for path in tried) if tried else "  - <no search roots>"
    raise FileNotFoundError(
        f'Could not resolve model "{model_spec}" (label={model_label}). Tried:\n{tried_text}'
    )


def detect_available_snapshots(model_dir: Path):
    snapshot_numbers = set(list_snapshot_numbers(model_dir))
    for subdir in sorted(path for path in model_dir.iterdir() if path.is_dir()):
        snapshot_numbers.update(list_snapshot_numbers(subdir))
    return sorted(snapshot_numbers)


def list_snapshot_numbers(folder: Path):
    numbers = []
    for file_path in Path(folder).glob("snapshot_*.hdf5"):
        match = re.match(r"snapshot_(\d+)\.hdf5", file_path.name)
        if match:
            numbers.append(int(match.group(1)))
    return sorted(numbers)


def select_snapshot(model_dir: Path, requested_snapshot: Optional[int]):
    snapshots = detect_available_snapshots(model_dir)
    if not snapshots:
        raise FileNotFoundError(f"No snapshots found under {model_dir}")
    if requested_snapshot is None:
        return snapshots[-1]
    if requested_snapshot not in snapshots:
        raise FileNotFoundError(
            f"Snapshot {requested_snapshot:03d} not found under {model_dir}. "
            f"Available snapshots: {snapshots}"
        )
    return requested_snapshot


def load_density_history_from_candidates(model_dir: Path):
    candidates = [
        model_dir,
        model_dir / "output",
        model_dir / "density_output",
    ]
    for candidate in candidates:
        if candidate.exists():
            history = load_density_history(candidate)
            if history is not None:
                return history
    return None


def parse_softening_values(param_path: Path):
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


def summarize_softening(param_path: Path):
    values = parse_softening_values(param_path)
    gas_value = values.get("SofteningGas")
    stars_value = values.get("SofteningStars")

    if gas_value is None and stars_value is None:
        return "n/a"
    if gas_value is None or stars_value is None:
        return "incomplete"
    if math.isclose(gas_value, stars_value, rel_tol=0.0, abs_tol=1e-12):
        return f"{gas_value:g}"
    return f"Gas={gas_value:g},Stars={stars_value:g}"


def load_softening_summary(model_dir: Path):
    candidate_paths = [
        model_dir / "GZWJL.PARAM",
        model_dir.parent / "GZWJL.PARAM",
    ]
    for param_path in candidate_paths:
        if param_path.exists():
            return summarize_softening(param_path)
    return "n/a"


def build_density_history(model_dir: Path, radius_max: float, nbin: int, sigma_radius: float):
    history = load_density_history_from_candidates(model_dir)
    if history is not None:
        return history

    records = []
    for snapshot in detect_available_snapshots(model_dir):
        _, _, tsnap, dwarf_df = load_centered_dwarf_dataframe(model_dir, snapshot)
        summary = compute_summary(dwarf_df, radius_max=radius_max, nbin=nbin, sigma_radius=sigma_radius)
        records.append((tsnap, summary["sigma_h_inner_gas"], summary["sigma_h_inner_star"]))

    if not records:
        return None

    records_array = np.array(sorted(records, key=lambda item: item[0]), dtype=float)
    return {
        "time": records_array[:, 0],
        "gas_density": records_array[:, 1],
        "stellar_density": records_array[:, 2],
    }


def build_model_result(model_request, args):
    model_spec = model_request["model_spec"]
    requested_snapshot = model_request["snapshot"] if model_request["snapshot"] is not None else args.snapshot

    model_dir = resolve_model_dir(model_spec, args.dwarf_name, args.model_root)
    snapshot = select_snapshot(model_dir, requested_snapshot)
    model_label = build_model_label(model_spec, args.dwarf_name)

    _, _, tsnap, dwarf_df = load_centered_dwarf_dataframe(model_dir, snapshot)
    summary = compute_summary(
        dwarf_df,
        radius_max=args.radius_max,
        nbin=args.nbin,
        sigma_radius=args.sigma_radius,
    )
    density_history = build_density_history(
        model_dir,
        radius_max=args.radius_max,
        nbin=args.nbin,
        sigma_radius=args.sigma_radius,
    )
    softening = load_softening_summary(model_dir)
    return {
        "raw_spec": model_request["raw_spec"],
        "model_spec": model_spec,
        "model_label": model_label,
        "model_dir": model_dir,
        "snapshot": snapshot,
        "tsnap": tsnap,
        "summary": summary,
        "density_history": density_history,
        "softening": softening,
    }


def make_output_dir(output_dir_arg: Optional[str]):
    if output_dir_arg:
        return Path(output_dir_arg).resolve()
    sandbox_root = WORKSPACE_ROOT / "sandbox_runs"
    if sandbox_root.exists():
        return (sandbox_root / "compare_ic").resolve()
    return (Path.cwd().resolve() / "compare_ic_output").resolve()


def make_output_name(results, custom_name: Optional[str]):
    if custom_name:
        return custom_name
    joined = "_vs_".join(f"{result['model_label']}-{result['snapshot']:03d}" for result in results)
    return f"compareIC_{joined}.png"


def apply_log_density_axis(ax, x, y, **plot_kwargs):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    positive = y[y > 0]
    if positive.size > 0:
        floor = positive.min() * 1e-3
        ax.plot(x, np.where(y <= 0, floor, y), **plot_kwargs)
        return True
    ax.plot(x, y, **plot_kwargs)
    return False


def make_comparison_figure(results, wlm_r=None, wlm_f=None):
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), dpi=150)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    radial_ax = axes[0]
    for idx, result in enumerate(results):
        color = colors[idx % len(colors)]
        summary = result["summary"]
        label = f"{result['model_label']}-{result['snapshot']:03d}"
        radial_ax.plot(
            summary["r_cent"],
            summary["sigma_h_cold"],
            lw=2.0,
            color=color,
            ls="-",
            label=label,
        )
    if wlm_r is not None and wlm_f is not None:
        radial_ax.plot(wlm_r, wlm_f, color="black", lw=1.8, ls=":", label="WLM (Neel et al.)")
    radial_ax.set_yscale("log")
    radial_ax.set_xlim(0, 10)
    radial_ax.set_ylim(0.01, 55)
    radial_ax.set_xlabel(r"$R\ [{\rm kpc}]$")
    radial_ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
    radial_ax.set_title("Gas Radial Surface Density Profile")
    radial_ax.grid(True, ls="--", alpha=0.4)
    radial_ax.legend(fontsize=9)
    radial_ax.tick_params(axis="both", which="both", direction="in")

    history_ax = axes[1]
    use_log = False
    for idx, result in enumerate(results):
        color = colors[idx % len(colors)]
        label = result["model_label"]
        history = result["density_history"]
        if history is None:
            continue
        use_log |= apply_log_density_axis(
            history_ax,
            history["time"],
            history["gas_density"],
            lw=2.0,
            color=color,
            ls="-",
            label=label,
        )
    if use_log:
        history_ax.set_yscale("log")
    history_ax.set_xlim(0, 8)
    history_ax.set_ylim(0.01, 55)
    history_ax.set_xlabel("Time [Gyr]")
    history_ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
    history_ax.set_title("Central Gas Surface Density Evolution")
    history_ax.grid(True, ls="--", alpha=0.4)
    handles, labels = history_ax.get_legend_handles_labels()
    if handles:
        history_ax.legend(fontsize=9)
    history_ax.tick_params(axis="both", which="both", direction="in")

    summary_lines = [
        f"{result['model_label']}: profile snap={result['snapshot']:03d}, T={result['tsnap']:.2f} Gyr, "
        f"sigma_y={result['summary']['sigma_y']:.2f} km/s, softening={result['softening']}"
        for result in results
    ]
    fig.suptitle("\n".join(summary_lines), fontsize=10)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.12, top=0.84, wspace=0.22)
    return fig


def main():
    args = parse_args()
    if not args.models:
        raise SystemExit("Please provide at least one model, e.g. python3 compareIC.py 119-500")

    model_requests = [parse_model_spec(model_spec) for model_spec in args.models]
    results = [build_model_result(model_request, args) for model_request in model_requests]
    wlm_r, wlm_f = load_wlm_reference(results[0]["model_dir"])

    output_dir = make_output_dir(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / make_output_name(results, args.output_name)

    fig = make_comparison_figure(results, wlm_r=wlm_r, wlm_f=wlm_f)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved comparison figure to: {output_path}")
    for result in results:
        print(
            f"{result['model_label']}: dir={result['model_dir']}, "
            f"profile_snapshot={result['snapshot']:03d}, T={result['tsnap']:.3f} Gyr"
        )


if __name__ == "__main__":
    main()
