import argparse
import os
import re
from pathlib import Path
from typing import Optional

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare gas X-Z and X-Y distributions for two models by reading "
            "existing Fornax.ini files."
        )
    )
    parser.add_argument(
        "model_a",
        help='First model identifier, e.g. "34" or "Fornax34"',
    )
    parser.add_argument(
        "model_b",
        help='Second model identifier, e.g. "139" or "Fornax139"',
    )
    parser.add_argument(
        "--limit",
        type=float,
        default=10.0,
        help="Axis half-width in kpc. Default: 10, matching funcini.py.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output PNG path. Default: compareIC_<model_a>_vs_<model_b>.png",
    )
    return parser.parse_args()


def resolve_model_name(model_ref: str, base_name: str) -> str:
    model_ref = model_ref.strip()
    if re.fullmatch(r"\d+", model_ref):
        return f"{base_name}{model_ref}"
    if re.fullmatch(r"[A-Za-z_]+\d+", model_ref):
        return model_ref
    raise ValueError(f"Invalid model identifier: {model_ref}")


def get_ini_path(model_name: str, root_dir: Path) -> Path:
    ini_path = root_dir / model_name / "Fornax.ini"
    if not ini_path.exists():
        raise FileNotFoundError(f"Fornax.ini not found: {ini_path}")
    return ini_path


def read_gas_positions(ini_path: Path) -> np.ndarray:
    with ini_path.open("r", encoding="utf-8") as fh:
        header = fh.readline().split()

    if len(header) < 1:
        raise ValueError(f"Invalid header in {ini_path}")

    gas_count = int(header[0])
    if gas_count <= 0:
        raise ValueError(f"No gas particles found in {ini_path}")

    data = np.loadtxt(
        ini_path,
        skiprows=1,
        max_rows=gas_count,
        usecols=(1, 2, 3),
    )
    return np.atleast_2d(data)


def get_plot_limit(gas_a: np.ndarray, gas_b: np.ndarray, user_limit: Optional[float]) -> float:
    del gas_a
    del gas_b
    if user_limit is None:
        return 10.0
    return user_limit


def plot_model_row(ax_xz, ax_xy, gas_xyz: np.ndarray, model_name: str, limit: float) -> None:
    x = gas_xyz[:, 0]
    y = gas_xyz[:, 1]
    z = gas_xyz[:, 2]
    label = f"IC_{model_name}"

    ax_xz.scatter(x, z, s=0.4, lw=0)
    ax_xz.set_xlim(-limit, limit)
    ax_xz.set_ylim(-limit, limit)
    ax_xz.set_xlabel("X")
    ax_xz.set_ylabel("Z")
    ax_xz.set_title(f"{label} Gas X-Z")

    ax_xy.scatter(x, y, s=0.4, lw=0)
    ax_xy.set_xlim(-limit, limit)
    ax_xy.set_ylim(-limit, limit)
    ax_xy.set_xlabel("X")
    ax_xy.set_ylabel("Y")
    ax_xy.set_title(f"{label} Gas X-Y")


def main():
    args = parse_args()
    root_dir = Path.cwd()
    base_name = root_dir.name

    model_a = resolve_model_name(args.model_a, base_name)
    model_b = resolve_model_name(args.model_b, base_name)

    ini_a = get_ini_path(model_a, root_dir)
    ini_b = get_ini_path(model_b, root_dir)

    gas_a = read_gas_positions(ini_a)
    gas_b = read_gas_positions(ini_b)

    limit = get_plot_limit(gas_a, gas_b, args.limit)

    fig, axes = plt.subplots(2, 2, figsize=(10, 10), dpi=150)
    plot_model_row(axes[0, 0], axes[0, 1], gas_a, model_a, limit)
    plot_model_row(axes[1, 0], axes[1, 1], gas_b, model_b, limit)

    for ax in axes.flat:
        ax.tick_params(labelsize=8)

    plt.tight_layout()

    output_name = args.output or f"compareIC_{model_a}_vs_{model_b}.png"
    output_path = root_dir / output_name
    plt.savefig(output_path)
    plt.close()

    print(f"Saved comparison figure: {output_path}")


if __name__ == "__main__":
    main()
