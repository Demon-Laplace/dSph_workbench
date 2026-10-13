#!/usr/bin/env python3

import argparse
import os
import re
from pathlib import Path
from model_paths import model_id

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from basefunc import Analysis, GalaxySimulation


def parse_target(text):
    match = re.fullmatch(r"(\d+)-(\d+)", text.strip())
    if not match:
        raise ValueError("argument must look like '144-200'")
    return int(match.group(1)), int(match.group(2))


def resolve_model_dir(root, modelnum, prefix=None):
    if prefix:
        candidate = root / f"{prefix}{modelnum}"
        if candidate.is_dir():
            return candidate
        raise FileNotFoundError(f"model folder not found: {candidate}")

    candidates = []
    for path in root.iterdir():
        if not path.is_dir():
            continue
        match = re.fullmatch(r"([A-Za-z]+)(\d+)", path.name)
        if match and int(match.group(2)) == modelnum:
            candidates.append(path)

    if not candidates:
        raise FileNotFoundError(f"cannot find a model folder ending with {modelnum}")

    candidates.sort(key=lambda p: (0 if p.name.startswith("Fornax") else 1, p.name))
    return candidates[0]


def resolve_snapshot_folder(model_dir, numsp):
    local_snapshot = model_dir / f"snapshot_{numsp:03d}.hdf5"
    output_snapshot = model_dir / "output" / f"snapshot_{numsp:03d}.hdf5"

    if local_snapshot.exists():
        return model_dir
    if output_snapshot.exists():
        return model_dir / "output"

    raise FileNotFoundError(
        f"snapshot_{numsp:03d}.hdf5 not found in {model_dir} or {model_dir / 'output'}"
    )


def resolve_ini_file(model_dir):
    dwarf_name = re.match(r"([A-Za-z]+)", model_id(model_dir)).group(1)
    preferred = model_dir / f"{dwarf_name}.ini"
    if preferred.exists():
        return preferred

    ini_files = sorted(model_dir.glob("*.ini"))
    for ini_file in ini_files:
        if not ini_file.name.startswith("IC_"):
            return ini_file
    if ini_files:
        return ini_files[0]

    raise FileNotFoundError(f"no ini file found in {model_dir}")


def classify_for_conversion(simulation, df, dw_particles):
    dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)

    mw_star_mask = (~dw_mask) & ((df.tp == 2) | (df.tp == 3))
    dw_star_mask = dw_mask & ((df.tp == 2) | (df.tp == 3))
    dw_gas_mask = dw_mask & (df.tp == 0)
    mw_gas_mask = (~dw_mask) & (df.tp == 0)
    dw_hot_gas_mask = dw_gas_mask & (df.temp >= 20000)
    dw_cold_gas_mask = dw_gas_mask & (df.temp < 20000)

    new_dw_star_mask, new_mw_star_mask = simulation.classify_new_dwarf_stars(dw_star_mask)
    total_dw_star_mask = dw_star_mask | new_dw_star_mask
    total_mw_star_mask = mw_star_mask | new_mw_star_mask
    total_dw_mask = dw_mask | new_dw_star_mask
    df_mw = df.loc[~total_dw_mask].copy()

    return {
        "df_mw": df_mw,
        "dw_star_mask": total_dw_star_mask,
        "dw_hot_gas_mask": dw_hot_gas_mask,
        "dw_cold_gas_mask": dw_cold_gas_mask,
        "mw_star_mask": total_mw_star_mask,
        "mw_gas_mask": mw_gas_mask,
    }


def build_output_dataframe(df_mw):
    cols = ["id", "x", "y", "z", "vx", "vy", "vz", "m", "ugas", "tp"]
    df_ini = df_mw.loc[:, cols].copy()
    df_ini.loc[df_ini["tp"] == 4, "tp"] = 2
    df_ini["m"] = df_ini["m"] / 1e10
    df_ini["ugas"] = df_ini["ugas"].fillna(0)
    df_ini["tp"] = df_ini["tp"].astype(int)
    df_ini = df_ini.sort_values("tp", kind="stable").reset_index(drop=True)
    return df_ini


def save_converted_files(df_ini, output_dir, model_name, numsp, save_ini=True):
    type_0_count = int((df_ini["tp"] == 0).sum())
    type_1_count = int((df_ini["tp"] == 1).sum())
    type_2_count = int((df_ini["tp"] == 2).sum())
    type_3_count = int((df_ini["tp"] == 3).sum())
    header_line = f"{type_0_count}  {type_1_count}  {type_2_count}  {type_3_count}"

    csv_path = output_dir / f"{model_name}_snapshot_{numsp:03d}_converted.csv"
    ini_path = output_dir / f"IC_{model_name}_snapshot_{numsp:03d}.ini"

    df_ini.to_csv(csv_path, index=False, float_format="%.8g")
    if save_ini:
        with ini_path.open("w") as f:
            f.write(header_line + "\n")
            df_ini.to_csv(f, sep=" ", index=False, header=False, float_format="%.8g")
    else:
        ini_path = None

    return csv_path, ini_path


def sample_for_plot(df, max_points=150000):
    if len(df) <= max_points:
        return df
    return df.sample(n=max_points, random_state=42)


def scatter_stars_and_gas(ax, df, xcol, ycol, marker_size=1):
    component_map = [
        ((df["tp"] == 0), "Gas", "tab:blue"),
        ((df["tp"] >= 2), "Stars", "black"),
    ]

    for mask, label, color in component_map:
        sub = sample_for_plot(df.loc[mask, [xcol, ycol]])
        if sub.empty:
            continue
        ax.scatter(sub[xcol], sub[ycol], s=marker_size, lw=0, alpha=0.5, color=color, label=label)


def save_projection_figure(df, output_dir, model_name, numsp):

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=180, sharey=True)

    scatter_stars_and_gas(axes[0], df, "x", "y", marker_size=1)
    axes[0].set_xlabel("x (kpc)")
    axes[0].set_ylabel("y (kpc)")
    axes[0].set_title("x-y")

    scatter_stars_and_gas(axes[1], df, "x", "z", marker_size=1)
    axes[1].set_xlabel("x (kpc)")
    axes[1].set_title("x-z")
    axes[1].tick_params(labelleft=False)

    for ax in axes:
        ax.minorticks_on()
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)
        ax.set_xlim(-40, 40)
        ax.set_ylim(-40, 40)

    fig.tight_layout()

    output_path = output_dir / f"{model_name}_snapshot_{numsp:03d}_xy_xz.png"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def save_surface_brightness_figure(df, dw_star_mask, output_dir, model_name, numsp):
    star_tp = 2
    star_mask = df["tp"] >= star_tp
    if int(star_mask.sum()) == 0:
        return None

    ellipticity, _ = Analysis.calculate_ellipticity(
        df.loc[star_mask, "x"].to_numpy(),
        df.loc[star_mask, "y"].to_numpy(),
        mass=df.loc[star_mask, "m"].to_numpy(),
        n_neighbors=30,
    )
    r_centers, mag_profile = Analysis.radial_magnitude_profile(
        df.loc[star_mask, "x"].to_numpy(),
        df.loc[star_mask, "y"].to_numpy(),
        mass=df.loc[star_mask, "m"].to_numpy(),
        d_kpc=139,
        ep=ellipticity,
        bins=50,
        m_l_relation=2.1,
    )

    fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
    ax.scatter(r_centers, mag_profile, label="Simulation data", c="black", s=10, zorder=3)

    density_candidates = [
        Path("Sample1_density.csv"),
        output_dir / "Sample1_density.csv",
        output_dir / f"{model_name}_density.csv",
        output_dir / "density_profile.csv",
    ]
    fndens = next((path for path in density_candidates if path.exists()), None)
    if fndens is not None:
        dprof = pd.read_csv(fndens)
        factor = 2.4885305633697762e-11
        I_obs = dprof["density"] * factor
        I_obs_err = dprof["density_error"] * factor
        mag_obs = -2.5 * np.log10(I_obs)
        mag_obs_err = 1.0857 * (I_obs_err / I_obs)
        d_today = 139.6
        ax.errorbar(
            np.radians(dprof["xr"]) * d_today,
            mag_obs,
            yerr=mag_obs_err,
            fmt="o",
            capsize=4,
            color="darkorchid",
            ecolor="gray",
            markersize=4,
            elinewidth=1,
            alpha=0.5,
            zorder=1,
            label="Yang 2022",
        )

    ax.scatter(0.791, 24.77, facecolors="none", edgecolors="orange", s=50, zorder=4, lw=1.5)
    ax.scatter(
        0.05,
        23.59,
        facecolors="none",
        edgecolors="orange",
        s=50,
        label="Munoz 2018",
        zorder=4,
        lw=1.5,
    )

    ax.set_ylim(34, 22.5)
    ax.set_xlim(0, 0.5)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=15)
    ax.grid()
    ax.legend(fontsize=12)
    ax.set_xlabel("Radius (kpc)", fontsize=16)
    ax.set_ylabel(r"Surface Brightness (mag arcsec$^{-2}$)", fontsize=16)
    fig.tight_layout()

    output_path = output_dir / f"{model_name}_snapshot_{numsp:03d}_surface_brightness.png"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Convert one model snapshot and save diagnostic figures."
    )
    parser.add_argument("target", help="format: modelnum-numsp, e.g. 1047-200")
    parser.add_argument(
        "--prefix",
        help="optional model prefix such as Fornax or MW; if omitted, auto-detect",
    )
    parser.add_argument(
        "--no-save-ini",
        action="store_false",
        dest="save_ini",
        help="do not save the old space-separated IC text file",
    )
    parser.set_defaults(save_ini=True)
    parser.add_argument(
        "--dw-particles",
        type=int,
        help="manually set dwarf particle count; if omitted, read it from the model ini file",
    )
    args = parser.parse_args()

    root = Path.cwd()
    modelnum, numsp = parse_target(args.target)
    model_dir = resolve_model_dir(root, modelnum, prefix=args.prefix)
    snapshot_dir = resolve_snapshot_folder(model_dir, numsp)
    model_name = model_id(model_dir)

    if args.dw_particles is not None:
        dw_particles = args.dw_particles
    else:
        ini_file = resolve_ini_file(model_dir)
        dw_particles = GalaxySimulation.get_dw_num(str(ini_file))
    if dw_particles < 0:
        raise ValueError("invalid dwarf particle count")

    simulation = GalaxySimulation(folder_path=str(snapshot_dir), snapshot_num=numsp)
    df, tsnap = simulation.load_snapshot()
    df = simulation.center_on_MW()
    df = simulation.convert_coordinates(df=df)

    masks = classify_for_conversion(simulation, df, dw_particles)
    df_ini = build_output_dataframe(masks["df_mw"])

    csv_path, ini_path = save_converted_files(
        df_ini,
        model_dir,
        model_name,
        numsp,
        save_ini=args.save_ini,
    )

    proj_path = save_projection_figure(df, model_dir, model_name, numsp)

    # --- Surface brightness profile (radial) ---
    try:
        sb_profile_path = save_surface_brightness_figure(
            df,
            masks["dw_star_mask"],
            model_dir,
            model_name,
            numsp,
        )
        if sb_profile_path is not None:
            print(f"Surface brightness profile saved: {sb_profile_path}")
    except Exception as e:
        print(f"Surface brightness profile failed: {e}")

    print(f"Model: {model_name}")
    print(f"Snapshot: {numsp:03d}")
    print(f"Snapshot time: {tsnap:.6f}")
    print(f"Converted CSV: {csv_path}")
    print(f"Converted INI: {ini_path}")
    print(f"Projection figure: {proj_path}")


if __name__ == "__main__":
    main()
