import argparse
import contextlib
import io
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import pandas as pd

from basefunc import Analysis, GalaxySimulation
from variable import d_today, fndens

import matplotlib.pyplot as plt


plt.switch_backend("Agg")

SURFACE_BRIGHTNESS_FACTOR = 2.6885305633697762e-11
DEFAULT_PROFILE_DISTANCE_KPC = 139.6
DEFAULT_PLOT_LIMIT_KPC = 20.0
DEFAULT_ML_RELATION = 2.1


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a snapshot to IC ini format and save diagnostic figures.",
    )
    parser.add_argument(
        "target",
        nargs="?",
        help="Compact target in the form <model-num>-<snapshot-num>, for example 152-800.",
    )
    parser.add_argument("--model-num", type=int, default=152, help="Model number, e.g. 152 for Fornax152.")
    parser.add_argument("--snapshot-num", type=int, default=800, help="Snapshot number to export.")
    parser.add_argument("--ntp", type=int, default=2, help="Minimum stellar particle type used for centering.")
    parser.add_argument("--star-tp", type=int, default=2, help="Minimum stellar particle type used for surface brightness.")
    parser.add_argument("--profile-distance-kpc", type=float, default=DEFAULT_PROFILE_DISTANCE_KPC)
    parser.add_argument("--surface-brightness-bins", type=int, default=50)
    parser.add_argument("--plot-limit-kpc", type=float, default=DEFAULT_PLOT_LIMIT_KPC)
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    if args.target:
        match = re.fullmatch(r"(\d+)-(\d+)", args.target)
        if match is None:
            parser.error("target must be in the form <model-num>-<snapshot-num>, for example 152-800")
        args.model_num = int(match.group(1))
        args.snapshot_num = int(match.group(2))

    return args


def get_paths(model_num):
    script_dir = Path.cwd()
    dwarf_name = script_dir.name
    model_name = f"{dwarf_name}{model_num}"
    model_dir = script_dir / model_name
    inifile = model_dir / f"{dwarf_name}.ini"
    density_profile_path = script_dir / Path(fndens).name
    return script_dir, dwarf_name, model_name, model_dir, inifile, density_profile_path


def load_and_center_snapshot(model_dir, snapshot_num, dw_particles, ntp):
    simulation = GalaxySimulation(folder_path=f"{model_dir}/", snapshot_num=snapshot_num)
    with contextlib.redirect_stdout(io.StringIO()):
        df, tsnap = simulation.load_snapshot(box=False)
    dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)
    df = df.loc[dw_mask].reset_index(drop=True)

    star_mask = df["tp"] >= ntp
    if not np.any(star_mask):
        raise ValueError(f"No particles with tp >= {ntp} were found after dwarf selection.")

    xc, _ = Analysis.find_center(df.loc[star_mask, "x"], df.loc[star_mask, "y"], units="kpc")
    yc, zc = Analysis.find_center(df.loc[star_mask, "y"], df.loc[star_mask, "z"], units="kpc")

    df = df.copy()
    df["x"] -= xc
    df["y"] -= yc
    df["z"] -= zc
    df["r"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2 + df["z"] ** 2)
    return df, tsnap


def summarize_dwarf_components(df, ntp):
    star_mask = df["tp"] >= ntp
    gas_mask = df["tp"] == 0

    sigma_x = float(np.std(df.loc[star_mask, "vx"])) if np.any(star_mask) else np.nan
    sigma_y = float(np.std(df.loc[star_mask, "vy"])) if np.any(star_mask) else np.nan
    sigma_z = float(np.std(df.loc[star_mask, "vz"])) if np.any(star_mask) else np.nan

    return {
        "star_count": int(np.count_nonzero(star_mask)),
        "gas_count": int(np.count_nonzero(gas_mask)),
        "sigma_x": sigma_x,
        "sigma_y": sigma_y,
        "sigma_z": sigma_z,
    }


def save_components_figure(df, ntp, output_path, model_name, snapshot_num, limit_kpc, dpi):
    star_mask = df["tp"] >= ntp

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=dpi)
    projections = [
        ("x", "y", "XY"),
        ("x", "z", "XZ"),
    ]

    for ax, (x_col, y_col, title) in zip(axes, projections):
        ax.scatter(df[x_col], df[y_col], s=1, alpha=0.25, color="gray", label="All dwarf particles")
        ax.scatter(df.loc[star_mask, x_col], df.loc[star_mask, y_col], s=1, alpha=0.5, color="royalblue", label="Stars")
        ax.scatter(0, 0, s=20, c="red", marker="x", zorder=5, label="Center")
        ax.set_xlim(-limit_kpc, limit_kpc)
        ax.set_ylim(-limit_kpc, limit_kpc)
        ax.set_xlabel(f"{x_col} (kpc)")
        ax.set_ylabel(f"{y_col} (kpc)")
        ax.set_title(f"{title} Projection")
        ax.grid(True, ls="--", alpha=0.3)
        ax.minorticks_on()
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)

    handles, labels = axes[0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    axes[0].legend(unique.values(), unique.keys(), loc="upper right", fontsize=9)
    fig.suptitle(f"{model_name} snapshot {snapshot_num}", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_surface_brightness_figure(
    df,
    star_tp,
    density_profile_path,
    output_path,
    model_name,
    snapshot_num,
    profile_distance_kpc,
    bins,
    dpi,
):
    star_mask = df["tp"] >= star_tp
    if not np.any(star_mask):
        raise ValueError(f"No particles with tp >= {star_tp} were found for surface brightness.")

    stellar_mass = df.loc[star_mask, "m"]
    ellipticity, _ = Analysis.calculate_ellipticity(
        df.loc[star_mask, "x"],
        df.loc[star_mask, "y"],
        mass=stellar_mass,
        n_neighbors=30,
    )
    r_centers, mag_profile = Analysis.radial_magnitude_profile(
        df.loc[star_mask, "x"],
        df.loc[star_mask, "y"],
        mass=stellar_mass,
        d_kpc=profile_distance_kpc,
        ep=ellipticity,
        bins=bins,
        m_l_relation=DEFAULT_ML_RELATION,
    )

    dprof = pd.read_csv(density_profile_path)
    i_obs = dprof["density"] * SURFACE_BRIGHTNESS_FACTOR
    i_obs_err = dprof["density_error"] * SURFACE_BRIGHTNESS_FACTOR
    mag_obs = -2.5 * np.log10(i_obs)
    mag_obs_err = 1.0857 * (i_obs_err / i_obs)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    ax.scatter(r_centers, mag_profile, label="Simulation data", c="black", s=10, zorder=3)
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
    ax.scatter(0.05, 23.59, facecolors="none", edgecolors="orange", s=50, label="Munoz 2018", zorder=4, lw=1.5)
    ax.set_ylim(34, 22.5)
    ax.set_xlim(0, 0.5)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=12)
    ax.grid(True, ls="--", alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_xlabel("Radius (kpc)", fontsize=14)
    ax.set_ylabel(r"Surface Brightness (mag arcsec$^{-2}$)", fontsize=14)
    ax.set_title(f"{model_name} snapshot {snapshot_num}", fontsize=13)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def build_export_dataframe(df):
    export_df = df.copy()
    export_df.loc[export_df["tp"] == 4, "tp"] = 2
    export_df["m"] = export_df["m"] / 1e10
    export_df["temp"] = export_df["temp"].fillna(0)
    export_df["tp"] = export_df["tp"].astype(int)
    export_df = export_df.loc[export_df["tp"].isin([0, 2])].copy()
    export_df = export_df.sort_values("tp").reset_index(drop=True)
    return export_df


def write_ini(export_df, output_path):
    type_0_count = int((export_df["tp"] == 0).sum())
    type_1_count = int((export_df["tp"] == 1).sum())
    type_2_count = int((export_df["tp"] == 2).sum())
    type_3_count = int((export_df["tp"] == 3).sum())
    header_line = f"{type_0_count}  {type_1_count}  {type_2_count}  {type_3_count}"

    columns = ["id", "x", "y", "z", "vx", "vy", "vz", "m", "ugas", "tp"]
    export_block = export_df.loc[:, columns].copy()
    export_block["id"] = export_block["id"].astype(int)
    export_block["tp"] = export_block["tp"].astype(int)

    with open(output_path, "w") as fh:
        fh.write(header_line + "\n")
        export_block.to_csv(fh, sep=" ", header=False, index=False)


def main():
    args = parse_args()
    script_dir, dwarf_name, model_name, model_dir, inifile, density_profile_path = get_paths(args.model_num)

    if not model_dir.exists():
        raise FileNotFoundError(f"Model directory not found: {model_dir}")
    if not inifile.exists():
        raise FileNotFoundError(f"Model ini file not found: {inifile}")
    if not density_profile_path.exists():
        raise FileNotFoundError(f"Observed density profile not found: {density_profile_path}")
    if not re.match(r"[A-Za-z]+\d+", model_name):
        raise ValueError(f"Invalid model name: {model_name}")

    dw_particles = GalaxySimulation.get_dw_num(inifile)
    df, tsnap = load_and_center_snapshot(model_dir, args.snapshot_num, dw_particles, args.ntp)
    summary = summarize_dwarf_components(df, args.ntp)

    components_path = model_dir / f"IC_snapshot{args.snapshot_num}_{model_name}_components.png"
    surface_brightness_path = model_dir / f"IC_snapshot{args.snapshot_num}_{model_name}_surface_brightness.png"
    export_path = model_dir / f"IC_snapshot{args.snapshot_num}_{model_name}.ini"

    save_components_figure(
        df=df,
        ntp=args.ntp,
        output_path=components_path,
        model_name=model_name,
        snapshot_num=args.snapshot_num,
        limit_kpc=args.plot_limit_kpc,
        dpi=args.dpi,
    )
    save_surface_brightness_figure(
        df=df,
        star_tp=args.star_tp,
        density_profile_path=density_profile_path,
        output_path=surface_brightness_path,
        model_name=model_name,
        snapshot_num=args.snapshot_num,
        profile_distance_kpc=args.profile_distance_kpc,
        bins=args.surface_brightness_bins,
        dpi=args.dpi,
    )

    export_df = build_export_dataframe(df)
    write_ini(export_df, export_path)

    print(f"Loaded snapshot {args.snapshot_num} at t={tsnap:.4f} Gyr")
    print(f"Identified dwarf stars: {summary['star_count']}")
    print(f"Identified dwarf gas: {summary['gas_count']}")
    print(f"Sigma x: {summary['sigma_x']:.4f} km/s")
    print(f"Sigma y: {summary['sigma_y']:.4f} km/s")
    print(f"Sigma z: {summary['sigma_z']:.4f} km/s")
    print(f"Saved components figure: {components_path}")
    print(f"Saved surface brightness figure: {surface_brightness_path}")
    print(f"Saved ini export: {export_path}")


if __name__ == "__main__":
    main()
