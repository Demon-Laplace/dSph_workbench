import argparse
import math
import os

import PlotFig
import matplotlib.pyplot as plt
import numpy as np
from snapshot_context import prepare_snapshot_context


PANEL_ORDER = [
    ((0, 0), "stellar_sky"),
    ((0, 1), "vlos_sky"),
    ((0, 2), "sigma_history"),
    ((0, 3), "scene_yz"),
    ((1, 0), "velocity_dispersion"),
    ((1, 1), "gas_sky"),
    ((1, 2), "gas_fraction"),
    ((1, 3), "orbit_history"),
    ((2, 0), "surface_brightness"),
    ((2, 1), "sfr_history"),
    ((2, 2), "gas_loss"),
    ((2, 3), "info"),
]

PANEL_INDEX_BY_COORD = {coord: index for index, (coord, _) in enumerate(PANEL_ORDER)}
PANEL_NAME_BY_COORD = {coord: name for coord, name in PANEL_ORDER}


def parse_panel_token(token):
    token = token.strip()
    if not token:
        raise ValueError("empty panel token")

    if token.lower() == "all":
        return list(range(len(PANEL_ORDER)))

    if "," in token:
        parts = token.split(",")
        if len(parts) != 2:
            raise ValueError(f"invalid panel coordinate '{token}'")
        row, col = (int(part.strip()) for part in parts)
        coord = (row, col)
        if coord not in PANEL_INDEX_BY_COORD:
            raise ValueError(f"panel coordinate out of range: '{token}'")
        return [PANEL_INDEX_BY_COORD[coord]]

    index = int(token)
    if index < 0 or index >= len(PANEL_ORDER):
        raise ValueError(f"panel index out of range: '{token}'")
    return [index]


def parse_panel_selection(panel_tokens):
    if not panel_tokens:
        return list(range(len(PANEL_ORDER)))

    selected = []
    seen = set()
    for token in panel_tokens:
        for index in parse_panel_token(token):
            if index not in seen:
                selected.append(index)
                seen.add(index)
    return selected


def print_panel_table():
    print("Available panels:")
    for index, (coord, name) in enumerate(PANEL_ORDER):
        print(f"  {index:2d}  {coord[0]},{coord[1]}  {name}")


def make_output_path(numsp, modelname, output_dir, output_path=None):
    if output_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
        return os.path.abspath(output_path)

    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, f"QuickCheck_{modelname}_{numsp}.png")


def draw_info_panel(ax, context, metrics):
    star_mass = metrics["star_mass"]
    gas_frac_contour = metrics["cold_gas_mass_contour"] / (metrics["cold_gas_mass_contour"] + star_mass)

    info_orbit_text = (
        f"Orbit : {context['ini_x']:.2f}, {context['ini_y']:.2f}, {context['ini_z']:.2f}\n"
        f"        {context['ini_vx']:.2f}, {context['ini_vy']:.2f}, {context['ini_vz']:.2f}\n"
    )

    ax.axis("off")
    ax.text(
        0.02,
        0.95,
        info_orbit_text,
        fontsize=15,
        ha="left",
        va="top",
        family="monospace",
    )
    ax.text(
        0.02,
        0.70,
        f"$R_{{half}}$: {metrics['r_half']:.2f} kpc\n"
        f"$D_{{peri}}$: {context['d_peri']:.2f} kpc\n"
        f"$Stellar\\ Mass$: {PlotFig.DataProcessor.sci_notation(metrics['star_mass'])} $M_{{\\odot}}$\n"
        f"$Stellar\\ Mass_{{Rhalf}}$: {PlotFig.DataProcessor.sci_notation(metrics['star_half_mass'])} $M_{{\\odot}}$\n"
        f"$Cold\\ Gas\\ Mass$: {PlotFig.DataProcessor.sci_notation(metrics['coldgas_mass'])} $M_{{\\odot}}$\n"
        f"$Cold\\ Gas\\ Mass_{{Rhalf}}$: {PlotFig.DataProcessor.sci_notation(metrics['coldgas_half_mass'])} $M_{{\\odot}}$\n"
        f"Gas Fraction in Contour: {gas_frac_contour:.3f}\n"
        f"Gas Fraction in {PlotFig.stellar_region_description}: {metrics['gas_frac_stellar_region']:.3f}\n"
        f"Gas Fraction in $R_{{half}}$: {metrics['gas_half_frac']:.3f}\n"
        f"Feedback: {context['fd_factor_str']} median feedback",
        fontsize=14,
        ha="left",
        va="top",
    )


def remove_sigma_xyz_from_axis(ax):
    sigma_xyz_label = r'$\sigma_{xyz}$'
    filtered_handles = []
    filtered_labels = []

    for line in list(ax.lines):
        if line.get_label() == sigma_xyz_label:
            line.remove()

    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if label != sigma_xyz_label:
            filtered_handles.append(handle)
            filtered_labels.append(label)

    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    if filtered_handles:
        ax.legend(filtered_handles, filtered_labels, loc='upper left', fontsize=15)


def remove_velocity_dispersion_markers(ax):
    red_rgba = np.array([1.0, 0.0, 0.0, 1.0])
    gray_rgba = np.array([0.50196078, 0.50196078, 0.50196078, 1.0])

    for collection in list(ax.collections):
        if hasattr(collection, "get_offsets"):
            offsets = collection.get_offsets()
            facecolors = collection.get_facecolors()
            if len(offsets) == 1 and len(facecolors) > 0 and np.allclose(facecolors[0], red_rgba):
                collection.remove()
                continue

        if hasattr(collection, "get_segments"):
            segments = collection.get_segments()
            colors = collection.get_colors()
            if len(segments) == 1 and len(colors) > 0 and np.allclose(colors[0], gray_rgba):
                collection.remove()


def add_munoz_2018_points(ax):
    ax.scatter(0.791, 24.77, facecolors="none", edgecolors="orange", s=50, zorder=4, lw=1.5)
    ax.scatter(0.05, 23.59, facecolors="none", edgecolors="orange", s=50, label="Munoz 2018", zorder=4, lw=1.5)
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    ax.legend(handles, labels, fontsize=12, loc='upper right')


def draw_panel(fig, ax, panel_index, context, snapshot, metrics, row, numsp):
    coord = PANEL_ORDER[panel_index][0]
    tsnap = snapshot["tsnap"]

    if coord == (0, 0):
        PlotFig._plot_stellar_sky_panel(
            fig,
            ax,
            metrics["img"],
            metrics["extent"],
            metrics["cold_gas_center_ra"],
            metrics["cold_gas_center_dec"],
            metrics["x_contour"],
            metrics["y_contour"],
            metrics["cd_smooth"],
            PlotFig.DataProcessor.sci_notation(PlotFig.cd_threshold),
            metrics["mean_pmra"],
            metrics["mean_pmdec"],
            metrics["r_half_deg"],
            metrics["eps"],
            metrics["pa"],
            metrics["shape_center_ra"],
            metrics["shape_center_dec"],
        )
    elif coord == (0, 1):
        PlotFig._plot_vlos_panel(fig, ax, metrics["vf"], metrics["extent_vf"])
    elif coord == (0, 2):
        PlotFig._plot_sigma_history_panel(
            ax,
            context["dw_elinfo"],
            context["sigma_mw_theoretical"],
            context["sigma_mw"],
            tsnap,
            row["sigma"],
            context["t_peri"],
            metrics["sigma_x"],
            metrics["sigma_y"],
            metrics["sigma_z"],
        )
        remove_sigma_xyz_from_axis(ax)
    elif coord == (0, 3):
        PlotFig._plot_scene_panel(
            ax,
            metrics["cold_gas_y"],
            metrics["cold_gas_z"],
            metrics["hot_gas_y"],
            metrics["hot_gas_z"],
            metrics["star_y"],
            metrics["star_z"],
            metrics["cy"],
            metrics["cz"],
        )
    elif coord == (1, 0):
        PlotFig._plot_velocity_dispersion_panel(
            ax,
            metrics["bin_centers"],
            metrics["vlos_dispersion"],
            row["tsigma"],
            metrics["r_half_circularized"],
        )
        remove_velocity_dispersion_markers(ax)
    elif coord == (1, 1):
        PlotFig._plot_gas_sky_panel(
            ax,
            snapshot["cold_gas_x_kpc"],
            snapshot["cold_gas_y_kpc"],
            snapshot["hot_gas_x_kpc"],
            snapshot["hot_gas_y_kpc"],
        )
    elif coord == (1, 2):
        PlotFig._plot_gas_fraction_panel(
            ax,
            context["dw_elinfo"],
            tsnap,
            metrics["gas_half_ratio"],
            metrics["gas_ratio_stellar_region"],
            context["t_peri"],
        )
    elif coord == (1, 3):
        PlotFig._plot_orbit_panel(
            ax,
            context["dw_elinfo"],
            tsnap,
            row["distance"],
            context["t_peri"],
            snapshot["d_mean_gal"],
            metrics["vr"],
            metrics["vtheta"],
            metrics["vphi"],
            metrics["vgsr"],
        )
    elif coord == (2, 0):
        PlotFig._plot_surface_brightness_profile(
            ax,
            context["dprof"],
            context["mag_obs"],
            context["mag_obs_err"],
            metrics["r_centers"],
            metrics["mag_profile"],
        )
        add_munoz_2018_points(ax)
    elif coord == (2, 1):
        PlotFig._plot_sfr_panel(
            ax,
            context["lookback_sim"],
            context["sfr_sim"],
            context["lookback_obs"],
            context["sfr_obs"],
            context["t_peri"],
        )
        ax.set_ylim(0, 100)
    elif coord == (2, 2):
        PlotFig._plot_gas_loss_panel(
            ax,
            context["dw_elinfo"],
            context["gas_loss_rate"],
            tsnap,
            context["t_peri"],
        )
    elif coord == (2, 3):
        draw_info_panel(ax, context, metrics)
    else:
        raise ValueError(f"unsupported panel coordinate: {coord}")

    ax.set_title(PANEL_NAME_BY_COORD[coord], fontsize=13)


def render_selected_panels(numsp, selected_indices, dpi=150, output_dir=None, output_path=None):
    context = PlotFig.get_plot_context()
    output_dir = os.getcwd() if output_dir is None else os.path.abspath(output_dir)

    snapshot = prepare_snapshot_context(
        folder_path=PlotFig.folder_path,
        snapshot_num=numsp,
        core_radius=context["core_radius"],
    )
    metrics = PlotFig._compute_plot_snapshot_metrics(snapshot, context, numsp)
    row = context["elinfo_by_numsp"].loc[numsp]

    panel_count = len(selected_indices)
    ncols = min(4, max(1, math.ceil(math.sqrt(panel_count))))
    nrows = math.ceil(panel_count / ncols)

    fig = plt.figure(figsize=(6.2 * ncols, 5.2 * nrows), dpi=dpi)
    grid = fig.add_gridspec(nrows, ncols)
    fig.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.12, wspace=0.38, hspace=0.34)
    fig.suptitle(f"{context['modelname']}  numsp={numsp}", fontsize=18)

    for slot, panel_index in enumerate(selected_indices):
        ax = fig.add_subplot(grid[slot // ncols, slot % ncols])
        draw_panel(fig, ax, panel_index, context, snapshot, metrics, row, numsp)

    output_file = make_output_path(
        numsp=numsp,
        modelname=context["modelname"],
        output_dir=output_dir,
        output_path=output_path,
    )
    fig.savefig(output_file, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Render selected PlotFig panels for one snapshot."
    )
    parser.add_argument("numsp", type=int, nargs="?", help="snapshot number")
    parser.add_argument(
        "panels",
        nargs="*",
        help="panel selectors such as 0,0 0,1 or flat indices 0 1 2 ...; default is all panels",
    )
    parser.add_argument("--dpi", type=int, default=150, help="output image DPI")
    parser.add_argument("--output-dir", type=str, default=None, help="directory for the generated image; default is the current directory")
    parser.add_argument("--output", type=str, default=None, help="exact output image path")
    parser.add_argument("--list-panels", action="store_true", help="print available panels and exit")
    args = parser.parse_args()

    if args.list_panels:
        print_panel_table()
        return

    if args.numsp is None:
        parser.error("numsp is required unless --list-panels is used")

    selected_indices = parse_panel_selection(args.panels)
    output_file = render_selected_panels(
        numsp=args.numsp,
        selected_indices=selected_indices,
        dpi=args.dpi,
        output_dir=args.output_dir,
        output_path=args.output,
    )
    print(f"Saved figure: {output_file}")


if __name__ == "__main__":
    main()
