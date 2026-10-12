import os
import re
import numpy as np
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.stats import zscore
from configparser import RawConfigParser

try:
    import agama  # noqa: F401
except ImportError:
    agama = None

GAS_DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "gas_particle_density",
    [
        (0.00, "#1d234f"),
        (0.10, "#383071"),
        (0.22, "#6b3c86"),
        (0.36, "#a24a78"),
        (0.50, "#d06055"),
        (0.64, "#ee9148"),
        (0.78, "#f6c85f"),
        (0.90, "#fff0a8"),
        (0.97, "#fffbe5"),
        (1.00, "#ffffff"),
    ],
)

STELLAR_DENSITY_CMAP = ListedColormap(
    [
        "#070707",
        "#0d1018",
        "#111a2d",
        "#18345c",
        "#2b5f92",
        "#5d8fb6",
        "#99bfd3",
        "#c5d8dd",
        "#dbe2df",
        "#ede8d8",
        "#faf4e8",
        "#ffffff",
    ],
    "stellar_particle_density",
)

GAS_DENSITY_CMAP.set_bad((0.0, 0.0, 0.0, 0.0))
STELLAR_DENSITY_CMAP.set_bad((0.0, 0.0, 0.0, 0.0))


def smooth_histogram2d(hist, passes=2):
    """Lightly smooth particle-count maps without changing the plotted extent."""
    smoothed = np.asarray(hist, dtype=float)
    if passes <= 0 or smoothed.size == 0:
        return smoothed

    kernel = np.array([1.0, 4.0, 6.0, 4.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    for _ in range(passes):
        padded_x = np.pad(smoothed, ((0, 0), (2, 2)), mode="edge")
        smoothed = (
            kernel[0] * padded_x[:, :-4]
            + kernel[1] * padded_x[:, 1:-3]
            + kernel[2] * padded_x[:, 2:-2]
            + kernel[3] * padded_x[:, 3:-1]
            + kernel[4] * padded_x[:, 4:]
        )
        padded_y = np.pad(smoothed, ((2, 2), (0, 0)), mode="edge")
        smoothed = (
            kernel[0] * padded_y[:-4, :]
            + kernel[1] * padded_y[1:-3, :]
            + kernel[2] * padded_y[2:-2, :]
            + kernel[3] * padded_y[3:-1, :]
            + kernel[4] * padded_y[4:, :]
        )
    return smoothed


def make_particle_density_image(x, y, extent, bins, smooth_passes=2):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    view = (np.abs(x) <= extent) & (np.abs(y) <= extent)
    x_view = x[view]
    y_view = y[view]
    if x_view.size == 0:
        return None

    x_edges = np.linspace(-extent, extent, bins + 1)
    y_edges = np.linspace(-extent, extent, bins + 1)
    hist, _, _ = np.histogram2d(x_view, y_view, bins=[x_edges, y_edges])
    hist = smooth_histogram2d(hist.T, passes=smooth_passes)
    if not np.any(hist > 0):
        return None

    image = np.full_like(hist, np.nan, dtype=float)
    image[hist > 0] = np.log10(hist[hist > 0])
    return image


def plot_particle_density_panel(ax, x, y, title, xlabel, ylabel, cmap, extent=10, bins=450, smooth_passes=2):
    image = make_particle_density_image(x, y, extent=extent, bins=bins, smooth_passes=smooth_passes)

    ax.set_facecolor("black")
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=10, pad=5)
    ax.tick_params(axis="both", which="both", direction="in", color="white", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("white")

    if image is None:
        ax.text(
            0.5,
            0.5,
            "No particles in view",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            bbox=dict(facecolor="black", alpha=0.8, edgecolor="white", linewidth=0.4),
        )
        return None

    values = image[np.isfinite(image)]
    vmin = float(np.percentile(values, 5.0))
    vmax = float(np.percentile(values, 99.5))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmin = float(np.nanmin(values))
        vmax = float(np.nanmax(values))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmax = float(values[0])
        vmin = vmax - 1.0

    return ax.imshow(
        image,
        origin="lower",
        extent=[-extent, extent, -extent, extent],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="bicubic",
    )


def overlay_particle_density_contours(ax, x, y, extent=20, bins=350, smooth_passes=2, color="white"):
    image = make_particle_density_image(x, y, extent=extent, bins=bins, smooth_passes=smooth_passes)
    if image is None:
        return

    values = image[np.isfinite(image)]
    if values.size < 4:
        return

    levels = np.unique(np.percentile(values, [65.0, 78.0, 90.0, 97.0]))
    if levels.size < 2:
        return

    centers = np.linspace(-extent + extent / bins, extent - extent / bins, bins)
    masked = np.ma.masked_invalid(image)
    ax.contour(centers, centers, masked, levels=levels, colors=color, linewidths=0.6, alpha=0.65)

def get_density_components_from_ini(ini_path):
    """Parse ini file to extract all [Density_...] component names"""
    ini = RawConfigParser()
    ini.read(ini_path)
    return [section[8:] for section in ini.sections() if section.lower().startswith("density_")]

def classify_component(name):
    """Classify component type based on name"""
    name = name.lower()
    if "gas" in name:
        return {"type": 0, "temp": 0}
    elif "star" in name or "stellar" in name:
        return {"type": 2, "temp": 0}
    elif "dm" in name or "dark" in name:
        return {"type": 1, "temp": 1}
    else:
        raise ValueError("please check the component name")

def load_and_merge_components():
    """Load and merge all component data
    
    IMPORTANT: Components are merged in type order (0, 1, 2, ...):
    - All type=0 (gas) components first
    - Then all type=1 (dark matter) components  
    - Then all type=2 (stellar) components
    This ensures consistent particle ordering in the merged dataframe.
    """
    frames_by_type = {0: [], 1: [], 2: []}  # type -> list of dataframes
    
    for filename in os.listdir():
        match = re.match(r'model_component_(.+)\.nbody', filename)
        if match:
            comp_name = match.group(1)
            tags = classify_component(comp_name)
            df = pd.read_csv(filename, delimiter='\t')
            df.rename(columns=lambda x: x.lstrip('#') if isinstance(x, str) and x.startswith('#') else x, inplace=True)
            df = df.assign(temp=tags["temp"], type=tags["type"])
            frames_by_type[tags["type"]].append(df)
            print(f"✅ Loaded: {filename} → type={tags['type']}, temp={tags['temp']}")
    
    # Merge in type order: 0 → 1 → 2
    frames = []
    for type_id in sorted(frames_by_type.keys()):
        frames.extend(frames_by_type[type_id])
    
    return pd.concat(frames, ignore_index=True) if frames else None

def process_components(model_name, ini_path):
    """Process model components and generate plots"""
    # Load components
    components = get_density_components_from_ini(ini_path)
    if not components:
        print(f"No components found in {ini_path}")
        return None
        
    # Load and merge data
    merged_data = load_and_merge_components()
    if merged_data is None:
        return None
        
    # Generate plots
    gas_data = merged_data[merged_data['type'] == 0]
    stellar_data = merged_data[merged_data['type'] == 2]
    
    generate_component_plots(model_name, gas_data, stellar_data)
    generate_mass_relation_plot(model_name, gas_data, stellar_data)
    generate_density_profile(model_name, gas_data, stellar_data)
    
    return merged_data

def save_ini_file(merged_data, current_folder):

    """Save final ini file with particle counts"""
    type_0_count = (merged_data['type'] == 0).sum()
    type_1_count = (merged_data['type'] == 1).sum()
    type_2_count = (merged_data['type'] == 2).sum()
    
    header_line = f"{type_0_count} {type_1_count} {type_2_count} 0"
    
    with open(f'./{current_folder}.ini', 'w') as f:
        f.write(header_line + '\n')
        merged_data.to_csv(f, sep=' ', index=True, header=False, mode='a')

def remove_outer_shell(data, keep_fraction=0.95):

    center = np.mean(data, axis=0)
    distances = np.linalg.norm(data - center, axis=1)
    threshold = np.quantile(distances, keep_fraction)
    return distances <= threshold

def generate_component_plots(model_name, gas_data, stellar_data, sight_field=10):
    """Generate projected particle-density maps."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8), dpi=160)
    plt.subplots_adjust(wspace=0.15, hspace=0.18)

    # Gas X-Z plot
    plot_particle_density_panel(
        ax1,
        gas_data["x"],
        gas_data["z"],
        "Gas Particle Density (x-z)",
        "x [kpc]",
        "z [kpc]",
        GAS_DENSITY_CMAP,
        extent=sight_field,
        bins=450,
        smooth_passes=2,
    )

    # Gas X-Y plot
    plot_particle_density_panel(
        ax2,
        gas_data["x"],
        gas_data["y"],
        "Gas Particle Density (x-y)",
        "x [kpc]",
        "y [kpc]",
        GAS_DENSITY_CMAP,
        extent=sight_field,
        bins=450,
        smooth_passes=2,
    )

    # Stellar X-Z plot
    plot_particle_density_panel(
        ax3,
        stellar_data["x"],
        stellar_data["z"],
        "Stellar Particle Density (x-z)",
        "x [kpc]",
        "z [kpc]",
        STELLAR_DENSITY_CMAP,
        extent=sight_field,
        bins=700,
        smooth_passes=2,
    )

    # Stellar X-Y plot
    plot_particle_density_panel(
        ax4,
        stellar_data["x"],
        stellar_data["y"],
        "Stellar Particle Density (x-y)",
        "x [kpc]",
        "y [kpc]",
        STELLAR_DENSITY_CMAP,
        extent=sight_field,
        bins=700,
        smooth_passes=2,
    )

    plt.tight_layout(pad=0.25, w_pad=0.2, h_pad=0.25)
    plt.savefig(f"{model_name}_components.png")
    plt.close()

def generate_mass_relation_plot(model_name, gas_data, stellar_data):
    """Generate mass relation plot"""
    gas_mass = gas_data['mass'].sum() * 1e10 / 1e8
    stellar_mass = stellar_data['mass'].sum() * 1e10 / 1e6

    log_M_star = np.linspace(7.5, 11.5, 8)
    log_M_HI_values = np.where(
        log_M_star <= 9,
        0.712 * log_M_star + 3.117,
        0.276 * log_M_star + 7.042
    )

    plt.figure(figsize=(8, 8))
    plt.scatter(log_M_star, log_M_HI_values, marker="s", color='black', s=80)
    plt.plot(log_M_star, log_M_HI_values, color='black', lw=2)
    plt.scatter(np.log10(1*1e6), np.log10(2*1e8), s=300, color="blue", label="Crater II", marker='*')
    plt.scatter(np.log10(7*1e6), np.log10(4.5*1e8), s=300, color="green", label="Antlia II", marker='*')
    plt.scatter(np.log10(20*1e6), np.log10(6*1e8), s=300, color="red", label="Sculptor", marker='*')
    plt.scatter(np.log10(stellar_mass*1e6), np.log10(gas_mass*1e8), s=300, color="blueviolet", label="Fornax", marker='*')

    plt.xlabel(r'$\log M_*$ (Stellar Mass)', fontsize=14)
    plt.ylabel(r'$\langle \log M_{\mathrm{HI}} \rangle$ (Hydrogen Mass)', fontsize=14)
    plt.title('Piecewise Relation Between Stellar Mass and Neutral Hydrogen Mass')
    plt.legend()
    plt.grid(True)
    plt.minorticks_on()
    plt.xlim(5.8, 12.5)
    plt.ylim(6.8, 11.5)
    plt.tight_layout()
    plt.savefig(f"{model_name}_Mstar_Mgas_relation.png")
    plt.close()

def compute_half_mass_radius(R_proj_cm, mass_g):
    """Compute half-mass radius"""
    sorted_indices = np.argsort(R_proj_cm)
    R_sorted = R_proj_cm[sorted_indices]
    m_sorted = mass_g[sorted_indices]
    m_cumsum = np.cumsum(m_sorted)
    half_mass = m_cumsum[-1] / 2
    idx_half = np.searchsorted(m_cumsum, half_mass)
    R_half_cm = R_sorted[idx_half]
    return R_half_cm / 3.086e21

def generate_density_profile(model_name, gas_data, stellar_data):
    """Generate density profile plot"""

    wlm_r = np.load('../radiusbin.npy')
    wlm_f = np.load('../flux.npy')

    x, y, z, m = gas_data[["x", "y", "z", "mass"]].values.T
    r3d = np.sqrt(x**2 + y**2 + z**2)
    mask = r3d <= 10.0
    x, y, z, m = x[mask], y[mask], z[mask], m[mask] * 1e10
    r = np.sqrt(x**2 + y**2)

    xs, ys, zs, ms = stellar_data[["x", "y", "z", "mass"]].values.T
    r3d_star = np.sqrt(xs**2 + ys**2 + zs**2)
    mask_star = r3d_star <= 10.0
    xs, ys, zs, ms = xs[mask_star], ys[mask_star], zs[mask_star], ms[mask_star] * 1e10
    r_star = np.sqrt(xs**2 + ys**2)

    # Use same bins as gas

    r_max = 10.0
    nbin = 100
    bins = np.linspace(0, r_max, nbin + 1)
    r_cent = 0.5 * (bins[:-1] + bins[1:])
    mass_in_bin, _ = np.histogram(r, bins=bins, weights=m)
    area = np.pi * (bins[1:]**2 - bins[:-1]**2)
    Sigma = mass_in_bin / area

    Sigma_H = Sigma * 1.25e-6 

    mass_in_bin_star, _ = np.histogram(r_star, bins=bins, weights=ms)
    Sigma_star = mass_in_bin_star / area
    Sigma_H_star = Sigma_star * 1.25e-6  # Same conversion factor for comparison

    # ======================================================
    # --- Plot & Save Figure ---
    # ======================================================
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=150)

    # Left: spatial density map
    ax = axes[0]
    plot_particle_density_panel(
        ax,
        gas_data.x,
        gas_data.y,
        "Gas Particle Density (x-y)",
        "x [kpc]",
        "y [kpc]",
        GAS_DENSITY_CMAP,
        extent=20,
        bins=650,
        smooth_passes=2,
    )
    overlay_particle_density_contours(ax, stellar_data.x, stellar_data.y, extent=20, bins=500, color="white")
    
    # Calculate gas fraction within 1kpc
    r3d_inner = np.sqrt(x**2 + y**2 + z**2)
    mask_inner = r3d_inner <= 1.0
    inner_gas_mass = m[mask_inner].sum()
    
    r3d_inner_star = np.sqrt(xs**2 + ys**2 + zs**2)
    mask_inner_star = r3d_inner_star <= 1.0
    inner_star_mass = ms[mask_inner_star].sum()
    
    gas_fraction = inner_gas_mass / (inner_gas_mass + inner_star_mass) * 100
    
    ax.text(0.05, 0.95, 
            f'Gas fraction (<1kpc):\n{gas_fraction:.1f}%',
            transform=ax.transAxes,
            verticalalignment='top',
            fontsize=8,
            color="white",
            bbox=dict(facecolor='black', alpha=0.75, edgecolor='white', linewidth=0.4))
            
    ax.tick_params(axis='both', which='both', direction='in')

    # Right: Σ(R)
    ax = axes[1]
    ax.plot(r_cent, Sigma_H, lw=2., c="cyan", label=f"Gas")
    ax.plot(r_cent, Sigma_H_star, lw=2., c="black", label="Stars")
    ax.plot(wlm_r, wlm_f, label="WLM (Neel et al.)", color="orange", lw=2, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R\ [{\rm kpc}]$")
    ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
    ax.set_xlim(0, 10)
    ax.set_ylim(0.01, 55)
    ax.set_title("Radial Surface Density Profile")
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend()
    ax.tick_params(axis='both', which='both', direction='in')

    fig.suptitle(f"{model_name} Initial Density Profile")

    plt.tight_layout()

    plt.savefig(f"{model_name}_column_density_profile.png")
    plt.close()
