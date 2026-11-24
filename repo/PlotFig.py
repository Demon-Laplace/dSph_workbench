import gc
import json
import os
import ctypes
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import numpy as np
import pandas as pd
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord
# from scipy.interpolate import UnivariateSpline
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LogNorm
import warnings

from basefunc import Analysis, DataProcessor, GalaxySimulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
# from multiprocessing import Pool
from multiprocessing import get_context
from snapshot_context import prepare_snapshot_context
from snapshot_metrics import old_star_projected_kinematics, member_star_projected_kinematics, hi_particle_masses, hi_fraction, require_hi_fraction_csv
from canonical_observer import require_canonical_csv, production_output_root
from wind_tunnel_frame import WIND_CSV_COLUMNS, require_summary_frame, wind_tunnel_frame
from canonical_display import (DISPLAY_RELATIVE_DIR, tangent_coordinates, tangent_velocity_grid40,
    tangent_motion_arrow, tangent_kpc_to_display, raw_offsets_to_tangent, tangent_pixel_area_arcsec2)
from display_outputs import (display_identity, frame_is_current, record_frame,
    record_display_manifest, ordinary_media_reuse, prepare_display_directory)

from matplotlib.patches import Ellipse
from variable import (
    cd_threshold,
    core_radius_dict,
    cube_path,
    d_today,
    err_sigma,
    fndens,
    folder_path,
    obs_plot_mass_to_light_ratio,
    obs_surface_brightness_factor_base,
    output_folder,
    r_pc,
    sim_plot_mass_to_light_ratio,
    sigma,
    stellar_region_description,
    stellar_region_rhalf_multiplier,
    stellar_region_subscript_label,
    stellar_region_tex_label,
    xylabel_size,
)

warnings.filterwarnings("ignore")
plt.switch_backend('Agg')

PLOT_CONTEXT = None
# New contract outputs are isolated from old CSV, PNG and videos.
CANONICAL_ANALYSIS_ROOT = production_output_root(folder_path).resolve()
DISPLAY_OUTPUT_ROOT = CANONICAL_ANALYSIS_ROOT / DISPLAY_RELATIVE_DIR
output_folder = str(DISPLAY_OUTPUT_ROOT / 'output_img')
SKY_FIELD_HALF_WIDTH_DEG = 2.1
GAS_CONTOUR_NBINS = 50
GAS_CONTOUR_SMOOTH_SIGMA_PIXELS = 2.0
GAS_PARTICLE_CONTOUR_PERCENTILES = (35.0, 55.0, 75.0, 90.0, 94.0, 97.0, 99.0)
GAS_PARTICLE_CONTOUR_COLOR = '#303030'
WIND_HI_HALF_WIDTH_KPC = 5.
WIND_HI_GRID_HALF_WIDTH_KPC = 15.
WIND_HI_COLOR_MAX = 100.  # 10^19 cm^-2
WIND_HI_CONTOURS = np.array([1., 3., 6., 10., 30., 100.]) * 1.e19
WIND_SURFACE_DENSITY_LIMITS = (1.e4, 1.e8)  # Msun kpc^-2
WIND_GAS_HISTORY = (
    ('wind_dwarf_gas_r2_msun', r'Dwarf-origin gas, $r<2$ kpc'),
    ('wind_dwarf_gas_r5_msun', r'Dwarf-origin gas, $r<5$ kpc'),
    ('wind_ambient_gas_r5_msun', r'CGM-origin gas, $r<5$ kpc'),
)
WIND_PRESSURE_HISTORY = (
    ('wind_prescribed_Pram_dyn_cm2', 'Input at the box centre'),
    ('wind_upstream12_20_Pram_dyn_cm2', 'Upstream, 12–20 kpc'),
    ('wind_upstream20_40_Pram_dyn_cm2', 'Upstream, 20–40 kpc'),
)
OBS_BACKGROUND_MIN_RADIUS_DEG = 5.0
OBS_BACKGROUND_MIN_BINS = 3
OBS_BACKGROUND_FALLBACK_RADIUS_QUANTILE = 0.9
OBS_SURFACE_BRIGHTNESS_PROFILE_TEMPLATE = "{dwarf_name}_surface_brightness_profile.csv"
ADAPTIVE_ML_SEGMENT_WEIGHTS = (5.0, 3.0, 1.0)
ADAPTIVE_ML_FIT_MAX_RADIUS_KPC = 5.0
ADAPTIVE_ML_MIN = 0.1
ADAPTIVE_ML_MAX = 50.0
SIM_SB_MARKER_SPACING_KPC = 0.2
SFH_SFR_UNIT_MSUN_PER_YEAR = 1.0e-4
PLOT_VIDEO_VERSION = 5

try:
    LIBC = ctypes.CDLL("libc.so.6")
except OSError:
    LIBC = None

try:
    from threadpoolctl import threadpool_limits
except ImportError:
    threadpool_limits = None

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

THREADPOOL_CONTROLLER = None
PNG_SIGNATURE = b'\x89PNG\r\n\x1a\n'
PNG_IEND_CHUNK = b'\x00\x00\x00\x00IEND\xaeB`\x82'


def make_plot_output_path(output_dir, numsp):
    return os.path.join(output_dir, f'check_nature_{numsp:03d}.png')


def is_complete_png(path):
    try:
        if not os.path.isfile(path):
            return False
        if os.path.getsize(path) < len(PNG_SIGNATURE) + len(PNG_IEND_CHUNK):
            return False
        with open(path, 'rb') as handle:
            if handle.read(len(PNG_SIGNATURE)) != PNG_SIGNATURE:
                return False
            handle.seek(-len(PNG_IEND_CHUNK), os.SEEK_END)
            return handle.read(len(PNG_IEND_CHUNK)) == PNG_IEND_CHUNK
    except OSError:
        return False


def reusable_frame(path, identity):
    if ordinary_media_reuse():
        return frame_is_current(path, identity)
    return is_complete_png(path) and frame_is_current(path, identity)


def save_figure_atomic(fig, output_file):
    output_dir = os.path.dirname(os.path.abspath(output_file))
    os.makedirs(output_dir, exist_ok=True)
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix='.png',
            prefix='.plotfig-',
            dir=output_dir,
            delete=False,
        ) as handle:
            tmp_path = handle.name

        fig.savefig(tmp_path)
        if not ordinary_media_reuse() and not is_complete_png(tmp_path):
            raise RuntimeError(f"Failed to write a complete PNG: {tmp_path}")
        os.replace(tmp_path, output_file)
    finally:
        if tmp_path is not None and os.path.exists(tmp_path):
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


def render_video_from_frames(output_dir, modelname, start_numsp, frame_count,
                             *, snapshots=None, video_name=None):
    if frame_count <= 0:
        print("[PlotFigVideo] Skip video generation: no frames requested")
        return None

    ffmpeg_path = shutil.which("ffmpeg")
    if ffmpeg_path is None:
        if ordinary_media_reuse():
            raise FileNotFoundError('ffmpeg is required to encode the requested video')
        print("[PlotFigVideo] Skip video generation: ffmpeg not found in PATH")
        return None

    frames_dir = os.path.abspath(output_dir)
    video_dir = os.path.dirname(os.path.normpath(frames_dir))
    output_name = video_name or f"{modelname}_tan_xi_eta_2p1_v1_v{PLOT_VIDEO_VERSION}.mp4"
    if Path(output_name).name != output_name or not output_name.endswith('.mp4'):
        raise ValueError('--video-name must be a plain .mp4 filename')
    output_path = os.path.join(video_dir, output_name)
    with tempfile.NamedTemporaryFile(
        suffix='.mp4',
        prefix=f".{output_name}.tmp.",
        dir=video_dir,
        delete=False,
    ) as handle:
        tmp_output_path = handle.name
    sequence_file = None
    if snapshots is not None:
        with tempfile.NamedTemporaryFile('w', encoding='utf-8', suffix='.txt',
                                        prefix='.plotfig-sequence-', dir=video_dir,
                                        delete=False) as handle:
            sequence_file = handle.name
            for number in snapshots:
                name = Path(make_plot_output_path(frames_dir, number)).as_posix()
                handle.write("file '" + name.replace("'", "'\\''") + "'\n")
                handle.write('duration 0.04\n')
        input_args = ['-f', 'concat', '-safe', '0', '-i', sequence_file]
    else:
        input_args = ['-start_number', str(start_numsp), '-i', 'check_nature_%03d.png']
    cmd = [ffmpeg_path] + input_args + [
        "-vf",
        "scale=3416:trunc(ow/a/2)*2,setsar=1",
        "-frames:v",
        str(frame_count),
        "-r",
        "25",
        "-crf",
        "18",
        "-pix_fmt",
        "yuv420p",
        "-y",
        tmp_output_path,
    ]

    try:
        subprocess.run(
            cmd,
            cwd=frames_dir,
            check=True,
            stdout=None if ordinary_media_reuse() else subprocess.DEVNULL,
            stderr=None if ordinary_media_reuse() else subprocess.STDOUT,
        )
        os.replace(tmp_output_path, output_path)
    except subprocess.CalledProcessError as exc:
        if ordinary_media_reuse():
            raise
        print(f"[PlotFigVideo] ffmpeg failed with exit code {exc.returncode}")
        return None
    finally:
        if sequence_file is not None:
            Path(sequence_file).unlink(missing_ok=True)
        if os.path.exists(tmp_output_path):
            try:
                os.unlink(tmp_output_path)
            except OSError:
                pass

    print(f"[PlotFigVideo] Saved video: {output_path}")
    return output_path


def finalize_display_outputs(output_dir, modelname, snapshots, *, identity=None,
                             video_name=None, encode_video=True):
    if not snapshots:
        raise ValueError('No snapshots selected for video generation')
    if identity is None:
        identity = get_plot_context()['display_identity']
    if video_name is None and PLOT_CONTEXT is not None:
        video_name = PLOT_CONTEXT.get('video_name')
    videos = []
    if encode_video and (snapshots_are_consecutive(snapshots) or ordinary_media_reuse()):
        video = render_video_from_frames(
            output_dir, modelname, snapshots[0], len(snapshots),
            snapshots=snapshots if ordinary_media_reuse() else None,
            video_name=video_name)
        if video is not None:
            videos.append(video)
    elif encode_video:
        print("[PlotFigVideo] Skip video generation: requested snapshots are not consecutive")
    record_display_manifest(Path(output_dir).parent, identity, snapshots, videos)


def build_consecutive_blocks(items, block_size):
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    return [items[i:i + block_size] for i in range(0, len(items), block_size)]


def filter_available_snapshots(start, end, available_snapshots):
    return [numsp for numsp in available_snapshots if start <= numsp <= end]


def format_limited_ints(values, limit=20):
    values = list(values)
    prefix = ",".join(str(value) for value in values[:limit])
    return prefix + ("..." if len(values) > limit else "")


def snapshots_are_consecutive(snapshot_numbers):
    if len(snapshot_numbers) < 2:
        return True
    return all(
        current - previous == 1
        for previous, current in zip(snapshot_numbers, snapshot_numbers[1:])
    )


def validate_plot_outputs(output_dir, requested_numsp):
    missing_numsp = []
    invalid_numsp = []
    for numsp in requested_numsp:
        output_file = make_plot_output_path(output_dir, numsp)
        if is_complete_png(output_file) and frame_is_current(output_file, get_plot_context()['display_identity']):
            continue
        if os.path.exists(output_file):
            invalid_numsp.append(numsp)
        else:
            missing_numsp.append(numsp)

    if invalid_numsp:
        print(
            f"PlotFig validation failed: {len(invalid_numsp)} incomplete PNG file(s): "
            + format_limited_ints(invalid_numsp)
        )
    if missing_numsp:
        print(
            f"PlotFig validation failed: {len(missing_numsp)} missing PNG file(s): "
            + format_limited_ints(missing_numsp)
        )

    if invalid_numsp or missing_numsp:
        return False

    print(f"PlotFig validation: complete PNG frames present ({len(requested_numsp)} requested).")
    return True


class ProgressReporter:
    def __init__(self, total, desc, enabled=True, unit='snap', mininterval=0.5):
        self.total = total
        self.desc = desc
        self.enabled = enabled and total > 0
        self.unit = unit
        self.mininterval = mininterval
        self.count = 0
        self._last_print = 0.0
        self._use_tqdm = self.enabled and tqdm is not None and sys.stdout.isatty()
        self._progress = None
        if self._use_tqdm:
            self._progress = tqdm(
                total=total,
                desc=desc,
                unit=unit,
                dynamic_ncols=True,
                mininterval=mininterval,
                file=sys.stdout,
            )
        elif self.enabled:
            print(f"[{self.desc}Progress] 0/{self.total} {self.unit}", flush=True)

    def update(self, n=1):
        if not self.enabled:
            return
        self.count += n
        if self._progress is not None:
            self._progress.update(n)
            return

        now = time.perf_counter()
        if self.count >= self.total or now - self._last_print >= self.mininterval:
            print(f"[{self.desc}Progress] {self.count}/{self.total} {self.unit}", flush=True)
            self._last_print = now

    def close(self):
        if self._progress is not None:
            self._progress.close()


def release_process_memory():
    gc.collect()
    if LIBC is not None:
        try:
            LIBC.malloc_trim(0)
        except Exception:
            pass


def get_rss_mb():
    try:
        with open("/proc/self/statm", "r") as fh:
            rss_pages = int(fh.readline().split()[1])
        page_size = os.sysconf("SC_PAGE_SIZE")
        return rss_pages * page_size / (1024.0 * 1024.0)
    except Exception:
        if os.name == "nt":
            return float('nan')
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def init_worker(worker_threads):
    global THREADPOOL_CONTROLLER
    if threadpool_limits is not None and worker_threads is not None and worker_threads > 0:
        THREADPOOL_CONTROLLER = threadpool_limits(limits=worker_threads)


def _estimate_observed_density_background(dprof):
    radius_deg = dprof['xr'].to_numpy(dtype=float)
    density = dprof['density'].to_numpy(dtype=float)
    good = np.isfinite(radius_deg) & np.isfinite(density)

    background_mask = good & (radius_deg >= OBS_BACKGROUND_MIN_RADIUS_DEG)
    if np.sum(background_mask) < OBS_BACKGROUND_MIN_BINS and np.any(good):
        fallback_min_radius = np.nanquantile(
            radius_deg[good],
            OBS_BACKGROUND_FALLBACK_RADIUS_QUANTILE,
        )
        background_mask = good & (radius_deg >= fallback_min_radius)

    if not np.any(background_mask):
        raise ValueError("Cannot estimate observed density background: no finite outer bins.")

    if 'area' in dprof.columns:
        area = dprof['area'].to_numpy(dtype=float)
        weighted = background_mask & np.isfinite(area) & (area > 0)
        if np.any(weighted):
            background_density = np.average(density[weighted], weights=area[weighted])
            return float(background_density), int(np.sum(weighted))

    return float(np.nanmedian(density[background_mask])), int(np.sum(background_mask))


def _observed_surface_brightness_from_counts(dprof, factor):
    density = dprof['density'].to_numpy(dtype=float)
    density_error = dprof['density_error'].to_numpy(dtype=float)
    background_density, background_bins = _estimate_observed_density_background(dprof)
    net_density = density - background_density

    mag_obs = np.full_like(net_density, np.nan, dtype=float)
    mag_obs_err = np.full_like(net_density, np.nan, dtype=float)

    valid = (
        np.isfinite(net_density)
        & np.isfinite(density_error)
        & (net_density > 0)
        & (density_error >= 0)
    )
    intensity = net_density[valid] * factor
    intensity_error = density_error[valid] * factor
    mag_obs[valid] = -2.5 * np.log10(intensity)
    mag_obs_err[valid] = 1.0857 * (intensity_error / intensity)

    return {
        'mag_obs': mag_obs,
        'mag_obs_err': mag_obs_err,
        'density_background': background_density,
        'background_bins': background_bins,
        'positive_net_bins': int(np.sum(valid)),
    }


def _observed_profile_candidate_paths(dwarf_name):
    filename = OBS_SURFACE_BRIGHTNESS_PROFILE_TEMPLATE.format(dwarf_name=dwarf_name)
    explicit = os.environ.get("FORNAX_PREPARED_OBSERVATION_PROFILE")
    if explicit:
        return [Path(explicit).expanduser()]
    code = next((p for p in Path(__file__).resolve().parents if p.name == "Code"), None)
    if code is not None:
        if str(code) not in sys.path:
            sys.path.insert(0, str(code))
        from project_paths import OBSERVATION_DATA
        return [OBSERVATION_DATA / "Yang_data" / filename]
    declared = os.environ.get("FORNAX_OBSERVATION_ROOT")
    if declared:
        return [Path(declared).expanduser() / "Yang_data" / filename]
    # Compatibility for explicitly retained standalone server packages.
    return [Path(__file__).resolve().parent.parent / filename]


def _observed_radius_kpc(dprof):
    if 'r_mid_kpc' in dprof.columns:
        return dprof['r_mid_kpc'].to_numpy(dtype=float)
    return np.radians(dprof['xr'].to_numpy(dtype=float)) * d_today


def _observed_profile_intervals_kpc(dprof):
    if {'r_inner_kpc', 'r_outer_kpc'}.issubset(dprof.columns):
        intervals = dprof[['r_inner_kpc', 'r_outer_kpc']].to_numpy(dtype=float)
    elif {'xr1', 'xr2'}.issubset(dprof.columns):
        intervals = np.radians(dprof[['xr1', 'xr2']].to_numpy(dtype=float)) * d_today
    else:
        return None

    good = np.isfinite(intervals).all(axis=1) & (intervals[:, 1] > intervals[:, 0])
    if np.sum(good) < 2:
        return None
    return intervals[good]


def _load_surface_brightness_observation(dwarf_name, factor):
    required = {'r_mid_kpc', 'mu_v_mag_arcsec2', 'mu_v_err_total'}
    for candidate in _observed_profile_candidate_paths(dwarf_name):
        if not candidate.exists():
            continue
        dprof = pd.read_csv(candidate)
        if not required.issubset(dprof.columns):
            raise ValueError(
                f"Preprocessed observation profile {candidate} is missing required columns: "
                f"{sorted(required - set(dprof.columns))}"
            )

        mag_obs = dprof['mu_v_mag_arcsec2'].to_numpy(dtype=float)
        mag_obs_err = dprof['mu_v_err_total'].to_numpy(dtype=float)
        positive = np.isfinite(mag_obs) & np.isfinite(mag_obs_err)
        return {
            'dprof': dprof,
            'mag_obs': mag_obs,
            'mag_obs_err': mag_obs_err,
            'density_background': (
                float(dprof['background_density_stars_arcmin2'].dropna().iloc[0])
                if 'background_density_stars_arcmin2' in dprof.columns
                and dprof['background_density_stars_arcmin2'].notna().any()
                else np.nan
            ),
            'background_bins': (
                int(dprof['background_bin_count'].dropna().iloc[0])
                if 'background_bin_count' in dprof.columns
                and dprof['background_bin_count'].notna().any()
                else 0
            ),
            'positive_net_bins': int(np.sum(positive)),
            'source_path': str(candidate),
            'preprocessed': True,
        }

    dprof = pd.read_csv(fndens)
    obs_surface_brightness = _observed_surface_brightness_from_counts(dprof, factor)
    obs_surface_brightness.update(
        {
            'dprof': dprof,
            'source_path': fndens,
            'preprocessed': False,
        }
    )
    return obs_surface_brightness


def _normalize_elinfo_schema(dw_elinfo):
    """Supply plotting-only defaults for elinfo files made by older GetInfo versions."""
    required_hi = {'hi_mass', 'hi_half_mass', 'hi_fraction', 'hi_half_fraction'}
    if not required_hi.issubset(dw_elinfo.columns):
        raise ValueError("Re-run GetInfo: this elinfo predates corrected H I mass statistics")
    for fraction, hi, stars in (
        ('hi_fraction', 'hi_mass', 'star_mass'),
        ('hi_half_fraction', 'hi_half_mass', 'star_half_mass'),
    ):
        expected = np.array([hi_fraction(h, s) for h, s in zip(dw_elinfo[hi], dw_elinfo[stars])])
        if not np.allclose(dw_elinfo[fraction], expected, rtol=2e-6, atol=1e-8, equal_nan=True):
            raise ValueError("Re-run GetInfo: stored H I fractions disagree with M_HI/(M_star+M_HI)")
    legacy_fallbacks = []

    if 'rhalf_circularized' not in dw_elinfo.columns:
        if 'rhalf' not in dw_elinfo.columns:
            raise KeyError("elinfo is missing both 'rhalf_circularized' and 'rhalf'")
        dw_elinfo['rhalf_circularized'] = dw_elinfo['rhalf']
        legacy_fallbacks.append("rhalf_circularized=rhalf")

    for column in ('shape_center_x_kpc', 'shape_center_y_kpc'):
        if column not in dw_elinfo.columns:
            dw_elinfo[column] = 0.0
            legacy_fallbacks.append(f"{column}=0")

    if legacy_fallbacks:
        print(
            "[PlotFig] Legacy elinfo compatibility: "
            + ", ".join(legacy_fallbacks)
        )

    return dw_elinfo


def build_plot_context():
    wind_tunnel = wind_tunnel_frame(folder_path)
    model_info = DataProcessor.parse_model_name_details(wind_tunnel.model_name if wind_tunnel else None)
    modelname = model_info['modelname']
    base_modelname = model_info['base_modelname']
    dwarf_name = model_info['dwarf_name']
    num = model_info['model_num']
    description = model_info['description']
    if description:
        print(f"[PlotFig] Model: {base_modelname}, description: {description}")
    # Retain compatibility context fields without requiring unused IC counts.
    inifile, dw_particles = None, None
    elinfo_path = str(CANONICAL_ANALYSIS_ROOT / f'elinfo_{modelname}.csv')

    core_radius = core_radius_dict.get(dwarf_name, None)
    if core_radius is None:
        raise ValueError(f"please input '{dwarf_name}''s radius")

    require_summary_frame(require_canonical_csv(elinfo_path), wind_tunnel)
    require_hi_fraction_csv(elinfo_path)
    dw_elinfo = DataProcessor.read_csv_with_comments(elinfo_path, dtype=DataProcessor.ELINFO_DTYPE_MAP)
    dw_elinfo = _normalize_elinfo_schema(dw_elinfo)
    if wind_tunnel:
        missing = set(WIND_CSV_COLUMNS) - set(dw_elinfo)
        if missing:
            raise ValueError('Run GetInfo for this wind tunnel summary. Missing columns: ' + ', '.join(sorted(missing)))
        # This copy is for display. The CSV and stellar-age calculations keep
        # native GIZMO time so formation times retain their original clock.
        dw_elinfo = dw_elinfo.assign(age=dw_elinfo['wind_elapsed_Gyr'])
    if 'numsp' not in dw_elinfo.columns:
        raise ValueError("The elinfo file does not contain a numsp column")
    elinfo_numsp_values = pd.to_numeric(dw_elinfo['numsp'], errors='coerce').dropna()
    elinfo_snapshots = sorted(set(elinfo_numsp_values.astype(int).tolist()))
    if not elinfo_snapshots:
        raise ValueError("The elinfo file does not contain any valid snapshot numbers")
    elinfo_by_numsp = dw_elinfo.set_index('numsp', drop=False)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
        print(f"Created output folder: {output_folder}")

    ini_x, ini_y, ini_z, ini_vx, ini_vy, ini_vz, fd_value = DataProcessor.get_ini_data(elinfo_path=elinfo_path)
    fd_factor_raw = fd_value / 8.3
    if fd_factor_raw >= 0.90:
        fd_factor_str = f"{round(fd_factor_raw):.0f}"
    else:
        fd_factor_str = f"{fd_factor_raw:.1f}"

    factor = obs_surface_brightness_factor_base
    obs_surface_brightness = _load_surface_brightness_observation(dwarf_name, factor)
    dprof = obs_surface_brightness['dprof']
    mag_obs = obs_surface_brightness['mag_obs']
    mag_obs_err = obs_surface_brightness['mag_obs_err']

    distance_series = dw_elinfo['distance_gal' if wind_tunnel else 'distance']
    sigma_mw = np.where(
        dw_elinfo.sigma >= dw_elinfo.tsigma,
        np.sqrt(dw_elinfo.sigma**2 - dw_elinfo.tsigma**2),
        -np.sqrt(dw_elinfo.tsigma**2 - dw_elinfo.sigma**2),
    )
    sigma_mw_theoretical = Analysis.compute_sigma_mw(
        dw_elinfo['mw_mass_r'].to_numpy(),
        distance_series.to_numpy(),
        dw_elinfo['rhalf_circularized'].to_numpy(dtype=float),
    )
    history_context = {}
    wind_sfh_snapshot = None
    sfh_path = os.environ.get('DSPH_SFH_OBSERVATION', f'../SFH_{dwarf_name}.csv')
    if wind_tunnel:
        history_limit = wind_tunnel.orbit[-1, 0] / 1000.
        matching = set(DataProcessor.list_snapshot_numbers(folder_path)).intersection(elinfo_snapshots)
        if not matching:
            raise FileNotFoundError('No snapshots are present in both the output directory and the elinfo file')
        wind_sfh_snapshot = max(matching)
        history_context.update(
            wind_history_limits=_wind_history_limits(wind_tunnel, dw_elinfo),
            wind_sfh=_wind_reference_sfh(folder_path, wind_sfh_snapshot, core_radius),
        )
        reference_epoch = wind_tunnel.sfh_reference_elapsed_gyr
        if reference_epoch is not None:
            observed_sfh = pd.read_csv(sfh_path)
            history_context['wind_observed_sfh'] = _observed_sfh_on_wind_clock(observed_sfh, reference_epoch)
    else:
        index_peri = distance_series.idxmin()
        t_peri = dw_elinfo.loc[index_peri, 'age']
        history_limit = t_peri
        Tage = dw_elinfo.loc[distance_series.loc[:index_peri].sub(d_today).abs().idxmin(), 'age']
        sfh_dw = pd.read_csv(sfh_path)
        dw_sfrmass = Analysis.GetSFH_sim(folder_path=folder_path, Tage=Tage, snapshot_num=elinfo_snapshots[-1])
        sfr_t_peri_gyr = t_peri * dw_sfrmass.attrs.get('time_code_to_gyr', 1.)
        lookback_sim, sfr_sim, lookback_obs, sfr_obs = _compute_sfr_series(sfh_dw, dw_sfrmass, sfr_t_peri_gyr)
        history_context.update(
            gas_loss_rate=Analysis.calculate_gas_loss_rate(dw_elinfo['age'].to_numpy(), dw_elinfo['coldgas_half_mass'].to_numpy()),
            d_peri=distance_series.loc[index_peri], t_peri=t_peri, Tage=Tage, sfh_dw=sfh_dw,
            dw_sfrmass=dw_sfrmass, sfr_t_peri_gyr=sfr_t_peri_gyr,
            lookback_sim=lookback_sim, sfr_sim=sfr_sim, lookback_obs=lookback_obs, sfr_obs=sfr_obs,
        )

    return {
        'modelname': modelname,
        'base_modelname': base_modelname,
        'dwarf_name': dwarf_name,
        'num': num,
        'description': description,
        'inifile': inifile,
        'dw_particles': dw_particles,
        'elinfo_path': elinfo_path,
        'display_identity': display_identity(elinfo_path, wind_sfh_snapshot=wind_sfh_snapshot),
        'core_radius': core_radius,
        'dw_elinfo': dw_elinfo,
        'elinfo_by_numsp': elinfo_by_numsp,
        'elinfo_snapshots': elinfo_snapshots,
        'ini_x': ini_x,
        'ini_y': ini_y,
        'ini_z': ini_z,
        'ini_vx': ini_vx,
        'ini_vy': ini_vy,
        'ini_vz': ini_vz,
        'fd_value': fd_value,
        'fd_factor_str': fd_factor_str,
        'dprof': dprof,
        'mag_obs': mag_obs,
        'mag_obs_err': mag_obs_err,
        'obs_density_background': obs_surface_brightness['density_background'],
        'obs_background_bins': obs_surface_brightness['background_bins'],
        'obs_positive_net_bins': obs_surface_brightness['positive_net_bins'],
        'obs_surface_brightness_path': obs_surface_brightness['source_path'],
        'obs_surface_brightness_preprocessed': obs_surface_brightness['preprocessed'],
        'distance_series': distance_series,
        'sigma_mw': sigma_mw,
        'sigma_mw_theoretical': sigma_mw_theoretical,
        'wind_tunnel': wind_tunnel,
        'history_limit': history_limit,
        'output_folder': output_folder,
        **history_context,
    }


def get_plot_context():
    global PLOT_CONTEXT
    if PLOT_CONTEXT is None:
        PLOT_CONTEXT = build_plot_context()
    return PLOT_CONTEXT


def safe_plot(numsp, show_mem=False, profile=False, profile_memory=False, verbose=False):
    try:
        return plot(
            numsp,
            show_mem=show_mem,
            profile=profile,
            profile_memory=profile_memory,
            verbose=verbose,
        )
    except Exception as e:
        print(f"Error in numsp={numsp}: {e}")
        raise


def _style_axis(ax, xlabel, ylabel, xlim=None, ylim=None, invert_x=False):
    ax.set_xlabel(xlabel, fontsize=xylabel_size)
    ax.set_ylabel(ylabel, fontsize=xylabel_size)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if invert_x:
        ax.invert_xaxis()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    ax.grid(True)


def _append_colorbar(fig, ax, mappable, label, **kwargs):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax, **kwargs)
    cbar.set_label(label, fontsize=14)
    cax.tick_params(axis='both', which='both', direction='in')
    cax.minorticks_on()
    return cbar


def _adaptive_ml_segment_weights(radius_kpc, r_half):
    radius_kpc = np.asarray(radius_kpc, dtype=float)
    weights = np.full(radius_kpc.shape, ADAPTIVE_ML_SEGMENT_WEIGHTS[2], dtype=float)
    labels = np.full(radius_kpc.shape, '2Rhalf-5kpc', dtype=object)

    first_segment = radius_kpc <= r_half
    second_segment = (radius_kpc > r_half) & (radius_kpc <= 2.0 * r_half)
    outside_fit = radius_kpc > ADAPTIVE_ML_FIT_MAX_RADIUS_KPC

    weights[first_segment] = ADAPTIVE_ML_SEGMENT_WEIGHTS[0]
    weights[second_segment] = ADAPTIVE_ML_SEGMENT_WEIGHTS[1]
    weights[outside_fit] = np.nan

    labels[first_segment] = '0-Rhalf'
    labels[second_segment] = 'Rhalf-2Rhalf'
    labels[outside_fit] = 'outside-fit'
    return weights, labels


def _fit_adaptive_scalar_mass_to_light(r_centers, mag_profile, dprof, mag_obs, mag_obs_err, r_half):
    r_centers = np.asarray(r_centers, dtype=float)
    mag_profile = np.asarray(mag_profile, dtype=float)
    obs_radius_kpc = _observed_radius_kpc(dprof)
    mag_obs = np.asarray(mag_obs, dtype=float)
    mag_obs_err = np.asarray(mag_obs_err, dtype=float)
    r_half = float(r_half)

    if not np.isfinite(r_half) or r_half <= 0:
        return None

    obs_good = np.isfinite(obs_radius_kpc) & np.isfinite(mag_obs) & np.isfinite(mag_obs_err)
    obs_good &= mag_obs_err > 0
    if np.sum(obs_good) < 2:
        return None

    order = np.argsort(obs_radius_kpc[obs_good])
    obs_r = obs_radius_kpc[obs_good][order]
    obs_m = mag_obs[obs_good][order]
    obs_e = mag_obs_err[obs_good][order]

    fit_good = np.isfinite(r_centers) & np.isfinite(mag_profile)
    fit_good &= r_centers >= obs_r.min()
    fit_good &= r_centers <= obs_r.max()
    fit_good &= r_centers <= ADAPTIVE_ML_FIT_MAX_RADIUS_KPC
    fit_r = r_centers[fit_good]
    if fit_r.size < 2:
        return None

    target_mag = np.interp(fit_r, obs_r, obs_m)
    target_err = np.interp(fit_r, obs_r, obs_e)
    segment_multiplier, segment_label = _adaptive_ml_segment_weights(fit_r, r_half)
    weight_good = np.isfinite(segment_multiplier) & (segment_multiplier > 0)
    if np.sum(weight_good) < 2:
        return None

    fit_r = fit_r[weight_good]
    target_mag = target_mag[weight_good]
    target_err = target_err[weight_good]
    fit_mag = mag_profile[fit_good][weight_good]
    segment_multiplier = segment_multiplier[weight_good]
    segment_label = segment_label[weight_good]

    weights = (1.0 / np.square(target_err)) * segment_multiplier
    target_offset_mag = target_mag - fit_mag
    offset_mag = float(np.sum(weights * target_offset_mag) / np.sum(weights))
    adaptive_ml = float(sim_plot_mass_to_light_ratio * np.power(10.0, offset_mag / 2.5))
    adaptive_ml = float(np.clip(adaptive_ml, ADAPTIVE_ML_MIN, ADAPTIVE_ML_MAX))
    offset_mag = float(2.5 * np.log10(adaptive_ml / sim_plot_mass_to_light_ratio))

    adaptive_mag_profile = mag_profile + offset_mag
    residual = fit_mag + offset_mag - target_mag
    within_2rhalf = fit_r <= 2.0 * r_half

    return {
        'mass_to_light': adaptive_ml,
        'mag_offset': offset_mag,
        'mag_profile': adaptive_mag_profile,
        'fit_radius_kpc': fit_r,
        'fit_weight': weights,
        'fit_segment': segment_label,
        'target_mag': target_mag,
        'residual_mag': residual,
        'r_half': r_half,
        'segment_weights': ADAPTIVE_ML_SEGMENT_WEIGHTS,
        'fit_max_radius_kpc': ADAPTIVE_ML_FIT_MAX_RADIUS_KPC,
        'rms_mag': float(np.sqrt(np.nanmean(np.square(residual)))),
        'within_2rhalf_rms_mag': (
            float(np.sqrt(np.nanmean(np.square(residual[within_2rhalf]))))
            if np.any(within_2rhalf)
            else np.nan
        ),
    }


def _gas_contour_levels(cd_smooth, mode="threshold"):
    if cd_smooth.size == 0:
        return None
    finite_density = cd_smooth[np.isfinite(cd_smooth)]
    if finite_density.size == 0:
        return None
    if mode == "threshold":
        density_max = finite_density.max()
        if density_max <= cd_threshold:
            return None
        return np.linspace(cd_threshold, density_max, 5)
    if mode == "morphology":
        positive_density = finite_density[finite_density > 0]
        if positive_density.size == 0:
            return None
        levels = np.percentile(positive_density, GAS_PARTICLE_CONTOUR_PERCENTILES)
        levels = np.unique(levels[np.isfinite(levels) & (levels > 0)])
        return levels if levels.size > 0 else None
    raise ValueError(f"Unsupported gas contour level mode: {mode}")


def _wind_hi_levels(column, minimum=0.):
    finite = column[np.isfinite(column)]
    if not finite.size:
        return np.array([])
    return WIND_HI_CONTOURS[(WIND_HI_CONTOURS >= minimum)
                           & (WIND_HI_CONTOURS > finite.min()) & (WIND_HI_CONTOURS < finite.max())]


def _plot_hi_sky_panel(
    fig,
    ax,
    x_contour_kpc,
    y_contour_kpc,
    cd_smooth,
    field_half_width_kpc,
    *,
    wind_tunnel=False,
):
    column = np.ma.masked_less_equal(np.ma.masked_invalid(cd_smooth.T), 0.)
    if wind_tunnel or column.count():
        image = ax.pcolormesh(x_contour_kpc, y_contour_kpc, column / 1.e19,
                              cmap='BuGn', vmin=0., vmax=WIND_HI_COLOR_MAX if wind_tunnel else None,
                              shading='auto')
        _append_colorbar(fig, ax, image, r'$N_{\mathrm{H\,I}}$ ($10^{19}\,\mathrm{cm}^{-2}$)',
                         extend='max' if wind_tunnel else 'neither')
    if wind_tunnel:
        levels = _wind_hi_levels(cd_smooth)
    else:
        levels = _gas_contour_levels(cd_smooth, mode="morphology")
    if levels is not None and len(levels):
        ax.contour(
            x_contour_kpc,
            y_contour_kpc,
            cd_smooth.T,
            levels=levels,
            colors=GAS_PARTICLE_CONTOUR_COLOR,
            linewidths=.8,
        )
    _style_axis(
        ax,
        'East (kpc)' if wind_tunnel else 'RA (kpc)',
        'North (kpc)' if wind_tunnel else 'DEC (kpc)',
        xlim=(-field_half_width_kpc, field_half_width_kpc),
        ylim=(-field_half_width_kpc, field_half_width_kpc),
        invert_x=True,
    )
    ax.set_aspect('equal', adjustable='box')


def _plot_stellar_sky_panel(
    fig,
    ax,
    img,
    extent,
    cold_gas_center_ra,
    cold_gas_center_dec,
    x_contour,
    y_contour,
    cd_smooth,
    cd_threshold_str,
    mean_pmra,
    mean_pmdec,
    r_half_deg,
    eps,
    pa,
    shape_center_ra=0.0,
    shape_center_dec=0.0,
    arrow_tip=None,
    shape_outline=None,
    gas_center_label='HI Center',
    fixed_gas_levels=False,
):
    arrow_length = 1.5
    norm = np.hypot(mean_pmra, mean_pmdec)
    dx_scaled = 0.0 if norm == 0 else mean_pmra / norm * arrow_length
    dy_scaled = 0.0 if norm == 0 else mean_pmdec / norm * arrow_length
    if arrow_tip is not None:
        dx_scaled,dy_scaled=np.asarray(arrow_tip,dtype=float)

    ax.arrow(
        0,
        0,
        dx_scaled,
        dy_scaled,
        head_width=0.1,
        head_length=0.1,
        fc='black',
        ec='black',
        alpha=0.7,
        lw=2,
        label='Heliocentric proper motion',
    )

    im = ax.imshow(img.T, origin='lower', extent=extent, cmap='bone_r', vmin=22, vmax=34)
    ax.scatter(cold_gas_center_ra, cold_gas_center_dec, s=50, marker='x', color='red', label=gas_center_label)

    levels = _wind_hi_levels(cd_smooth, cd_threshold) if fixed_gas_levels else _gas_contour_levels(cd_smooth, mode="threshold")
    if levels is not None and len(levels):
        ax.contour(x_contour, y_contour, cd_smooth.T, levels=levels, colors='teal')
        ax.plot([], [], color='teal', label=f'$N_{{HI}} > {cd_threshold_str} atoms/cm^2$')

    _style_axis(
        ax,
        r'$\xi$ (deg)',
        r'$\eta$ (deg)',
        xlim=(-SKY_FIELD_HALF_WIDTH_DEG, SKY_FIELD_HALF_WIDTH_DEG),
        ylim=(-SKY_FIELD_HALF_WIDTH_DEG, SKY_FIELD_HALF_WIDTH_DEG),
        invert_x=True,
    )
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='lower right', fontsize=10)

    if shape_outline is not None:
        ax.plot(shape_outline[0],shape_outline[1],color='dimgray',lw=2,zorder=10)
    _append_colorbar(fig, ax, im, r'$\mu_V$ (mag arcsec$^{-2}$)')


def _plot_vlos_panel(fig, ax, vf, extent_vf, arrow_tip=None):
    sc = ax.imshow(vf.T, origin='lower', extent=extent_vf, cmap='coolwarm', vmin=-12, vmax=12)
    _style_axis(
        ax,
        r'$\xi$ (deg)',
        r'$\eta$ (deg)',
        xlim=(extent_vf[0], extent_vf[1]),
        ylim=(extent_vf[2], extent_vf[3]),
        invert_x=True,
    )
    ax.set_aspect('equal', adjustable='box')
    if arrow_tip is not None:
        ax.arrow(0.,0.,*arrow_tip,head_width=.08,head_length=.1,color='black',lw=2,label='Heliocentric proper motion')
        ax.legend(loc='lower right',fontsize=8)
    _append_colorbar(fig, ax, sc, r'$\langle v_{los}\rangle-v_0$ (km/s)')


def _plot_scene_panel(fig, ax, snapshot, star_y, star_z):
    df = snapshot['df']
    cy, cz = snapshot['cy'], snapshot['cz']
    wind_tunnel = snapshot.get('wind_tunnel') is not None
    y_limits = z_limits = (-20., 20.)
    arrow_length = 7.
    star_y_offset, star_z_offset = star_y - cy, star_z - cz
    if wind_tunnel:
        host = snapshot['observer_host_state']
        width = host['box_size_code'] * host['length_to_kpc']
        y_limits = z_limits = (-width / 2, width / 2)
        arrow_length *= width / 40.

        # Centre each periodic image on the stars while retaining the full box.
        gas = df.loc[df['tp'] == 0, ['y', 'z', 'm', 'temp']]
        gas_y, gas_z = _periodic_scene_offsets(gas['y'] - cy, gas['z'] - cz, width)
        star_y_offset, star_z_offset = _periodic_scene_offsets(star_y_offset, star_z_offset, width)
        bins = max(1, round(width / 2.))
        mass, y_edges, z_edges = np.histogram2d(
            gas_y, gas_z, bins=bins, range=(y_limits, z_limits), weights=gas['m'])
        area = np.outer(np.diff(y_edges), np.diff(z_edges))
        surface_density = np.ma.masked_less_equal(mass.T, 0.) / area.T
        background = ax.pcolormesh(y_edges, z_edges, surface_density, shading='flat', cmap='YlOrBr',
                                   norm=LogNorm(*WIND_SURFACE_DENSITY_LIMITS), zorder=0)
        _append_colorbar(fig, ax, background, r'$\Sigma_{\mathrm{gas}}$ ($M_\odot\,\mathrm{kpc}^{-2}$)',
                         extend='both')
        cold = gas['temp'].to_numpy() < snapshot['membership_state']['gas_temperature_split_K']
        cold_y, cold_z = gas_y[cold], gas_z[cold]
        cold_label = 'Cold gas'
        ax.set_title('Full wind tunnel', fontsize=13)
        ax.set_aspect('equal', adjustable='box')
    else:
        cold_gas = df.loc[snapshot['dw_cold_gas_mask']]
        cold_y, cold_z = cold_gas['y'] - cy, cold_gas['z'] - cz
        cold_label = 'HI Gas'

    total_distance = np.hypot(cy, cz)
    if total_distance > 0:
        normalized_cz = (-cy / total_distance) * arrow_length
        normalized_cy = (-cz / total_distance) * arrow_length
    else:
        normalized_cz = 0
        normalized_cy = 0

    ax.arrow(
        0,
        0,
        normalized_cz,
        normalized_cy,
        head_width=1.2,
        head_length=1.2,
        fc='black',
        ec='black',
        zorder=6,
    )
    ax.scatter(cold_y, cold_z, c='aquamarine', s=4, lw=0, label=cold_label)
    if not wind_tunnel:
        hot_gas = df.loc[snapshot['dw_hot_gas_mask']]
        ax.scatter(hot_gas['y'] - cy, hot_gas['z'] - cz, c='orange', s=2, lw=0, label='HII Gas')
    ax.scatter(star_y_offset, star_z_offset, c='gray', s=2, lw=0, label='Stars')
    _style_axis(ax, 'Y (kpc)', 'Z (kpc)', xlim=y_limits, ylim=z_limits)
    ax.legend(loc='lower left', fontsize=14, markerscale=4.)


def _periodic_scene_offsets(y_offset, z_offset, width):
    return tuple(np.remainder(np.asarray(offset) + width / 2, width) - width / 2
                 for offset in (y_offset, z_offset))


def _plot_velocity_dispersion_panel(ax, bin_centers, vlos_dispersion, tsigma, r_half):
    ax.plot(
        bin_centers,
        vlos_dispersion,
        marker='o',
        ms=3,
        lw=1.4,
        color='black',
        label='Member stars, detrended',
        zorder=3,
    )
    ax.errorbar(
        r_pc / 1000,
        sigma,
        yerr=err_sigma,
        fmt='o',
        color='darkorchid',
        ecolor='gray',
        elinewidth=1,
        label='Walker 2009',
        capsize=4,
        alpha=0.6,
        markersize=4,
        zorder=2,
    )
    ax.hlines(tsigma, 0, 3, lw=2, ls='--', color='gray', zorder=0)
    ax.vlines(r_half, 0, 15, lw=2, ls='--', color='gray', zorder=1)
    ax.scatter(r_half, tsigma, c='red', s=80, zorder=10)
    _style_axis(ax, 'Radius (kpc)', 'Velocity Dispersion (km/s)', xlim=(0, 2.2), ylim=(4, 13.5))
    ax.legend(fontsize=12, loc='upper right')


def _surface_brightness_profile_masks(radius_kpc, profile):
    radius_kpc = np.asarray(radius_kpc, dtype=float)
    profile = np.asarray(profile, dtype=float)
    finite_profile = np.isfinite(radius_kpc) & np.isfinite(profile)
    marker_mask = np.zeros(radius_kpc.shape, dtype=bool)
    finite_indices = np.flatnonzero(finite_profile)
    if finite_indices.size == 0:
        return finite_profile, marker_mask

    ordered_indices = finite_indices[np.argsort(radius_kpc[finite_indices])]
    next_marker_radius = -np.inf
    for idx in ordered_indices:
        if radius_kpc[idx] >= next_marker_radius:
            marker_mask[idx] = True
            next_marker_radius = radius_kpc[idx] + SIM_SB_MARKER_SPACING_KPC
    marker_mask[ordered_indices[-1]] = True
    return finite_profile, marker_mask


def _plot_surface_brightness_profile(
    ax,
    dprof,
    mag_obs,
    mag_obs_err,
    r_centers,
    mag_profile,
    adaptive_ml_profile=None,
    adaptive_ml_summary=None,
    *,
    fixed_view_profile=None,
):
    r_centers_plot = np.asarray(r_centers, dtype=float)
    mag_profile_plot = np.asarray(mag_profile, dtype=float)
    finite_profile, marker_mask = _surface_brightness_profile_masks(r_centers_plot, mag_profile_plot)

    ax.errorbar(
        _observed_radius_kpc(dprof),
        mag_obs,
        yerr=mag_obs_err,
        fmt='o',
        capsize=4,
        color='darkorchid',
        ecolor='gray',
        markersize=4,
        elinewidth=1,
        alpha=0.5,
        zorder=1,
        label='Yang 2022',
    )
    ax.plot(
        r_centers_plot[finite_profile],
        mag_profile_plot[finite_profile],
        label=(f'Orbit view, M/L={sim_plot_mass_to_light_ratio:g}' if fixed_view_profile is not None
               else f'Simulation M/L={sim_plot_mass_to_light_ratio:g}'),
        color='black',
        lw=1.2,
        zorder=3,
    )
    ax.plot(
        r_centers_plot[marker_mask],
        mag_profile_plot[marker_mask],
        linestyle='None',
        label='_nolegend_',
        color='black',
        marker='o',
        ms=3,
        zorder=3.1,
    )
    if fixed_view_profile is not None:
        finite, _ = _surface_brightness_profile_masks(r_centers_plot, fixed_view_profile)
        ax.plot(r_centers_plot[finite], fixed_view_profile[finite], color='tab:blue', ls='--',
                lw=1.4, label=f'Fixed view, circular, M/L={sim_plot_mass_to_light_ratio:g}')
    elif adaptive_ml_profile is not None and adaptive_ml_summary is not None:
        adaptive_ml_profile_plot = np.asarray(adaptive_ml_profile, dtype=float)
        adaptive_finite_profile, adaptive_marker_mask = _surface_brightness_profile_masks(
            r_centers_plot,
            adaptive_ml_profile_plot,
        )
        ax.plot(
            r_centers_plot[adaptive_finite_profile],
            adaptive_ml_profile_plot[adaptive_finite_profile],
            color='tab:blue',
            lw=1.4,
            label=f"Adaptive fixed M/L={adaptive_ml_summary['mass_to_light']:.3g}",
            zorder=4,
        )
        ax.plot(
            r_centers_plot[adaptive_marker_mask],
            adaptive_ml_profile_plot[adaptive_marker_mask],
            linestyle='None',
            label='_nolegend_',
            color='tab:blue',
            marker='s',
            ms=3,
            zorder=4.1,
        )
        ax.text(
            0.03,
            0.04,
            f"fixed M/L={sim_plot_mass_to_light_ratio:g}\n"
            f"adaptive M/L={adaptive_ml_summary['mass_to_light']:.3g}",
            transform=ax.transAxes,
            ha='left',
            va='bottom',
            fontsize=9,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.75, 'pad': 3},
        )
    _style_axis(ax, 'Radius (kpc)', 'Surface Brightness (mag arcsec$^{-2}$)', xlim=(0, 5.1), ylim=(36, 22))
    ax.legend(fontsize=12, loc='upper right')


def _plot_gas_fraction_panel(ax, dw_elinfo, tsnap, gas_half_ratio, gas_ratio_stellar_region, t_peri,
                             *, time_label='Simulation Time (Gyr)'):
    gas_half_series = np.asarray(
        dw_elinfo.hi_half_fraction.to_numpy(),
        dtype=float,
    )
    gas_full_series = np.asarray(
        dw_elinfo.hi_fraction.to_numpy(),
        dtype=float,
    )
    ax.plot(
        dw_elinfo.age.to_numpy(),
        gas_half_series,
        lw=3,
        c='black',
        label=r'within $R_{\mathrm{half}}$',
    )
    ax.plot(
        dw_elinfo.age.to_numpy(),
        gas_full_series,
        lw=3,
        c='gray',
        label=stellar_region_description,
    )
    ax.scatter(tsnap, gas_half_ratio, c='red', marker='o', s=100, zorder=5)
    ax.scatter(
        tsnap,
        gas_ratio_stellar_region,
        c='salmon',
        marker='o',
        s=100,
        zorder=5,
    )
    ax.vlines(tsnap, 0, 1.05, color='gray', lw=2, ls='--', zorder=0)
    _style_axis(ax, time_label, r'$f_{\mathrm{HI}} = M_{\mathrm{HI}}/(M_\star+M_{\mathrm{HI}})$', xlim=(0, t_peri), ylim=(1e-4, 1.05))
    ax.set_yscale('log')
    ax.legend(loc='lower left', fontsize=15)


def _plot_sigma_history_panel(ax, dw_elinfo, sigma_mw_theoretical, sigma_mw, tsnap, sigma_value, t_peri, sigma_x, sigma_y, sigma_z,
                              *, time_label='Simulation Time (Gyr)', time_end=None):
    ax.plot(dw_elinfo.age.to_numpy(), dw_elinfo.sigma.to_numpy(), lw=2, color='black', label=r'$\sigma_{los}$')
    ax.plot(dw_elinfo.age.to_numpy(), dw_elinfo.tsigma.to_numpy(), lw=2, color='gray', label=r'$\sigma_{baryon}$')
    if 'tsigma_xyz' in dw_elinfo.columns:
        ax.plot(dw_elinfo.age.to_numpy(), dw_elinfo.tsigma_xyz.to_numpy(), lw=2, color='orange', label=r'$\sigma_{xyz}$')
    ax.plot(dw_elinfo.age.to_numpy(), sigma_mw_theoretical, lw=2, color='lightcoral', label=r'$\sigma_{MW}$')
    ax.plot(dw_elinfo.age.to_numpy(), sigma_mw, lw=2, color='c', label=r'$\sigma_{est}$')
    ax.scatter(tsnap, sigma_value, s=70, zorder=5, color='red')
    _style_axis(ax, time_label, 'velocity dispersion (km/s)', xlim=(0, t_peri * 1.2 if time_end is None else time_end),
                ylim=(-0.08, 20.9))
    ax.vlines(tsnap, 0, 25, color='gray', lw=2, ls='--', zorder=0)
    ax.legend(loc='upper left', fontsize=15)
    ax.text(
        0.98,
        0.95,
        (r"$\sigma_x = %.2f$ km/s" % sigma_x) + "\n"
        + (r"$\sigma_y = %.2f$ km/s" % sigma_y) + "\n"
        + (r"$\sigma_z = %.2f$ km/s" % sigma_z),
        transform=ax.transAxes,
        ha='right',
        va='top',
        fontsize=14,
    )


def _fixed_view_light_profile(star_position, mass, center, los, distance, **profile_kwargs):
    offsets = star_position - center
    radius_squared = np.einsum('ij,ij->i', offsets, offsets) - (offsets @ los) ** 2
    radius = np.sqrt(np.maximum(radius_squared, 0.))
    _, profile = Analysis.radial_magnitude_profile(
        radius, np.zeros_like(radius), mass, distance, ep=0.,
        m_l_relation=sim_plot_mass_to_light_ratio, **profile_kwargs,
    )
    return profile


def _compute_plot_snapshot_metrics(snapshot, context, numsp):
    df = snapshot['df']
    dw_cold_gas_mask = snapshot['dw_cold_gas_mask']
    total_dw_star_mask = snapshot['total_dw_star_mask']
    old_kinematics = old_star_projected_kinematics(snapshot)
    member_kinematics = member_star_projected_kinematics(snapshot)
    old_star_local_mask = old_kinematics['old_local_mask']
    xi_all, eta_all = tangent_coordinates(snapshot['star_coords']['_sky'])
    sky_xi_star = xi_all[old_star_local_mask]
    sky_eta_star = eta_all[old_star_local_mask]
    photometry_front = np.asarray(snapshot['star_coords']['_sky']['front_mask'])[old_star_local_mask]
    rotra_dw_cold_gas = snapshot['rotra_dw_cold_gas']
    rotdec_dw_cold_gas = snapshot['rotdec_dw_cold_gas']
    d_mean = snapshot['d_mean']
    x_kpc = old_kinematics['x_kpc']
    y_kpc = old_kinematics['y_kpc']
    vlos_array = member_kinematics['vlos']
    vlos_detrended = member_kinematics['vlos_detrended']
    cold_gas_x_kpc = snapshot['cold_gas_x_kpc']
    cold_gas_y_kpc = snapshot['cold_gas_y_kpc']
    dw_elinfo = context['dw_elinfo']
    row = context['elinfo_by_numsp'].loc[numsp]
    star_mass = row['star_mass']
    star_half_mass = row['star_half_mass']
    coldgas_mass = row['hi_mass']
    coldgas_half_mass = row['hi_half_mass']
    vr = row['vr']
    vtheta = row['vtheta']
    vphi = row['vphi']

    star_df = df.loc[
        total_dw_star_mask,
        ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz'],
    ].iloc[old_star_local_mask]
    # The old-star cut applies to photometry; kinematics include all current members.
    star_mass_array = star_df['m'].to_numpy()
    x_star = star_df['x'].to_numpy()
    y_star = star_df['y'].to_numpy()
    z_star = star_df['z'].to_numpy()
    vx_star = star_df['vx'].to_numpy()
    vy_star = star_df['vy'].to_numpy()
    vz_star = star_df['vz'].to_numpy()
    cold_gas_df = df.loc[dw_cold_gas_mask, ['nh', 'm']]
    hot_gas_df = df.loc[snapshot['dw_hot_gas_mask'], ['nh', 'm']]

    gas_ratio_stellar_region = hi_fraction(coldgas_mass, star_mass)
    gas_half_ratio = hi_fraction(coldgas_half_mass, star_half_mass)
    gas_frac_stellar_region = hi_fraction(coldgas_mass, star_mass)
    gas_half_frac = hi_fraction(coldgas_half_mass, star_half_mass)

    eps = row['eps']
    pa = row['pa']
    r_half = row['rhalf']
    r_half_circularized = row['rhalf_circularized']
    shape_center_x_kpc = row['shape_center_x_kpc']
    shape_center_y_kpc = row['shape_center_y_kpc']

    r_half_deg = np.degrees(r_half / d_mean)
    shape_center_ra,shape_center_dec=tangent_kpc_to_display(snapshot['star_coords']['_sky'],shape_center_x_kpc,shape_center_y_kpc)
    phi=np.linspace(0,2*np.pi,181)
    major=stellar_region_rhalf_multiplier*r_half
    minor=major*(1-eps)
    xp=shape_center_x_kpc+major*np.cos(phi)*np.cos(pa)-minor*np.sin(phi)*np.sin(pa)
    yp=shape_center_y_kpc+major*np.cos(phi)*np.sin(pa)+minor*np.sin(phi)*np.cos(pa)
    shape_outline=tangent_kpc_to_display(snapshot['star_coords']['_sky'],xp,yp)

    if 'cold_gas_center_ra' in row.index and 'cold_gas_center_dec' in row.index:
        cold_gas_center_ra = row['cold_gas_center_ra']
        cold_gas_center_dec = row['cold_gas_center_dec']
    else:
        cold_gas_r = np.asarray(snapshot['cold_gas_coords']['_sky']['rho_deg'])
        cold_gas_within_mask = cold_gas_r <= 2.1
        rotra_dw_cold_gas_within = np.array(rotra_dw_cold_gas)[cold_gas_within_mask]
        rotdec_dw_cold_gas_within = np.array(rotdec_dw_cold_gas)[cold_gas_within_mask]
        if len(rotra_dw_cold_gas_within) > 0:
            cold_gas_center_ra, cold_gas_center_dec = Analysis.find_center_2d(
                rotra_dw_cold_gas_within,
                rotdec_dw_cold_gas_within,
                units='offset',
            )
        else:
            cold_gas_center_ra, cold_gas_center_dec = np.nan, np.nan

    # CSV centre fields retain their historical raw offsets; transform this same
    # sky point for display, without changing the centre estimator or the CSV.
    cold_gas_center_ra, cold_gas_center_dec = raw_offsets_to_tangent(
        snapshot['star_coords']['_sky'], cold_gas_center_ra, cold_gas_center_dec)

    # Current snapshot canonical six-vectors, never CSV/legacy PM precedence.
    mean_pmra,mean_pmdec=snapshot['motion_R5']['HRF_EN_masyr']
    arrow_tip=tangent_motion_arrow(snapshot['motion_R5']['HRF_unit_EN'],length_deg=.65)['tip_xi_eta_deg']

    obs_r_intervals = _observed_profile_intervals_kpc(context['dprof'])
    profile_kwargs = {'r_intervals': obs_r_intervals} if obs_r_intervals is not None else {'bins': 40}
    r_centers, mag_profile = Analysis.radial_magnitude_profile(
        x_kpc,
        y_kpc,
        mass=star_mass_array,
        d_kpc=d_mean,
        ep=eps,
        m_l_relation=sim_plot_mass_to_light_ratio,
        pa=pa,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
        **profile_kwargs,
    )
    wind_tunnel = context.get('wind_tunnel')
    fixed_view_profile = adaptive_ml_summary = None
    if wind_tunnel:
        native = snapshot['star_coords']['_native']
        sun = native['observer']['sun_position_kpc']
        los = wind_tunnel.orbit[0, 1:4] - sun
        los /= np.linalg.norm(los)
        pointing = snapshot['star_coords']['_sky']['center_h_native_kpc'] + sun
        fixed_view_profile = _fixed_view_light_profile(
            np.column_stack((x_star, y_star, z_star)), star_mass_array, pointing,
            los, d_mean, **profile_kwargs,
        )
    else:
        adaptive_ml_summary = _fit_adaptive_scalar_mass_to_light(
            r_centers, mag_profile, context['dprof'], context['mag_obs'], context['mag_obs_err'], r_half,
        )
    adaptive_ml_profile = (
        None if adaptive_ml_summary is None else adaptive_ml_summary['mag_profile']
    )

    cy = snapshot['cy']
    cz = snapshot['cz']
    cx = snapshot['cx']
    dw_xc = snapshot['dw_xc']
    dw_yc = snapshot['dw_yc']
    dw_zc = snapshot['dw_zc']
    
    if 'sigma_x' in row.index and 'sigma_y' in row.index and 'sigma_z' in row.index:
        sigma_x = row['sigma_x']
        sigma_y = row['sigma_y']
        sigma_z = row['sigma_z']
    else:
        members = df.loc[total_dw_star_mask]
        r = np.sqrt((members['x'].to_numpy() - dw_xc) ** 2
                    + (members['y'].to_numpy() - dw_yc) ** 2
                    + (members['z'].to_numpy() - dw_zc) ** 2)
        within_500pc = r <= 0.5
        sigma_x = np.nan if not np.any(within_500pc) else np.std(members['vx'].to_numpy()[within_500pc])
        sigma_y = np.nan if not np.any(within_500pc) else np.std(members['vy'].to_numpy()[within_500pc])
        sigma_z = np.nan if not np.any(within_500pc) else np.std(members['vz'].to_numpy()[within_500pc])

    bin_centers, vlos_dispersion = Analysis.get_vlos_dispersion(
        member_kinematics['x_kpc'], member_kinematics['y_kpc'], vlos_detrended)
    img, extent = Analysis.surface_brightness_map(
        sky_xi_star[photometry_front],
        sky_eta_star[photometry_front],
        star_mass_array[photometry_front],
        d_mean,
        npix=120,
        coordinate_size=2.0 * SKY_FIELD_HALF_WIDTH_DEG,
        mag_sys='vega',
        m_l_relation=sim_plot_mass_to_light_ratio,
        coordinate_unit='deg',
    )
    # Replace the precise planar area used inside surface_brightness_map by
    # each TAN pixel's integrated solid angle, preserving M/L and zero point.
    actual_area=tangent_pixel_area_arcsec2(120,half_width_deg=SKY_FIELD_HALF_WIDTH_DEG)
    old_area=(np.radians(2.*SKY_FIELD_HALF_WIDTH_DEG/120)*206265.)**2
    img=img+2.5*np.log10(actual_area/old_area)
    sky_field_half_width_kpc=np.radians(SKY_FIELD_HALF_WIDTH_DEG)*d_mean
    display_grid=tangent_velocity_grid40(snapshot['star_coords']['_sky'],vlos_array,
        display_v0_kms=snapshot['walker_R5']['v0_kms'])
    vf=display_grid['mean_HRF_minus_v0']
    extent_vf=(-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG,
               -SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG)


    vgsr = np.sqrt(vr**2 + vtheta**2 + vphi**2)

    gas_hi_mass = hi_particle_masses(
        np.concatenate((cold_gas_df['m'].to_numpy(), hot_gas_df['m'].to_numpy())),
        np.concatenate((cold_gas_df['nh'].to_numpy(), hot_gas_df['nh'].to_numpy())))
    contour_range_kpc = [
        [-sky_field_half_width_kpc, sky_field_half_width_kpc],
        [-sky_field_half_width_kpc, sky_field_half_width_kpc],
    ]
    gas_front = np.concatenate((snapshot['cold_gas_coords']['_sky']['front_mask'],
                                snapshot['hot_gas_coords']['_sky']['front_mask']))
    gas_east = np.concatenate((cold_gas_x_kpc, snapshot['hot_gas_x_kpc']))[gas_front]
    gas_north = np.concatenate((cold_gas_y_kpc, snapshot['hot_gas_y_kpc']))[gas_front]
    gas_mass = gas_hi_mass[gas_front]
    cold_gas_mass_contour, x_contour_kpc, y_contour_kpc, cd_smooth = Analysis.est_gas_contour(
        gas_east, gas_north, gas_mass,
        box=None,
        nbins=GAS_CONTOUR_NBINS,
        smooth_sigma=GAS_CONTOUR_SMOOTH_SIGMA_PIXELS,
        histogram_range=contour_range_kpc,
        smooth_mode="constant",
        unit="kpc",
        d_mean=d_mean,
    )
    star_front = np.asarray(snapshot['star_coords']['_sky']['front_mask'])
    star_hist, _, _ = np.histogram2d(np.asarray(snapshot['x_kpc'])[star_front], np.asarray(snapshot['y_kpc'])[star_front], bins=GAS_CONTOUR_NBINS,
        range=contour_range_kpc, weights=df.loc[total_dw_star_mask, 'm'].to_numpy()[star_front])
    star_mass_contour = star_hist[cd_smooth > cd_threshold].sum(dtype=np.float64)
    hi_map = (x_contour_kpc, y_contour_kpc, cd_smooth, sky_field_half_width_kpc)
    overlay_east, overlay_north, overlay_column = x_contour_kpc, y_contour_kpc, cd_smooth
    if wind_tunnel:
        half_width = WIND_HI_GRID_HALF_WIDTH_KPC
        limits = (-half_width, half_width)
        bins = round(GAS_CONTOUR_NBINS * half_width / WIND_HI_HALF_WIDTH_KPC)
        _, overlay_east, overlay_north, overlay_column = Analysis.est_gas_contour(
            gas_east, gas_north, gas_mass, box=None, nbins=bins,
            smooth_sigma=GAS_CONTOUR_SMOOTH_SIGMA_PIXELS, histogram_range=(limits, limits),
            smooth_mode='constant',
        )
        hi_map = (overlay_east, overlay_north, overlay_column, WIND_HI_HALF_WIDTH_KPC)
    x_contour, y_contour = tangent_kpc_to_display(snapshot['star_coords']['_sky'], overlay_east, overlay_north)

    return {
        'display_grid': display_grid,
        'pixel_area_arcsec2': actual_area,
        'row': row,
        'star_mass': star_mass,
        'star_half_mass': star_half_mass,
        'coldgas_mass': coldgas_mass,
        'coldgas_half_mass': coldgas_half_mass,
        'gas_ratio_stellar_region': gas_ratio_stellar_region,
        'gas_half_ratio': gas_half_ratio,
        'gas_frac_stellar_region': gas_frac_stellar_region,
        'gas_half_frac': gas_half_frac,
        'star_y': y_star,
        'star_z': z_star,
        'cold_gas_center_ra': cold_gas_center_ra,
        'cold_gas_center_dec': cold_gas_center_dec,
        'arrow_tip': arrow_tip,
        'shape_outline': shape_outline,
        'mean_pmra': mean_pmra,
        'mean_pmdec': mean_pmdec,
        'eps': eps,
        'pa': pa,
        'r_half': r_half,
        'r_half_circularized': r_half_circularized,
        'r_half_deg': r_half_deg,
        'shape_center_x_kpc': shape_center_x_kpc,
        'shape_center_y_kpc': shape_center_y_kpc,
        'shape_center_ra': shape_center_ra,
        'shape_center_dec': shape_center_dec,
        'r_centers': r_centers,
        'mag_profile': mag_profile,
        'adaptive_ml_profile': adaptive_ml_profile,
        'adaptive_ml_summary': adaptive_ml_summary,
        'fixed_view_profile': fixed_view_profile,
        'vlos_array': vlos_array,
        'vlos_detrended': vlos_detrended,
        'velocity_gradient': member_kinematics['velocity_gradient'],
        'cy': cy,
        'cz': cz,
        'sigma_x': sigma_x,
        'sigma_y': sigma_y,
        'sigma_z': sigma_z,
        'bin_centers': bin_centers,
        'vlos_dispersion': vlos_dispersion,
        'img': img,
        'extent': extent,
        'vf': vf,
        'extent_vf': extent_vf,
        'vr': vr,
        'vtheta': vtheta,
        'vphi': vphi,
        'vgsr': vgsr,
        'cold_gas_mass_contour': cold_gas_mass_contour,
        'hi_fraction_contour': hi_fraction(cold_gas_mass_contour, star_mass_contour),
        'x_contour_kpc': x_contour_kpc,
        'y_contour_kpc': y_contour_kpc,
        'x_contour': x_contour,
        'y_contour': y_contour,
        'cd_smooth': cd_smooth,
        'overlay_column': overlay_column,
        'sky_field_half_width_kpc': sky_field_half_width_kpc,
        'hi_map': hi_map,
    }


def _binned_simulated_sfr(formation_gyr, mass_msun, edges_gyr):
    mass_by_bin, _ = np.histogram(formation_gyr, bins=edges_gyr, weights=mass_msun)
    return mass_by_bin / (np.diff(edges_gyr) * 1.e9) / SFH_SFR_UNIT_MSUN_PER_YEAR


def _compute_sfr_series(sfh_dw, dw_sfrmass, t_peri, bin_width=0.5):
    if dw_sfrmass is None or len(dw_sfrmass) == 0:
        bins = np.arange(0 - bin_width, t_peri * 1.2, bin_width)
        lookback_obs, sfr_obs = DataProcessor.GetSFR(sfh_dw, bins=bins)
        return np.array([]), np.array([]), lookback_obs, sfr_obs

    filtered_sfrmass = dw_sfrmass[dw_sfrmass['m'] < 10000]
    if len(filtered_sfrmass) == 0:
        bins = np.arange(0 - bin_width, t_peri * 1.2, bin_width)
        lookback_obs, sfr_obs = DataProcessor.GetSFR(sfh_dw, bins=bins)
        return np.array([]), np.array([]), lookback_obs, sfr_obs

    lookback = filtered_sfrmass['age'].to_numpy()
    sfh_mass = filtered_sfrmass['m'].to_numpy()
    bins = np.arange(0 - bin_width, t_peri * 1.2, bin_width)

    sfr_sim = _binned_simulated_sfr(lookback, sfh_mass, bins)
    lookback_sim = bins[:-1]
    lookback_obs, sfr_obs = DataProcessor.GetSFR(sfh_dw, bins=bins)

    return lookback_sim, sfr_sim, lookback_obs, sfr_obs


def _wind_sfr_series(snapshot, elapsed_gyr, bin_width_gyr=.05):
    """Bin the current member stars formed in the simulation by formation time."""
    df = snapshot['df']
    stars = df.loc[snapshot['total_dw_star_mask'] & (df['tp'] == 4), ['birth', 'stellar_age_gyr', 'm']]
    frame = snapshot['simulation'].wind_tunnel
    time_to_gyr = frame.time_to_myr / 1000.
    birth = (stars['birth'].to_numpy(dtype=float) - frame.time_begin_code) * time_to_gyr
    age = stars['stellar_age_gyr'].to_numpy(dtype=float)
    birth = np.where(np.isfinite(birth), birth, elapsed_gyr - age)
    mass = stars['m'].to_numpy(dtype=float)
    # The stored stellar times can be rounded to float32 at the snapshot boundary.
    time_slop = abs(np.spacing(np.float32(snapshot['tsnap']))) * time_to_gyr
    valid = np.isfinite(birth) & (birth <= elapsed_gyr + time_slop) & np.isfinite(mass) & (mass > 0.)
    if not valid.any():
        return np.array([]), np.array([])
    birth = np.minimum(birth[valid], elapsed_gyr)
    first_bin = int(np.floor(min(0., birth.min()) / bin_width_gyr))
    last_bin = int(np.ceil(elapsed_gyr / bin_width_gyr))
    edges = np.arange(first_bin, last_bin) * bin_width_gyr
    edges = np.append(edges[edges < elapsed_gyr - time_slop], elapsed_gyr)
    return edges, _binned_simulated_sfr(birth, mass[valid], edges)


def _wind_reference_sfh(source_directory, snapshot_number, core_radius):
    """Use one reference output for the SFH shown throughout the sequence."""
    simulation = GalaxySimulation(folder_path=source_directory, snapshot_num=snapshot_number)
    _, time = simulation.load_snapshot(
        dftype=3, box=False, gas_fields=(), star_fields=('birth',), particle_types=(2, 3, 4),
    )
    stars = simulation.center_on_MW()
    _, members = simulation.classify_mw_dwarf(df=stars, r_exclude=5, dwarf_radius=3 * core_radius,
                                              k_density=16)
    elapsed = simulation.wind_tunnel.at_time(time)['elapsed_since_insertion_Myr'] / 1000.
    return _wind_sfr_series(dict(simulation=simulation, tsnap=time, df=stars,
                                total_dw_star_mask=members), elapsed)


def _history_upper_limit(values, step):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return max(step, np.ceil(1.05 * finite.max(initial=0.) / step) * step)


def _wind_history_limits(frame, history):
    limits = {}
    for name, series, units, step in (
        ('gas_mass', WIND_GAS_HISTORY, 1.e6, 20.),
        ('ram_pressure', WIND_PRESSURE_HISTORY, 1.e-15, 10.),
    ):
        data = history[[column for column, _ in series]].to_numpy(dtype=float)
        limits[name] = (0., _history_upper_limit(data / units, step))
    distances = np.r_[np.linalg.norm(frame.orbit[:, 1:4], axis=1),
                      history[['distance', 'distance_gal']].to_numpy().ravel(), d_today]
    finite = distances[np.isfinite(distances)]
    limits['distance'] = (max(0., 20. * np.floor((finite.min() - 5.) / 20.)),
                          _history_upper_limit(finite, 20.))
    return limits


def _observed_sfh_on_wind_clock(sfh, reference_elapsed_gyr):
    """Transform observed lookback bins to the wind tunnel elapsed-time axis."""
    observed = sfh.loc[np.isfinite(sfh['age']) & np.isfinite(sfh['sfr']) & (sfh['age'] >= 0.)]
    if len(observed) < 2:
        raise ValueError('The observed SFH needs at least two finite lookback and SFR measurements')
    n_bins = max(1, int(np.ceil(observed['age'].max() / .5)))
    lookback_edges = .5 * np.arange(n_bins + 1)
    _, sfr = DataProcessor.GetSFR(observed, bins=lookback_edges)
    counts, _ = np.histogram(observed['age'], bins=lookback_edges)
    sfr[counts == 0] = np.nan
    return reference_elapsed_gyr - lookback_edges[::-1], sfr[::-1]


def _plot_wind_sfh_panel(ax, series, elapsed_gyr, clock_range, observed=None):
    edges, sfr = series
    first, last = clock_range
    values = list(sfr)
    if len(sfr):
        first = min(first, edges[0])
        ax.stairs(sfr, edges, baseline=None, color='black', lw=2, label='Simulated SFH')
    if observed is not None:
        obs_edges, obs_sfr = observed
        overlaps = (obs_edges[:-1] < last) & (obs_edges[1:] > first)
        if np.any(overlaps & np.isfinite(obs_sfr)):
            ax.stairs(np.where(overlaps, obs_sfr, np.nan), obs_edges,
                      baseline=None, color='darkorange', lw=2, label='Observed SFH')
            values.extend(obs_sfr[overlaps])
    ax.axvline(0, color='gray', lw=1.4, ls='--', label='Insertion')
    ax.axvline(elapsed_gyr, color='firebrick', lw=1.4, label='Current time')
    _style_axis(ax, 'Time since insertion (Gyr)', r'SFR ($10^{-4}\,M_\odot\,\mathrm{yr}^{-1}$)',
                xlim=(first, last), ylim=(0., max(25., _history_upper_limit(values, 5.))))
    ax.legend(loc='upper right', fontsize=11)


def _plot_gas_loss_panel(ax, dw_elinfo, gas_loss_rate, tsnap, t_peri):
    ax.plot(dw_elinfo.age.to_numpy(), gas_loss_rate, lw=2.5, c='black')
    ax.vlines(tsnap, -1e10, 1e10, color='gray', lw=2, ls='--')
    _style_axis(ax, 'Simulation Time (Gyr)', 'Gas Loss Rate (Msun / Gyr)', xlim=(0, t_peri * 1.2), ylim=(0.5e7, -2.5e7))


def _plot_sfr_panel(ax, lookback_sim, sfr_sim, lookback_obs, sfr_obs, t_peri):
    ax.step(lookback_sim, sfr_sim, where='post', label='Simulated SFR', color='black', lw=3, alpha=0.7)
    ax.step(lookback_obs, sfr_obs, where='post', label='Observed SFR', color='orange', lw=3, alpha=0.7)
    _style_axis(ax, 'Lookback Time (Gyr)', r'SFR ($10^{-4} M_{\odot}$ yr$^{-1}$)', xlim=(0, t_peri), ylim=(0, 25))
    ax.legend(loc='upper left')


def _plot_wind_history(ax, history, tsnap, history_limit, series, ylabel, scale, ylim):
    for column, label in series:
        available = history[column].notna()
        if available.any():
            ax.scatter(history.loc[available, 'age'], history.loc[available, column] / scale, label=label)
    ax.axvline(tsnap, color='gray', lw=1.5, ls='--', zorder=0)
    _style_axis(ax, 'Time since insertion (Gyr)', ylabel, xlim=(0, history_limit), ylim=ylim)
    if ax.collections:
        ax.legend(loc='upper right', fontsize=10)


def _plot_orbit_panel(ax, dw_elinfo, tsnap, current_distance, t_peri, vr, vtheta, vphi, vgsr,
                      *, wind_tunnel=None, wind_distance_limits=None):
    if wind_tunnel:
        elapsed = wind_tunnel.orbit[:, 0] / 1000.
        box_radius = np.linalg.norm(wind_tunnel.orbit[:, 1:4], axis=1)
        ax.plot(elapsed, box_radius, color='gray', ls='--', lw=2, label='Box Galactocentric')
        ax.plot(dw_elinfo.age, dw_elinfo.distance_gal, color='tab:blue', lw=1.4, label='Stellar Galactocentric')
        ax.plot(dw_elinfo.age, dw_elinfo.distance, color='black', lw=2, label='Stellar heliocentric')
        ax.axhline(d_today, color='gray', lw=1., ls=':', zorder=0)
        xlabel, ylabel, ylim = 'Time since insertion (Gyr)', 'Distance (kpc)', wind_distance_limits
    else:
        ax.plot(dw_elinfo.age.to_numpy(), dw_elinfo.distance.to_numpy(), color='black', lw=3, zorder=1)
        ax.hlines(d_today, -0.1, t_peri * 3, color='gray', lw=2, ls='--', zorder=0)
        xlabel, ylabel, ylim = 'Simulation Time (Gyr)', 'Distance (kpc)', (30, 820)
        distance_label = '$D$'
    ax.scatter(
        tsnap,
        current_distance,
        color='red',
        s=80,
        zorder=3,
        label=f'{tsnap:.2f} Gyr',
    )
    _style_axis(ax, xlabel, ylabel, xlim=(0, t_peri if wind_tunnel else t_peri * 1.2), ylim=ylim)
    if wind_tunnel:
        ax.legend(loc='upper right', fontsize=10)
        return
    ax.legend(loc='lower left', fontsize=16)
    ax.text(
        0.05,
        0.95,
        f'{distance_label} : {current_distance:.2f} kpc\n'
        f'$v_{{r}}$ : {vr:.2f} km/s\n'
        f'$v_{{\\theta}}$ : {vtheta:.2f} km/s\n'
        f'$v_{{\\phi}}$ : {vphi:.2f} km/s\n'
        f'$v_{{gsr}}$ : {vgsr:.2f} km/s',
        transform=ax.transAxes,
        fontsize=20,
        verticalalignment='top',
        horizontalalignment='left',
    )


def plot(numsp, show_mem=False, context=None, profile=False, profile_memory=False, verbose=False):
    start_time = time.perf_counter()
    start_rss = get_rss_mb() if profile_memory else None
    context = get_plot_context() if context is None else context
    modelname = context['modelname']
    core_radius = context['core_radius']
    dw_elinfo = context['dw_elinfo']
    elinfo_by_numsp = context['elinfo_by_numsp']
    ini_x = context['ini_x']
    ini_y = context['ini_y']
    ini_z = context['ini_z']
    ini_vx = context['ini_vx']
    ini_vy = context['ini_vy']
    ini_vz = context['ini_vz']
    fd_factor_str = context['fd_factor_str']
    dprof = context['dprof']
    mag_obs = context['mag_obs']
    mag_obs_err = context['mag_obs_err']
    sigma_mw = context['sigma_mw']
    sigma_mw_theoretical = context['sigma_mw_theoretical']
    wind_tunnel = context.get('wind_tunnel')
    history_limit = context['history_limit']
    time_label = 'Time since insertion (Gyr)' if wind_tunnel else 'Simulation Time (Gyr)'
    output_folder = context['output_folder']

    if show_mem:
        try:
            import psutil
            process = psutil.Process(os.getpid())
            mem_mb = process.memory_info().rss / 1024 / 1024
            print(f"[numsp={numsp}] Memory usage at start: {mem_mb:.2f} MB")
        except Exception as e:
            print(f"[numsp={numsp}] Memory usage check failed: {e}")

    snapshot = prepare_snapshot_context(
        folder_path=folder_path,
        snapshot_num=numsp,
        core_radius=core_radius,
        include_star_birth=True,
    )
    after_prepare = time.perf_counter()
    after_prepare_rss = get_rss_mb() if profile_memory else None
    simulation = snapshot['simulation']
    df = snapshot['df']
    tsnap = (snapshot['observer_host_state']['elapsed_since_insertion_Myr'] / 1000.
             if wind_tunnel else snapshot['tsnap'])
    metrics = _compute_plot_snapshot_metrics(snapshot, context, numsp)
    after_metrics = time.perf_counter()
    after_metrics_rss = get_rss_mb() if profile_memory else None
    row = elinfo_by_numsp.loc[numsp]

    star_mass = metrics['star_mass']
    star_half_mass = metrics['star_half_mass']
    coldgas_mass = metrics['coldgas_mass']
    coldgas_half_mass = metrics['coldgas_half_mass']
    gas_ratio_stellar_region = metrics['gas_ratio_stellar_region']
    gas_half_ratio = metrics['gas_half_ratio']
    gas_frac_stellar_region = metrics['gas_frac_stellar_region']
    gas_half_frac = metrics['gas_half_frac']

    sm_str = DataProcessor.sci_notation(star_mass)
    shm_str = DataProcessor.sci_notation(star_half_mass)
    cg_str = DataProcessor.sci_notation(coldgas_mass)
    cghm_str = DataProcessor.sci_notation(coldgas_half_mass)
    cd_threshold_str = DataProcessor.sci_notation(cd_threshold)

    ### --- ###

    fig = plt.figure(figsize=(21, 15), dpi=150)
    gs = gridspec.GridSpec(3, 4)

    # The formal-orbit frame uses its final cell for figure annotations.
    axes = [[None]*4 for _ in range(3)]
    for i in range(3):
        for j in range(4):
            if wind_tunnel or (i, j) != (2, 3):
                axes[i][j] = fig.add_subplot(gs[i, j])

    axes = np.array(axes) 

    plt.subplots_adjust(
        left=0.05,   
        right=0.95, 
        top=0.95,   
        bottom=0.08,
        wspace=0.48 if wind_tunnel else 0.35,
        hspace=0.22
    )

    #--- Title ---#

    title = f'{modelname} wind tunnel, snapshot {numsp:03d}' if wind_tunnel else modelname
    fig.suptitle(title, fontsize=22)

    if verbose:
        print('Begin plotting...')

    #--- H I sky map (kpc) ---#

    _plot_hi_sky_panel(
        fig,
        axes[1, 1],
        *metrics['hi_map'],
        wind_tunnel=wind_tunnel,
    )

    if verbose:
        print('Gas sky map plotted.')

    #--- Stellar sky map (deg) ---#
    _plot_stellar_sky_panel(
        fig,
        axes[0, 0],
        metrics['img'],
        metrics['extent'],
        metrics['cold_gas_center_ra'],
        metrics['cold_gas_center_dec'],
        metrics['x_contour'],
        metrics['y_contour'],
        metrics['overlay_column'],
        cd_threshold_str,
        metrics['mean_pmra'],
        metrics['mean_pmdec'],
        metrics['r_half_deg'],
        metrics['eps'],
        metrics['pa'],
        metrics['shape_center_ra'],
        metrics['shape_center_dec'],
        arrow_tip=metrics['arrow_tip'],
        shape_outline=metrics['shape_outline'],
        gas_center_label='Cold-gas centre' if wind_tunnel else 'HI Center',
        fixed_gas_levels=bool(wind_tunnel),
    )

    if verbose:
        print('Stellar sky map plotted.')

    #--- vlos sky map (deg) ---#
    _plot_vlos_panel(fig, axes[0, 1], metrics['vf'], metrics['extent_vf'],arrow_tip=metrics['arrow_tip'])

    if verbose:
        print('vlos sky map plotted.')

    #--- 3D Scene in Y-Z ---#
    _plot_scene_panel(
        fig,
        axes[0, 3],
        snapshot,
        metrics['star_y'],
        metrics['star_z'],
    )

    if verbose:
        print('3D scene plotted.')

    #--- velocity dispersion ---#
    _plot_velocity_dispersion_panel(
        axes[1, 0],
        metrics['bin_centers'],
        metrics['vlos_dispersion'],
        row['tsigma'],
        metrics['r_half_circularized'],
    )

    if verbose:
        print('Velocity dispersion plotted.')

    #--- Surface luminosity profile ---#
    _plot_surface_brightness_profile(
        axes[2, 0],
        dprof,
        mag_obs,
        mag_obs_err,
        metrics['r_centers'],
        metrics['mag_profile'],
        metrics['adaptive_ml_profile'],
        metrics['adaptive_ml_summary'],
        fixed_view_profile=metrics['fixed_view_profile'],
    )

    if verbose:
        if metrics['adaptive_ml_summary'] is not None:
            print(
                "Surface luminosity profile plotted. "
                f"Adaptive M/L={metrics['adaptive_ml_summary']['mass_to_light']:.6g}, "
                f"RMS={metrics['adaptive_ml_summary']['rms_mag']:.4f} mag."
            )
        elif not wind_tunnel:
            print('Surface luminosity profile plotted. Adaptive M/L unavailable.')
        else:
            print('Surface luminosity profiles plotted at a fixed M/L for both views.')

    #--- Gas-to-stellar mass ratio ---#
    _plot_gas_fraction_panel(
        axes[1, 2],
        dw_elinfo,
        tsnap,
        gas_half_frac,
        gas_frac_stellar_region,
        history_limit,
        time_label=time_label,
    )

    if verbose:
        print('Gas-to-stellar mass ratio plotted.')

    #--- Sigma history ---#
    _plot_sigma_history_panel(
        axes[0, 2],
        dw_elinfo,
        sigma_mw_theoretical,
        sigma_mw,
        tsnap,
        row['sigma'],
        history_limit,
        metrics['sigma_x'],
        metrics['sigma_y'],
        metrics['sigma_z'],
        time_label=time_label,
        time_end=history_limit if wind_tunnel else None,
    )
    
    if verbose:
        print('Sigma history plotted.')

    #------#

    #--- Central Gas Density ---#

    if wind_tunnel:
        wind_limits = context['wind_history_limits']
        _plot_wind_history(axes[2, 2], dw_elinfo, tsnap, history_limit, WIND_GAS_HISTORY,
                           r'Gas mass ($10^6\,M_\odot$)', 1.e6, wind_limits['gas_mass'])
    else:
        _plot_gas_loss_panel(axes[2, 2], dw_elinfo, context['gas_loss_rate'], tsnap, history_limit)

    #--- Evolution Histroy ---#

    # axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.coldgas_mass.to_numpy(), lw=2.5, c='aquamarine', label=f'HI mass {stellar_region_description}')
    # axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.coldgas_half_mass.to_numpy(), lw=2.5, c='teal', label='HI mass within $R_{{half}}$')
    # axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.star_mass.to_numpy(), lw=2, c='lightblue', alpha=0.9, label='Stellar mass')
    # axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.star_half_mass.to_numpy(), lw=2, c='darkblue', label='Stellar mass within $R_{{half}}$')

    # axes[2, 2].step(sfh_dw['cosmic_time_shifted'].to_numpy(), 1e9*sfh_dw['cumulative_mass'].to_numpy(), 
    #                 lw=2.5, color='maroon', ls='--', label='SFH in Observation')
    # axes[2, 2].step(sfh_dw['cosmic_time_shifted'].to_numpy(), 1e9*sfh_dw['pure_sfh'].to_numpy(), 
    #                 lw=2.5, color='coral', ls='--', label='Pure SFH in Observation')
    # if dw_sfrmass is not None:
    #     axes[2, 2].step(dw_sfrmass['birth'].to_numpy(), dw_sfrmass['cumulative_mass'].to_numpy(), lw=2, color='indigo', label='SFH in Simulation')

    # axes[2, 2].vlines(tsnap, 0, 1e11, color='gray', lw=2, ls='--')
    # axes[2, 2].hlines(dw_cold_gas_lower_limit, -0.1, T_peri*2, lw=2, ls='--', color='paleturquoise', label='lower limit of HI mass')

    if wind_tunnel:
        time_begin = -max(0., wind_tunnel.time_begin_code) * wind_tunnel.time_to_myr / 1000.
        _plot_wind_sfh_panel(axes[2, 1], context['wind_sfh'], tsnap, (time_begin, history_limit),
                             context.get('wind_observed_sfh'))
        _plot_wind_history(axes[2, 3], dw_elinfo, tsnap, history_limit, WIND_PRESSURE_HISTORY,
                           r'Ram pressure ($10^{-15}$ dyn cm$^{-2}$)', 1.e-15, wind_limits['ram_pressure'])
    else:
        _plot_sfr_panel(axes[2, 1], context['lookback_sim'], context['sfr_sim'],
                        context['lookback_obs'], context['sfr_obs'], context['sfr_t_peri_gyr'])

    #--- oribt ---#

    _plot_orbit_panel(
        axes[1, 3],
        dw_elinfo,
        tsnap,
        row['distance'],
        history_limit,
        metrics['vr'],
        metrics['vtheta'],
        metrics['vphi'],
        metrics['vgsr'],
        wind_tunnel=wind_tunnel,
        wind_distance_limits=wind_limits['distance'] if wind_tunnel else None,
    )

    #--- Tidal Effect ---#

    # te = axes[2, 1].scatter(np.log10(tidal_eff), np.log10(dw_energy), c=dw_elinfo['age'].to_numpy(), s=10, vmin=0, vmax=max(dw_elinfo['age']), cmap='PuBu')
    # axes[2, 1].scatter(np.log10(tidal_eff_numsp), np.log10(dw_energy_nusmp), c='orangered', s=260, marker='+', lw=3)
    # axes[2, 1].minorticks_on()
    # axes[2, 1].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    # axes[2, 1].grid(True)
    # axes[2, 1].set_xlim(-33.9, -29.1)
    # axes[2, 1].set_ylim(-33.9, -29.1)
    # axes[2, 1].set_ylabel(r'$log_{10} (\frac{\sigma^2_{los}}{r^2_{half}}) ~ (s^{-2})$', fontsize=xylabel_size)
    # axes[2, 1].set_xlabel(r'$log_{10} (\frac{128}{27} \frac{g^{2}_{MW}}{v^2}) ~ (s^{-2})$', fontsize=xylabel_size)

    # x_tidal = np.linspace(-100, 0, 20)
    # y_tidal = x_tidal

    # axes[2, 1].plot(x_tidal, y_tidal, lw=2, ls='--', c='orchid')

    # divider = make_axes_locatable(axes[2, 1])
    # cax = divider.append_axes("right", size="4%", pad=0.05)
    # cbar = fig.colorbar(te, cax=cax)
    # cbar.set_label(r'Simulation Time (Gyr)', fontsize=12)
    # cax.tick_params(axis='both', which='both', direction='in')
    # cax.minorticks_on()

    if wind_tunnel:
        stellar_offset = np.linalg.norm(snapshot['wind_tunnel']['stellar_center_box_kpc'])
        stellar_drift = np.linalg.norm(snapshot['wind_tunnel']['stellar_bulk_r5_box_kms'])
        fig.text(
            .05, .023,
            f'Stellar offset from box centre  {stellar_offset:.3f} kpc    Drift speed  {stellar_drift:.2f} km/s    '
            f'$R_{{\\mathrm{{half}}}}$  {metrics["r_half"]:.2f} kpc    '
            f'$M_\\star(<{stellar_region_tex_label})$  {sm_str} $M_\\odot$    '
            f'$M_{{\\mathrm{{HI}}}}(<{stellar_region_tex_label})$  {cg_str} $M_\\odot$    '
            f'Feedback decay time  {context["fd_value"]:g} Myr',
            fontsize=12, va='bottom',
        )
    else:
        d_peri = context["d_peri"]
        info_orbit_text = (
            f"Orbit : {ini_x:.2f}, {ini_y:.2f}, {ini_z:.2f}\n"
            f"        {ini_vx:.2f}, {ini_vy:.2f}, {ini_vz:.2f}\n"
        )


        fig.text(
            0.935, 0.031,
            info_orbit_text,
            fontsize=20,
            verticalalignment='bottom',
            horizontalalignment='right'
        )

        gas_frac_contour = metrics['hi_fraction_contour']

        fig.text(
            0.74, 0.32,
            f'$R_{{half}}$: {metrics["r_half"]:.2f} kpc \n'
            f'$D_{{peri}}$: {d_peri:.2f} kpc \n'
            f'Stellar Mass within ${stellar_region_tex_label}$: {sm_str} $M_{{\\odot}}$ \n'
            f'Stellar Mass within $R_{{half}}$: {shm_str} $M_{{\\odot}}$ \n'
            f'$HI Mass_{{{stellar_region_subscript_label}}}$: {cg_str} $M_{{\\odot}}$ \n'
            f'$HI Mass_{{Rhalf}}$: {cghm_str} $M_{{\\odot}}$ \n'
            f'HI Fraction in Contour: {gas_frac_contour:.3f} \n'
            f'HI Fraction in ${stellar_region_tex_label}$: {gas_frac_stellar_region:.3f} \n'
            f'HI Fraction in $R_{{half}}$: {gas_half_frac:.3f} \n'
            f'Feedback: {fd_factor_str} median feedback \n',
            fontsize=20,
            ha='left', va='top'
        )

    if verbose:
        print('Information added.')

    output_file = make_plot_output_path(output_folder, numsp)
  
    try:
        save_figure_atomic(fig, output_file)
        record_frame(output_file, context['display_identity'], numsp)
    finally:
        plt.close(fig)
    after_render = time.perf_counter()
    after_render_rss = get_rss_mb() if profile_memory else None
    
    if verbose:
        print(f"Figure saved to {output_file}")

    df._msg = None
    
    del df
    del simulation
    del fig
    del dw_elinfo
    del metrics
    del sm_str, shm_str, cg_str, cghm_str, gas_ratio_stellar_region, gas_half_ratio, gas_frac_stellar_region, gas_half_frac
    del snapshot

    release_process_memory()
    end_time = time.perf_counter()

    if show_mem:
        try:
            mem_mb = process.memory_info().rss / 1024 / 1024
            print(f"[numsp={numsp}] Memory usage at end: {mem_mb:.2f} MB")
        except Exception as e:
            print(f"[numsp={numsp}] Memory usage check failed: {e}")

    if profile:
        parts = [
            f"pid={os.getpid()}",
            f"numsp={numsp}",
            f"prepare={after_prepare - start_time:.2f}s",
            f"metrics={after_metrics - after_prepare:.2f}s",
            f"render={after_render - after_metrics:.2f}s",
            f"cleanup={end_time - after_render:.2f}s",
            f"total={end_time - start_time:.2f}s",
        ]
        if profile_memory:
            parts.extend([
                f"rss_start={start_rss:.0f}MB",
                f"rss_prepare={after_prepare_rss:.0f}MB",
                f"rss_metrics={after_metrics_rss:.0f}MB",
                f"rss_render={after_render_rss:.0f}MB",
                f"rss_end={get_rss_mb():.0f}MB",
            ])
        print("[PlotFigProfile] " + " ".join(parts), flush=True)

    return {'numsp': numsp, 'status': 'plotted', 'output_file': output_file}

def main():
    global PLOT_CONTEXT
    os.environ.setdefault("DSPH_RUN_ID", str(os.getpid()))
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--range', type=str, default=None, help='the range of snapshots, e.g., "1,20"')
    parser.add_argument('--snapshot-numbers', type=str, help='explicit comma-separated snapshot numbers, including actual gaps')
    parser.add_argument('--video-only', action='store_true', help='encode current frames without loading particle snapshots or the plotting context')
    parser.add_argument('--frames-only', action='store_true', help='render the selected snapshots without encoding a video')
    parser.add_argument('--video-name', type=str, help='plain output .mp4 filename within the current display directory')
    parser.add_argument('--processes', type=int, default=1, help='number of processes to use')
    parser.add_argument('--chunksize', type=int, default=1, help='chunksize for pool.imap_unordered')
    parser.add_argument('--dispatch', choices=('blocks', 'dynamic'), default='blocks', help='task dispatch strategy: contiguous snapshot blocks or per-snapshot dynamic scheduling')
    parser.add_argument('--block-size', type=int, default=4, help='number of consecutive snapshots per worker task when --dispatch=blocks')
    parser.add_argument('--maxtasksperchild', type=int, default=4, help='restart each worker after this many pool tasks')
    parser.add_argument('--worker-threads', type=int, default=1, help='limit native BLAS/OpenMP threads inside each worker')
    parser.add_argument('--profile', action='store_true', help='print per-snapshot timing information from each worker')
    parser.add_argument('--profile-memory', action='store_true', help='include per-snapshot RSS information in profile logs')
    parser.add_argument('--no-progress', action='store_true', help='disable tqdm progress bar output')
    parser.add_argument('--verbose-plot', action='store_true', help='print per-panel plotting messages')
    parser.add_argument('--show-mem', action='store_true', help='show memory usage for each snapshot')
    parser.add_argument('--replace', action='store_true', help='replace existing .png files')
    args = parser.parse_args()
    if args.video_only and args.frames_only:
        parser.error('--video-only and --frames-only are mutually exclusive')

    explicit_numbers = None
    if args.snapshot_numbers is not None:
        explicit_numbers = sorted(set(int(n) for n in args.snapshot_numbers.split(',') if n.strip()))
        if not explicit_numbers:
            raise ValueError('--snapshot-numbers must contain at least one snapshot')

    if args.video_only:
        wind_tunnel = wind_tunnel_frame(folder_path)
        modelname = DataProcessor.parse_model_name_details(wind_tunnel.model_name if wind_tunnel else None)['modelname']
        summary = CANONICAL_ANALYSIS_ROOT / f'elinfo_{modelname}.csv'
        require_summary_frame(require_canonical_csv(summary), wind_tunnel)
        summary_rows = pd.read_csv(summary, comment='#', usecols=['numsp'])
        summary_numbers = sorted(set(summary_rows['numsp'].astype(int).tolist()))
        if not summary_numbers:
            raise ValueError('The summary contains no snapshot rows')
        wind_sfh_snapshot = None
        if wind_tunnel:
            manifest_path = DISPLAY_OUTPUT_ROOT / 'display_manifest.json'
            manifest = json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.is_file() else {}
            wind_sfh_snapshot = manifest.get('wind_sfh_reference_snapshot', summary_numbers[-1])
        identity = display_identity(summary, wind_sfh_snapshot=wind_sfh_snapshot)
        start, end = map(int, (args.range or f'{summary_numbers[0]},{summary_numbers[-1]}').split(','))
        selected = filter_available_snapshots(start, end, summary_numbers)
        if explicit_numbers is not None:
            missing = set(explicit_numbers) - set(selected)
            if missing:
                raise ValueError('Requested snapshots absent from the current summary: ' + format_limited_ints(sorted(missing)))
            selected = explicit_numbers
        prepare_display_directory(DISPLAY_OUTPUT_ROOT, identity)
        video_name = args.video_name
        if video_name is None and ordinary_media_reuse():
            video_name = f"{modelname}_{identity['member_selection_id']}_tan_xi_eta_2p1_v1_v{PLOT_VIDEO_VERSION}.mp4"
        finalize_display_outputs(output_folder, modelname, selected,
                                 identity=identity, video_name=video_name)
        print(f'PlotFig video-only completed: {len(selected)} existing frames encoded')
        return

    PLOT_CONTEXT = build_plot_context()
    modelname = PLOT_CONTEXT['modelname']
    PLOT_CONTEXT['video_name'] = args.video_name
    if args.video_name is None and ordinary_media_reuse():
        member_id = PLOT_CONTEXT['display_identity']['member_selection_id']
        PLOT_CONTEXT['video_name'] = f'{modelname}_{member_id}_tan_xi_eta_2p1_v1_v{PLOT_VIDEO_VERSION}.mp4'
    print('CubePath: ' + str(cube_path))
    available_snapshots = DataProcessor.list_snapshot_numbers(folder_path)
    if not available_snapshots:
        raise FileNotFoundError(f"No snapshots found under {folder_path}")
    elinfo_snapshots = PLOT_CONTEXT['elinfo_snapshots']
    drawable_snapshots = sorted(set(available_snapshots) & set(elinfo_snapshots))
    if not drawable_snapshots:
        raise FileNotFoundError('No snapshots are present in both the output directory and the elinfo file')
    start, end = map(int, (args.range or f'{elinfo_snapshots[0]},{elinfo_snapshots[-1]}').split(','))
    requested_numsp = filter_available_snapshots(start, end, drawable_snapshots)
    if explicit_numbers is not None:
        unavailable = set(explicit_numbers) - set(requested_numsp)
        if unavailable:
            raise ValueError('Requested snapshots absent from available data: ' + format_limited_ints(sorted(unavailable)))
        requested_numsp = explicit_numbers
    missing_numsp = [numsp for numsp in range(start, end + 1) if numsp not in set(requested_numsp)]
    if not requested_numsp:
        raise FileNotFoundError(
            f"No snapshot numbers within requested range {start},{end} are present in both "
            f"the output directory and the elinfo file"
        )
    if missing_numsp:
        print(
            f"PlotFig: skipping {len(missing_numsp)} snapshot ids without both snapshot and elinfo data: "
            + ",".join(str(numsp) for numsp in missing_numsp[:10])
            + ("..." if len(missing_numsp) > 10 else "")
        )
    output_dir = get_plot_context()['output_folder']
    identity = get_plot_context()['display_identity']
    prepare_display_directory(Path(output_dir).parent, identity)
    if not ordinary_media_reuse():
        for numsp in requested_numsp:
            path = make_plot_output_path(output_dir, numsp)
            if is_complete_png(path) and not frame_is_current(path, identity) and not args.replace:
                raise ValueError('Existing frame has a different or missing display identity; '
                                 'use a fresh DSPH_CANONICAL_OUTPUT_ROOT: ' + path)
    existing_numsp = [
        numsp for numsp in requested_numsp
        if reusable_frame(make_plot_output_path(output_dir, numsp), identity)
    ]
    incomplete_existing_numsp = [
        numsp for numsp in requested_numsp
        if (
            not ordinary_media_reuse()
            and os.path.exists(make_plot_output_path(output_dir, numsp))
            and not is_complete_png(make_plot_output_path(output_dir, numsp))
        )
    ]
    if incomplete_existing_numsp:
        print(
            f"PlotFig: replotting {len(incomplete_existing_numsp)} incomplete PNG file(s): "
            + format_limited_ints(incomplete_existing_numsp)
        )
    if args.replace:
        pending_numsp = requested_numsp
    else:
        existing_set = set(existing_numsp)
        pending_numsp = [numsp for numsp in requested_numsp if numsp not in existing_set]

    print(
        f"PlotFig resume: total={len(requested_numsp)}, "
        f"existing={len(existing_numsp)}, pending={len(pending_numsp)}, replace={args.replace}"
    )

    if not pending_numsp:
        if not ordinary_media_reuse() and not validate_plot_outputs(output_dir, requested_numsp):
            sys.exit(2)
        finalize_display_outputs(output_dir, modelname, requested_numsp, encode_video=not args.frames_only)
        print("PlotFig completed: plotted=0, skipped=0")
        return

    tasks = [
        (
            numsp,
            args.show_mem,
            args.replace,
            args.profile or args.profile_memory,
            args.profile_memory,
            args.verbose_plot,
        )
        for numsp in pending_numsp
    ]

    with get_context("spawn" if os.name == "nt" else "fork").Pool(
        processes=args.processes,
        maxtasksperchild=args.maxtasksperchild,
        initializer=init_worker,
        initargs=(args.worker_threads,),
    ) as pool:
        plotted = 0
        skipped = 0
        progress = ProgressReporter(len(tasks), 'PlotFig', enabled=not args.no_progress, unit='snap', mininterval=0.5)

        try:
            if args.dispatch == 'blocks':
                task_blocks = build_consecutive_blocks(tasks, args.block_size)
                result_iter = pool.imap_unordered(plot_if_needed_block, task_blocks, chunksize=1)
                for block_results in result_iter:
                    for result in block_results:
                        if result['status'] == 'plotted':
                            plotted += 1
                        elif result['status'] == 'skipped':
                            skipped += 1
                    if progress is not None:
                        progress.update(len(block_results))
            else:
                result_iter = pool.imap_unordered(plot_if_needed_task, tasks, chunksize=args.chunksize)
                for result in result_iter:
                    if result['status'] == 'plotted':
                        plotted += 1
                    elif result['status'] == 'skipped':
                        skipped += 1
                    if progress is not None:
                        progress.update(1)
        finally:
            progress.close()

    if not ordinary_media_reuse() and not validate_plot_outputs(output_dir, requested_numsp):
        sys.exit(2)

    finalize_display_outputs(output_dir, modelname, requested_numsp, encode_video=not args.frames_only)
    print(f"PlotFig completed: plotted={plotted}, skipped={skipped}")


def plot_if_needed(numsp, show_mem, replace, profile=False, profile_memory=False, verbose=False):
    output_folder = get_plot_context()['output_folder']
    output_file = make_plot_output_path(output_folder, numsp)
    if (not replace) and reusable_frame(output_file, get_plot_context()['display_identity']):
        return {'numsp': numsp, 'status': 'skipped', 'output_file': output_file}
    return safe_plot(
        numsp,
        show_mem=show_mem,
        profile=profile,
        profile_memory=profile_memory,
        verbose=verbose,
    )


def plot_if_needed_task(task):
    return plot_if_needed(*task)


def plot_if_needed_block(task_block):
    return [plot_if_needed_task(task) for task in task_block]

if __name__ == "__main__":
    main() 
