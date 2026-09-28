from pathlib import Path
import re
import warnings
import numpy as np

MODEL_DIRNAME = 'output'
OUTPUT_IMG_DIRNAME = 'output_img'
LEGACY_CUBE_PATH = '../VubeL2000d4.0e-06m1.4e+05.ini'


def resolve_runtime_paths(base_dir=None):
    runtime_root = Path.cwd() if base_dir is None else Path(base_dir).resolve()
    model_dir = runtime_root / MODEL_DIRNAME
    return {
        'path': f"{runtime_root}/",
        'pathout': f"{model_dir}/",
        'folder_path': f"{model_dir}/",
        'output_folder': f"{runtime_root / OUTPUT_IMG_DIRNAME}/",
    }


def get_cube_path(base_dir=None, fallback=LEGACY_CUBE_PATH):
    runtime_root = Path.cwd() if base_dir is None else Path(base_dir).resolve()
    file_path = runtime_root / 'mergecube.par'
    pattern = re.compile(r'\S*VubeL\S*\.ini')

    if not file_path.exists():
        warnings.warn(
            f"mergecube.par not found under {runtime_root}; falling back to legacy cube_path={fallback}",
            RuntimeWarning,
            stacklevel=2,
        )
        return fallback

    with file_path.open('r', encoding='utf-8') as file:
        first_line = file.readline().strip()

    match = pattern.search(first_line)
    if match:
        return match.group(0)

    warnings.warn(
        f"Failed to parse cube path from {file_path}; falling back to legacy cube_path={fallback}",
        RuntimeWarning,
        stacklevel=2,
    )
    return fallback


_runtime_paths = resolve_runtime_paths()

path = _runtime_paths['path']
pathout = _runtime_paths['pathout']
folder_path = _runtime_paths['folder_path']
output_folder = _runtime_paths['output_folder']

orbit_file = 'mergegalaxy.par'
# dw_elinfo = './SFH_Fornax.csv'

gizmo_config = './GZWJL.PARAM'

T_min = 0.001
dw_filed = 7.5
sculptor_core_radius = 1
d_today = 139.4
dw_filed_deg = 2.1
fornax_core_radius = np.radians(dw_filed_deg) * d_today
dw_cold_gas_lower_limit = 1.46e5
xylabel_size = 16
cd_threshold = 6e19
r_3d_cut = 0.3
stellar_region_rhalf_multiplier = 8.0
obs_surface_brightness_factor_base = 2.8085305633697762e-11
obs_plot_mass_to_light_ratio = 2.0
sim_plot_mass_to_light_ratio = 2.6


def _format_multiplier_label(value):
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:g}"


stellar_region_multiplier_label = _format_multiplier_label(stellar_region_rhalf_multiplier)
stellar_region_subscript_label = f"{stellar_region_multiplier_label}Rhalf"
stellar_region_tex_label = rf"{stellar_region_multiplier_label}R_{{half}}"
stellar_region_description = rf"within ${stellar_region_tex_label}$"

mag_lower_limit, mag_upper_limit = 23, 36
sigma_lower_limit, sigma_upper_limit = 4, 13
filed_sight_lower_limit, filed_sight_upper_limit = -6, 6
dw_sight_lower_limit, dw_sight_upper_limit = -40, 40

hot_gas_color = 'orangered'
cold_gas_color = 'aquamarine'
stars_color = 'royalblue' 

Sname = './Sample1_'
fcat =Sname+'candidates.csv'
fndens = Sname+'density.csv'
cube_path = get_cube_path()

r_pc = np.array([
    56.8182, 96.5909, 126.420, 149.148, 174.716, 194.602, 214.489, 234.375, 254.261, 272.727,
    294.034, 308.239, 325.284, 342.330, 366.477, 384.943, 406.250, 424.716, 446.023, 470.170,
    491.477, 512.784, 535.511, 558.239, 579.545, 600.852, 622.159, 647.727, 666.193, 690.341,
    708.807, 735.795, 755.682, 775.568, 801.136, 828.125, 856.534, 884.943, 913.352, 943.182,
    978.693, 1018.47, 1061.08, 1109.38, 1170.45, 1284.09, 1420.45, 1680.40
])

sigma = np.array([
    11.0062, 9.77554, 9.77554, 10.4721, 10.8204, 10.4721, 10.4954, 11.7957, 10.4954, 10.7740,
    10.4025, 8.56811, 10.2864, 10.4954, 10.0077, 9.45046, 11.9814, 8.40557, 9.47368, 9.77554,
    8.38235, 9.58978, 10.8669, 9.70588, 9.12539, 9.38081, 11.1919, 10.1703, 10.2167, 10.0774,
    9.77554, 11.1687, 9.49690, 9.96130, 9.79876, 9.70588, 9.40402, 10.8437, 8.93963, 9.07895,
    9.28792, 9.77554, 8.66099, 9.56656, 7.22136, 8.19659, 7.40712, 7.15170
])

err_sigma = np.array([
    1.04489, 1.11455, 1.02167, 1.11455, 1.06811, 1.11455, 1.06811, 1.23065, 1.06811, 1.13777,
    1.09133, 0.928792, 1.06811, 1.11455, 1.04489, 1.02167, 1.27709, 0.928793, 1.11455, 1.04489,
    0.882353, 0.998450, 1.23065, 0.975233, 0.952013, 0.998452, 1.20743, 1.13777, 1.06811, 1.09133,
    1.02167, 1.23065, 0.998453, 1.02167, 0.998452, 0.998452, 0.975233, 1.13777, 0.952010, 1.02167,
    1.02167, 1.02167, 0.928794, 0.998452, 0.882354, 0.882354, 0.789474, 0.859133
])

core_radius_dict = {
    'Fornax': fornax_core_radius,
    'Sculptor': sculptor_core_radius,
}
