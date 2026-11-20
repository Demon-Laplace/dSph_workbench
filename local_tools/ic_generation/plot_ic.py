#!/usr/bin/env python3
"""Default initial-model figures, shared by model_tool/workflow builds and backfill.

Uses the in-memory exported particles during a build. Standalone --build only
reads that explicitly selected IC's component files, never simulation snapshots.
"""
import argparse
import ast
import configparser
import csv
import json
import os
from pathlib import Path
import tempfile

import numpy as np

SUFFIXES = ('components', 'central_components', 'column_density_profile',
            'Mstar_Mgas_relation', 'kinematics')
CODE_VERSION = 'ic-figures-v4-munoz-published-profile'


def default_directory(build):
    build = Path(build)
    receipt = build / 'build.json'
    if receipt.is_file():
        destination = json.loads(receipt.read_text()).get('local_output_directory')
        if destination:
            return Path(destination)
    return build.parent.parent/'isolation' if build.parent.name == 'ic' else build


def project_settings(build=None):
    if build is not None:
        build = Path(build).resolve()
        receipt = build / 'build.json'
        model = json.loads(receipt.read_text()).get('local_model_directory') if receipt.is_file() else None
        if model:
            project = Path(model).parent
            path = project/'workflow.json'
            return project, json.loads(path.read_text()) if path.is_file() else {}
        if build.parent.name == 'ic':
            project = build.parent.parent.parent
            path = project/'workflow.json'
            return project, json.loads(path.read_text()) if path.is_file() else {}
    return None, {}


def photometric_settings(build=None, mass_to_light=None, distance_kpc=None):
    """Use the existing analysis photometry without importing its runtime paths."""
    project, config = project_settings(build)
    candidates = []
    if project is not None and config.get('analysis_code'):
        candidates.append((project/config['analysis_code']).resolve()/'variable.py')
    candidates.append(Path(__file__).resolve().parents[2]/'repo/variable.py')
    value, source = 2.6, 'existing project default'
    distance, distance_source = 139.4, 'existing project distance default'
    for path in candidates:
        if not path.is_file():
            continue
        settings = {}
        for node in ast.parse(path.read_text()).body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id in (
                            'sim_plot_mass_to_light_ratio', 'd_today'):
                        settings[target.id] = float(ast.literal_eval(node.value))
        if 'sim_plot_mass_to_light_ratio' not in settings:
            continue
        value, source = settings['sim_plot_mass_to_light_ratio'], str(path)
        if 'd_today' in settings:
            distance, distance_source = settings['d_today'], str(path)
        break
    if mass_to_light is not None:
        value, source = float(mass_to_light), 'explicit --mass-to-light'
    if not np.isfinite(value) or value <= 0:
        raise ValueError('V-band mass-to-light ratio must be positive and finite')
    if distance_kpc is not None:
        distance, distance_source = float(distance_kpc), 'explicit --distance-kpc'
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError('Reference display distance must be positive and finite')
    return {'mass_to_light_V': value, 'mass_to_light_source': source,
            'reference_display_distance_kpc': distance, 'distance_source': distance_source,
            'band': 'V', 'magnitude_system': 'Vega, existing project convention',
            'solar_magnitude_input': 4.83, 'AB_to_Vega_offset': -.02,
            'arcsec_per_radian': 206265.,
            'normalization_fit_to_observations': False,
            'projection': 'Initial IC face-on (native XY); circular model annuli, no observer-pose fit.'}


def surface_mass_to_mu_v(surface_Msun_pc2, mass_to_light):
    """Exactly the existing L_to_ABmag + angular-area convention; D cancels."""
    surface = np.asarray(surface_Msun_pc2, dtype=float)
    result = np.full_like(surface, np.nan)
    positive = np.isfinite(surface) & (surface > 0)
    result[positive] = (4.83-.02+5*np.log10(206265.)-5
                        -2.5*np.log10(surface[positive]/mass_to_light))
    return result


def stellar_brightness_profile(x, mass, inner, outer, mass_to_light):
    """Face-on IC stellar light in supplied radial intervals, no gas or HI."""
    R = np.hypot(x[:, 0], x[:, 1])
    order = np.argsort(R)
    accumulated = np.r_[0., np.cumsum(mass[order])*1e10]
    lower = np.searchsorted(R[order], inner, side='left')
    upper = np.searchsorted(R[order], outer, side='left')
    annular_mass = accumulated[upper]-accumulated[lower]
    area_pc2 = np.pi*(outer**2-inner**2)*1e6
    surface = np.divide(annular_mass, area_pc2,
                        out=np.full_like(area_pc2, np.nan), where=area_pc2 > 0)
    return .5*(inner+outer), surface_mass_to_mu_v(surface, mass_to_light)


def load_surface_brightness_observation(build=None, observation_profile=None):
    """Use the project's prepared Fornax stellar mu_V values/errors unchanged."""
    project, config = project_settings(build)
    explicit = observation_profile or os.environ.get('FORNAX_PREPARED_OBSERVATION_PROFILE')
    paths = []
    if explicit:
        paths.append(Path(explicit).expanduser().resolve())
    elif project is not None and config.get('ic_plotting', {}).get('observation_profile'):
        paths.append((project/config['ic_plotting']['observation_profile']).resolve())
    else:
        if os.environ.get('FORNAX_OBSERVATION_ROOT'):
            paths.append(Path(os.environ['FORNAX_OBSERVATION_ROOT']).expanduser()/
                         'Yang_data/Fornax_surface_brightness_profile.csv')
        if project is not None:
            paths.extend([project/'workspace/inputs/observations/Fornax_surface_brightness_profile.csv',
                          project/'Fornax_surface_brightness_profile.csv'])
    for path in paths:
        if not path.is_file():
            continue
        with path.open(newline='') as stream:
            reader = csv.DictReader(stream)
            required = {'r_mid_kpc', 'mu_v_mag_arcsec2', 'mu_v_err_total'}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError('Prepared Fornax observation lacks radius, mu_V or uncertainty columns')
            rows = list(reader)
        def values(key):
            return np.array([float(row.get(key) or 'nan') for row in rows])
        data = {'radius_kpc': values('r_mid_kpc'), 'mu_v': values('mu_v_mag_arcsec2'),
                'mu_v_error': values('mu_v_err_total'),
                'inner_kpc': values('r_inner_kpc'), 'outer_kpc': values('r_outer_kpc'),
                'radius_deg': values('r_mid_deg') if 'r_mid_deg' in rows[0] else values('xr'),
                'inner_deg': values('r_inner_deg'), 'outer_deg': values('r_outer_deg')}
        finite = (np.isfinite(data['radius_kpc']) & np.isfinite(data['mu_v']) &
                  np.isfinite(data['mu_v_error']) & (data['mu_v_error'] >= 0))
        return data, {
            'state': 'loaded', 'observable': 'Fornax stellar V-band surface brightness',
            'label': 'Fornax observed: Yang 2022, project Muñoz 2018 calibration',
            'source_file': str(path), 'radius_unit': 'kpc, prepared elliptical radii as supplied',
            'surface_brightness_unit': 'mag arcsec^-2', 'uncertainty_column': 'mu_v_err_total',
            'reference_points': len(rows), 'finite_points_with_uncertainty': int(finite.sum()),
            'calibration_names': sorted({row.get('calibration_name', '') for row in rows}),
            'prepared_versions': sorted({row.get('code_version', '') for row in rows}),
            'observation_recalibrated': False, 'include_WLM': False,
            'interpretation': 'Present-day Fornax reference versus an initial face-on IC, not an evolved sky-view fit.',
        }
    return None, {'state': 'unavailable', 'observable': 'Fornax stellar V-band surface brightness',
                  'searched_paths': [str(path) for path in paths], 'include_WLM': False}


def observation_at_distance(data, distance_kpc):
    """Display angular observations at a common distance; never alter their CSV."""
    result = dict(data)
    for physical, angular in [('radius_kpc', 'radius_deg'), ('inner_kpc', 'inner_deg'),
                              ('outer_kpc', 'outer_deg')]:
        if angular in data:
            good = np.isfinite(data[angular])
            result[physical] = np.where(good, np.deg2rad(data[angular])*distance_kpc, data[physical])
    return result


def surface_image(x, y, mass, extent, bins=400):
    edges = np.linspace(-extent, extent, bins+1)
    histogram = np.histogram2d(x, y, bins=(edges, edges), weights=mass*1e10)[0]
    return histogram.T/((2*extent/bins)**2*1e6)


def curve(x, v, mass, edges):
    R = np.hypot(x[:, 0], x[:, 1])
    vR = np.sum(x[:, :2]*v[:, :2], axis=1)/np.maximum(R, 1e-15)
    vp = (x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])/np.maximum(R, 1e-15)
    result = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        use = (R >= lo) & (R < hi)
        if use.sum() < 20:
            continue
        w = mass[use]
        values = np.column_stack((vR[use], vp[use], v[use, 2]))
        mean = np.average(values, weights=w, axis=0)
        sigma = np.sqrt(np.average((values-mean)**2, weights=w, axis=0))
        result.append([np.sqrt(lo*hi), mean[1], *sigma,
                       w[vp[use] > 0].sum()/w.sum()])
    return np.asarray(result).reshape(-1, 6)


def cumulative(radius, mass, radii):
    order = np.argsort(radius)
    accumulated = np.r_[0., np.cumsum(mass[order])*1e10]
    return accumulated[np.searchsorted(radius[order], radii)]


def generate(data, model, directory, build=None, prefix=None,
             observation_profile=None, mass_to_light=None, distance_kpc=None, profile_only=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    # Temporary display-only files stay with the explicitly selected output,
    # including on servers where all writes must remain in the data workspace.
    cache = tempfile.TemporaryDirectory(prefix='.ic-matplotlib-', dir=directory)
    os.environ.setdefault('MPLCONFIGDIR', cache.name)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from scipy.ndimage import gaussian_filter
    import munoz2018
    observation, reference = load_surface_brightness_observation(build, observation_profile)
    photometry = photometric_settings(build, mass_to_light, distance_kpc)
    distance = photometry['reference_display_distance_kpc']
    if observation is not None:
        observation = observation_at_distance(observation, distance)
        reference['display_radius_convention'] = (
            'Use supplied angular radii at common D='+str(distance)+' kpc when available; '
            'mu_V/errors and source CSV unchanged.')
    munoz = munoz2018.metadata(distance)
    if reference['state'] != 'loaded':
        print('IC figures: Fornax stellar-light reference unavailable; configure --observation-profile.')

    prefix = prefix or model
    selected_suffixes = ('column_density_profile',) if profile_only else SUFFIXES
    if any((directory/(prefix+'_'+suffix+'.png')).exists() for suffix in selected_suffixes):
        if build is None:
            raise FileExistsError('Existing IC figures are preserved; choose another output prefix.')
        prefix = model+'_'+Path(build).name
    paths = [directory/(prefix+'_'+suffix+'.png') for suffix in SUFFIXES]
    selected_paths = [paths[2]] if profile_only else paths
    if any(path.exists() for path in selected_paths):
        raise FileExistsError('Figures for this build already exist; use --prefix for a new rendering.')

    populations = {}
    for name, kind in [('Gas', 0), ('Stars', 2)]:
        rows = data[data[:, 8] == kind]
        populations[name] = [rows[:, :3], rows[:, 3:6], rows[:, 6]]
    base = populations['Stars'] if len(populations['Stars'][2]) else [
        data[:, :3], data[:, 3:6], data[:, 6]]
    x, v, mass = base
    center = np.average(x, weights=mass, axis=0)
    for radius in (5., 3., 2., 1.5):
        use = np.sum((x-center)**2, axis=1) < radius**2
        if use.sum() >= 30:
            center = np.average(x[use], weights=mass[use], axis=0)
    use = np.sum((x-center)**2, axis=1) < 25.
    bulk = np.average(v[use], weights=mass[use], axis=0)
    for values in populations.values():
        values[0] = values[0]-center
        values[1] = values[1]-bulk

    potential, densities = None, {}
    if build is not None and (Path(build)/'model_potential_gas').is_file():
        try:
            import agama
            agama.setUnits(length=1, velocity=1, mass=1e10)
            potential = agama.Potential(*[
                agama.Potential(file=str(Path(build)/('model_potential_'+name)))
                for name in ('gas', 'stellar')])
            config = configparser.RawConfigParser()
            config.read(Path(build)/'parameters.ini')
            for label, section in [('Gas', 'Density_gas'), ('Stars', 'Density_stellar')]:
                densities[label] = agama.Density(**dict(config.items(section)))
        except (ImportError, ValueError, RuntimeError) as error:
            print('IC figures: analytic overlays unavailable:', error)

    def save(fig, path):
        fig.savefig(path, dpi=170, bbox_inches='tight')
        plt.close(fig)

    def result():
        return {'state': 'complete', 'directory': str(directory), 'prefix': prefix,
                'files': [path.name for path in selected_paths], 'render_recipe': CODE_VERSION,
                'observational_reference': reference, 'published_munoz_reference': munoz,
                'photometric_convention': photometry,
                'spatial_views': {'components': 'overview, half-width 5 kpc',
                                  'central_components': 'same IC, central zoom, half-width 1.5 kpc',
                                  'panels': 'rows: gas/stars; columns: XY/XZ; equal spatial scales'},
                'note': 'Initial IC diagnostics. Muñoz line is a published fit, not newly measured bins. No WLM.'}

    with plt.rc_context({'font.size': 10, 'axes.grid': False}):
        for path, extent in (() if profile_only else zip(paths[:2], (5., 1.5))):
            fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
            for row, (name, (x, _, mass)) in enumerate(populations.items()):
                for col, (a, b, plane) in enumerate([(0, 1, 'XY / face-on'), (0, 2, 'XZ / edge-on')]):
                    ax = axes[row, col]
                    image = gaussian_filter(surface_image(x[:, a], x[:, b], mass, extent), .7)
                    positive = image[image > 0]
                    if len(positive):
                        vmin = max(float(np.percentile(positive, 2)), 1e-5)
                        vmax = max(float(np.percentile(positive, 99.9)), vmin*1.01)
                        artist = ax.imshow(np.ma.masked_less_equal(image, 0), origin='lower',
                            extent=(-extent, extent, -extent, extent),
                            norm=LogNorm(vmin, vmax), cmap='magma' if name == 'Gas' else 'cividis',
                            interpolation='nearest', aspect='equal')
                        fig.colorbar(artist, ax=ax, shrink=.8, label=r'$\Sigma\ [M_\odot\,{\rm pc}^{-2}]$')
                    ax.set(xlim=(-extent, extent), ylim=(-extent, extent),
                           xlabel='XYZ'[a]+' [kpc]', ylabel='XYZ'[b]+' [kpc]',
                           title=name+' — '+plane)
                    ax.set_facecolor('#080c16')
            view = 'overview: +/-5 kpc' if extent == 5. else 'central zoom: +/-1.5 kpc'
            fig.suptitle(model+' initial IC — '+view+'\n'
                         'Rows: gas / stars; columns: XY / XZ; equal spatial scales')
            save(fig, path)

        colors = {'Gas': '#008fbd', 'Stars': '#b23b38'}
        edges = np.geomspace(.025, 10., 48)
        radii = np.sqrt(edges[:-1]*edges[1:])
        fig, (ax, light_ax) = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
        sampled_surface_values = []
        for name, (x, _, mass) in populations.items():
            R = np.hypot(x[:, 0], x[:, 1])
            binned = np.histogram(R, bins=edges, weights=mass*1e10)[0]
            surface = binned/(np.pi*np.diff(edges**2)*1e6)
            use = surface > 0
            sampled_surface_values.extend(surface[use].tolist())
            ax.plot(radii[use], surface[use], '.-', color=colors[name], label=name+' particles')
            if name in densities:
                z = np.r_[0., np.geomspace(1e-5, 10., 400)]
                rr, zz = np.meshgrid(radii, z, indexing='ij')
                points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
                rho = densities[name].density(points).reshape(rr.shape)
                target = 2*np.trapz(rho, z, axis=1)*1e4
                ax.plot(radii, target, '--', color=colors[name], alpha=.7, label=name+' analytic target')
        ax.set(xscale='log', yscale='log', xlim=(.025, 10),
               xlabel='Cylindrical R [kpc]', ylabel=r'$\Sigma\ [M_\odot\,{\rm pc}^{-2}]$',
               title='Mass surface density (total gas / stars)')
        if sampled_surface_values:
            # Do not let an analytic Sersic tail at 1e-30 hide every resolved
            # particle bin. These are display limits, not a data selection.
            maximum = max(sampled_surface_values)
            ax.set_ylim(max(1e-4, maximum*1e-5), maximum*2)
        ax.grid(alpha=.2); ax.legend(loc='lower left', fontsize=8)

        intervals = np.geomspace(.025, 5.1, 45)
        inner, outer = intervals[:-1], intervals[1:]
        if observation is not None:
            good_intervals = (np.isfinite(observation['inner_kpc']) &
                              np.isfinite(observation['outer_kpc']) &
                              (observation['inner_kpc'] >= 0) &
                              (observation['outer_kpc'] > observation['inner_kpc']))
            if good_intervals.any():
                inner = observation['inner_kpc'][good_intervals]
                outer = observation['outer_kpc'][good_intervals]
            good = (np.isfinite(observation['radius_kpc']) & np.isfinite(observation['mu_v']) &
                    np.isfinite(observation['mu_v_error']) & (observation['mu_v_error'] >= 0))
            light_ax.errorbar(observation['radius_kpc'][good], observation['mu_v'][good],
                yerr=observation['mu_v_error'][good], fmt='o', ms=3, capsize=2,
                lw=.7, color='darkorchid', ecolor='gray', alpha=.5,
                label='Yang 2022 (project Muñoz calibration)', zorder=2)
        else:
            light_ax.text(.03, .05, 'Fornax observation profile unavailable',
                          transform=light_ax.transAxes, fontsize=8, color='#a86600')
        x, _, mass = populations['Stars']
        light_radius, light_mu = stellar_brightness_profile(
            x, mass, inner, outer, photometry['mass_to_light_V'])
        light_ax.plot(light_radius, light_mu, '.-', color='black', ms=3, lw=1.2,
                      label=f"IC stars, face-on; M/L={photometry['mass_to_light_V']:g}", zorder=3)
        published_R = np.linspace(0., 5.1, 400)
        published_mu = munoz2018.sersic_mu_v(published_R, distance)
        light_ax.fill_between(published_R, published_mu-.16, published_mu+.16,
                              color='#df8a00', alpha=.14, linewidth=0)
        light_ax.plot(published_R, published_mu, '--', color='#d78300', lw=1.8,
                      label="Muñoz 2018 Sérsic fit (n=0.71, Re=18.4')", zorder=4)
        angular = light_ax.secondary_xaxis('top', functions=(
            lambda r: r/distance*180/np.pi*60, lambda theta: np.deg2rad(theta/60)*distance))
        angular.set_xlabel('Angular semi-major radius [arcmin]')
        light_ax.set(xlim=(0, 5.1), ylim=(36, 22),
                     xlabel='Radius [kpc]', ylabel=r'$\mu_V\ [{\rm mag\ arcsec}^{-2}]$',
                     title='Stellar V-band surface brightness (Vega)')
        light_ax.grid(alpha=.2); light_ax.legend(loc='lower left', fontsize=8)
        fig.suptitle(model+' initial IC — mass profiles and Fornax stellar-light reference')
        fig.supxlabel(f'Common reference distance {distance:g} kpc; IC is face-on, fixed M/L (not a sky-view fit).\n'
                      'Muñoz: published fit, not new bins; outer continuation is extrapolation. Band: ±0.16 mag normalization only.\n'
                      'Yang already uses Muñoz calibration; no independent joint-photometry fit, no WLM.', fontsize=8)
        save(fig, paths[2])
        if profile_only:
            cache.cleanup()
            return result()

        radii = np.geomspace(.03, 15., 160)
        enclosed = {}
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
        for name, (x, _, mass) in populations.items():
            enclosed[name] = cumulative(np.linalg.norm(x, axis=1), mass, radii)
            axes[0].plot(radii, enclosed[name], color=colors[name], label=name)
        denominator = enclosed['Gas']+enclosed['Stars']
        fraction = np.divide(enclosed['Gas'], denominator,
                             out=np.full_like(denominator, np.nan), where=denominator > 0)
        axes[1].plot(radii, fraction, color=colors['Gas'])
        axes[0].set(xscale='log', yscale='log', xlabel='3D radius [kpc]',
                    ylabel=r'$M(<r)\ [M_\odot]$', title='Enclosed component masses')
        axes[0].legend()
        axes[1].set(xscale='log', xlabel='3D aperture radius [kpc]', ylim=(0, 1),
                    ylabel=r'$M_{\rm gas}/(M_{\rm gas}+M_*)$', title='Total gas fraction — same aperture')
        for ax in axes:
            for radius in (1., 5., 15.):
                ax.axvline(radius, color='gray', alpha=.25, lw=.8)
            ax.grid(alpha=.2)
        fig.suptitle(model+' initial mass budget (total gas, not HI)')
        save(fig, paths[3])

        fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
        curves = {}
        for col, name in enumerate(('Stars', 'Gas')):
            x, v, mass = populations[name]
            values = curve(x, v, mass, np.geomspace(.035, 5., 28))
            curves[name] = values
            ax = axes[0, col]
            if len(values):
                for k, label, color in [(1, r'$\overline{v_\phi}$', '#202020'),
                                        (2, r'$\sigma_R$', '#c34a36'),
                                        (3, r'$\sigma_\phi$', '#336da4'),
                                        (4, r'$\sigma_z$', '#4e9563')]:
                    ax.plot(values[:, 0], values[:, k], '.-', label=label, color=color)
            if potential is not None:
                R = np.geomspace(.035, 5., 100)
                xyz = np.column_stack((R, np.zeros_like(R), np.zeros_like(R)))
                vc = np.sqrt(np.maximum(-R*potential.force(xyz)[:, 0], 0))
                ax.plot(R, vc, '--', color='#9e78b5', label=r'$V_c$ (smooth IC model)')
            ax.set(xlim=(0, 3 if name == 'Stars' else 5), xlabel='R [kpc]',
                   ylabel='Velocity [km/s]', title=name+' — annular mean and dispersions')
            ax.legend(fontsize=8); ax.grid(alpha=.2)
        values = curves['Stars']
        if len(values):
            for sigma, label, color in [
                (values[:, 2], r'$|V|/\sigma_R$', '#c34a36'),
                (values[:, 3], r'$|V|/\sigma_\phi$', '#336da4'),
                (np.sqrt(np.mean(values[:, 2:5]**2, axis=1)), r'$|V|/\sigma_{\rm 1d,cyl}$', '#202020')]:
                ratio = np.divide(abs(values[:, 1]), sigma,
                                  out=np.full_like(sigma, np.nan), where=sigma > 1e-9)
                axes[1, 0].plot(values[:, 0], ratio, '.-', color=color, label=label)
        axes[1, 0].axhline(1., ls=':', color='gray')
        axes[1, 0].axhline(2., ls='--', color='gray', alpha=.5)
        axes[1, 0].set(xlim=(0, 3), xlabel='R [kpc]',
                       ylabel='Intrinsic rotation / dispersion', title='Stellar rotation dominance')
        axes[1, 0].set_ylim(bottom=0)
        axes[1, 0].legend(fontsize=8)
        for name, values in curves.items():
            if len(values):
                axes[1, 1].plot(values[:, 0], values[:, 5], '.-', color=colors[name], label=name)
        axes[1, 1].axhline(.5, ls=':', color='gray')
        axes[1, 1].set(xlim=(0, 5), ylim=(0, 1.02), xlabel='R [kpc]',
                       ylabel='Prograde mass fraction', title='Direction occupancy is not disk coldness')
        axes[1, 1].legend()
        fig.suptitle(model+' initial kinematics — intrinsic, not LOS; plotted bins have N>=20')
        save(fig, paths[4])
    cache.cleanup()
    return result()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--prefix')
    parser.add_argument('--observation-profile', type=Path,
                        help='prepared Fornax stellar mu_V CSV; default: configured ic_plotting.observation_profile')
    parser.add_argument('--mass-to-light', type=float,
                        help='explicit stellar V-band M/L; default: existing analysis sim_plot_mass_to_light_ratio')
    parser.add_argument('--distance-kpc', type=float,
                        help='common reference-radius distance; default: existing analysis d_today')
    parser.add_argument('--profile-only', action='store_true',
                        help='refresh only the mass/surface-brightness profile, not the other IC figures')
    args = parser.parse_args()
    receipt = json.loads((args.build/'build.json').read_text())
    arrays = []
    for kind, name in [(0, 'gas'), (2, 'stellar')]:
        data = np.loadtxt(args.build/('model_component_'+name+'.nbody'), ndmin=2)
        arrays.append(np.column_stack((data, np.zeros(len(data)), np.full(len(data), kind))))
    model = Path(receipt['source_ini']).stem.removeprefix('IC_')
    result = generate(np.vstack(arrays), model, args.output or default_directory(args.build),
                       build=args.build, prefix=args.prefix,
                      observation_profile=args.observation_profile, mass_to_light=args.mass_to_light,
                      distance_kpc=args.distance_kpc, profile_only=args.profile_only)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
