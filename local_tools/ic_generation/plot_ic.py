#!/usr/bin/env python3
"""One initial-model PDF report, shared by model_tool/workflow builds and refresh.

Uses the in-memory exported particles during a build. Standalone --build only
reads that explicitly selected IC's component files, never simulation snapshots.
"""
import argparse
import ast
import configparser
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile

import numpy as np

LEGACY_PNG_SUFFIXES = ('components', 'central_components', 'column_density_profile',
            'Mstar_Mgas_relation', 'kinematics', 'gas_density_profile')
CODE_VERSION = 'ic-report-v8-particle-scatter'


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
    """Annular R, mean v_phi, cylindrical sigmas, prograde fraction, XYZ sigmas.

    Each dispersion subtracts its own mass-weighted annular component mean.
    Cartesian scatter over a whole annulus includes azimuthally varying rotation.
    """
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
        xyz_mean = np.average(v[use], weights=w, axis=0)
        xyz_sigma = np.sqrt(np.average((v[use]-xyz_mean)**2, weights=w, axis=0))
        result.append([np.sqrt(lo*hi), mean[1], *sigma,
                       w[vp[use] > 0].sum()/w.sum(), *xyz_sigma])
    return np.asarray(result).reshape(-1, 9)


def cumulative(radius, mass, radii):
    order = np.argsort(radius)
    accumulated = np.r_[0., np.cumsum(mass[order])*1e10]
    return accumulated[np.searchsorted(radius[order], radii)]


def spherical_density_profile(x, mass, edges):
    """Mass density in spherical shells, in Msun/pc^3; empty shells are undefined.

    Input masses use the IC unit of 1e10 Msun and positions/edges use kpc.
    The error is a particle shot-noise estimate sqrt(sum(m_i^2))/shell_volume.
    """
    edges = np.asarray(edges, dtype=float)
    radius = np.linalg.norm(x, axis=1)
    mass_msun = np.asarray(mass) * 1e10
    counts = np.histogram(radius, bins=edges)[0]
    shell_mass = np.histogram(radius, bins=edges, weights=mass_msun)[0]
    shell_mass2 = np.histogram(radius, bins=edges, weights=mass_msun**2)[0]
    volume_pc3 = 4*np.pi/3*np.diff(np.asarray(edges)**3)*1e9
    rho = np.divide(shell_mass, volume_pc3,
                    out=np.full_like(shell_mass, np.nan), where=counts > 0)
    error = np.divide(np.sqrt(shell_mass2), volume_pc3,
                      out=np.full_like(shell_mass2, np.nan), where=counts > 0)
    return .5*(edges[:-1]+edges[1:]), rho, error, counts


def kinematic_summary(x, v, mass):
    """Mass-weighted intrinsic moments about native +Z in an already centered frame."""
    if not len(mass):
        return {'N': 0, 'mass_Msun': 0., 'mean_speed_3d_kms': None,
                'mean_vphi_kms': None, 'sigma_xyz_kms': None,
                'prograde_mass_fraction': None, 'retrograde_mass_fraction': None,
                'zero_Lz_mass_fraction': None}
    R = np.hypot(x[:, 0], x[:, 1])
    Lz = x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0]
    vp = np.divide(Lz, R, out=np.zeros_like(R), where=R > 0)
    mean = np.average(v, weights=mass, axis=0)
    rotating = R > 0
    total = mass.sum()
    return {'N': int(len(mass)), 'mass_Msun': float(total*1e10),
            'mean_speed_3d_kms': float(np.average(np.linalg.norm(v, axis=1), weights=mass)),
            'mean_vphi_kms': float(np.average(vp[rotating], weights=mass[rotating]))
                if rotating.any() else None,
            'sigma_xyz_kms': np.sqrt(np.average((v-mean)**2, weights=mass, axis=0)).tolist(),
            'prograde_mass_fraction': float(mass[Lz > 0].sum()/total),
            'retrograde_mass_fraction': float(mass[Lz < 0].sum()/total),
            'zero_Lz_mass_fraction': float(mass[Lz == 0].sum()/total)}


def radial_kinematics(x, v, mass, edges):
    """Spherical-shell r, mean speed, mean v_phi, sigma_xyz, prograde fraction."""
    radius = np.linalg.norm(x, axis=1)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        use = (radius >= lo) & (radius < hi)
        if use.sum() < 20:
            continue
        values = kinematic_summary(x[use], v[use], mass[use])
        rows.append([np.sqrt(lo*hi), values['mean_speed_3d_kms'],
                     np.nan if values['mean_vphi_kms'] is None else values['mean_vphi_kms'],
                     *values['sigma_xyz_kms'], values['prograde_mass_fraction']])
    return np.asarray(rows).reshape(-1, 7)


def _structure_page(populations, model, plt):
    from matplotlib.colors import LogNorm
    from scipy.ndimage import gaussian_filter
    fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27), constrained_layout=True)
    extent = 3.
    for row, name in enumerate(('Stars', 'Gas')):
        x, _, mass = populations[name]
        for col, (a, b, plane) in enumerate(((0, 1, 'XY'), (0, 2, 'XZ'))):
            ax = axes[row, col]
            values = gaussian_filter(surface_image(x[:, a], x[:, b], mass, extent), .7)
            positive = values[values > 0]
            if len(positive):
                low = max(float(np.percentile(positive, 2)), 1e-5)
                high = max(float(np.percentile(positive, 99.9)), low*1.01)
                artist = ax.imshow(np.ma.masked_less_equal(values, 0), origin='lower',
                                   extent=(-extent, extent, -extent, extent),
                                   norm=LogNorm(low, high), aspect='equal',
                                   cmap='cividis' if name == 'Stars' else 'magma')
                fig.colorbar(artist, ax=ax, shrink=.8, label=r'$\Sigma\ [M_\odot\,{\rm pc}^{-2}]$')
            ax.set(xlim=(-extent, extent), ylim=(-extent, extent),
                   xlabel='XYZ'[a]+' [kpc]', ylabel='XYZ'[b]+' [kpc]', title=name+' — '+plane)
            ax.set_facecolor('#080c16')
    fig.suptitle(model+' initial conditions — spatial structure', fontsize=15)
    fig.supxlabel('One common +/-3 kpc window; native coordinates. Statistics use the stated apertures, not this display crop.', fontsize=9)
    return fig


def _profiles_page(populations, model, densities, observation, photometry, plt):
    import munoz2018
    colors = {'Gas': '#008fbd', 'Stars': '#b23b38'}
    fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27), constrained_layout=True)
    edges = np.geomspace(.025, 10., 48)
    radii = np.sqrt(edges[:-1]*edges[1:])
    ax, light_ax = axes[0]
    surfaces = []
    for name, (x, _, mass) in populations.items():
        R = np.hypot(x[:, 0], x[:, 1])
        surface = np.histogram(R, bins=edges, weights=mass*1e10)[0]/(np.pi*np.diff(edges**2)*1e6)
        surface[surface <= 0] = np.nan
        surfaces.extend(surface[np.isfinite(surface)])
        ax.plot(radii, surface, '.-', color=colors[name], ms=3, label=name)
        if name in densities:
            z = np.r_[0., np.geomspace(1e-5, 10., 400)]
            rr, zz = np.meshgrid(radii, z, indexing='ij')
            points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
            target = 2*np.trapz(densities[name].density(points).reshape(rr.shape), z, axis=1)*1e4
            ax.plot(radii, target, '--', color=colors[name], alpha=.65)
    ax.set(xscale='log', yscale='log', xlim=(.025, 10), xlabel='Cylindrical R [kpc]',
           ylabel=r'$\Sigma\ [M_\odot\,{\rm pc}^{-2}]$', title='Mass surface density')
    if surfaces:
        ax.set_ylim(max(1e-4, max(surfaces)*1e-5), max(surfaces)*2)
    ax.legend(fontsize=8)
    if densities:
        ax.text(.03, .04, 'Dashed: analytic recipe', transform=ax.transAxes, fontsize=8)

    intervals = np.geomspace(.025, 5.1, 45)
    inner, outer = intervals[:-1], intervals[1:]
    if observation is not None:
        good = (np.isfinite(observation['inner_kpc']) & np.isfinite(observation['outer_kpc']) &
                (observation['inner_kpc'] >= 0) & (observation['outer_kpc'] > observation['inner_kpc']))
        if good.any():
            inner, outer = observation['inner_kpc'][good], observation['outer_kpc'][good]
        good = (np.isfinite(observation['radius_kpc']) & np.isfinite(observation['mu_v']) &
                np.isfinite(observation['mu_v_error']) & (observation['mu_v_error'] >= 0))
        light_ax.errorbar(observation['radius_kpc'][good], observation['mu_v'][good],
                          yerr=observation['mu_v_error'][good], fmt='o', ms=2, lw=.6,
                          color='darkorchid', ecolor='gray', alpha=.45, label='Yang: project calibration')
    else:
        light_ax.text(.03, .05, 'Observed profile unavailable', transform=light_ax.transAxes, fontsize=8)
    x, _, mass = populations['Stars']
    light_radius, mu = stellar_brightness_profile(x, mass, inner, outer, photometry['mass_to_light_V'])
    light_ax.plot(light_radius, mu, '.-', ms=3, color='black', label='IC stars, face-on')
    distance = photometry['reference_display_distance_kpc']
    published_r = np.linspace(0., 5.1, 400)
    published_mu = munoz2018.sersic_mu_v(published_r, distance)
    light_ax.fill_between(published_r, published_mu-.16, published_mu+.16, color='#df8a00', alpha=.14)
    light_ax.plot(published_r, published_mu, '--', color='#d78300', label='Muñoz 2018 Sérsic fit')
    light_ax.set(xlim=(0, 5.1), ylim=(36, 22), xlabel='Radius [kpc]',
                 ylabel=r'$\mu_V\ [{\rm mag\ arcsec}^{-2}]$', title='Stellar V-band surface brightness')
    angular = light_ax.secondary_xaxis('top', functions=(
        lambda r: r/distance*180/np.pi*60, lambda theta: np.deg2rad(theta/60)*distance))
    angular.set_xlabel('Angular radius [arcmin]', fontsize=8)
    light_ax.legend(fontsize=7, loc='lower left')

    ax = axes[1, 0]
    x, _, mass = populations['Gas']
    radius, rho, error, count = spherical_density_profile(x, mass, np.geomspace(.025, 10., 40))
    ax.plot(radius, rho, '-', color=colors['Gas'], lw=1)
    for use, filled in ((count >= 20, True), ((count > 0) & (count < 20), False)):
        if use.any():
            ax.errorbar(radius[use], rho[use], yerr=error[use], fmt='o', ms=3, lw=.6,
                        color=colors['Gas'], mfc=colors['Gas'] if filled else 'white')
    ax.set(xscale='log', yscale='log', xlim=(.025, 10), xlabel='3D radius r [kpc]',
           ylabel=r'$\rho_{\rm gas}\ [M_\odot\,{\rm pc}^{-3}]$', title='Gas volume density: spherical-shell mean')
    positive = np.isfinite(rho) & (rho > 0)
    if positive.any():
        ax.set_ylim(float(rho[positive].min())*.3, float((rho+error)[positive].max())*2)
    ax.text(.03, .04, 'Open: N < 20; bars: particle shot noise', transform=ax.transAxes, fontsize=8)

    ax = axes[1, 1]
    radii = np.geomspace(.03, 15., 160)
    for name, (x, _, mass) in populations.items():
        values = cumulative(np.linalg.norm(x, axis=1), mass, radii)
        values[values <= 0] = np.nan
        ax.plot(radii, values, color=colors[name], label=name)
    ax.set(xscale='log', yscale='log', xlabel='3D aperture radius r [kpc]',
           ylabel=r'$M(<r)\ [M_\odot]$', title='Enclosed mass — total gas / stars')
    ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.18)
    fig.suptitle(model+' initial conditions — density, light and mass', fontsize=15)
    fig.supxlabel(f'Light: native face-on IC, fixed M/L={photometry["mass_to_light_V"]:g}, D={distance:g} kpc; not an evolved sky-view fit.\n'
                  'Muñoz fit is extrapolated outside its footprint; band is ±0.16 mag normalization only. Yang already uses Muñoz calibration.\n'
                  'Gas density is a spherical-shell average, not the disk mid-plane density; total gas is not HI.', fontsize=8)
    return fig


def _kinematics_page(populations, model, summaries, plt):
    fig = plt.figure(figsize=(11.69, 8.27), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(.75, 1., 1.))
    table_ax = fig.add_subplot(grid[0, :])
    table_ax.axis('off')
    rows = []
    for name in ('Stars', 'Gas'):
        for label, values in summaries[name].items():
            def number(value, pattern='.2f'):
                return '—' if value is None else format(value, pattern)
            rows.append([name, label, str(values['N']), number(values['mass_Msun'], '.3g'),
                         number(values['mean_speed_3d_kms']), number(values['mean_vphi_kms']),
                         number(None if values['prograde_mass_fraction'] is None else 100*values['prograde_mass_fraction'])])
    table = table_ax.table(cellText=rows,
                          colLabels=['Component', '3D aperture', 'N', 'Mass [Msun]',
                                     '<v3D> [km/s]', 'Vrot [km/s]', 'Prograde [% mass]'],
                          cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.25)
    colors = {'Gas': '#008fbd', 'Stars': '#b23b38'}
    axes = np.array([[fig.add_subplot(grid[row, col]) for col in range(2)] for row in (1, 2)])
    maximum_radius = max((float(np.linalg.norm(values[0], axis=1).max())
                          for values in populations.values() if len(values[0])), default=0.)
    for col, name in enumerate(('Stars', 'Gas')):
        x, v, mass = populations[name]
        radius = np.linalg.norm(x, axis=1)
        speed = np.linalg.norm(v, axis=1)
        R = np.hypot(x[:, 0], x[:, 1])
        defined = R > 0
        vp = (x[defined, 0]*v[defined, 1]-x[defined, 1]*v[defined, 0])/R[defined]
        style = {'s': .5, 'alpha': .12 if name == 'Stars' else .3,
                 'color': colors[name], 'edgecolors': 'none', 'rasterized': True, 'label': name}
        axes[0, 0].scatter(radius, speed, **style)
        axes[0, 1].scatter(radius[defined], vp, **style)
        values = radial_kinematics(x, v, mass, np.geomspace(.025, 5., 28))
        if len(values):
            for index, label, color in ((3, r'$\sigma_x$', '#c34a36'),
                                        (4, r'$\sigma_y$', '#336da4'), (5, r'$\sigma_z$', '#4e9563')):
                axes[1, col].plot(values[:, 0], values[:, index], '.-', ms=3, color=color, label=label)
        axes[1, col].set(title=name+' — three-axis velocity dispersion', ylabel=r'$\sigma\ [\mathrm{km/s}]$', ylim=(0, None))
    axes[0, 0].set(title='Particle 3D speed vs radius', ylabel=r'$|\mathbf{v}-\mathbf{v}_{\rm bulk}|\ [\mathrm{km/s}]$', ylim=(0, None))
    axes[0, 1].set(title='Particle rotation velocity vs radius', ylabel=r'$v_\phi\ [\mathrm{km/s}]$')
    axes[0, 1].axhline(0, color='gray', lw=.7, alpha=.5)
    for ax in axes.flat:
        ax.set(xlabel='3D radius r [kpc]')
        ax.grid(alpha=.18)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8, markerscale=5 if ax in axes[0] else 1)
    for ax in axes[0]:
        ax.set_xlim(0, max(.1, maximum_radius*1.03))
    for ax in axes[1]:
        ax.set_xlim(0, 5)
    fig.suptitle(model+' initial conditions — speed, rotation and dispersion', fontsize=15)
    fig.supxlabel('Top panels: one point per particle, no subsampling; v_phi excludes R=0. Table statistics are mass weighted.\n'
                  'Prograde: Lz > 0 about native +Z. XYZ dispersions use N>=20 shells and include rotation varying across azimuth.\n'
                  'Velocities are relative to the stellar-core bulk frame (gas frame if no stars); no LOS projection.', fontsize=8)
    return fig


def generate(data, model, directory, build=None, prefix=None,
             observation_profile=None, mass_to_light=None, distance_kpc=None, profile_only=False,
             overwrite=False):
    """Write one atomic PDF report. Derived figures never modify particles."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    prefix = prefix or model
    suffix = 'initial_profiles' if profile_only else 'initial_diagnostics'
    destination = directory/(prefix+'_'+suffix+'.pdf')
    if destination.exists() and not overwrite:
        if build is None:
            raise FileExistsError('Existing IC report is preserved; use --overwrite or another prefix.')
        prefix = model+'_'+Path(build).name
        destination = directory/(prefix+'_'+suffix+'.pdf')
        if destination.exists():
            raise FileExistsError('Report for this build exists; use --overwrite or another prefix.')
    with tempfile.TemporaryDirectory(prefix='.ic-matplotlib-', dir=directory) as cache:
        os.environ.setdefault('MPLCONFIGDIR', cache)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import munoz2018
        observation, reference = load_surface_brightness_observation(build, observation_profile)
        photometry = photometric_settings(build, mass_to_light, distance_kpc)
        distance = photometry['reference_display_distance_kpc']
        if observation is not None:
            observation = observation_at_distance(observation, distance)
            reference['display_radius_convention'] = 'Supplied angular radii at common D='+str(distance)+' kpc; mu_V/errors unchanged.'
        populations = {}
        for name, kind in (('Gas', 0), ('Stars', 2)):
            rows = data[data[:, 8] == kind]
            populations[name] = [rows[:, :3].copy(), rows[:, 3:6].copy(), rows[:, 6]]
        base = populations['Stars'] if len(populations['Stars'][2]) else populations['Gas']
        if not len(base[2]):
            base = [data[:, :3], data[:, 3:6], data[:, 6]]
        x, v, mass = base
        center = np.average(x, weights=mass, axis=0)
        for radius in (5., 3., 2., 1.5):
            use = np.sum((x-center)**2, axis=1) < radius**2
            if use.sum() >= 30:
                center = np.average(x[use], weights=mass[use], axis=0)
        use = np.sum((x-center)**2, axis=1) < 25.
        if not use.any():
            use = np.ones(len(x), dtype=bool)
        bulk = np.average(v[use], weights=mass[use], axis=0)
        summaries = {}
        for name, values in populations.items():
            values[0] -= center
            values[1] -= bulk
            x, v, mass = values
            r = np.linalg.norm(x, axis=1)
            summaries[name] = {label: kinematic_summary(x[use], v[use], mass[use]) for label, use in (
                ('All', np.ones(len(mass), dtype=bool)), ('r < 1 kpc', r < 1.), ('r < 5 kpc', r < 5.))}
        densities = {}
        if build is not None and (Path(build)/'parameters.ini').is_file():
            try:
                import agama
                agama.setUnits(length=1, velocity=1, mass=1e10)
                config = configparser.RawConfigParser()
                config.read(Path(build)/'parameters.ini')
                for label, section in (('Gas', 'Density_gas'), ('Stars', 'Density_stellar')):
                    if config.has_section(section):
                        densities[label] = agama.Density(**dict(config.items(section)))
            except (ImportError, ValueError, RuntimeError) as error:
                print('IC report: analytic surface-density overlays unavailable:', error)
        pages = ['Density, light and mass'] if profile_only else ['Spatial structure', 'Density, light and mass', 'Kinematics and aperture statistics']
        descriptor, temporary = tempfile.mkstemp(prefix='.'+destination.stem+'-', suffix='.pdf', dir=directory)
        os.close(descriptor)
        try:
            with PdfPages(temporary, metadata={'Title': model+' initial-condition diagnostics',
                                              'Subject': 'Intrinsic IC structure, density and kinematics'}) as pdf, \
                 plt.rc_context({'font.size': 9, 'axes.titlesize': 11}):
                renderers = [] if profile_only else [lambda: _structure_page(populations, model, plt)]
                renderers.append(lambda: _profiles_page(populations, model, densities, observation, photometry, plt))
                if not profile_only:
                    renderers.append(lambda: _kinematics_page(populations, model, summaries, plt))
                for renderer in renderers:
                    fig = renderer()
                    try:
                        pdf.savefig(fig, dpi=170)
                    finally:
                        plt.close(fig)
            os.replace(temporary, destination)
        finally:
            Path(temporary).unlink(missing_ok=True)
        return {'state': 'complete', 'directory': str(directory), 'prefix': prefix,
                'files': [destination.name], 'render_recipe': CODE_VERSION,
                'format': 'PDF', 'page_count': len(pages), 'page_titles': pages,
                'observational_reference': reference, 'published_munoz_reference': munoz2018.metadata(distance),
                'photometric_convention': photometry,
                'frame': {'center_kpc': center.tolist(), 'bulk_velocity_kms': bulk.tolist(),
                          'rotation_axis': 'native +Z; not fitted or flipped'},
                'kinematic_definitions': {
                    'weighting': 'particle mass', 'radial_bins': '3D spherical shells, N >= 20',
                    'speed_3d': 'individual |v - v_bulk| vs 3D radius; table reports mass mean, not net rotation',
                    'system_rotation': 'mass mean v_phi about native +Z; R=0 particles excluded from v_phi mean',
                    'particle_scatter': 'all particles, no subsampling or bin averaging; rasterized points, full radial extent',
                    'scatter_particle_counts': {name: {'speed_3d': int(len(values[2])),
                        'v_phi': int(np.count_nonzero(np.hypot(values[0][:, 0], values[0][:, 1]) > 0))}
                        for name, values in populations.items()},
                    'prograde': 'mass with Lz > 0 / total mass in the stated 3D aperture',
                    'dispersion_xyz': 'sqrt(<(v_i-<v_i>_shell)^2>_mass); includes azimuthally varying ordered rotation'},
                'aperture_statistics': summaries,
                'gas_density_profile': {'geometry': '3D spherical-shell average; not disk mid-plane density',
                                        'density_unit': 'Msun/pc^3', 'component': 'total gas, not HI',
                                        'error': 'particle shot noise only; open markers indicate N < 20'},
                'note': 'Initial IC diagnostics; no evolution, sky-view fitting or stability claim.'}


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
                        help='write a separate one-page density/light/mass PDF instead of the full report')
    parser.add_argument('--overwrite', action='store_true',
                        help='replace PDF in place, remove superseded PNGs and update initial_figures.json; never rebuild particles')
    args = parser.parse_args()
    receipt = json.loads((args.build/'build.json').read_text())
    arrays = []
    for kind, name in [(0, 'gas'), (2, 'stellar')]:
        data = np.loadtxt(args.build/('model_component_'+name+'.nbody'), ndmin=2)
        arrays.append(np.column_stack((data, np.zeros(len(data)), np.full(len(data), kind))))
    model = Path(receipt['source_ini']).stem.removeprefix('IC_')
    prefix = args.prefix
    if args.overwrite and prefix is None:
        prefix = receipt.get('plots', {}).get('prefix')
    result = generate(np.vstack(arrays), model, args.output or default_directory(args.build),
                       build=args.build, prefix=prefix,
                      observation_profile=args.observation_profile, mass_to_light=args.mass_to_light,
                      distance_kpc=args.distance_kpc, profile_only=args.profile_only,
                      overwrite=args.overwrite)
    if args.overwrite:
        removed = []
        suffixes = ('column_density_profile',) if args.profile_only else LEGACY_PNG_SUFFIXES
        for suffix in suffixes:
            legacy = Path(result['directory'])/(result['prefix']+'_'+suffix+'.png')
            if legacy.is_file():
                legacy.unlink()
                removed.append(legacy.name)
        refresh = {**result, 'source_build_id': receipt['build_id'],
                   'refreshed_utc': datetime.now(timezone.utc).isoformat(),
                   'particles_regenerated': False, 'superseded_pngs_removed': removed}
        (Path(result['directory'])/'initial_figures.json').write_text(json.dumps(refresh, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
