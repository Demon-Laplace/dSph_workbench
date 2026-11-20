#!/usr/bin/env python3
"""Default initial-model figures, shared by model_tool/workflow builds and backfill.

Uses the in-memory exported particles during a build. Standalone --build only
reads that explicitly selected IC's component files, never simulation snapshots.
"""
import argparse
import configparser
import json
import os
from pathlib import Path

import numpy as np

SUFFIXES = ('components', 'central_components', 'column_density_profile',
            'Mstar_Mgas_relation', 'kinematics')


def default_directory(build):
    build = Path(build)
    return build.parent.parent/'isolation' if build.parent.name == 'ic' else build


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


def generate(data, model, directory, build=None, prefix=None):
    # Config/cache writes stay inside the caller's approved temporary directory.
    cache = Path(os.environ.get('TMPDIR', '/tmp/opencode'))/'ic-matplotlib'
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(cache))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from scipy.ndimage import gaussian_filter

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    prefix = prefix or model
    if any((directory/(prefix+'_'+suffix+'.png')).exists() for suffix in SUFFIXES):
        if build is None:
            raise FileExistsError('Existing IC figures are preserved; choose another output prefix.')
        prefix = model+'_'+Path(build).name
    paths = [directory/(prefix+'_'+suffix+'.png') for suffix in SUFFIXES]
    if any(path.exists() for path in paths):
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

    with plt.rc_context({'font.size': 10, 'axes.grid': False}):
        for path, extent in zip(paths[:2], (5., 1.5)):
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
            fig.suptitle(model+' initial IC — mass-weighted projections; equal spatial scales')
            save(fig, path)

        colors = {'Gas': '#008fbd', 'Stars': '#b23b38'}
        edges = np.geomspace(.025, 10., 48)
        radii = np.sqrt(edges[:-1]*edges[1:])
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
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
               title=model+' initial face-on surface density')
        if sampled_surface_values:
            # Do not let an analytic Sersic tail at 1e-30 hide every resolved
            # particle bin. These are display limits, not a data selection.
            maximum = max(sampled_surface_values)
            ax.set_ylim(max(1e-4, maximum*1e-5), maximum*2)
        ax.grid(alpha=.2); ax.legend()
        fig.supxlabel('Gas means total IC gas, not an assumed neutral-HI mass.', fontsize=9)
        save(fig, paths[2])

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
    return {'state': 'complete', 'directory': str(directory), 'prefix': prefix,
            'files': [path.name for path in paths],
            'render_recipe': 'mass_weighted_ic_figures_v2_resolved_profile_range',
            'note': 'Initial IC scientific figures; no snapshot evolution or HI ionization assumption.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--prefix')
    args = parser.parse_args()
    receipt = json.loads((args.build/'build.json').read_text())
    arrays = []
    for kind, name in [(0, 'gas'), (2, 'stellar')]:
        data = np.loadtxt(args.build/('model_component_'+name+'.nbody'), ndmin=2)
        arrays.append(np.column_stack((data, np.zeros(len(data)), np.full(len(data), kind))))
    model = Path(receipt['source_ini']).stem.removeprefix('IC_')
    result = generate(np.vstack(arrays), model, args.output or default_directory(args.build),
                      build=args.build, prefix=args.prefix)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
