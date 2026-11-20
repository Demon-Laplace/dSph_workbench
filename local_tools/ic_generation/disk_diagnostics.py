"""Cold stellar-disk measurements in the initial smooth axisymmetric potential.

No particle IDs or provenance gates. Q is a local single-component thin-disk
indicator, not a promise of coupled gas/star or global non-axisymmetric stability.
"""
import configparser
from pathlib import Path

import numpy as np


def average(values, mass):
    return np.average(values, weights=mass, axis=0)


def quantiles(values, mass):
    order = np.argsort(values)
    cumulative = (np.cumsum(mass[order])-.5*mass[order])/mass.sum()
    return np.interp([.1, .5, .9], cumulative, values[order]).tolist()


def toomre_q(sigma_R, epicycle, surface_density, G):
    return sigma_R*epicycle/(3.36*G*surface_density)


def circularity(potential, x, v):
    energy = potential.potential(x)+.5*np.sum(v*v, axis=1)
    valid = np.isfinite(energy) & (energy < 0)
    epsilon = np.full(len(x), np.nan)
    if valid.any():
        rc = np.asarray(potential.Rcirc(E=energy[valid]))
        xyz = np.column_stack((rc, np.zeros_like(rc), np.zeros_like(rc)))
        circular_speed = np.sqrt(np.maximum(-rc*potential.force(xyz)[:, 0], 0.))
        Lc = rc*circular_speed
        Lz = x[valid, 0]*v[valid, 1]-x[valid, 1]*v[valid, 0]
        epsilon[valid] = np.divide(Lz, Lc, out=np.full(len(Lc), np.nan), where=Lc > 0)
    return epsilon, energy


def kinematics(x, v, mass, use):
    x, v, mass = x[use], v[use], mass[use]
    if not len(mass):
        return None
    R = np.hypot(x[:, 0], x[:, 1])
    vR = np.sum(x[:, :2]*v[:, :2], axis=1)/np.maximum(R, 1e-15)
    vp = (x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])/np.maximum(R, 1e-15)
    velocities = np.column_stack((vR, vp, v[:, 2]))
    mean = average(velocities, mass)
    sigma = np.sqrt(average((velocities-mean)**2, mass))
    sigma1d = np.sqrt(np.mean(sigma*sigma))
    return {'N': len(mass), 'mass_Msun': float(mass.sum()*1e10),
            'mean_vphi_kms': float(mean[1]), 'sigma_R_phi_z_kms': sigma.tolist(),
            'V_over_sigma_R': float(abs(mean[1])/sigma[0]) if sigma[0] else None,
            'V_over_sigma_phi': float(abs(mean[1])/sigma[1]) if sigma[1] else None,
            'V_over_sigma_1d_cyl': float(abs(mean[1])/sigma1d) if sigma1d else None,
            'prograde_mass_fraction': float(mass[vp > 0].sum()/mass.sum()),
            'z_rms_kpc': float(np.sqrt(average(x[:, 2]**2, mass)))}


def epsilon_summary(epsilon, mass, use):
    good = use & np.isfinite(epsilon)
    denominator = mass[use].sum()
    return {
        'epsilon10_50_90': quantiles(epsilon[good], mass[good]) if good.any() else None,
        'mass_fraction_epsilon_gt_0p7': float(mass[use & (epsilon > .7)].sum()/denominator),
        'mass_fraction_epsilon_gt_0p8': float(mass[use & (epsilon > .8)].sum()/denominator),
        'mass_fraction_epsilon_gt_0p9': float(mass[use & (epsilon > .9)].sum()/denominator),
        'mass_fraction_without_bound_circular_reference': float(mass[use & ~np.isfinite(epsilon)].sum()/denominator),
        'mass_fraction_epsilon_gt_1p02': float(mass[use & (epsilon > 1.02)].sum()/denominator),
    }


def disk_summary(population, directory, center, bulk):
    import agama
    from scipy.integrate import quad

    directory = Path(directory)
    agama.setUnits(length=1, velocity=1, mass=1e10)
    config = configparser.RawConfigParser()
    config.read(directory/'parameters.ini')
    stellar_density = agama.Density(**dict(config.items('Density_stellar')))
    potential = agama.Potential(*[
        agama.Potential(file=str(directory/('model_potential_'+name)))
        for name in ('gas', 'stellar')])
    x, v, mass = population
    x, v = x-center, v-bulk
    R = np.hypot(x[:, 0], x[:, 1])
    core = np.sum(x*x, axis=1) < 25.
    main = core & (R >= .25) & (R < 1.5)
    epsilon, energy = circularity(potential, x, v)
    main_stats = kinematics(x, v, mass, main)
    main_epsilon = epsilon_summary(epsilon, mass, main)

    def vc2(radii):
        points = np.column_stack((radii, np.zeros_like(radii), np.zeros_like(radii)))
        return -radii*potential.force(points)[:, 0]

    edges = np.array([.15, .25, .4, .6, .85, 1.2, 1.7, 2.5])
    radii = np.sqrt(edges[:-1]*edges[1:])
    circular2 = vc2(radii)
    step = radii*.01
    derivative = (vc2(radii+step)-vc2(radii-step))/(2*step)
    kappa2 = derivative/radii+2*circular2/radii**2
    annuli = []
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        use = core & (R >= lo) & (R < hi)
        stats = kinematics(x, v, mass, use)
        if stats is None:
            continue
        sigma_target = 2*quad(
            lambda z: stellar_density.density([float(radii[i]), 0., z]),
            0., 5., epsabs=1e-10, epsrel=1e-5)[0]
        sigma_particles = mass[use].sum()/(np.pi*(hi*hi-lo*lo))
        epicycle = float(np.sqrt(max(kappa2[i], 0.)))
        stats.update(
            Rmin_kpc=float(lo), Rmax_kpc=float(hi), R_reference_kpc=float(radii[i]),
            circular_speed_kms=float(np.sqrt(max(circular2[i], 0.))),
            epicycle_frequency_kms_kpc=epicycle,
            Sigma_target_Msun_kpc2=float(sigma_target*1e10),
            Sigma_particles_annular_Msun_kpc2=float(sigma_particles*1e10),
            Q_star_target_Sigma=float(toomre_q(
                stats['sigma_R_phi_z_kms'][0], epicycle, sigma_target, agama.G)),
            Q_star_particle_Sigma=float(toomre_q(
                stats['sigma_R_phi_z_kms'][0], epicycle, sigma_particles, agama.G)),
            circularity=epsilon_summary(epsilon, mass, use))
        annuli.append(stats)
    dense_annuli = [a for a in annuli if a['Rmin_kpc'] >= .25 and a['Rmax_kpc'] <= 1.2]
    qmin = min(a['Q_star_target_Sigma'] for a in dense_annuli)
    criteria = {
        'main_V_over_sigma_R_ge_2': main_stats['V_over_sigma_R'] >= 2.,
        'main_V_over_sigma_1d_ge_2': main_stats['V_over_sigma_1d_cyl'] >= 2.,
        'main_prograde_ge_0p95': main_stats['prograde_mass_fraction'] >= .95,
        'main_circularity_median_ge_0p8': main_epsilon['epsilon10_50_90'][1] >= .8,
        'initial_Q_thin_target_floor_ge_1p2': qmin >= 1.2,
    }
    return {
        'reference': 'Initial gas+stellar smooth axisymmetric potential; native +Z symmetry axis.',
        'membership': 'Current particles within 3D5; no ID or boundness gate.',
        'main_region': '0.25 <= cylindrical R < 1.5 kpc, within 3D5.',
        'main_mass_fraction_of_stars_3D5': float(mass[main].sum()/mass[core].sum()),
        'main': main_stats, 'main_circularity': main_epsilon,
        'all_stars_3D5_circularity': epsilon_summary(epsilon, mass, core),
        'mass_fraction_positive_energy_in_initial_reference': float(mass[energy >= 0].sum()/mass.sum()),
        'annuli': annuli, 'min_Q_target_in_R025_12': qmin,
        'declared_exploratory_criteria': criteria,
        'all_declared_indicators_met': all(criteria.values()),
        'Q_definition': 'Q_star=sigma_R*kappa_epicycle/(3.36*G*Sigma_star); total initial potential sets kappa.',
        'Q_limit': 'Approximate finite-bin, single-component, local thin-disk indicator. Does not prove '
                   'coupled gas/star stability, bending or global bar stability.',
        'cold_limit': 'Finite-dispersion rotation-dominated candidate, not sigma=0 or long-term validation.',
    }
