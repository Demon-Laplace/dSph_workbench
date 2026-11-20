"""Published Fornax Sersic-fit parameters, not newly measured radial bins.

Muñoz et al. 2018, ApJ 860, 66, doi:10.3847/1538-4357/aac16b,
Tables 1, 5, 6; Sections 3.1--3.2 and 4.4.4. The original Fornax footprint
does not reach the galaxy edge. Large-radius model continuation is extrapolation.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

CODE_VERSION = 'munoz2018-fornax-published-sersic-v1'
DOI = '10.3847/1538-4357/aac16b'
PUBLISHED = {
    'n': .71, 'n_error': .01,
    'theta_e_arcmin': 18.4, 'theta_e_error_arcmin': .2,
    'mu0_V': 23.59, 'mu0_V_error': .16,
    'mu_e_V': 24.77, 'mu_e_V_error': .16,
    'ellipticity': .28, 'ellipticity_error': .01,
    'position_angle_deg': 46., 'position_angle_error_deg': 1.,
    'r_e_at_paper_distance_pc': 787., 'r_e_error_pc': 9.,
    'paper_distance_kpc': 147.,
    'M_V': -13.46, 'M_V_error': .14,
    'log10_L_V_Lsun': 7.317, 'log10_L_V_error': .056,
}


def sersic_mu_v(radius_kpc, distance_kpc):
    """Evaluate the published rounded fit using its angular effective radius."""
    if not np.isfinite(distance_kpc) or distance_kpc <= 0:
        raise ValueError('Reference display distance must be positive and finite')
    radius = np.asarray(radius_kpc, dtype=float)
    theta_arcmin = radius/distance_kpc*180/np.pi*60
    bn = 1.999*PUBLISHED['n']-.327  # the approximation stated in the paper
    return PUBLISHED['mu0_V']+(2.5/np.log(10))*bn*(
        theta_arcmin/PUBLISHED['theta_e_arcmin'])**(1/PUBLISHED['n'])


def metadata(distance_kpc):
    return {
        'data_kind': 'published_sersic_fit_reconstruction',
        'new_independent_radial_measurements': 0,
        'doi': DOI, 'article_url': 'https://iopscience.iop.org/article/'+DOI,
        'parameter_tables': [5, 6], 'published_parameters': dict(PUBLISHED),
        'formula': 'mu_V=mu0_V+(2.5/ln(10))*(1.999*n-0.327)*(theta/theta_e)^(1/n)',
        'display_distance_kpc': float(distance_kpc),
        'display_r_e_kpc': float(np.deg2rad(PUBLISHED['theta_e_arcmin']/60)*distance_kpc),
        'display_r_e_error_kpc_at_fixed_distance': float(
            np.deg2rad(PUBLISHED['theta_e_error_arcmin']/60)*distance_kpc),
        'radius_definition': 'elliptical semi-major radius; angular theta_e is the published input',
        'uncertainty_shown': 'mu0 +/-0.16 mag normalization only; not a full profile confidence band',
        'parameter_covariance_used': False,
        'coverage_limit': 'Paper Section 4.4.4: approximately 0.6 deg2, does not reach the Fornax edge; '
                          'outer Sersic continuation is extrapolation.',
        'independent_absolute_calibration_of_existing_yang_profile': False,
        'calibration_note': 'The project Yang profile was already normalized using Muñoz summaries; '
                            'these are not independent absolute-photometry constraints.',
        'mu_e_definition': 'local surface brightness at r_e, not the mean inside r_e',
        'rounding_note': 'Reconstruction uses rounded table values and the paper b_n approximation; '
                         'the value at r_e agrees with tabulated mu_e to published rounding precision.',
    }


def export_curve(path, distance_kpc, radius_max=5.1, points=400):
    radii = np.linspace(0., radius_max, points)
    theta = radii/distance_kpc*180/np.pi*60
    mu = sersic_mu_v(radii, distance_kpc)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['radius_kpc', 'theta_arcmin', 'mu_V_published_fit',
                         'mu0_normalization_lower', 'mu0_normalization_upper', 'data_kind'])
        for r, t, value in zip(radii, theta, mu):
            writer.writerow([r, t, value, value-.16, value+.16, 'published_fit_evaluation_not_observed_bin'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--distance-kpc', type=float, default=139.4)
    args = parser.parse_args()
    export_curve(args.output, args.distance_kpc)
