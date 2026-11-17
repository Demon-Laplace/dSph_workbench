"""Conservative, spherical baryonic Jeans conditional reference.

Only NumPy is required.  Lengths are kpc, masses Msun, velocities km/s.
The primary model fixes beta=0; passing beta changes the declared condition,
not the measurements.  The active equal-star pipeline supplies stellar counts as the tracer.
Physical stellar and gas masses remain the gravitational source.

The supplied histogram is interpreted as uniform density in each spherical
volume bin.  Cumulative stellar and gas mass are interpolated linearly in r^3
inside the same bins.  Thus M(<r)=A+B*r^3 in a bin, and Jeans pressure is solved
analytically there, rather than by a midpoint circular-speed approximation.
The gas input must represent the declared full-phase baryonic source sample;
this module cannot decide particle membership or host/background attribution.

Outer pressure means nu*sigma_r^2 at the supplied *outer* radius, in
Msun/kpc^3*(km/s)^2.  Zero is a finite outer boundary condition, never assumed
at the last plotted 2-kpc radius.  Call boundary_sensitivity for this choice.

Empty tracer bins are kept empty (no density floor).  Formal pressure can be
nonzero across a gap, but there is no finite point dispersion there.  Such
points are NaN and empty bins make no tracer-weighted prediction contribution.
This histogram-support limitation is explicitly flagged.
"""

from __future__ import annotations

import numpy as np


# Actual production GIZMO allvars.h G=6.672e-8 (cgs), converted with this
# project's adopted kpc=3.08568e21 cm and Msun=1.989e33 g.
G_KPC_KMS2_MSUN = 4.300710573683437e-6
METHOD = "conservative_piecewise_constant_density_and_volume_linear_mass"


def _float_array(value, name):
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or np.any(~np.isfinite(result)):
        raise ValueError(name + " must be a finite one-dimensional array")
    return result


def _power_difference(log_ratio, exponent):
    """(exp(exponent*log_ratio)-1)/exponent, including log limit."""
    if abs(exponent) < 1e-10:
        return log_ratio
    return np.expm1(exponent * log_ratio) / exponent


def _pressure_inside(model, indices, radius):
    """Analytic inward Jeans solution from each bin's outer pressure."""
    indices = np.asarray(indices, dtype=int)
    radius = np.asarray(radius, dtype=float)
    outer = model["edges_kpc"][indices + 1]
    beta = model["beta"][indices]
    exponent = 2.0 * beta
    density = model["stellar_density_Msun_kpc3"][indices]
    aa = model["mass_A_Msun"][indices]
    bb = model["mass_B_Msun_kpc3"][indices]
    log_ratio = np.log(outer / radius)
    answer = model["pressure_at_edges"][indices + 1] * np.exp(exponent * log_ratio)
    # Group equal beta values: normally this loop has exactly one member.
    forcing = np.zeros_like(radius)
    for value in np.unique(exponent):
        selected = exponent == value
        forcing[selected] = density[selected] * model["G"] * (
            aa[selected] / radius[selected]
            * _power_difference(log_ratio[selected], value - 1.0)
            + bb[selected] * radius[selected] ** 2
            * _power_difference(log_ratio[selected], value + 2.0)
        )
    return answer + forcing


def solve_jeans(edges_kpc, stellar_bin_mass_Msun,
                gas_cumulative_at_edges_Msun, G=G_KPC_KMS2_MSUN,
                beta=None, outer_pressure=0.0, quadrature_order=24, tracer_bin_counts=None):
    """Return a model dictionary accepted by the evaluation helpers.

    ``edges_kpc`` must start at exactly zero and extend beyond 2 kpc.
    ``stellar_bin_mass_Msun`` has len(edges)-1 entries; gas cumulative has
    len(edges), starts at zero, and is nondecreasing.  beta is None (zero),
    a scalar, or a per-bin array.  It is piecewise constant and must be <1.
    No beta clipping, fitting or matching to measured dispersion occurs.

    Result contains original inputs, density and cumulative masses, exact
    pressure at grid edges, bin-centre evaluation, and metadata/quality flags.
    Predictions are computed by shell_predictions and
    projected_aperture_predictions; neither replaces a shell by its midpoint.
    """
    edges = _float_array(edges_kpc, "edges_kpc")
    masses = _float_array(stellar_bin_mass_Msun, "stellar_bin_mass_Msun")
    gas = _float_array(gas_cumulative_at_edges_Msun,
                       "gas_cumulative_at_edges_Msun")
    if len(edges) < 2 or edges[0] != 0 or np.any(np.diff(edges) <= 0):
        raise ValueError("edges must strictly increase from zero")
    if len(masses) != len(edges) - 1 or len(gas) != len(edges):
        raise ValueError("mass array lengths do not match radial edges")
    if np.any(masses < 0) or gas[0] != 0 or np.any(np.diff(gas) < 0):
        raise ValueError("negative mass, nonzero central mass, or decreasing gas mass")
    if not np.isfinite(G) or G <= 0 or not np.isfinite(outer_pressure) or outer_pressure < 0:
        raise ValueError("G must be positive and outer pressure nonnegative")
    if beta is None:
        beta_array = np.zeros_like(masses)
    elif np.ndim(beta) == 0:
        beta_array = np.full_like(masses, float(beta))
    else:
        beta_array = _float_array(beta, "beta")
    if beta_array.shape != masses.shape or np.any(~np.isfinite(beta_array)) or np.any(beta_array >= 1):
        raise ValueError("beta must be finite, <1, and scalar or one value per bin")
    if not isinstance(quadrature_order, (int, np.integer)) or quadrature_order < 4:
        raise ValueError("quadrature_order must be an integer >=4")
    volume = 4 * np.pi / 3 * np.diff(edges ** 3)
    tracer = masses if tracer_bin_counts is None else _float_array(tracer_bin_counts, "tracer_bin_counts")
    if tracer.shape != masses.shape or np.any(tracer < 0):
        raise ValueError("invalid tracer counts")
    density = tracer / volume
    star_cumulative = np.r_[0.0, np.cumsum(masses)]
    baryonic = star_cumulative + gas
    bb = np.diff(baryonic) / np.diff(edges ** 3)
    aa = baryonic[:-1] - bb * edges[:-1] ** 3
    flags = []
    positive = np.flatnonzero(masses > 0)
    empty = np.flatnonzero(masses == 0)
    if len(positive) == 0:
        flags.append("no_tracer_mass")
    elif np.any(masses[:positive[-1] + 1] == 0):
        flags.append("empty_density_gaps_no_finite_point_dispersion")
    if len(empty):
        flags.append("zero_density_bins_preserved_no_positive_floor")
    if edges[-1] <= 2.0:
        flags.append("outer_radius_not_extended_beyond_main_survey")
    if beta_array[0] > 0:
        flags.append("radial_anisotropy_constant_density_core_central_pressure_diverges")
    model = dict(
        edges_kpc=edges.copy(), stellar_bin_mass_Msun=masses.copy(),
        gas_cumulative_at_edges_Msun=gas.copy(),
        stellar_cumulative_at_edges_Msun=star_cumulative,
        baryonic_cumulative_at_edges_Msun=baryonic,
        stellar_density_Msun_kpc3=density, tracer_bin_counts=tracer_bin_counts, beta=beta_array,
        mass_A_Msun=aa, mass_B_Msun_kpc3=bb,
        pressure_at_edges=np.zeros_like(edges), G=float(G),
        outer_pressure=float(outer_pressure),
        quadrature_order=int(quadrature_order), quality_flags=flags,
        metadata=dict(method=METHOD,
                      tracer="stellar number density" if tracer_bin_counts is not None else "stellar mass density",
                      beta_condition="fixed_beta_zero" if np.all(beta_array == 0)
                      else "supplied_piecewise_constant_beta_condition",
                      outer_radius_kpc=float(edges[-1]),
                      outer_pressure_Msun_kpc3_kms2=float(outer_pressure),
                      interpolation="uniform density per volume bin; M linear in r^3",
                      pressure_solver="analytic inward solution in each radial bin",
                      empty_bins="density zero; sigma NaN; no tracer prediction contribution",
                      gravitational_kernel="unsoftened baryonic Newtonian monopole",
                      units="kpc, Msun, km/s; tracer pressure count/kpc^3*(km/s)^2" if tracer_bin_counts is not None else "kpc, Msun, km/s; pressure Msun/kpc^3*(km/s)^2"),
    )
    model["pressure_at_edges"][-1] = outer_pressure
    for index in range(len(masses) - 1, 0, -1):
        model["pressure_at_edges"][index] = _pressure_inside(
            model, np.array([index]), np.array([edges[index]]))[0]
    # At r=0 the isotropic core has A=0 and the pressure is regular.
    if beta_array[0] == 0:
        model["pressure_at_edges"][0] = (
            model["pressure_at_edges"][1] + density[0] * G * bb[0] * edges[1] ** 2 / 2)
    elif beta_array[0] > 0:
        model["pressure_at_edges"][0] = np.inf
    else:
        model["pressure_at_edges"][0] = 0.0
    mid = ((edges[:-1] ** 3 + edges[1:] ** 3) / 2) ** (1 / 3)
    model["bin_centres"] = evaluate_jeans(model, mid)
    return model


def evaluate_jeans(model, radii_kpc):
    """Evaluate exact within-bin density, M, g, pressure and Jeans variance.

    Outside the supplied radial coverage or in empty tracer bins, point
    dispersion is NaN.  Density zero is never reported as a cold population.
    At r=0 g is its regular monopole limit zero; optional beta>0 variance is
    divergent and flagged in the model metadata.
    """
    radius = np.atleast_1d(np.asarray(radii_kpc, dtype=float))
    shape = radius.shape
    flat = radius.ravel()
    edges = model["edges_kpc"]
    valid = np.isfinite(flat) & (flat >= 0) & (flat <= edges[-1])
    index = np.clip(np.searchsorted(edges, flat, side="right") - 1, 0, len(edges) - 2)
    out = {name: np.full_like(flat, np.nan) for name in [
        "density_Msun_kpc3", "Mbar_Msun", "Mstar_Msun", "Mgas_Msun",
        "g_sph_kms2_kpc", "pressure_Msun_kpc3_kms2", "sigma_r2_kms2",
        "s1d_equiv2_kms2", "beta"]}
    if np.any(valid):
        rr, ii = flat[valid], index[valid]
        frac = (rr ** 3 - edges[ii] ** 3) / (edges[ii + 1] ** 3 - edges[ii] ** 3)
        stellar = model["stellar_cumulative_at_edges_Msun"][ii] + frac * model["stellar_bin_mass_Msun"][ii]
        gas = model["gas_cumulative_at_edges_Msun"][ii] + frac * np.diff(model["gas_cumulative_at_edges_Msun"])[ii]
        density, beta = model["stellar_density_Msun_kpc3"][ii], model["beta"][ii]
        pressure = np.empty_like(rr)
        nonzero = rr > 0
        pressure[nonzero] = _pressure_inside(model, ii[nonzero], rr[nonzero])
        pressure[~nonzero] = model["pressure_at_edges"][0]
        grav = np.zeros_like(rr)
        grav[nonzero] = model["G"] * (stellar[nonzero] + gas[nonzero]) / rr[nonzero] ** 2
        dispersion = np.full_like(rr, np.nan)
        np.divide(pressure, density, out=dispersion, where=density > 0)
        values = [density, stellar + gas, stellar, gas, grav, pressure,
                  dispersion, (3 - 2 * beta) * dispersion / 3, beta]
        for key, value in zip(out, values):
            out[key][valid] = value
    answer = {key: value.reshape(shape) for key, value in out.items()}
    answer["r_kpc"] = radius
    answer["supported_tracer"] = np.isfinite(answer["sigma_r2_kms2"])
    return answer


def _quadrature_intervals(model, boundaries, transform_lower=False):
    cuts = np.unique(np.r_[model["edges_kpc"], boundaries])
    cuts = cuts[(cuts >= 0) & (cuts <= model["edges_kpc"][-1])]
    lo, hi = cuts[:-1], cuts[1:]
    nodes, weights = np.polynomial.legendre.leggauss(model["quadrature_order"])
    if transform_lower:
        upper = np.sqrt(hi - lo)
        uu = (nodes[None, :] + 1) * upper[:, None] / 2
        rr = lo[:, None] + uu ** 2
        ww = weights[None, :] * upper[:, None] * uu
    else:
        rr = (lo[:, None] + hi[:, None]) / 2 + (hi - lo)[:, None] * nodes[None, :] / 2
        ww = (hi - lo)[:, None] * weights[None, :] / 2
    return rr, ww, lo, hi


def shell_predictions(model, shell_edges_kpc=(0.0, 0.25, 0.5, 1.0, 2.0)):
    """Whole-shell s_bar^2 averaged over the supplied tracer density.

    Integrals split at all density/mass bin edges and all requested shell
    edges.  Gauss-Legendre integration is exact (to rounding) for the beta=0
    piecewise-volume model's shell numerator.  ``M_tracer_Msun`` is the
    conservative volume integral, not a midpoint density estimate.
    """
    shells = _float_array(shell_edges_kpc, "shell_edges_kpc")
    if np.any(np.diff(shells) <= 0) or shells[0] < 0 or shells[-1] > model["edges_kpc"][-1]:
        raise ValueError("shell bounds outside model or not increasing")
    rr, ww, lo, hi = _quadrature_intervals(model, shells)
    evaluated = evaluate_jeans(model, rr)
    support = evaluated["density_Msun_kpc3"] > 0
    mass_integral = np.sum(4 * np.pi * rr ** 2 * evaluated["density_Msun_kpc3"] * ww, axis=1)
    pressure_integral = np.sum(np.where(support,
        4 * np.pi * rr ** 2 * evaluated["pressure_Msun_kpc3_kms2"]
        * (3 - 2 * evaluated["beta"]) / 3 * ww, 0), axis=1)
    masses, variances, flags = [], [], []
    for inner, outer in zip(shells[:-1], shells[1:]):
        selected = (lo >= inner) & (hi <= outer)
        mass = float(np.sum(mass_integral[selected]))
        value = float(np.sum(pressure_integral[selected]) / mass) if mass > 0 else np.nan
        masses.append(mass); variances.append(value)
        flags.append("ok" if mass > 0 and np.isfinite(value) and value >= 0 else "no_valid_tracer_prediction")
    return dict(r_inner_kpc=shells[:-1], r_outer_kpc=shells[1:],
                M_tracer_Msun=np.array(masses), s_bar2_kms2=np.array(variances),
                sigma_bar_1d_equiv_kms=np.sqrt(variances), quality_flags=flags,
                model_quality_flags=list(model["quality_flags"]))


def _projected_mass(model, cylinder_radius, minimum_3d_radius=0.0):
    edges = model["edges_kpc"]
    lower = np.maximum(edges[:-1], minimum_3d_radius)
    upper = edges[1:]
    good = upper > lower
    volume_function = lambda r: 4 * np.pi / 3 * (r ** 3 - np.maximum(r ** 2 - cylinder_radius ** 2, 0) ** 1.5)
    return float(np.sum(model["stellar_density_Msun_kpc3"][good]
        * (volume_function(upper[good]) - volume_function(lower[good]))))


def projected_aperture_predictions(model, aperture_edges_kpc=(0, 0.25, 0.5, 1, 2),
                                   tail_radius_kpc=2.0):
    """LOS variance in circular projected annuli, including all r>2 tails.

    This is a uniform complete-sampling, stellar-mass-weighted reference.
    Other selection functions/physical tracers require a separate matching
    model or selection integrator; no observational completeness is assumed.
    For annulus a<R<b, angle integration is done analytically and radius
    integration uses Gauss-Legendre with r=lower+u^2 at split boundaries.
    The result exactly corresponds to the projected Jeans formula integrated
    over that aperture; it is not the full 3D shell covariance projection.
    """
    apertures = _float_array(aperture_edges_kpc, "aperture_edges_kpc")
    if apertures[0] < 0 or np.any(np.diff(apertures) <= 0) or apertures[-1] > model["edges_kpc"][-1]:
        raise ValueError("invalid projected aperture bounds")
    rr, ww, lo, hi = _quadrature_intervals(model, np.r_[apertures, tail_radius_kpc], True)
    evaluated = evaluate_jeans(model, rr)
    support = evaluated["density_Msun_kpc3"] > 0
    pressure_weight = np.where(support, 4 * np.pi * rr ** 2
        * evaluated["pressure_Msun_kpc3_kms2"] * ww, 0)
    beta = evaluated["beta"]
    masses, values, tails, numerator_tails, flags = [], [], [], [], []
    for inner, outer in zip(apertures[:-1], apertures[1:]):
        mu_hi = np.sqrt(np.maximum(1 - (inner / rr) ** 2, 0))
        mu_lo = np.sqrt(np.maximum(1 - (outer / rr) ** 2, 0))
        kernel = (1 - beta) * (mu_hi - mu_lo) + beta * (mu_hi ** 3 - mu_lo ** 3) / 3
        mass = _projected_mass(model, outer) - _projected_mass(model, inner)
        numerator = float(np.sum(pressure_weight * kernel))
        mass_tail = _projected_mass(model, outer, tail_radius_kpc) - _projected_mass(model, inner, tail_radius_kpc)
        numerator_tail = float(np.sum(np.where(rr >= tail_radius_kpc, pressure_weight * kernel, 0)))
        value = numerator / mass if mass > 0 else np.nan
        masses.append(mass); values.append(value)
        tails.append(mass_tail / mass if mass > 0 else np.nan)
        numerator_tails.append(numerator_tail / numerator if numerator > 0 else np.nan)
        flags.append("ok" if mass > 0 and np.isfinite(value) and value >= 0 else "no_valid_projected_tracer_prediction")
    return dict(R_inner_kpc=apertures[:-1], R_outer_kpc=apertures[1:],
                M_projected_tracer_Msun=np.array(masses),
                sigma_bar_los2_kms2=np.array(values),
                sigma_bar_los_kms=np.sqrt(values),
                tail_mass_fraction_r_ge_2=np.array(tails),
                tail_variance_numerator_fraction_r_ge_2=np.array(numerator_tails),
                tail_radius_kpc=float(tail_radius_kpc), quality_flags=flags,
                model_quality_flags=list(model["quality_flags"]))


def projected_profile(model, radii_kpc):
    """LOS surface density and point-R variance using nonsingular z integrals."""
    radius = _float_array(radii_kpc, "radii_kpc")
    values, surface = [], []
    nodes, weights = np.polynomial.legendre.leggauss(model["quadrature_order"])
    for projected in radius:
        if projected < 0 or projected >= model["edges_kpc"][-1]:
            surface.append(np.nan); values.append(np.nan); continue
        edges = model["edges_kpc"]
        radial_edges = np.r_[projected, edges[edges > projected]]
        z = np.sqrt(np.maximum(radial_edges ** 2 - projected ** 2, 0))
        zz = (z[:-1, None] + z[1:, None]) / 2 + (z[1:] - z[:-1])[:, None] * nodes[None, :] / 2
        ww = (z[1:] - z[:-1])[:, None] * weights[None, :] / 2
        rr = np.sqrt(projected ** 2 + zz ** 2)
        evaluated = evaluate_jeans(model, rr)
        support = evaluated["density_Msun_kpc3"] > 0
        denominator = float(np.sum(2 * evaluated["density_Msun_kpc3"] * ww))
        numerator = float(np.sum(np.where(support,
            2 * evaluated["pressure_Msun_kpc3_kms2"]
            * (1 - evaluated["beta"] * projected ** 2 / rr ** 2) * ww, 0)))
        surface.append(denominator); values.append(numerator / denominator if denominator > 0 else np.nan)
    return dict(R_kpc=radius, surface_density_Msun_kpc2=np.array(surface),
                sigma_bar_los2_kms2=np.array(values))


def regrid_model(model, new_edges_kpc, outer_pressure=0.0):
    """Conservatively rebin/truncate the declared density model for sensitivity.

    This changes the piecewise density approximation if bins are coarsened.
    Refining an existing uniform-volume model preserves it, and is a solver
    consistency check, not a substitute for particle-based fine/coarse grids.
    beta is transferred at volume midpoints; measured-beta tests should pass
    their independently measured grid instead when it is available.
    """
    edges = _float_array(new_edges_kpc, "new_edges_kpc")
    if edges[0] != 0 or edges[-1] > model["edges_kpc"][-1]:
        raise ValueError("new coverage must start at zero and stay within original")
    evaluated = evaluate_jeans(model, edges)
    mid = ((edges[:-1] ** 3 + edges[1:] ** 3) / 2) ** (1 / 3)
    return solve_jeans(edges, np.diff(evaluated["Mstar_Msun"]),
                       evaluated["Mgas_Msun"], G=model["G"],
                       beta=evaluate_jeans(model, mid)["beta"],
                       outer_pressure=outer_pressure,
                       quadrature_order=model["quadrature_order"])


def boundary_sensitivity(model, radii_kpc=(32, 64, 128),
                         shell_edges_kpc=(0, 0.25, 0.5, 1, 2),
                         aperture_edges_kpc=(0, 0.25, 0.5, 1, 2)):
    """Compare declared finite zero-pressure boundaries to the full model.

    Output preserves signed relative differences in variance.  The pressure
    inherited from the full model at each shorter radius is also tested: this
    must reproduce the full model's interior and is a pressure-closure QA.
    """
    base_shell = shell_predictions(model, shell_edges_kpc)["s_bar2_kms2"]
    base_los = projected_aperture_predictions(model, aperture_edges_kpc)["sigma_bar_los2_kms2"]
    answer = []
    for radius in radii_kpc:
        if radius > model["edges_kpc"][-1] or radius < max(shell_edges_kpc[-1], aperture_edges_kpc[-1]):
            continue
        new_edges = np.r_[model["edges_kpc"][model["edges_kpc"] < radius], float(radius)]
        zero = regrid_model(model, new_edges, 0.0)
        inherited_pressure = float(evaluate_jeans(model, [radius])["pressure_Msun_kpc3_kms2"][0])
        inherited = regrid_model(model, new_edges, inherited_pressure)
        sv = shell_predictions(zero, shell_edges_kpc)["s_bar2_kms2"]
        lv = projected_aperture_predictions(zero, aperture_edges_kpc)["sigma_bar_los2_kms2"]
        iv = shell_predictions(inherited, shell_edges_kpc)["s_bar2_kms2"]
        answer.append(dict(outer_radius_kpc=float(radius),
            zero_pressure_shell_variance=sv,
            shell_relative_change=(sv - base_shell) / base_shell,
            zero_pressure_los_variance=lv,
            los_relative_change=(lv - base_los) / base_los,
            inherited_outer_pressure=inherited_pressure,
            inherited_pressure_shell_relative_error=(iv - base_shell) / base_shell))
    return answer


def grid_sensitivity(model, alternative_edges, shell_edges_kpc=(0, 0.25, 0.5, 1, 2),
                     aperture_edges_kpc=(0, 0.25, 0.5, 1, 2)):
    """Compare conservative alternative histogram approximations to baseline.

    alternative_edges is a dict(label -> radial edge array).  The same outer
    coverage and pressure are required to isolate grid approximation changes.
    For the strongest validation, construct models directly from actual
    particle histograms on each grid, rather than rebinned baseline masses.
    """
    base_shell = shell_predictions(model, shell_edges_kpc)["s_bar2_kms2"]
    base_los = projected_aperture_predictions(model, aperture_edges_kpc)["sigma_bar_los2_kms2"]
    answer = []
    for label, edges in alternative_edges.items():
        if float(edges[-1]) != float(model["edges_kpc"][-1]):
            raise ValueError("grid sensitivity requires same outer radius")
        changed = regrid_model(model, edges, model["outer_pressure"])
        sv = shell_predictions(changed, shell_edges_kpc)["s_bar2_kms2"]
        lv = projected_aperture_predictions(changed, aperture_edges_kpc)["sigma_bar_los2_kms2"]
        answer.append(dict(grid_label=str(label), n_radial_bins=len(edges) - 1,
            shell_relative_change=(sv - base_shell) / base_shell,
            los_relative_change=(lv - base_los) / base_los,
            model_quality_flags=changed["quality_flags"]))
    return answer


def limited_linear_density_model(model, subdivisions=4):
    """Predeclared conservative positive linear-density sensitivity model.

    The primary uniform-bin tracer is replaced by a slope-limited linear
    density in each original bin.  Slopes use minmod of adjacent density
    differences and are reduced, only if necessary, to remain nonnegative
    at both bin edges.  The expansion is about the volume-weighted mean
    radius, so each original bin's mass is preserved.  Empty bins and their
    support remain exactly empty; no density floor or gap filling occurs.

    The linear density is integrated exactly into ``subdivisions`` bins of
    equal shell volume, then passed through the same analytic solver.
    Gas cumulative mass and beta retain their original interpolation and
    conditions.  This is an interpolation-sensitivity condition, not fitted
    to velocities.  It requires >=2 bins and never replaces primary results.
    """
    if not isinstance(subdivisions, (int, np.integer)) or subdivisions < 2:
        raise ValueError("subdivisions must be integer >=2")
    edges = model["edges_kpc"]
    lower, upper = edges[:-1], edges[1:]
    density = model["stellar_density_Msun_kpc3"]
    mean_radius = 0.75 * (upper ** 4 - lower ** 4) / (upper ** 3 - lower ** 3)
    slope = np.zeros_like(density)
    if len(density) >= 3:
        left = np.diff(density)[:-1] / np.diff(mean_radius)[:-1]
        right = np.diff(density)[1:] / np.diff(mean_radius)[1:]
        same = left * right > 0
        limited = np.where(same, np.sign(left) * np.minimum(np.abs(left), np.abs(right)), 0)
        slope[1:-1] = limited
    for bound in (lower, upper):
        increment = slope * (bound - mean_radius)
        too_low = increment < -density
        factor = np.ones_like(density)
        factor[too_low] = density[too_low] / (-increment[too_low])
        slope *= factor
    slope[density == 0] = 0
    fractions = np.linspace(0, 1, subdivisions + 1)
    subedges = (lower[:, None] ** 3 + (upper ** 3 - lower ** 3)[:, None] * fractions[None, :]) ** (1 / 3)
    # Preserve supplied boundaries bit-for-bit, including 0.25,0.5,1,2.
    subedges[:, 0], subedges[:, -1] = lower, upper
    ll, uu = subedges[:, :-1], subedges[:, 1:]
    dv3, dv4 = uu ** 3 - ll ** 3, uu ** 4 - ll ** 4
    submass = 4 * np.pi * (density[:, None] * dv3 / 3
        + slope[:, None] * (dv4 / 4 - mean_radius[:, None] * dv3 / 3))
    if np.any(submass < 0):
        raise FloatingPointError("limited linear density produced negative submass")
    subtotal = np.sum(submass, axis=1)
    correction = np.ones_like(subtotal)
    np.divide(model["stellar_bin_mass_Msun"], subtotal, out=correction, where=subtotal > 0)
    submass *= correction[:, None]
    new_edges = np.r_[subedges[:, :-1].ravel(), edges[-1]]
    gas = evaluate_jeans(model, new_edges)["Mgas_Msun"]
    changed = solve_jeans(new_edges, submass.ravel(), gas, G=model["G"],
                         beta=np.repeat(model["beta"], subdivisions),
                         outer_pressure=model["outer_pressure"],
                         quadrature_order=model["quadrature_order"])
    changed["metadata"]["density_variant"] = "mass_conservative_positive_minmod_linear_in_r"
    changed["metadata"]["density_variant_subdivisions"] = int(subdivisions)
    return changed


def density_interpolation_sensitivity(model, subdivisions=4,
                                     shell_edges_kpc=(0, 0.25, 0.5, 1, 2),
                                     aperture_edges_kpc=(0, 0.25, 0.5, 1, 2)):
    """Compare predeclared limited-linear density to primary uniform density."""
    changed = limited_linear_density_model(model, subdivisions)
    base_shell = shell_predictions(model, shell_edges_kpc)["s_bar2_kms2"]
    base_los = projected_aperture_predictions(model, aperture_edges_kpc)["sigma_bar_los2_kms2"]
    sv = shell_predictions(changed, shell_edges_kpc)["s_bar2_kms2"]
    lv = projected_aperture_predictions(changed, aperture_edges_kpc)["sigma_bar_los2_kms2"]
    return dict(density_variant=changed["metadata"]["density_variant"],
                subdivisions=int(subdivisions),
                shell_variance=sv, shell_relative_change=(sv - base_shell) / base_shell,
                los_variance=lv, los_relative_change=(lv - base_los) / base_los,
                mass_relative_error=(np.sum(changed["stellar_bin_mass_Msun"])
                                     / np.sum(model["stellar_bin_mass_Msun"]) - 1)
                if np.sum(model["stellar_bin_mass_Msun"]) > 0 else np.nan,
                model_quality_flags=changed["quality_flags"])
