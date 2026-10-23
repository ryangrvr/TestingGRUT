#!/usr/bin/env python3
"""Post-scan witnesses at the ORIGINAL supplied net-charge yield.

Test charges remain supplied. This calibrates v to total late rest mass;
it does not predict a charge distribution, v, or a primordial state.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from widths import ROOT, Widths, population
from profiles import solve, save
from scan import Object, detector, row_for
from projection import original

YQ_OLD = 4.8294009076988175e-9
BDM = 4.4e-10  # GeV, the parent's explicit density calibration


class CSVSeed:
    def __init__(self, row):
        a = np.loadtxt(ROOT/row['profile_file'], delimiter=',', skiprows=1)
        self.x, self.p = a[:, 0], np.array([row['Omega']])
        j = cumulative_trapezoid(self.x**2*a[:, 1]**2, self.x, initial=0.)
        self.y = np.vstack([a[:, 1:].T, j])

    def sol(self, r):
        return np.array([np.interp(r, self.x, y) for y in self.y])


def load_toys(cache):
    path = Path(cache)/'xe_likelihood/data/xent_nr_20ty_projection/XENONnT_binwise_NR_LL_toyMC.json.gz'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == '9874ff4f0d5b2e3e0c48e497d1e685195baf7c16211645abf6e6dd921ab02885'
    with gzip.open(path, 'rt') as f:
        data = json.load(f)
    axis = data['Yedges']
    edges = np.round(np.linspace(axis['low_edge'], axis['high_edge'], axis['N']), 6)
    toys = np.array([data['Z']['Data'][str(i)] for i in range(len(data['Z']['Data']))])
    return (edges[1:]+edges[:-1])/2, toys, digest


def onset(v):
    """Controlled dilute NR limit; not nonlinear formation or final Q."""
    lam, g, lx, temp, degrees = .03, .1, .001, .02, 10.75
    entropy = 2*np.pi**2*degrees*temp**3/45
    n = YQ_OLD*entropy
    attraction = g**4/lam-lx
    f2 = n/(g*v)  # dilute approximation; quantified below
    k2 = attraction*f2
    gamma = k2/(2*g*v)
    hubble = np.sqrt(np.pi**2*degrees/90)*temp**2/2.435e18
    return dict(scope='Leading dilute linear onset only, at a supplied coherent state and 20 MeV visible bath.',
                attractive_quartic=attraction, charge_density_GeV3=n,
                f_squared_GeV2=f2, Gamma_over_H=gamma/hubble,
                Gamma_inverse_seconds=6.582119569e-25/gamma,
                k_peak_GeV=np.sqrt(k2),
                charge_per_peak_wavelength_cube_DIAGNOSTIC=n*(2*np.pi/np.sqrt(k2))**3,
                dilute_fractional_mass_correction_bound=k2/(g*v)**2,
                adiabatic_mediator_ratio=k2/(2*lam*v*v))


def main(cache):
    profile_record = json.loads((ROOT/'PROFILE_RESULTS.json').read_text())
    profiles = profile_record['profiles']
    start = next(x for x in profiles if x['Q'] == 1000000 and x['lambda_phi'] == .03)
    seed = CSVSeed(start)
    for q in np.geomspace(1e6, 2e5, 15)[1:]:
        seed = solve(float(q), .03, seed, length=100.)
    row = save(seed, 200000, .03)
    refined = solve(200000, .03, seed, length=120., tolerance=2e-9)
    refinement = save(refined, 200000, .03, True)
    refinement['tolerance'] = 2e-9
    ys, toys, digest = load_toys(cache)
    widths = Widths()
    old_q = next(x for x in profiles if x['Q'] == 39944 and x['lambda_phi'] == .021)
    results = []
    for profile in (old_q, row):
        q, lam = profile['Q'], profile['lambda_phi']
        obj = Object(profile)
        def budget(v):
            pop = population(q, obj.mass(v), v, lam, 1e-9,
                             widths.unit(np.sqrt(2*lam)*v))
            return YQ_OLD/q*pop['late_mass_per_object_GeV']/BDM-1
        v = brentq(budget, .5, 1.5, xtol=1e-14)
        r = row_for(obj, v, widths, shape=True)
        _, _, signals = detector(obj, v, r['object_halo_fraction'])
        strength, curves = original.curve_grid(signals, ys, toys, 1200)
        ul = original.limits(strength, curves, sum(signals))
        pct = np.quantile(1e-9*np.sqrt(ul), [.1586552539, .5, .8413447461])
        upper = min(r['kappa_LZ_same_mass_count_ESTIMATE'],
                    r['kappa_LZ_shape_matched_ESTIMATE'] or np.inf)
        lower = max(pct[1], r['kappa_for_lifetime_1s_central'])
        kval = float(np.sqrt(lower*upper)) if lower < upper else None
        pop = None if kval is None else population(q, obj.mass(v), v, lam,
                                                  kval, widths.unit(r['m_s_GeV']))
        r.update(original_fixed_charge_yield=YQ_OLD,
                 scale_status='v calibrated to late matter at fixed original YQ; not predicted.',
                 fixed_yield_matter_relative_residual=budget(v),
                 kappa_projected90_percentiles=pct.tolist(),
                 approximate_LZ_screen_upper=upper,
                 prompt_and_median_projection_screen_overlap=bool(lower < upper),
                 overlap_lower_kappa=float(lower), overlap_upper_kappa=upper,
                 joint_screen_kappa=kval, joint_screen_population=pop,
                 joint_screen_expected_signals=None if kval is None else float(sum(signals)*(kval/1e-9)**2))
        if lam == .03:
            r['dilute_linear_onset'] = onset(v)
        results.append(r)
        print('fixed-yield witness', q, lam, 'v', v, 'overlap', lower < upper,
              'kappa', kval, 'binding', 1-profile['energy_per_free_mass'], flush=True)
    output = dict(scope='Post-scan fixed-original-YQ witnesses; test Q and calibrated v, no population genesis or safe cosmology claimed.',
                  original_charge_yield=YQ_OLD, calibrated_total_DM_GeV_per_entropy=BDM,
                  toy_source_sha256=digest, public_background_toys=len(toys),
                  new_profile=row, refinement=refinement, rows=results)
    (ROOT/'FIXED_YIELD_RESULTS.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cache', type=Path)
    main(parser.parse_args().cache)
