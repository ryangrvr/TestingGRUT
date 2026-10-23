#!/usr/bin/env python3
"""Analytic coupling windows with declared detector and population screens."""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

from widths import ROOT, Widths, HBAR, hidden, visible_coefficient, population

PARENT = ROOT.parent/'grut_portal_empirical_test_20261022'
spec = importlib.util.spec_from_file_location('window_empirical_model', PARENT/'code/model.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
PUB = np.load(PARENT/'data/public_inputs.npz')
KREF = 1e-9


class Object(base.Object):
    def __init__(self, row):
        self.row = row
        self.charge = row['Q']
        self.lam = row['lambda_phi']
        data = np.loadtxt(ROOT/row['profile_file'], delimiter=',', skiprows=1)
        self.rho, self.s = data[:, 0], data[:, 3]
        self.weight = self.rho**2*(1-self.s**2)
        self.norm = simpson(self.weight, x=self.rho)
        qh = np.linspace(0, 10, 5001)
        vals = []
        for qs in np.array_split(qh, 50):
            vals.extend(simpson(self.weight[None, :]*np.sinc(qs[:, None]*self.rho[None, :]/np.pi),
                                x=self.rho, axis=1)/self.norm)
        self.form = CubicSpline(qh, vals, extrapolate=False)

    def mass(self, v):
        return self.row['energy_per_v_GeV']*v

    def sigma_n(self, v, kappa):
        mass = self.mass(v)
        mu = mass*base.MN/(mass+base.MN)
        c = kappa*base.FN*base.MN*mass/(2*self.lam*v*v*base.MH**2)
        return mu*mu*c*c/np.pi*base.GEV2CM2


def point_limit(mass):
    if not PUB['lz_si_x'][0] <= mass <= PUB['lz_si_x'][-1]:
        raise ValueError('No extrapolation of published LZ point limit')
    return float(np.exp(np.interp(np.log(mass), np.log(PUB['lz_si_x']), np.log(PUB['lz_si_y0']))))


def detector(obj, v, fraction, shape=False, points=1601):
    e = np.linspace(.001, 100, points)
    eff = np.interp(e, PUB['lz_efficiency_x'], PUB['lz_efficiency_y0'], left=0, right=0)
    r = fraction*obj.rate(e, v)*eff
    total = simpson(r, x=e)
    ref_mass = min(obj.mass(v), PUB['lz_si_x'][-1])
    ref = obj.rate(e, v, standard_mass=ref_mass, standard_sigma=point_limit(ref_mass))*eff
    count_upper = KREF*np.sqrt(simpson(ref, x=e)/total)
    true = 20*base.bin_integrals(lambda ee: fraction*obj.rate(ee, v), PUB['future_response_xedges'], 8)
    accepted = true@PUB['future_response_values']
    result = dict(kappa_LZ_same_mass_count_ESTIMATE=float(count_upper),
                  LZ_point_reference_mass_GeV=float(ref_mass),
                  LZ_count_at_point_limit_per_tonne_year=float(simpson(ref, x=e)),
                  xenonnt20ty_accepted_at_kappa_1e9=float(sum(accepted)),
                  kappa_three_accepted_signal_PROXY=float(KREF*np.sqrt(3/sum(accepted))))
    ma = 136*.9315
    mu = obj.mass(v)*ma/(obj.mass(v)+ma)
    result['Born_row_sum_bound_at_LZ_count_screen'] = float(mu*136*base.FN*base.MN*count_upper/(base.G**2*base.MH**2)
        *simpson(obj.rho*np.maximum(1-obj.s**2, 0), x=obj.rho))
    if shape:
        p = r/total

        def div(logm):
            mass = np.exp(logm)
            rp = obj.rate(e, v, standard_mass=mass, standard_sigma=1e-45)*eff
            norm = simpson(rp, x=e)
            if norm <= 0:
                return np.inf
            return float(.5*simpson((np.sqrt(p)-np.sqrt(rp/norm))**2, x=e))

        grid = np.linspace(np.log(3.), np.log(10000.), 100)
        scores = [div(x) for x in grid]
        j = int(np.argmin(scores))
        fit = minimize_scalar(div, bounds=(grid[max(0, j-1)], grid[min(len(grid)-1, j+1)]), method='bounded')
        fit_mass = float(np.exp(fit.x))
        result['shape_matched_point_mass_GeV'] = fit_mass
        result['shape_squared_Hellinger'] = float(fit.fun)
        if PUB['lz_si_x'][0] <= fit_mass <= PUB['lz_si_x'][-1]:
            rp = obj.rate(e, v, standard_mass=fit_mass, standard_sigma=point_limit(fit_mass))*eff
            result['kappa_LZ_shape_matched_ESTIMATE'] = float(KREF*np.sqrt(simpson(rp, x=e)/total))
        else:
            result['kappa_LZ_shape_matched_ESTIMATE'] = None
        result['shape_limit_scope'] = 'Named 2024 LZ heavy-WIMP table only; fit below 9 GeV needs the separate 2025 low-mass analysis.'
    return result, true, accepted


def row_for(obj, v, width, shape=False):
    lam = obj.lam
    ms = np.sqrt(2*lam)*v
    lo, central, hi = width.ensemble(ms)
    # In the open hidden-decay wedge this branching is effectively independent
    # of detector-sized kappa; the exact witness recalculates at its own kappa.
    pop = population(obj.charge, obj.mass(v), v, lam, KREF, central)
    det, _, _ = detector(obj, v, pop['object_halo_fraction'], shape)
    gx, gd = hidden(v, lam)
    coeff = visible_coefficient(v, lam, central)
    k_prompt = np.sqrt(max(HBAR-gx-gd, 0.)/coeff)
    k_eonly = np.sqrt(max(HBAR-gx-gd, 0.)/visible_coefficient(v, lam, width.lepton(ms)))
    k_fast = np.sqrt(max(HBAR-gx-gd, 0.)/visible_coefficient(v, lam, hi))
    k_slow = np.sqrt(max(HBAR-gx-gd, 0.)/visible_coefficient(v, lam, lo))
    k_max = det['kappa_LZ_same_mass_count_ESTIMATE']
    return dict(Q=obj.charge, lambda_phi=lam, v_phi_GeV=float(v), m_s_GeV=float(ms),
                mass_GeV=obj.mass(v), energy_per_free_mass=obj.row['energy_per_free_mass'],
                SM_unit_width_low_GeV=lo, SM_unit_width_central_GeV=central, SM_unit_width_high_GeV=hi,
                hidden_X_width_GeV=gx, hidden_D_width_GeV=gd,
                kappa_for_lifetime_1s_central=float(k_prompt),
                kappa_for_lifetime_1s_band_fast=float(k_fast),
                kappa_for_lifetime_1s_band_slow=float(k_slow),
                kappa_for_lifetime_01s_central=float(np.sqrt(max(10*HBAR-gx-gd, 0.)/coeff)),
                kappa_for_lifetime_1s_using_leptons_only=float(k_eonly),
                minimum_lifetime_at_LZ_count_screen_seconds=float(HBAR/(gx+gd+coeff*k_max*k_max)),
                prompt_and_count_reach_window_ESTIMATE=bool(max(k_prompt, det['kappa_three_accepted_signal_PROXY']) < k_max),
                **pop, **det)


def main():
    profiles = json.loads((ROOT/'PROFILE_RESULTS.json').read_text())
    widths = Widths()
    rows = []
    objects = {}
    for p in profiles['profiles']:
        if p['energy_per_free_mass'] >= 1:
            continue
        obj = Object(p)
        objects[p['Q'], p['lambda_phi']] = obj
        for v in np.unique(np.r_[np.geomspace(.5, 14, 61), 1., 2., 5., 7., 10., 14.]):
            if np.sqrt(2*obj.lam)*v > 2:
                continue
            rows.append(row_for(obj, float(v), widths))
        print('scanned Q', p['Q'], 'lambda', p['lambda_phi'], flush=True)
    selected = []
    for q, lam, v in [(39944,.01,1.),(39944,.01,7.),(39944,.01,10.),(39944,.01,14.),
                       (39944,.021,1.),(1000000,.03,1.),(1000000,.03,2.),(1000000,.03,5.)]:
        obj = objects[q, lam]
        result = row_for(obj, v, widths, True)
        shape_limit = result['kappa_LZ_shape_matched_ESTIMATE']
        upper = min(result['kappa_LZ_same_mass_count_ESTIMATE'], shape_limit or np.inf)
        # Pick a reviewable point below BOTH approximate LZ upper screens.
        kval = .8*upper
        central = widths.unit(result['m_s_GeV'])
        exact_pop = population(q, obj.mass(v), v, lam, kval, central)
        result['witness_kappa'] = kval
        result['witness_population'] = exact_pop
        result['witness_expected_XENONnT20ty_signals'] = result['xenonnt20ty_accepted_at_kappa_1e9']*(kval/KREF)**2
        result['witness_lifetime_seconds'] = exact_pop['scalar_lifetime_seconds']
        selected.append(result)
    original = next(r for r in selected if r['Q'] == 39944 and r['lambda_phi'] == .01 and r['v_phi_GeV'] == 1.)
    original['lifetime_at_kappa_1e10_seconds'] = HBAR/(visible_coefficient(1, .01, widths.unit(np.sqrt(.02)))*1e-20)
    original['kappa_for_1s_photon_C100_stress'] = float(np.sqrt(HBAR/visible_coefficient(1, .01, widths.unit(np.sqrt(.02), photon_coefficient=100))))
    result = dict(scope='Coupling/width/abundance screens, not complete cosmological viability or a formed population.',
                  parent_commit='ac96f9054b06a3d3c8cff1c1724b01e67fca1a4e',
                  population_assumptions='All net charge in identical test objects; full ground-state release as initially cold s; free prompt hidden daughters do not react or accrete; common halo fraction; charge yield derived for total cold matter, not predicted.',
                  lifetime_screen_seconds=1., signal_count_proxy=3,
                  kinematic_thresholds=dict(bulk_binding_lambda_upper=.1, XXstar_decay_lambda_lower=.02, DD_decay_lambda_lower=.0202),
                  total_grid_points=len(rows), selected=selected,
                  counts=dict(prompt_and_count_reach_window_points=sum(r['prompt_and_count_reach_window_ESTIMATE'] for r in rows),
                              hidden_decay_points=sum(r['hidden_X_width_GeV'] > 0 for r in rows)),
                  grid_file='data/scan.csv')
    (ROOT/'SCAN_RESULTS.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with (ROOT/'data/scan.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({k:v for k,v in result.items() if k != 'grid'}, indent=2))


if __name__ == '__main__':
    main()
