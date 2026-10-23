#!/usr/bin/env python3
"""Same-author numerical controls, including analytic and alternative routes."""
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy import fft
from scipy.integrate import trapezoid, solve_ivp

from fragment import L, initial, periodic_components
from transport import (ROOT, kernel, calculate, H0, OMEGA_R, MS, MD, MX,
                       PHOTON_ENERGY_PER_ENTROPY0, scattering_sigma_v, MPL)

checks = []


def check(name, value, limit, passed=None):
    passed = bool(abs(value) <= limit if passed is None else passed)
    checks.append(dict(name=name, value=float(value), limit=float(limit), passed=passed))


def analytic_mode(mode, dt):
    # Actual split algorithm in a 3D cube, compared with the linear PDE's
    # growing eigenvector, not with the same stepper or a desired clump count.
    n = 32
    coord = np.arange(n)*L/n
    phase = sum(mode[i]*.5*x for i, x in enumerate(np.meshgrid(coord, coord, coord, indexing='ij')))
    ksq = .25*sum(i*i for i in mode)
    growth = np.sqrt(ksq*(2-ksq))
    epsilon = 1e-6
    cos = np.cos(phase)
    psi = (1+epsilon*(1+1j*growth/ksq)*cos).astype(np.complex128)
    kval = 2*np.pi*fft.fftfreq(n, L/n)
    k2 = kval[:, None, None]**2+kval[None, :, None]**2+kval[None, None, :]**2
    steps = round(1/dt)
    for _ in range(steps):
        psi *= np.exp(.5j*dt*np.abs(psi)**2)
        psi = fft.ifftn(fft.fftn(psi)*np.exp(-1j*dt*k2))
        psi *= np.exp(.5j*dt*np.abs(psi)**2)
    measured = np.mean(2*np.real(psi*np.exp(-1j)-1)*cos)
    return float(measured/(epsilon*np.exp(growth))-1)


def main():
    data = ROOT/'data'
    sims = [json.loads(p.read_text()) for p in sorted(data.glob('n*.json'))]
    if len(sims) != 5:
        raise RuntimeError(f'Expected five completed controlled runs, got {len(sims)}')
    for sim in sims:
        tag = f"n{sim['n']}_seed{sim['seed']}_dt{sim['maximum_dt']:g}"
        check(tag+' charge', sim['max_norm_drift'], 1e-10)
        check(tag+' Hamiltonian', sim['max_relative_hamiltonian_drift'], 1e-3)
        check(tag+' spectral stop reached', 1, 1, sim['reason'] == 'spectral_stop')
        check(tag+' last safe high-k gradient', sim['last_safe']['high_k_gradient_fraction'], 1e-3)
        check(tag+' last safe high-k norm', sim['last_safe']['high_k_norm_fraction'], 1e-6)
    comparisons = []
    base = np.load(data/'n96_seed17_dt0.02.npz')
    for filename in ('n64_seed17_dt0.02.npz', 'n128_seed17_dt0.02.npz', 'n96_seed17_dt0.01.npz'):
        other = np.load(data/filename)
        for tau in (0, 6):
            a = base[f'tau_{tau}_psi32']
            b = other[f'tau_{tau}_psi32']
            error = float(np.linalg.norm(a-b)/np.linalg.norm(a))
            density_error = float(np.linalg.norm(np.abs(a)**2-np.abs(b)**2)/np.linalg.norm(np.abs(a)**2))
            comparisons.append(dict(file=filename, tau=tau, field_relative_L2=error, density_relative_L2=density_error))
            check(filename+f' tau{tau} field', error, 2e-4)
            check(filename+f' tau{tau} density', density_error, 3e-4)
    # Initial states are the same continuous Fourier polynomial, not unrelated
    # random grids. Even-grid 96 and 64 share every second/third point.
    p64, p96 = initial(64, 17), initial(96, 17)
    check('same continuous initial field shared grid', np.max(np.abs(p64[::2, ::2, ::2]-p96[::3, ::3, ::3])), 2e-15)
    del p64, p96
    rho = np.ones((8, 8, 8))
    rho[0, 3, 3] = rho[-1, 3, 3] = 3
    comps = periodic_components(rho, 2, 1.)
    check('periodic face components join', len(comps), 1, len(comps) == 1 and comps[0]['cells'] == 2)
    mode_results = []
    for mode in ((1, 0, 0), (2, 0, 0), (2, 1, 0)):
        coarse = analytic_mode(mode, .01)
        fine = analytic_mode(mode, .005)
        check(str(mode)+' linear growth fine', fine, 2e-5)
        check(str(mode)+' time refinement reduces error', abs(fine), abs(coarse), abs(fine) < abs(coarse))
        mode_results.append(dict(mode=mode, relative_error_dt_001=coarse, relative_error_dt_0005=fine))
    coarse = json.loads((ROOT/'TRANSPORT_RESULTS.json').read_text())
    fine, _ = calculate(32000)
    refine = []
    for i, (a, b) in enumerate(zip(coarse['decay_kernels'], fine['decay_kernels'])):
        error = a['pre_recombination_visible_heat_proxy']/b['pre_recombination_visible_heat_proxy']-1
        refine.append(error)
        check(f'kernel{i} heat grid refinement', error, 2e-5)
        check(f'kernel{i} decay probability', a['probability_integral_residual'], 1e-5)
        check(f'kernel{i} particle-energy continuity', a['energy_continuity_residual'], 1e-5)
    # Radiation-only closed form, nonrelativistic particle decay, all z.
    budget, aform, tau = 4.29400907698818e-11, 8.3e-12, 1.23796497995e11
    result, curve = kernel(budget, aform, 1., tau, 32000, True)
    z, inj, _ = curve
    actual = trapezoid(inj, -np.log1p(z))
    expected = budget/PHOTON_ENERGY_PER_ENTROPY0*np.sqrt(2*H0*np.sqrt(OMEGA_R)*tau)*math.gamma(1.5)
    check('radiation-era analytic decay energy kernel', actual/expected-1, 2e-6)
    # Generic COM angular integration reproduces the hand-derived threshold.
    sdd = coarse['threshold']['sigma_v_ss_DD_GeV_minus2']
    check('ss DD threshold from full angle integral', scattering_sigma_v(0, MD, .0001)/sdd-1, 1e-13)
    for p in (.01, .1, .3, 1):
        a = scattering_sigma_v(p, MX, .01, False, 64)
        b = scattering_sigma_v(p, MX, .01, False, 128)
        check(f'ss XXstar angular refinement p={p}', a/b-1, 1e-11)
    # Different integration route for the stated irreversible cold-channel model.
    channel = coarse['threshold']['channel_only'][1]
    yi = coarse['decay_kernels'][0]['scalar_yield']
    temp, deg = .02, 10.75
    ratio = (2*np.pi**2*deg*temp**3/45)*sdd/(np.sqrt(np.pi**2*deg/90)*temp**2/MPL)
    sol = solve_ivp(lambda n, y: [-ratio*np.exp(-n)*yi*y[0]**2], (0, 30), [1.], rtol=1e-10, atol=1e-12)
    check('irreversible channel analytic vs ODE', sol.y[0, -1]/channel['scalar_number_fraction_remaining']-1, 1e-8)
    # Full-field nonlinear amplitude is nowhere approached in archived onsets.
    peak = max(s['final']['peak_density'] for s in sims)
    f2 = 1.8218345666575245e-12
    check('NR onset max f squared over v squared', peak*f2, 1e-8)
    result = dict(author='Same author as construction; not independent review.',
                  passed=sum(c['passed'] for c in checks), total=len(checks), checks=checks,
                  field_comparisons_at_common_saved_states=comparisons,
                  linear_mode_controls=mode_results,
                  maximum_heat_refinement_relative_error=max(abs(x) for x in refine),
                  simulation_summaries=[{k:v for k,v in s.items() if k != 'component_diagnostics'} for s in sims])
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks', 'simulation_summaries')}, indent=2))
    if result['passed'] != result['total']:
        print(json.dumps([c for c in checks if not c['passed']], indent=2))
        raise SystemExit(1)


if __name__ == '__main__':
    main()
