#!/usr/bin/env python3
"""Item 4: is the frozen-snapshot eigen-split a valid bound/free charge separation? Own code.

Uses an exact stationary soliton of the packet's own two-field action (soliton.py), my
evolution code (evolve.py) and my re-implementation of the packet's frozen split.
 A  exact soliton: split is exact (sanity);
 B  soliton + neutral S-velocity kick of size eps (adds NO charge): compare the split's
    'free' charge Q_phys - Q_bound(t) with the charge that actually leaves r < 30;
 C  soliton with in-phase amplitude change (1+eps): same comparison;
 D  rotating-frame coupled linear fluctuation spectrum (X and S, with the 0.2 f^2 and
    2 s f couplings) versus the frozen operator's eigenvalues.
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from soliton import solve_fixed_norm
from evolve import Evolver, frozen_bound_split
from lattice import L_matrix

OUT = Path(__file__).resolve().parents[1] / 'results'


def embed(sol, n, h):
    U = np.zeros((3, n + 1)); V = np.zeros((3, n + 1))
    m = len(sol['U'])
    U[0, 1:m + 1] = sol['U']
    U[2, 1:m + 1] = sol['Z']
    V[1, 1:m + 1] = -sol['w0'] * sol['U']            # X = f exp(-i w0 t)
    return U, V


def run(sol, kind, eps, T=1500., dt=.025, L=200., R=150., h=.25, every=2.0):
    n = int(round(L / h))
    ev = Evolver(n, h, R)
    U, V = embed(sol, n, h)
    r = ev.r
    if kind == 'neutral_kick':
        V[2] += eps * r * np.exp(-r ** 2 / (2 * 5.0 ** 2))
    elif kind == 'amplitude':
        U[0] *= 1 + eps; V[1] *= 1 + eps
    Y = [U, V, np.zeros((3, n)), np.zeros((3, n + 1))]
    Q0 = ev.Q(U, V)
    rows = []
    steps = int(round(T / dt)); rec = int(round(every / dt))
    for j in range(steps + 1):
        if j % rec == 0:
            U, V = Y[0], Y[1]
            d = frozen_bound_split(ev, U, V)
            rows.append(dict(t=j * dt, Q=ev.Q(U, V), Q30=ev.Q_within(U, V, 30.), Q60=ev.Q_within(U, V, 60.),
                             Qbound=d['Qbound'], nb=len(d['omegas']), omegas=d['omegas']))
        if j < steps:
            Y = ev.step(Y, dt)
    t = np.array([x['t'] for x in rows]); Q = np.array([x['Q'] for x in rows])
    Q30 = np.array([x['Q30'] for x in rows]); Qb = np.array([x['Qbound'] for x in rows])
    free_split = Q - Qb                      # charge the frozen split calls 'not bound'
    emitted = Q0 - Q30                       # charge that has actually left r < 30
    late = t >= T - 300
    return dict(kind=kind, eps=eps, Q0=float(Q0), T=T,
                split_free_t0=float(free_split[0]),
                split_free_max=float(free_split.max()), split_free_min=float(free_split.min()),
                split_free_late_mean=float(free_split[late].mean()),
                split_free_late_range=float(np.ptp(free_split[late])),
                emitted_from_r30_final=float(emitted[-1]), emitted_from_r30_late_mean=float(emitted[late].mean()),
                Q30_late_range=float(np.ptp(Q30[late])),
                n_bound_modes_values=sorted(set(int(x['nb']) for x in rows)),
                series_every_10=[dict(t=x['t'], Q30=x['Q30'], Qbound=x['Qbound'], nb=x['nb']) for x in rows[::5]])


def bdg_spectrum(sol, h, Rb):
    """Rotating frame, X = e^{-i w0 t}(f + a + i b), S = s + zeta, U-form (A = r a etc.):
       A'' + 2 w0 B' = L A + w0^2 A - (s^2 + .3 f^2) A - 2 s f Zt
       B'' - 2 w0 A' = L B + w0^2 B - (s^2 + .1 f^2) B
       Zt''          = L Zt - (9 s^2 - 3 + f^2) Zt - 2 s f A
    Normal modes ~ e^{-i Omega t}; lab frequencies of the charged part are w0 +- Omega."""
    n = int(round(Rb / h)); r = np.arange(1, n) * h; L = L_matrix(n, h); N = n - 1
    f = sol['U'] / r; s = 1 + sol['Z'] / r; w0 = sol['w0']
    I = np.eye(N); Z0 = np.zeros((N, N))
    K = np.block([[L + np.diag(w0 ** 2 - s * s - .3 * f * f), Z0, -np.diag(2 * s * f)],
                  [Z0, L + np.diag(w0 ** 2 - s * s - .1 * f * f), Z0],
                  [-np.diag(2 * s * f), Z0, L - np.diag(9 * s * s - 3 + f * f)]])
    Gy = np.block([[Z0, -2 * w0 * I, Z0], [2 * w0 * I, Z0, Z0], [Z0, Z0, Z0]])   # x'' = K x + Gy x'
    M = np.block([[np.zeros((3 * N, 3 * N)), np.eye(3 * N)], [K, Gy]])
    lam, vec = np.linalg.eig(M)
    Om = 1j * lam                                  # e^{lam t} = e^{-i Omega t}
    # localisation: fraction of mode norm within r < 20
    loc = []
    for j in range(len(lam)):
        v = vec[:3 * N, j]
        a2 = abs(v[:N]) ** 2 + abs(v[N:2 * N]) ** 2 + abs(v[2 * N:]) ** 2
        loc.append(float(np.sum(a2[r < 20]) / np.sum(a2)))
    loc = np.array(loc)
    sel = (abs(Om.imag) < 1e-6) & (Om.real > 1e-4) & (loc > .99)
    internal = sorted(set(np.round(Om.real[sel], 6)))
    growing = float(np.max(lam.real))
    return dict(w0=w0, charged_continuum_threshold_Omega=1 - w0,
                localized_real_modes_Omega=[float(x) for x in internal[:12]],
                max_growth_rate=growing)


def main():
    out = {}
    t0 = time.time()
    sol = solve_fixed_norm(1029.)
    out['soliton'] = dict(w0=sol['w0'], charge=sol['charge'], residual=sol['residual'], fmax=sol['fmax'], smin=sol['smin'])
    print(out['soliton'], flush=True)
    # A: exact soliton split at t=0 and after evolution
    out['A_exact_soliton'] = run(sol, 'none', 0.0, T=300.)
    print('A', {k: v for k, v in out['A_exact_soliton'].items() if k != 'series_every_10'}, time.time() - t0, flush=True)
    # frozen operator eigenvalues for the exact soliton
    n = 800; ev = Evolver(n, .25, 150.); U, V = embed(sol, n, .25)
    out['A_frozen_eigs'] = frozen_bound_split(ev, U, V)
    for kind, eps in [('neutral_kick', .02), ('neutral_kick', .05), ('neutral_kick', .1), ('amplitude', .02), ('amplitude', .05)]:
        key = f'{kind}_{eps}'
        out[key] = run(sol, kind, eps)
        print(key, {k: v for k, v in out[key].items() if k != 'series_every_10'}, time.time() - t0, flush=True)
    try:
        out['D_coupled_fluctuation_spectrum'] = bdg_spectrum(sol, .25, 40.)
        print(out['D_coupled_fluctuation_spectrum'], flush=True)
    except Exception as e:  # keep the other results
        out['D_coupled_fluctuation_spectrum'] = dict(error=str(e))
    (OUT / 'bound_split_counterexample.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
