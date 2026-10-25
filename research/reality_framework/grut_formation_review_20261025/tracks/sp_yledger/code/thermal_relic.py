#!/usr/bin/env python3
"""sp_yledger (Claude, own code built on our gp_leftovers benchmark): symmetric relic per net charge when the
Y (and vector) debris energy thermalises in the hidden sector.  Item (b).

Reuses tracks/gp_leftovers/code (our earlier track): Bose-Einstein thermodynamics with a condensed phase,
4->2 relaxation, case H (homogeneous condensate) and case L (lumps + bb211d1 conversion).  The only change is
the hidden energy at the onset: E0 = rho_c (1 + R) + V0, where R = (deposited debris energy)/(X condensate
energy), both radiation-like after deposition, so the ratio at onset equals the ratio at deposition.  The
onset budget rho_c + V0 is the single X budget for remainder, binding pairs and quench debris (gp_critic); the
debris energy R rho_c is disjoint from it (E_Y,exc and E_sigma are separate pieces of the Hamiltonian; the
vector energy is taken out of the Y and sigma budgets).  No energy is added twice.
ymult: hidden charge density relative to Y_Q s during freeze-out (66 = dilution only after freeze-out).
Usage: python3 thermal_relic.py <point_ledger.json> <out.json>   (~5-10 min, 4 cores)
"""
import json, sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np
from scipy.optimize import brentq

GPL = Path(__file__).resolve().parents[2] / 'gp_leftovers' / 'code'
sys.dont_write_bytecode = True               # write nothing into the gp_leftovers track
sys.path.insert(0, str(GPL))
import equilibrium as EQ                      # noqa: E402  (our gp_leftovers code)
from common import YQ, BDM, MX, OMEGA_REF, E_PER_FREE_MASS, Q_REF, entropy, gstar_s, hubble  # noqa: E402


def run_inj(eps, case, R, ymult=1.0):
    tr = EQ.BUD['transition']
    rho_c, V0, n_c = tr['rho_c'], tr['V0_GeV4'], tr['n_c']
    m = MX
    nQ0 = eps * n_c
    Ye = YQ * ymult
    Tv0 = brentq(lambda T: Ye * entropy(T) - nQ0, 1e-3, 1e4)
    E0 = rho_c * (1 + R) + V0
    mu_c, eps_c = (m, m) if case == 'H' else (OMEGA_REF, E_PER_FREE_MASS * m)
    T0, mu0, nc0, th0 = EQ.solve_state(nQ0, mu_c, energy=E0, eps_c=eps_c)
    sigma = th0['s'] / nQ0
    Ym = th0['nmi'] / nQ0
    lna_prev, frozen_at, hist = None, None, []
    for Tv in np.geomspace(Tv0, 1e-6, 260):
        nQ = Ye * entropy(Tv)
        lna = -np.log(gstar_s(Tv)**(1 / 3) * Tv)
        T, mu, nc, th = EQ.solve_state(nQ, mu_c, sigma=sigma)
        p_eq = th['nmi'] / nQ
        nmi = Ym * nQ
        npl = nmi + nQ - nc if case == 'L' else nmi + nQ
        G = EQ.A['a'] * npl**2 * nmi / 4 + EQ.A['b'] * npl**3 / 6 + EQ.A['b'] * npl * nmi**2 / 6
        if case == 'L':
            G += nc / Q_REF * EQ.SV_CONV
        Hh = hubble(Tv)
        if lna_prev is not None:
            x = G / Hh * (lna - lna_prev)
            Ym = p_eq + (Ym - p_eq) * np.exp(-min(x, 700))
        lna_prev = lna
        if G / Hh < 1 and (not hist or hist[-1]['G_over_H'] >= 1):
            frozen_at = Tv
        hist.append(dict(T_vis=Tv, T_h_over_m=T / m, cond_frac=nc / nQ, p_eq=p_eq, p=Ym, G_over_H=G / Hh))
        if case == 'H' and frozen_at is not None and Tv < frozen_at / 30:
            break
    f = 2 * Ym * MX * YQ / BDM          # after the yield is brought to Y_Q
    return dict(eps=eps, case=case, R=R, ymult=ymult, T_vis_onset=Tv0, E_per_charge_over_mX=E0 / nQ0 / m,
                T_h0_over_m=T0 / m, cond_frac0=nc0 / nQ0, sigma=sigma, p_final=Ym, f_extra=f,
                f_over_1pct=f / 0.01, T_vis_freeze=frozen_at, max_cond_frac=max(h['cond_frac'] for h in hist),
                final_cond_frac=hist[-1]['cond_frac'])


def job(a):
    try:
        return run_inj(*a)
    except Exception as exc:
        return dict(args=list(a), status='solver failure %r' % exc)


def main():
    L = json.load(open(sys.argv[1]))
    eps = abs(L['windows']['2000-3000']['eps'])
    e = L['energies_over_EX']
    Rs = dict(X_only=0.0, codex_Y=274.85203825803666,
              Y_deposited=e['Y_deposited'], Y_plus_A_low=e['total_injected_range'][0],
              Y_plus_A_high=e['total_injected_range'][1])
    jobs = [(eps, c, R, 1.0) for R in Rs.values() for c in ('H', 'L')]
    jobs += [(eps, 'H', Rs['Y_plus_A_low'], 65.98), (eps, 'L', Rs['Y_plus_A_low'], 65.98),
             (eps, 'H', 0.0, 65.98)]
    with Pool(4) as p:
        rows = p.map(job, jobs)
    names = {v: k for k, v in Rs.items()}
    for r in rows:
        if 'R' in r:
            r['scenario'] = names.get(r['R'], '?')
    out = dict(eps_used=eps, R_scenarios=Rs, rows=rows,
               scope='Maximal-equilibrium benchmark of gp_leftovers (full thermalisation, 4->2, conversion on lumps). '
                     'It is a benchmark, not a kinetic simulation. f_extra is the symmetric X/Xbar mass relative '
                     'to calibrated dark matter once the net yield equals Y_Q.')
    json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)
    for r in rows:
        print({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == '__main__':
    main()
