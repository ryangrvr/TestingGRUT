#!/usr/bin/env python3
"""Elastic decoupling against source epoch (Claude review, track `conversion`), a sensitivity estimate.

Uses my own momentum-transfer cross section sigma_T(p) (results/elastic.json, frozen
Q = 200000 object).  For each escaped cohort the code finds a_dec, where n_obj sigma_T v/H = 1,
and the collisionless horizon from a_dec.  The diffusion length while coupled is
    lambda_D^2 = int 2 v^2 / (Gamma_T a^2) dt   (comoving).
The effective comoving momentum q' is the one whose collisionless horizon from the source
equals sqrt(lambda_dec^2 + lambda_D^2).  CLASS cannot hold the X tightly coupled to the
objects, so q' only mimics the shortened streaming.  The same treatment is applied to the
original positive daughters.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transport_mc as T

C_CM = 2.99792458e10
E = json.loads((T.RES / 'elastic.json').read_text())['rows']
LP = np.log([r['p_over_m'] for r in E]); LS = np.log([r['sigma_T_cm2'] for r in E])


def sigT(p):
    return np.exp(np.interp(np.log(np.clip(p, np.exp(LP[0]), np.exp(LP[-1]))), LP, LS))


def gam(a, q):
    p = q / a; v = p / np.sqrt(1 + p * p)
    return T.NOBJ0 / a ** 3 * sigT(p) * v * C_CM / T.hubble(a)


def cohort(a_start, q):
    """Segments with Gamma_T/H >= 1 are diffusive (random walk), those with < 1 free-streaming.
    lambda_eff^2 = sum(free segments' path^2) + sum(diffusive lambda_D^2); the last free
    segment runs to today."""
    g0 = gam(a_start, q)
    ys = np.linspace(np.log(a_start), 0.0, 4000)
    g = np.array([gam(np.exp(y), q) for y in ys])
    coupled = g >= 1
    lam2 = 0.0; lamD2 = 0.0; last_dec = a_start
    i = 0
    while i < len(ys) - 1:
        j = i
        while j < len(ys) - 1 and coupled[j + 1] == coupled[i]:
            j += 1
        a0, a1 = np.exp(ys[i]), np.exp(ys[j + 1]) if j + 1 < len(ys) else 1.0
        if coupled[i]:
            def f(y):
                a = np.exp(y); p = q / a; v = p / np.sqrt(1 + p * p)
                G = gam(a, q) * T.hubble(a)
                return 2 * (T.C_MPC_S * v) ** 2 / (G * a * a * T.hubble(a))
            lamD2 += quad(f, np.log(a0), np.log(a1), limit=200)[0]
            last_dec = a1
        else:
            lam2 += T.horizon(a0, q, a1) ** 2
        i = j + 1
    lam_dec = T.horizon(last_dec, q, 1.0)
    lam_eff = np.sqrt(lam2 + lamD2) if last_dec > a_start else T.horizon(a_start, q, 1.0)
    return dict(GammaT_over_H_start=float(g0), GammaT_over_H_max=float(g.max()),
                a_dec=float(last_dec), z_dec=float(1 / last_dec - 1),
                horizon_collisionless_Mpc=T.horizon(a_start, q, 1.0), horizon_from_dec_Mpc=lam_dec,
                diffusion_length_Mpc=float(np.sqrt(lamD2)), effective_Mpc=float(lam_eff))


def q_eff(a_start, lam_target, q0):
    return brentq(lambda q: T.horizon(a_start, q, 1.0) - lam_target, q0 * 1e-4, q0 * 1.5)


def main():
    rows = []
    files = [T.SCR / 'mc_aesc_T1keV_codex_1.npy'] + sorted(T.SCR.glob('mc_aesc_epoch_codex_*.npy'))
    for f in files:
        ae = np.load(f)
        meta = json.loads(f.with_suffix('.json').read_text()) if f.with_suffix('.json').exists() else {'a_form': 2.34865418e-07}
        af = meta['a_form']
        a_med = float(np.median(ae)); q = T.UESC * a_med
        c = cohort(a_med, q)
        qe = q_eff(a_med, c['effective_Mpc'], q)
        qd = af * T.U0
        d = cohort(af, qd)
        qde = q_eff(af, d['effective_Mpc'], qd)
        rows.append(dict(tag=f.stem.replace('mc_aesc_', ''), a_form=af, a_escape_median=a_med, q_over_m=q,
                         converted=c, q_eff_converted=qe, q_eff_ratio=qe / q,
                         daughters=d, q_daughters=qd, q_eff_daughters=qde))
        print(rows[-1]['tag'], c, qe / q, d['effective_Mpc'], qde / qd, flush=True)
    (T.RES / 'elastic_epochs.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main()
