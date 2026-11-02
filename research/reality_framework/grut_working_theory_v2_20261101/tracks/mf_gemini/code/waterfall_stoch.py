"""mf_gemini: stochastic (phase-space Langevin) first-passage delta N for the mild waterfall [C, approximate].

Cross-check of the delta N tail of waterfall_dn.py with an independent method:
  e-fold time; dphi_i = pi_i dN + (h hub/2 pi) dW_i (noise on both fields, every e-fold step);
  dpi_i = -[(3 - eps) pi_i + V_i/h^2] dN;  h^2 = V/(3 - eps).
Realisations start at sigma = 1.5 sigma_c (chi heavy, m^2 = 5.6 H^2), chi = 0, slow-roll pi_sigma.
Stop at first passage |chi| = chi_stop = 1e-3 chi_trap (classical: >> H/2pi, still linear); the remaining
evolution is deterministic and nearly the same for all realisations, so delta N = delta N_stop.
Limits: the white-noise amplitude h/2pi is the light-field value; chi is not light (|m^2| <= 4.5 H^2),
so the amplitude of the noise is uncertain by O(1) and this is a check of the tail SHAPE (exponential,
set by the growth rate), not of the spectrum.

Run: python3 -B waterfall_stoch.py -> results/waterfall_stoch.json
"""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfg_common import Model, MP_GEV, eps_star

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')


def run(D, kind='quad', H=1.3e11, nreal=20000, dN=0.004, Nmax=80.0, seed=1):
    rng = np.random.default_rng(seed)
    hub = H / MP_GEV
    pc = np.sqrt(2 * eps_star(H)) * np.exp(0.54)
    sigc = pc / D
    m = Model(sigc, dict(kind=kind, eta=D, p=pc), hub=hub)
    chi_trap = np.sqrt(2 * pc * sigc / 3)
    chi_stop = 1e-3 * chi_trap
    s = np.full(nreal, 1.5 * sigc)
    c = np.zeros(nreal)
    V = m.V(s, c)
    ps = -m.Vs(s, c) / V          # slow roll: pi = -V'/V (M_P = 1)
    pch = np.zeros(nreal)
    N = np.zeros(nreal)
    done = np.zeros(nreal, bool)
    Nstop = np.full(nreal, np.nan)
    sstop = np.full(nreal, np.nan)
    n = 0
    sq = np.sqrt(dN)
    while (~done).any() and n * dN < Nmax:
        a = ~done
        eps = 0.5 * (ps[a]**2 + pch[a]**2)
        Va = m.V(s[a], c[a])
        h2 = Va / (3 - eps)
        h = np.sqrt(h2)
        Fs = m.Vs(s[a], c[a]) / h2
        Fc = m.Vc(s[a], c[a]) / h2
        noise = h * hub / (2 * np.pi) * sq
        s[a] = s[a] + ps[a] * dN + noise * rng.standard_normal(a.sum())
        c[a] = c[a] + pch[a] * dN + noise * rng.standard_normal(a.sum())
        ps[a] = ps[a] - ((3 - eps) * ps[a] + Fs) * dN
        pch[a] = pch[a] - ((3 - eps) * pch[a] + Fc) * dN
        n += 1
        hit = a.copy()
        hit[a] = np.abs(c[a]) >= chi_stop
        Nstop[hit] = n * dN
        sstop[hit] = s[hit]
        done |= hit
    ok = np.isfinite(Nstop)
    x = Nstop[ok] - Nstop[ok].mean()
    tails = {f'{z}': float((x > z).mean()) for z in [0.5, 0.67, 1.0, 2.0, 3.0]}
    # exponential tail slope from the survival function between 1 and 3
    zz = np.linspace(1.0, 3.0, 9)
    surv = np.array([(x > z).mean() for z in zz])
    good = surv > 0
    slope = float(-np.polyfit(zz[good], np.log(surv[good]), 1)[0]) if good.sum() > 3 else None
    return dict(Delta=D, kind=kind, sigc=float(sigc), chi_stop=float(chi_stop), nreal=nreal, dN=dN,
                frac_finished=float(ok.mean()), mean_Nstop=float(Nstop[ok].mean()), varN=float(x.var()),
                tails=tails, tail_slope_1to3=slope,
                sigma_at_stop_over_sigc_mean=float(np.nanmean(sstop / sigc)))


def main():
    out = dict(method=__doc__.split('\n')[0], runs=[])
    for D in [1.0, 0.3, 0.1]:
        r = run(D)
        print(r, flush=True)
        out['runs'].append(r)
    r = run(1.0, dN=0.002, nreal=10000, seed=2)
    r['note'] = 'step-size check (dN halved)'
    print(r, flush=True)
    out['runs'].append(r)
    with open(os.path.join(RES, 'waterfall_stoch.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
