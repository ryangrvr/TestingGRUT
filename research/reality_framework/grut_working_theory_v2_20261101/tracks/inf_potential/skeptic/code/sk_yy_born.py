#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude, independent method): holder-induced Y-pair production phi_I -> YY in first-order
(Born) perturbation theory, computed in the TIME DOMAIN on our own post-inflation background (sk_cmb background +
Gamma_I = 0.5 H_i friction and radiation), instead of the track's per-period Fourier (harmonic) decomposition.
For a heavy real field with V = m^2 y^2/2 + dm2(t) y^2/2,  beta_k = (i/2) Int dt dm2(t)/w_k(t) exp(-2 i Int w_k dt);
for the off-diagonal  V = -c' phi y_L y_H,  beta_LH = (i/2) Int dt c' phi(t)/sqrt(w_L w_H) exp(-i Int (w_L + w_H) dt);
w_k = sqrt(k^2/a^2 + m^2) (m >> H, the 9H^2/4 terms dropped).  Comoving numbers N = (1/2 pi^2) Int k^2 |beta_k|^2 dk,
a = 1 at eps_H = 1.  Late-vacuum couplings (u = mu, Sigma = h/sqrt2) as in the track:
dm2_L = -sqrt2 h Re(G) phi, dm2_H = +sqrt2 h Re(G) phi, c' = sqrt2 h Im(G), |G phihat_I(S1)| = |10i - mu|.
This includes every harmonic, the amplitude decay and the expansion automatically.  A smooth switch-on window
sin^2(pi t/2 t_on) suppresses the spurious boundary term at eps_H = 1 (t_on varied).  Not included (as in the track):
second order in c (phi phi -> YY), Y/Sigma condensates in the first H_i^-1.
Output ../results/sk_yy_born.json.  Label [C] (Born) on our background.
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_cmb as cm
OUT = Path(__file__).resolve().parent.parent / 'results'
MU, H = 7.159456159, 40.0
ALPHA = math.atan2(10.0, -MU)
ML = math.sqrt(900 - math.sqrt(2) * MU * H); MH = math.sqrt(900 + math.sqrt(2) * MU * H)


def post(name, gam=0.5, tmax=12.0, Hi=1.29e11):
    d = json.load(open(cm.TRK / 'results' / ('bench_%s.json' % name)))['by_Hi'][0]
    P = cm.Pot(name, d['params'])
    sol, tend = cm.background(P, cm.f_start(P, 72.0))
    ts = np.linspace(0, tend, 200001); Y = sol.sol(ts); Nb = Y[2, -1] - Y[2]
    Hs = np.sqrt((0.5 * Y[1] ** 2 + P.V(Y[0])) / 3)
    hs = float(np.interp(-1.0, -Nb, Hs))               # H(S1) in code units -> t_Hi = t_code * hs
    fS1 = float(np.interp(-1.0, -Nb, Y[0]))
    pmin = P.phimin()
    f0, fd0 = Y[0, -1], Y[1, -1]
    G = gam * hs

    def rhs(t, y):
        f, fd, lna, rr = y
        Hh = math.sqrt((0.5 * fd * fd + P.V(f) + rr) / 3)
        return [fd, -(3 * Hh + G) * fd - P.dV(f), Hh, -4 * Hh * rr + G * fd * fd]
    s = solve_ivp(rhs, (0, tmax / hs), [f0, fd0, 0.0, 0.0], method='DOP853', rtol=1e-11, atol=1e-16, dense_output=True)
    return P, s, hs, (fS1 - pmin), pmin, d


def born(name, gam=0.5, tmax=12.0, t_on=0.25, nk=241, kmax=120.0, dt=0.001):
    P, s, hs, phiS1, pmin, d = post(name, gam, tmax)
    t = np.arange(0, tmax + dt / 2, dt)                 # H_i^-1
    Y = s.sol(t / hs)
    a = np.exp(Y[2]); sfrac = (Y[0] - pmin) / phiS1     # s(t) = phi_I/phi_I(S1)
    w = np.where(t < t_on, np.sin(0.5 * math.pi * t / t_on) ** 2, 1.0)
    gphi = abs(10j - MU) * sfrac * w                   # |G| phihat_I(t) in H_i
    dmL = math.sqrt(2) * H * abs(math.cos(ALPHA)) * gphi
    cLH = math.sqrt(2) * H * abs(math.sin(ALPHA)) * gphi
    ks = np.linspace(0, kmax, nk)
    out = {}
    for lab, m1, m2, amp in (('LL', ML, ML, dmL), ('HH', MH, MH, dmL), ('LH', ML, MH, cLH)):
        bk = np.zeros(nk)
        for i0 in range(0, nk, 30):
            kk = ks[i0:i0 + 30][:, None]
            w1 = np.sqrt((kk / a[None, :]) ** 2 + m1 ** 2); w2 = np.sqrt((kk / a[None, :]) ** 2 + m2 ** 2)
            th = np.concatenate([np.zeros((kk.shape[0], 1)), np.cumsum(0.5 * ((w1 + w2)[:, 1:] + (w1 + w2)[:, :-1]) * dt, axis=1)], axis=1)
            if lab == 'LH':
                integ = amp[None, :] / np.sqrt(w1 * w2) * np.exp(-1j * th)
            else:
                integ = amp[None, :] / w1 * np.exp(-1j * th)
            b = 0.5 * np.trapezoid(integ, dx=dt, axis=1)
            bk[i0:i0 + 30] = np.abs(b) ** 2
        Nc = float(np.trapezoid(ks ** 2 * bk, ks) / (2 * math.pi ** 2))
        out[lab] = dict(N=Nc, spectrum_peak_k=float(ks[np.argmax(ks ** 2 * bk)]))
    NL = out['LL']['N'] + out['LH']['N']; NH = out['HH']['N'] + out['LH']['N']
    return dict(benchmark=name, gamma=gam, t_max=tmax, t_on=t_on, channels=out, N_YL=NL, N_YH=NH,
                track=dict(N_YL=[x for x in json.load(open(cm.TRK / 'results' / 'yy_harmonics.json'))['runs']
                                 if x['benchmark'] == name and x['gamma_over_Hi'] == gam][0]['res']['S1']['N_YL'],
                           N_YH=[x for x in json.load(open(cm.TRK / 'results' / 'yy_harmonics.json'))['runs']
                                 if x['benchmark'] == name and x['gamma_over_Hi'] == gam][0]['res']['S1']['N_YH']))


if __name__ == '__main__':
    names = sys.argv[1:] or ['HT6t']
    runs = []
    for n in names:
        for t_on in (0.15, 0.3, 0.6):
            r = born(n, t_on=t_on)
            runs.append(r)
            print('%s t_on=%.2f: N_YL=%.4g N_YH=%.4g (LL %.3g, HH %.3g, LH %.3g; peak k %s) | track N_YL=%.4g N_YH=%.4g' % (
                n, t_on, r['N_YL'], r['N_YH'], r['channels']['LL']['N'], r['channels']['HH']['N'], r['channels']['LH']['N'],
                [r['channels'][c]['spectrum_peak_k'] for c in ('LL', 'HH', 'LH')], r['track']['N_YL'], r['track']['N_YH']), flush=True)
    json.dump(dict(method=__doc__, runs=runs), open(OUT / 'sk_yy_born.json', 'w'), indent=1)
