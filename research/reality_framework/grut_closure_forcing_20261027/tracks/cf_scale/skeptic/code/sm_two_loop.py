#!/usr/bin/env python3
"""cf_scale_skeptic: P4-SM re-check with an independent SM running code (own code).

1. Two-loop SM RGEs (g1 GUT-normalised, g2, g3 + 3-loop pure-QCD term, y_t, lambda with V = lambda |H|^4) from the
   NNLO MS-bar matching at mu = M_t of Buttazzo et al. arXiv:1307.3536 eqs. (55)-(60) [Q, re-fetched 27 Oct 2026].
   Compared with their NNLO fit lambda(M_Pl) eq. (61) [Q] (M_Pl = 1.22e19 GeV): central value and slopes.
2. PDG 2026 inputs (pdglive, fetched 27 Oct 2026): M_h = 125.13 +- 0.11, M_t = 172.60 +- 0.27 (direct; +-0.3 MC-pole
   added as cf_scale did) or 172.1 +- 0.6 (pole from cross sections); alpha_s = 0.1180 +- 0.0009 (PDG value as cited;
   pdglive QCD page returned HTTP 500 on 27 Oct 2026, not re-verified).
3. The comparator. cf_scale fixed lambda(M_P) = 0 "in advance". Under the stated postulate (beta_lambda += A_lam lambda
   above M_P, sharp switch) the prediction is the critical-surface value lambda_crit(M_P; A_lam), which equals 0 only
   when A_lam dominates the SM one-loop terms. Computed here for A_lam = 0.0237 ... 2.83 under three transplanckian
   treatments of the other SM couplings:
     (a) frozen: g1, g2, g3, y_t held at their M_P values (quasi-static fixed point);
     (b) g1, y_t at fixed points (frozen), g2, g3 on their Gaussian branch with beta += -f_g g (f_g = 9.83e-3);
     (c) as (b) but y_t also on its Gaussian branch, beta_yt += -f_y y_t, f_y = 0.01 or 0.03 [supplied, illustrative].
         Variant (c) is NOT stable against the UV start depth N at A_lam = 0.0237 (multi-branch); it is reported but not used.
   One-loop SM betas are used above M_P (the gravity parametrisation is one-loop-level anyway).
4. Absolute stability (the q-vacuum P_q reading) versus lambda(M_P) = 0: the M_t thresholds at PDG M_h, alpha_s.
Usage: python3 -B sm_two_loop.py <out.json>
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.dont_write_bytecode = True
K = 16 * np.pi**2
MPL, MPRED = 1.22e19, 2.435e18
FG = 9.83e-3


def beta_sm(y, loops=2):
    g1, g2, g3, yt, lam = y
    g12, g22, g32, yt2 = g1 * g1, g2 * g2, g3 * g3, yt * yt
    gp2 = 0.6 * g12                                        # hypercharge g'^2
    b1 = np.array([41 / 10 * g1**3, -19 / 6 * g2**3, -7 * g3**3,
                   yt * (4.5 * yt2 - 17 / 20 * g12 - 9 / 4 * g22 - 8 * g32),
                   24 * lam**2 - 6 * yt2**2 + 3 / 8 * (2 * g22**2 + (g22 + gp2)**2) + lam * (12 * yt2 - 9 / 5 * g12 - 9 * g22)]) / K
    if loops == 1:
        return b1
    b2g1 = g1**3 * (199 / 50 * g12 + 27 / 10 * g22 + 44 / 5 * g32 - 17 / 10 * yt2)
    b2g2 = g2**3 * (9 / 10 * g12 + 35 / 6 * g22 + 12 * g32 - 3 / 2 * yt2)
    b2g3 = g3**3 * (11 / 10 * g12 + 9 / 2 * g22 - 26 * g32 - 2 * yt2)
    b2yt = yt * (-12 * yt2**2 + yt2 * (393 / 80 * g12 + 225 / 16 * g22 + 36 * g32 - 12 * lam)
                 + 1187 / 600 * g12**2 - 9 / 20 * g12 * g22 + 19 / 15 * g12 * g32 - 23 / 4 * g22**2 + 9 * g22 * g32
                 - 108 * g32**2 + 6 * lam**2)
    g4, gp4 = g22**2, gp2**2
    b2lam = (-312 * lam**3 - 144 * lam**2 * yt2 + 36 * lam**2 * (3 * g22 + gp2) - 3 * lam * yt2**2 + 80 * lam * g32 * yt2
             + lam * yt2 * (45 / 2 * g22 + 85 / 6 * gp2) - 73 / 8 * lam * g4 + 39 / 4 * lam * g22 * gp2 + 629 / 24 * lam * gp4
             + 30 * yt2**3 - 32 * g32 * yt2**2 - 8 / 3 * gp2 * yt2**2 - 9 / 4 * g4 * yt2 + 21 / 2 * g22 * gp2 * yt2
             - 19 / 4 * gp4 * yt2 + 305 / 16 * g22**3 - 289 / 48 * g4 * gp2 - 559 / 48 * g22 * gp4 - 379 / 48 * gp2**3)
    b3g3 = 65 / 2 * g3**7                                   # leading 3-loop QCD
    return b1 + np.array([b2g1, b2g2, b2g3 + b3g3 / K, b2yt, b2lam]) / K**2


def match(Mt, Mh, als, MW=80.384):
    """Buttazzo et al. 1307.3536 eqs. (55)-(60): MS-bar couplings at mu = M_t [Q]."""
    da = (als - 0.1184) / 0.0007; dt = Mt - 173.34; dh = Mh - 125.15; dw = (MW - 80.384) / 0.014
    lam = 0.12604 + 0.00206 * dh - 0.00004 * dt
    yt = 0.93690 + 0.00556 * dt - 0.00042 * da
    g3 = 1.1666 + 0.00314 * da - 0.00046 * dt
    g2 = 0.64779 + 0.00004 * dt + 0.00011 * dw
    gY = 0.35830 + 0.00011 * dt - 0.00020 * dw
    return np.array([np.sqrt(5 / 3) * gY, g2, g3, yt, lam])


def run(Mt, Mh, als, mu_list, loops=2, dense=False):
    y0 = match(Mt, Mh, als)
    t1 = np.log(max(mu_list) / Mt)
    sol = solve_ivp(lambda t, y: beta_sm(y, loops), (0, t1), y0, method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
    return {mu: sol.sol(np.log(mu / Mt)) for mu in mu_list}, sol


def lam_formula(Mt, Mh, als):
    return -0.0143 - 0.0066 * (Mt - 173.34) + 0.0018 * (als - 0.1184) / 0.0007 + 0.0029 * (Mh - 125.15)


# ------------------------------------------------------------------ AS critical surface for lambda
def crit_lambda(yP, A, mode, f_y=0.0, N=60.0):
    """Critical-surface lambda(M_P) for beta_lambda = beta_SM^(1) + A lambda above M_P (t = ln k/M_P >= 0)."""
    g1P, g2P, g3P, ytP, _ = yP
    if mode == 'frozen':
        g1, g2, g3, yt = g1P, g2P, g3P, ytP
        a = 24 / K; b = A + (12 * yt**2 - 9 / 5 * g1**2 - 9 * g2**2) / K
        c = (-6 * yt**4 + 3 / 8 * (2 * g2**4 + (g2**2 + 0.6 * g1**2)**2)) / K
        disc = b * b - 4 * a * c
        if disc < 0:
            return None
        return (-b + np.sqrt(disc)) / (2 * a)            # the root continuously connected to 0 as c -> 0 (irrelevant one)
    # Gaussian-branch nonabelian couplings: exact solution of beta = b g^3/K - f_g g
    def gauss(gP, b, f, t):
        c = b / (K * f)
        return 1 / np.sqrt(c + (1 / gP**2 - c) * np.exp(2 * f * t))

    def rhs(t, y):
        lam, yt = y
        g2 = gauss(g2P, -19 / 6, FG, t); g3 = gauss(g3P, -7.0, FG, t); g1 = g1P
        bl = (24 * lam**2 - 6 * yt**4 + 3 / 8 * (2 * g2**4 + (g2**2 + 0.6 * g1**2)**2)
              + lam * (12 * yt**2 - 9 / 5 * g1**2 - 9 * g2**2)) / K + A * lam
        if mode == 'gauss_g23':
            byt = 0.0                                     # y_t at its fixed point (frozen)
        else:
            byt = yt * (4.5 * yt**2 - 17 / 20 * g1**2 - 9 / 4 * g2**2 - 8 * g3**2) / K - f_y * yt
        return [bl, byt]
    tUV = N / A
    # y_t trajectory: for the Gaussian-branch case integrate y_t upward first (relevant: fix its M_P value)
    if mode == 'gauss_g23':
        ytUV = ytP
    else:
        s = solve_ivp(lambda t, y: [rhs(t, [0.0, y[0]])[1]], (0, tUV), [ytP], method='LSODA', rtol=1e-10, atol=1e-14)
        ytUV = s.y[0, -1]
    out = []
    for lam0 in (0.0, 0.05):                              # two UV starts: irrelevance check
        s = solve_ivp(rhs, (tUV, 0.0), [lam0, ytUV], method='LSODA', rtol=1e-10, atol=1e-14)
        out.append(s.y[0, -1])
    return float(out[0]) if abs(out[0] - out[1]) < 1e-6 else None


def main(out):
    res = dict(method=__doc__)
    # ---- 1. reproduce Buttazzo central value and slopes
    ref = dict(Mt=173.34, Mh=125.15, als=0.1184)
    v, _ = run(ref['Mt'], ref['Mh'], ref['als'], [MPL, MPRED])
    lam_ref = v[MPL][4]
    slopes = {}
    for k, h in (('Mt', 0.5), ('Mh', 0.5), ('als', 0.0007)):
        p = dict(ref); m = dict(ref); p[k] += h; m[k] -= h
        lp = run(p['Mt'], p['Mh'], p['als'], [MPL])[0][MPL][4]; lm = run(m['Mt'], m['Mh'], m['als'], [MPL])[0][MPL][4]
        slopes[k] = (lp - lm) / (2 * h) * (0.0007 if k == 'als' else 1.0)
    res['reproduction'] = dict(lambda_MPl_2loop=lam_ref, lambda_MPl_Buttazzo_NNLO=-0.0143,
                               couplings_MPl_2loop=dict(zip(['g1', 'g2', 'g3', 'yt', 'lam'], v[MPL].tolist())),
                               Buttazzo_MPl=dict(g1=0.6154, g2=0.5055, g3=0.4873, yt=0.3825, lam=-0.0143),
                               slopes_2loop=dict(dlam_dMt=slopes['Mt'], dlam_dMh=slopes['Mh'], dlam_per_0p0007_als=slopes['als']),
                               slopes_Buttazzo=dict(dlam_dMt=-0.0066, dlam_dMh=0.0029, dlam_per_0p0007_als=0.0018),
                               lambda_MPred_minus_lambda_MPl=v[MPRED][4] - v[MPL][4])
    print('reference: 2-loop lambda(M_Pl) = %.5f (Buttazzo NNLO -0.0143); slopes %s; lambda(2.435e18)-lambda(1.22e19) = %.5f' % (
        lam_ref, {k: round(s, 5) for k, s in slopes.items()}, v[MPRED][4] - v[MPL][4]), flush=True)
    # ---- 2. PDG 2026
    Mh, sMh, als, sals = 125.13, 0.11, 0.1180, 0.0009
    pdg = {}
    for tag, Mt, sMt in (('mt_direct', 172.60, float(np.hypot(0.27, 0.3))), ('mt_xsec', 172.1, 0.6)):
        lf = lam_formula(Mt, Mh, als)
        sig = float(np.sqrt((0.0066 * sMt)**2 + (0.0018 * sals / 0.0007)**2 + (0.0029 * sMh)**2 + 0.0003**2))
        vv, _ = run(Mt, Mh, als, [MPL, MPRED])
        pdg[tag] = dict(Mt=Mt, sigma_Mt=sMt, lambda_MPl_formula=lf, sigma=sig, lambda_MPl_2loop=float(vv[MPL][4]),
                        lambda_MPred_2loop=float(vv[MPRED][4]),
                        lambda_MPred_formula_shifted=float(lf + (vv[MPRED][4] - vv[MPL][4])),
                        couplings_MPred_2loop=dict(zip(['g1', 'g2', 'g3', 'yt', 'lam'], vv[MPRED].tolist())),
                        pull_vs_0_formula=lf / sig)
        print('%s: lambda(M_Pl) formula %.4f +- %.4f (pull %.2f); 2-loop own %.4f; at reduced M_P: %.4f' % (
            tag, lf, sig, lf / sig, vv[MPL][4], vv[MPRED][4]), flush=True)
    res['pdg2026'] = pdg
    # ---- 3. comparator: critical-surface lambda(M_P) as a function of A_lam
    yP = np.array([pdg['mt_direct']['couplings_MPred_2loop'][k] for k in ['g1', 'g2', 'g3', 'yt', 'lam']])
    comp = []
    for A in (0.0237, 0.056, 0.1, 0.3, 1.0, 2.83):
        row = dict(A_lam=A)
        row['frozen'] = crit_lambda(yP, A, 'frozen')
        row['g23_gaussian_yt_FP'] = crit_lambda(yP, A, 'gauss_g23')
        row['g23_yt_gaussian_fy0.01'] = crit_lambda(yP, A, 'gauss_all', f_y=0.01)
        row['g23_yt_gaussian_fy0.03'] = crit_lambda(yP, A, 'gauss_all', f_y=0.03)
        # pulls of the measured value (shifted to the reduced Planck mass, formula central) against each comparator
        for tag in ('mt_direct', 'mt_xsec'):
            meas = pdg[tag]['lambda_MPred_formula_shifted']; sig = pdg[tag]['sigma']
            row['pull_' + tag] = {k: (None if row[k] is None else (meas - row[k]) / sig)
                                  for k in ('frozen', 'g23_gaussian_yt_FP', 'g23_yt_gaussian_fy0.01', 'g23_yt_gaussian_fy0.03')}
        comp.append(row)
        print('A=%.4g crit lambda(M_P): frozen %s | g2,g3 Gaussian %s | +yt Gauss f_y=.01 %s | .03 %s' % (
            A, *['%.4f' % row[k] if row[k] is not None else 'none' for k in ('frozen', 'g23_gaussian_yt_FP',
                                                                                 'g23_yt_gaussian_fy0.01', 'g23_yt_gaussian_fy0.03')]),
              ' pulls(direct):', {k: (None if p is None else round(p, 2)) for k, p in row['pull_mt_direct'].items()}, flush=True)
    res['comparator_scan'] = comp
    # ---- 4. absolute stability vs lambda(M_P) = 0 (own 2-loop code, PDG M_h, alpha_s)
    def min_lam(Mt):
        _, sol = run(Mt, Mh, als, [MPL])
        ts = np.linspace(np.log(1e3 / Mt), np.log(MPL / Mt), 3000)
        return float(np.min(sol.sol(ts)[4]))
    Mt_stab = brentq(min_lam, 168, 176, xtol=1e-4)
    Mt_sw_MPl = brentq(lambda m: run(m, Mh, als, [MPL])[0][MPL][4], 168, 176, xtol=1e-4)
    Mt_sw_MPred = brentq(lambda m: run(m, Mh, als, [MPRED])[0][MPRED][4], 168, 176, xtol=1e-4)
    res['stability_vs_SW'] = dict(Mt_absolute_stability_2loop=Mt_stab, Mt_lambda_MPl_zero_2loop=Mt_sw_MPl,
                                  Mt_lambda_MPred_zero_2loop=Mt_sw_MPred,
                                  Mt_lambda_MPl_zero_formula=173.34 + (-0.0143 + 0.0018 * (als - 0.1184) / 0.0007 + 0.0029 * (Mh - 125.15)) / 0.0066,
                                  Mt_stability_Buttazzo_at_PDG=173.34 + (Mh - 129.6 + 0.5 * (als - 0.1184) / 0.0007) / 2.0)
    print('M_t thresholds (2-loop own): absolute stability %.2f, lambda(M_Pl)=0 %.2f, lambda(reduced M_P)=0 %.2f GeV' % (
        Mt_stab, Mt_sw_MPl, Mt_sw_MPred), flush=True)
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
