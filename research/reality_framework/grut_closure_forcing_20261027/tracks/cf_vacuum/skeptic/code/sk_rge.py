"""Skeptic re-derivation of the SM absolute-stability line (track cf_vacuum, closure-forcing pass, 27 Oct 2026).

Independent of the track's code: nothing is imported from ../../code. RGEs transcribed again, by hand, from
Buttazzo, Degrassi, Giardino, Giudice, Sala, Salvio, Strumia, arXiv:1307.3536, appendix A (TeX artv2.tex,
lines 2656-2900): full 3-loop beta functions of g1^2 (GUT), g2^2, g3^2 (+4-loop pure QCD), y_t^2 and lambda,
written as d/d ln(mu^2), with V = -m^2|H|^2/2 + lambda |H|^4.  y_b, y_tau enter at 1-2 loops (their 1-loop
running is standard; effect on the line << 0.01 GeV).

Two matchings at the weak scale, neither of which is SMDR (the track's line used SMDR):
  'B13'  Buttazzo et al. 1307.3536 sec. 3 fits at mu = M_t (3-loop QCD in y_t; TeX lines 1155-1225).
  'B15'  Bednyakov, Kniehl, Pikelner, Veretin 1507.08833 Table 1 (prl.tex lines 306-340) at mu = 173.21 GeV:
         full 2-loop EW + O(alpha alpha_s) + 4-loop QCD (Marquard et al. 1502.01030); coefficients per 1 sigma
         of their PDG-2014 inputs (alpha_s 0.1185(6), M_W 80.385(15), M_H 125.7(4), M_t 173.21(87)).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

K1 = (4 * np.pi) ** 2
K2 = K1 ** 2
K3 = K1 ** 3
K4 = K1 ** 4
MPL = 1.2209e19          # non-reduced Planck mass (same cut as the track)
MPLR = 2.435e18          # reduced Planck mass


def betas(t, y, loops=3):
    """d/d ln mu^2 of (G1, G2, G3, YT, L, YB, YTAU) = (g1^2, g2^2, g3^2, yt^2, lambda, yb^2, ytau^2)."""
    G1, G2, G3, YT, L, YB, YTAU = y
    # ---- g1^2
    b1 = G1 ** 2 / K1 * (41 / 10)
    b2 = G1 ** 2 / K2 * (44 * G3 / 5 + 27 * G2 / 10 + 199 * G1 / 50 - 17 * YT / 10 - YB / 2 - 3 * YTAU / 2)
    b3 = G1 ** 2 / K3 * (YT * (189 * YT / 16 - 29 * G3 / 5 - 471 * G2 / 32 - 2827 * G1 / 800)
                         + L * (-9 * L / 5 + 9 * G2 / 10 + 27 * G1 / 50)
                         + 297 * G3 ** 2 / 5 + 789 * G2 ** 2 / 64 - 388613 * G1 ** 2 / 24000
                         - 3 * G3 * G2 / 5 - 137 * G3 * G1 / 75 + 123 * G2 * G1 / 160)
    dG1 = b1 + (b2 if loops >= 2 else 0) + (b3 if loops >= 3 else 0)
    # ---- g2^2
    b1 = G2 ** 2 / K1 * (-19 / 6)
    b2 = G2 ** 2 / K2 * (12 * G3 + 35 * G2 / 6 + 9 * G1 / 10 - 3 * YT / 2 - 3 * YB / 2 - YTAU / 2)
    b3 = G2 ** 2 / K3 * (YT * (147 * YT / 16 - 7 * G3 - 729 * G2 / 32 - 593 * G1 / 160)
                         + L * (-3 * L + 3 * G2 / 2 + 3 * G1 / 10)
                         + 81 * G3 ** 2 + 324953 * G2 ** 2 / 1728 - 5597 * G1 ** 2 / 1600
                         + 39 * G3 * G2 - G3 * G1 / 5 + 873 * G2 * G1 / 160)
    dG2 = b1 + (b2 if loops >= 2 else 0) + (b3 if loops >= 3 else 0)
    # ---- g3^2
    b1 = G3 ** 2 / K1 * (-7)
    b2 = G3 ** 2 / K2 * (-26 * G3 + 9 * G2 / 2 + 11 * G1 / 10 - 2 * YT - 2 * YB)
    b3 = G3 ** 2 / K3 * (YT * (15 * YT - 40 * G3 - 93 * G2 / 8 - 101 * G1 / 40)
                         + 65 * G3 ** 2 / 2 + 109 * G2 ** 2 / 8 - 523 * G1 ** 2 / 120
                         + 21 * G3 * G2 + 77 * G3 * G1 / 15 - 3 * G2 * G1 / 40)
    b4 = G3 ** 5 / K4 * (-2472.28)
    dG3 = b1 + (b2 if loops >= 2 else 0) + (b3 + b4 if loops >= 3 else 0)
    # ---- lambda
    b1 = (L * (12 * L + 6 * YT + 6 * YB + 2 * YTAU - 9 * G2 / 2 - 9 * G1 / 10) - 3 * YT ** 2 - 3 * YB ** 2
          - YTAU ** 2 + 9 * G2 ** 2 / 16 + 27 * G1 ** 2 / 400 + 9 * G2 * G1 / 40) / K1
    b2 = (L ** 2 * (-156 * L - 72 * YT - 72 * YB - 24 * YTAU + 54 * G2 + 54 * G1 / 5)
          + L * YT * (-3 * YT / 2 - 21 * YB + 40 * G3 + 45 * G2 / 4 + 17 * G1 / 4)
          + L * YB * (-3 * YB / 2 + 40 * G3 + 45 * G2 / 4 + 5 * G1 / 4)
          + L * YTAU * (-YTAU / 2 + 15 * G2 / 4 + 15 * G1 / 4)
          + L * (-73 * G2 ** 2 / 16 + 1887 * G1 ** 2 / 400 + 117 * G2 * G1 / 40)
          + YT ** 2 * (15 * YT - 3 * YB - 16 * G3 - 4 * G1 / 5)
          + YT * (-9 * G2 ** 2 / 8 - 171 * G1 ** 2 / 200 + 63 * G2 * G1 / 20)
          + YB ** 2 * (-3 * YT + 15 * YB - 16 * G3 + 2 * G1 / 5)
          + YB * (-9 * G2 ** 2 / 8 + 9 * G1 ** 2 / 40 + 27 * G2 * G1 / 20)
          + YTAU ** 2 * (5 * YTAU - 6 * G1 / 5)
          + YTAU * (-3 * G2 ** 2 / 8 - 9 * G1 ** 2 / 8 + 33 * G2 * G1 / 20)
          + 305 * G2 ** 3 / 32 - 3411 * G1 ** 3 / 4000 - 289 * G2 ** 2 * G1 / 160 - 1677 * G2 * G1 ** 2 / 800) / K2
    b3 = (L ** 3 * (6011.35 * L + 873 * YT - 387.452 * G2 - 77.490 * G1)
          + L ** 2 * YT * (1768.26 * YT + 160.77 * G3 - 359.539 * G2 - 63.869 * G1)
          + L ** 2 * (-790.28 * G2 ** 2 - 185.532 * G1 ** 2 - 316.64 * G2 * G1)
          + L * YT ** 2 * (-223.382 * YT - 662.866 * G3 - 5.470 * G2 - 21.015 * G1)
          + L * YT * (356.968 * G3 ** 2 - 319.664 * G2 ** 2 - 74.8599 * G1 ** 2 + 15.1443 * G3 * G2
                      + 17.454 * G3 * G1 + 5.615 * G2 * G1)
          + L * G2 ** 2 * (-57.144 * G3 + 865.483 * G2 + 79.638 * G1)
          + L * G1 ** 2 * (-8.381 * G3 + 61.753 * G2 + 28.168 * G1)
          + YT ** 3 * (-243.149 * YT + 250.494 * G3 + 74.138 * G2 + 33.930 * G1)
          + YT ** 2 * (-50.201 * G3 ** 2 + 15.884 * G2 ** 2 + 15.948 * G1 ** 2 + 13.349 * G3 * G2
                       + 17.570 * G3 * G1 - 70.356 * G2 * G1)
          + YT * G3 * (16.464 * G2 ** 2 + 1.016 * G1 ** 2 + 11.386 * G2 * G1)
          + YT * G2 ** 2 * (62.500 * G2 + 13.041 * G1)
          + YT * G1 ** 2 * (10.627 * G2 + 11.117 * G1)
          + G3 * (7.536 * G2 ** 3 + 0.663 * G1 ** 3 + 1.507 * G2 ** 2 * G1 + 1.105 * G2 * G1 ** 2)
          - 114.091 * G2 ** 4 - 1.508 * G1 ** 4 - 37.889 * G2 ** 3 * G1 + 6.500 * G2 ** 2 * G1 ** 2
          - 1.543 * G2 * G1 ** 3) / K3
    dL = b1 + (b2 if loops >= 2 else 0) + (b3 if loops >= 3 else 0)
    # ---- y_t^2
    b1 = YT / K1 * (9 * YT / 2 + 3 * YB / 2 + YTAU - 8 * G3 - 9 * G2 / 4 - 17 * G1 / 20)
    b2 = YT / K2 * (YT * (-12 * YT - 11 * YB / 4 - 9 * YTAU / 4 - 12 * L + 36 * G3 + 225 * G2 / 16 + 393 * G1 / 80)
                    + YB * (-YB / 4 + 5 * YTAU / 4 + 4 * G3 + 99 * G2 / 16 + 7 * G1 / 80)
                    + YTAU * (-9 * YTAU / 4 + 15 * G2 / 8 + 15 * G1 / 8)
                    + 6 * L ** 2 - 108 * G3 ** 2 - 23 * G2 ** 2 / 4 + 1187 * G1 ** 2 / 600
                    + 9 * G3 * G2 + 19 * G3 * G1 / 15 - 9 * G2 * G1 / 20)
    b3 = YT / K3 * (YT ** 2 * (58.6028 * YT + 198 * L - 157 * G3 - 1593 * G2 / 16 - 2437 * G1 / 80)
                    + L * YT * (15 * L / 4 + 16 * G3 - 135 * G2 / 2 - 127 * G1 / 10)
                    + YT * (363.764 * G3 ** 2 + 16.990 * G2 ** 2 - 24.422 * G1 ** 2 + 48.370 * G3 * G2
                            + 18.074 * G3 * G1 + 34.829 * G2 * G1)
                    + L ** 2 * (-36 * L + 45 * G2 + 9 * G1)
                    + L * (-171 * G2 ** 2 / 16 - 1089 * G1 ** 2 / 400 + 117 * G2 * G1 / 40)
                    - 619.35 * G3 ** 3 + 169.829 * G2 ** 3 + 16.099 * G1 ** 3 + 73.654 * G3 ** 2 * G2
                    - 15.096 * G3 ** 2 * G1 - 21.072 * G3 * G2 ** 2 - 22.319 * G3 * G1 ** 2
                    - 321 * G3 * G2 * G1 / 20 - 4.743 * G2 ** 2 * G1 - 4.442 * G2 * G1 ** 2)
    dYT = b1 + (b2 if loops >= 2 else 0) + (b3 if loops >= 3 else 0)
    # ---- small Yukawas, 1 loop
    dYB = YB / K1 * (9 * YB / 2 + 3 * YT / 2 + YTAU - 8 * G3 - 9 * G2 / 4 - G1 / 4)
    dYTAU = YTAU / K1 * (5 * YTAU / 2 + 3 * YT + 3 * YB - 9 * G2 / 4 - 9 * G1 / 4)
    return [dG1, dG2, dG3, dYT, dL, dYB, dYTAU]


def match_B13(Mt, Mh, als, MW=80.369):
    a = (als - 0.1184) / 0.0007
    w = (MW - 80.384) / 0.014
    lam = 0.12604 + 0.00206 * (Mh - 125.15) - 0.00004 * (Mt - 173.34)
    yt = 0.93690 + 0.00556 * (Mt - 173.34) - 0.00042 * a
    g2 = 0.64779 + 0.00004 * (Mt - 173.34) + 0.00011 * w
    gY = 0.35830 + 0.00011 * (Mt - 173.34) - 0.00020 * w
    g3 = 1.1666 + 0.00314 * a - 0.00046 * (Mt - 173.34)
    return Mt, dict(gY=gY, g2=g2, g3=g3, yt=yt, lam=lam, yb=0.0156, ytau=0.0100)


def match_B15(Mt, Mh, als, MW=80.369):
    a = (als - 0.1185) / 0.0006
    w = (MW - 80.385) / 0.015
    hh = (Mh - 125.7) / 0.4
    tt = (Mt - 173.21) / 0.87
    gY = 0.35838 - 3.8e-6 * a - 2.3e-4 * w - 2.5e-6 * hh + 7.1e-5 * tt
    g2 = 0.64812 + 8.5e-7 * a + 1.2e-4 * w - 6.6e-7 * hh - 9.8e-6 * tt
    g3 = 1.16540 + 2.7e-3 * a + 8.9e-8 * w + 7.8e-8 * hh - 4.0e-5 * tt
    yt = 0.93517 - 3.6e-4 * a - 1.3e-7 * w - 8.6e-6 * hh + 5.1e-3 * tt
    lam = 0.12714 - 6.2e-6 * a - 4.2e-7 * w + 8.2e-4 * hh + 6.4e-5 * tt
    yb = 0.01706 - 5.7e-5 * a
    return 173.21, dict(gY=gY, g2=g2, g3=g3, yt=yt, lam=lam, yb=yb, ytau=0.0100)


MATCH = {'B13': match_B13, 'B15': match_B15}


def run(Mt, Mh, als, scheme='B13', loops=3, mu_max=MPL, dlam=0.0, dyt=0.0):
    mu0, c = MATCH[scheme](Mt, Mh, als)
    y0 = [5 / 3 * c['gY'] ** 2, c['g2'] ** 2, c['g3'] ** 2, (c['yt'] + dyt) ** 2, c['lam'] + dlam,
          c['yb'] ** 2, c['ytau'] ** 2]
    t0, t1 = np.log(mu0 ** 2), np.log(mu_max ** 2)
    sol = solve_ivp(betas, (t0, t1), y0, method='DOP853', rtol=1e-10, atol=1e-13, dense_output=True,
                    args=(loops,))
    return sol


def lam_of_mu(sol, mu):
    return sol.sol(np.log(np.asarray(mu) ** 2))[4]


def min_lambda(sol, mu_lo=1e8, mu_hi=MPL, c6=0.0):
    """min over mu in [mu_lo, mu_hi] of lambda(mu) + c6 mu^2/(2 Mbar^2)  (c6 |H|^6/Mbar^2 operator)."""
    grid = np.logspace(np.log10(mu_lo), np.log10(mu_hi), 600)
    f = lambda m: lam_of_mu(sol, m) + c6 * m ** 2 / (2 * MPLR ** 2)
    vals = f(grid)
    i = int(np.argmin(vals))
    if i in (0, len(grid) - 1):
        return float(vals[i]), float(grid[i])
    r = minimize_scalar(lambda x: f(10 ** x), bracket=(np.log10(grid[i - 1]), np.log10(grid[i]),
                                                       np.log10(grid[i + 1])), tol=1e-10)
    return float(r.fun), float(10 ** r.x)


def mt_crit(Mh, als, scheme='B13', loops=3, mu_hi=MPL, c6=0.0, lo=166.0, hi=176.0, dlam=0.0, dyt=0.0):
    def g(Mt):
        sol = run(Mt, Mh, als, scheme, loops, mu_max=max(mu_hi, 1e19), dlam=dlam, dyt=dyt)
        return min_lambda(sol, mu_hi=mu_hi, c6=c6)[0]
    glo, ghi = g(lo), g(hi)
    if glo * ghi > 0:
        return None
    mt = brentq(g, lo, hi, xtol=1e-5)
    sol = run(mt, Mh, als, scheme, loops, mu_max=max(mu_hi, 1e19), dlam=dlam, dyt=dyt)
    return mt, min_lambda(sol, mu_hi=mu_hi, c6=c6)[1]


def lambda_zero(Mt, Mh, als, scheme='B13', loops=3):
    sol = run(Mt, Mh, als, scheme, loops)
    f = lambda x: lam_of_mu(sol, 10 ** x)
    xs = np.linspace(np.log10(Mt) + 0.1, np.log10(MPL), 400)
    v = f(xs)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(idx) == 0:
        return None
    return 10 ** brentq(f, xs[idx[0]], xs[idx[0] + 1])
