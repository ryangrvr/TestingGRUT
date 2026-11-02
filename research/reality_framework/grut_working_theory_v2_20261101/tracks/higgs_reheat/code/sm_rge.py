"""Standard Model MSbar matching at the top pole mass and RG evolution to the Planck scale.

Track `vacua`, q-theory construction stage 2 (issue #13).  Own code; formulas transcribed (not copied as text)
from the sources named below, with their equation locations.

Matching at mu = M_t (NNLO fits, valid near the reference point):
  Buttazzo, Degrassi, Giardino, Giudice, Sala, Salvio, Strumia, JHEP 12 (2013) 089, arXiv:1307.3536 (TeX v4),
  section 3 ("SM couplings at the electroweak scale"):
    lambda(M_t) = 0.12604 + 0.00206 (M_h - 125.15) - 0.00004 (M_t - 173.34)        +- 0.00030 (th)
    y_t(M_t)    = 0.93690 + 0.00556 (M_t - 173.34) - 0.00042 (a3 - 0.1184)/0.0007    +- 0.00050 (th)
    g3(M_t)     = 1.1666  + 0.00314 (a3 - 0.1184)/0.0007 - 0.00046 (M_t - 173.34)
    g2(M_t)     = 0.64779 + 0.00004 (M_t - 173.34) + 0.00011 (M_W - 80.384)/0.014
    gY(M_t)     = 0.35830 + 0.00011 (M_t - 173.34) - 0.00020 (M_W - 80.384)/0.014
  (masses in GeV; a3 = alpha_s^(5)(M_Z); y_t includes the 3-loop pure-QCD term, per the source.)

RG equations: the same paper, appendix A ("SM RGE equations up to three loops"): full 3-loop beta functions of
g1 (GUT normalised, g1^2 = 5 gY^2/3), g2, g3 (plus the 4-loop pure-QCD term), y_t and lambda; 2-loop y_b, y_tau.
Written there as d(coupling^2)/d ln(mu^2) (lambda: d lambda/d ln mu^2), with V = -m^2|H|^2/2 + lambda |H|^4.
Higgs-field anomalous dimension and the Landau-gauge effective potential lambda_eff at one and two loops (the
two-loop piece in the lambda = 0 approximation): appendix B of the same paper.

Original sources of the 3-loop RGEs (credited by Buttazzo et al.): Mihaila-Salomon-Steinhauser 1201.5868;
Chetyrkin-Zoller 1205.2892, 1303.2890; Bednyakov-Pikelner-Velizhanin 1210.6873, 1212.6829, 1303.4364.
"""
import numpy as np
from scipy.integrate import solve_ivp
import mpmath as mp

K16 = 16.0 * np.pi ** 2
REF = dict(Mt=173.34, Mh=125.15, as_=0.1184, MW=80.384, dMW=0.014, das=0.0007)

# small Yukawas at mu = M_t (only enter at 1-2 loops; their effect is quantified in sm_criticality.py)
YB_MT = 0.0156   # sqrt(2) m_b(M_t)/v with m_b(M_t) ~ 2.72 GeV  [E; Buttazzo use m_b(m_b) = 4.2 GeV]
YTAU_MT = 0.0100  # sqrt(2) m_tau/v                               [E]


def matching(Mt, Mh, as_, MW=80.384, dlam=0.0, dyt=0.0, yb=YB_MT, ytau=YTAU_MT):
    """MSbar couplings at mu = M_t from the Buttazzo et al. (2013) fits. dlam, dyt: theory shifts."""
    a = (as_ - REF['as_']) / REF['das']
    w = (MW - REF['MW']) / REF['dMW']
    lam = 0.12604 + 0.00206 * (Mh - REF['Mh']) - 0.00004 * (Mt - REF['Mt']) + dlam
    yt = 0.93690 + 0.00556 * (Mt - REF['Mt']) - 0.00042 * a + dyt
    g3 = 1.1666 + 0.00314 * a - 0.00046 * (Mt - REF['Mt'])
    g2 = 0.64779 + 0.00004 * (Mt - REF['Mt']) + 0.00011 * w
    gY = 0.35830 + 0.00011 * (Mt - REF['Mt']) - 0.00020 * w
    g1 = np.sqrt(5.0 / 3.0) * gY
    return dict(g1=g1, g2=g2, g3=g3, yt=yt, yb=yb, ytau=ytau, lam=lam)


def betas(y, loops=3):
    """d/d ln(mu^2) of (G1, G2, G3, T, B, L, lam), Gi = gi^2, T = yt^2, B = yb^2, L = ytau^2."""
    G1, G2, G3, T, B, L, lam = y
    k1, k2, k3, k4 = 1.0 / K16, 1.0 / K16 ** 2, 1.0 / K16 ** 3, 1.0 / K16 ** 4
    # gauge couplings
    bG1 = G1 ** 2 * k1 * (41.0 / 10)
    bG2 = G2 ** 2 * k1 * (-19.0 / 6)
    bG3 = G3 ** 2 * k1 * (-7.0)
    if loops >= 2:
        bG1 += G1 ** 2 * k2 * (44 / 5 * G3 + 27 / 10 * G2 + 199 / 50 * G1 - 17 / 10 * T - B / 2 - 3 / 2 * L)
        bG2 += G2 ** 2 * k2 * (12 * G3 + 35 / 6 * G2 + 9 / 10 * G1 - 3 / 2 * T - 3 / 2 * B - L / 2)
        bG3 += G3 ** 2 * k2 * (-26 * G3 + 9 / 2 * G2 + 11 / 10 * G1 - 2 * T - 2 * B)
    if loops >= 3:
        bG1 += G1 ** 2 * k3 * (T * (189 / 16 * T - 29 / 5 * G3 - 471 / 32 * G2 - 2827 / 800 * G1)
                               + lam * (-9 / 5 * lam + 9 / 10 * G2 + 27 / 50 * G1)
                               + 297 / 5 * G3 ** 2 + 789 / 64 * G2 ** 2 - 388613 / 24000 * G1 ** 2
                               - 3 / 5 * G3 * G2 - 137 / 75 * G3 * G1 + 123 / 160 * G2 * G1)
        bG2 += G2 ** 2 * k3 * (T * (147 / 16 * T - 7 * G3 - 729 / 32 * G2 - 593 / 160 * G1)
                               + lam * (-3 * lam + 3 / 2 * G2 + 3 / 10 * G1)
                               + 81 * G3 ** 2 + 324953 / 1728 * G2 ** 2 - 5597 / 1600 * G1 ** 2
                               + 39 * G3 * G2 - G3 * G1 / 5 + 873 / 160 * G2 * G1)
        bG3 += G3 ** 2 * k3 * (T * (15 * T - 40 * G3 - 93 / 8 * G2 - 101 / 40 * G1)
                               + 65 / 2 * G3 ** 2 + 109 / 8 * G2 ** 2 - 523 / 120 * G1 ** 2
                               + 21 * G3 * G2 + 77 / 15 * G3 * G1 - 3 / 40 * G2 * G1)
        bG3 += G3 ** 5 * k4 * (-2472.28)          # 4-loop pure QCD (as included by the source)
    # Higgs quartic
    bl = k1 * (lam * (12 * lam + 6 * T + 6 * B + 2 * L - 9 / 2 * G2 - 9 / 10 * G1)
               - 3 * T ** 2 - 3 * B ** 2 - L ** 2 + 9 / 16 * G2 ** 2 + 27 / 400 * G1 ** 2 + 9 / 40 * G2 * G1)
    if loops >= 2:
        bl += k2 * (lam ** 2 * (-156 * lam - 72 * T - 72 * B - 24 * L + 54 * G2 + 54 / 5 * G1)
                    + lam * T * (-3 / 2 * T - 21 * B + 40 * G3 + 45 / 4 * G2 + 17 / 4 * G1)
                    + lam * B * (-3 / 2 * B + 40 * G3 + 45 / 4 * G2 + 5 / 4 * G1)
                    + lam * L * (-L / 2 + 15 / 4 * G2 + 15 / 4 * G1)
                    + lam * (-73 / 16 * G2 ** 2 + 1887 / 400 * G1 ** 2 + 117 / 40 * G2 * G1)
                    + T ** 2 * (15 * T - 3 * B - 16 * G3 - 4 / 5 * G1)
                    + T * (-9 / 8 * G2 ** 2 - 171 / 200 * G1 ** 2 + 63 / 20 * G2 * G1)
                    + B ** 2 * (-3 * T + 15 * B - 16 * G3 + 2 / 5 * G1)
                    + B * (-9 / 8 * G2 ** 2 + 9 / 40 * G1 ** 2 + 27 / 20 * G2 * G1)
                    + L ** 2 * (5 * L - 6 / 5 * G1)
                    + L * (-3 / 8 * G2 ** 2 - 9 / 8 * G1 ** 2 + 33 / 20 * G2 * G1)
                    + 305 / 32 * G2 ** 3 - 3411 / 4000 * G1 ** 3 - 289 / 160 * G2 ** 2 * G1
                    - 1677 / 800 * G2 * G1 ** 2)
    if loops >= 3:
        bl += k3 * (lam ** 3 * (6011.35 * lam + 873 * T - 387.452 * G2 - 77.490 * G1)
                    + lam ** 2 * T * (1768.26 * T + 160.77 * G3 - 359.539 * G2 - 63.869 * G1)
                    + lam ** 2 * (-790.28 * G2 ** 2 - 185.532 * G1 ** 2 - 316.64 * G2 * G1)
                    + lam * T ** 2 * (-223.382 * T - 662.866 * G3 - 5.470 * G2 - 21.015 * G1)
                    + lam * T * (356.968 * G3 ** 2 - 319.664 * G2 ** 2 - 74.8599 * G1 ** 2 + 15.1443 * G3 * G2
                                 + 17.454 * G3 * G1 + 5.615 * G2 * G1)
                    + lam * G2 ** 2 * (-57.144 * G3 + 865.483 * G2 + 79.638 * G1)
                    + lam * G1 ** 2 * (-8.381 * G3 + 61.753 * G2 + 28.168 * G1)
                    + T ** 3 * (-243.149 * T + 250.494 * G3 + 74.138 * G2 + 33.930 * G1)
                    + T ** 2 * (-50.201 * G3 ** 2 + 15.884 * G2 ** 2 + 15.948 * G1 ** 2 + 13.349 * G3 * G2
                                + 17.570 * G3 * G1 - 70.356 * G2 * G1)
                    + T * G3 * (16.464 * G2 ** 2 + 1.016 * G1 ** 2 + 11.386 * G2 * G1)
                    + T * G2 ** 2 * (62.500 * G2 + 13.041 * G1)
                    + T * G1 ** 2 * (10.627 * G2 + 11.117 * G1)
                    + G3 * (7.536 * G2 ** 3 + 0.663 * G1 ** 3 + 1.507 * G2 ** 2 * G1 + 1.105 * G2 * G1 ** 2)
                    - 114.091 * G2 ** 4 - 1.508 * G1 ** 4 - 37.889 * G2 ** 3 * G1 + 6.500 * G2 ** 2 * G1 ** 2
                    - 1.543 * G2 * G1 ** 3)
    # top Yukawa
    bT = T * k1 * (9 / 2 * T + 3 / 2 * B + L - 8 * G3 - 9 / 4 * G2 - 17 / 20 * G1)
    if loops >= 2:
        bT += T * k2 * (T * (-12 * T - 11 / 4 * B - 9 / 4 * L - 12 * lam + 36 * G3 + 225 / 16 * G2 + 393 / 80 * G1)
                        + B * (-B / 4 + 5 / 4 * L + 4 * G3 + 99 / 16 * G2 + 7 / 80 * G1)
                        + L * (-9 / 4 * L + 15 / 8 * G2 + 15 / 8 * G1)
                        + 6 * lam ** 2 - 108 * G3 ** 2 - 23 / 4 * G2 ** 2 + 1187 / 600 * G1 ** 2
                        + 9 * G3 * G2 + 19 / 15 * G3 * G1 - 9 / 20 * G2 * G1)
    if loops >= 3:
        bT += T * k3 * (T ** 2 * (58.6028 * T + 198 * lam - 157 * G3 - 1593 / 16 * G2 - 2437 / 80 * G1)
                        + lam * T * (15 / 4 * lam + 16 * G3 - 135 / 2 * G2 - 127 / 10 * G1)
                        + T * (363.764 * G3 ** 2 + 16.990 * G2 ** 2 - 24.422 * G1 ** 2 + 48.370 * G3 * G2
                               + 18.074 * G3 * G1 + 34.829 * G2 * G1)
                        + lam ** 2 * (-36 * lam + 45 * G2 + 9 * G1)
                        + lam * (-171 / 16 * G2 ** 2 - 1089 / 400 * G1 ** 2 + 117 / 40 * G2 * G1)
                        - 619.35 * G3 ** 3 + 169.829 * G2 ** 3 + 16.099 * G1 ** 3 + 73.654 * G3 ** 2 * G2
                        - 15.096 * G3 ** 2 * G1 - 21.072 * G3 * G2 ** 2 - 22.319 * G3 * G1 ** 2
                        - 321 / 20 * G3 * G2 * G1 - 4.743 * G2 ** 2 * G1 - 4.442 * G2 * G1 ** 2)
    # bottom and tau Yukawas (two loops)
    bB = B * k1 * (3 / 2 * T + 9 / 2 * B + L - 8 * G3 - 9 / 4 * G2 - G1 / 4)
    bL = L * k1 * (3 * T + 3 * B + 5 / 2 * L - 9 / 4 * G2 - 9 / 4 * G1)
    if loops >= 2:
        bB += B * k2 * (T * (-T / 4 - 11 / 4 * B + 5 / 4 * L + 4 * G3 + 99 / 16 * G2 + 91 / 80 * G1)
                        + B * (-12 * B - 9 / 4 * L - 12 * lam + 36 * G3 + 225 / 16 * G2 + 237 / 80 * G1)
                        + L * (-9 / 4 * L + 15 / 8 * G2 + 15 / 8 * G1)
                        + 6 * lam ** 2 - 108 * G3 ** 2 - 23 / 4 * G2 ** 2 - 127 / 600 * G1 ** 2
                        + 9 * G3 * G2 + 31 / 15 * G3 * G1 - 27 / 20 * G2 * G1)
        bL += L * k2 * (6 * lam ** 2 - 23 / 4 * G2 ** 2 + 1371 / 200 * G1 ** 2 + 27 / 20 * G2 * G1
                        + T * (-27 / 4 * T + 3 / 2 * B - 27 / 4 * L + 20 * G3 + 45 / 8 * G2 + 17 / 8 * G1)
                        + B * (-27 / 4 * B - 27 / 4 * L + 20 * G3 + 45 / 8 * G2 + 5 / 8 * G1)
                        + L * (-3 * L - 12 * lam + 165 / 16 * G2 + 537 / 80 * G1))
    return np.array([bG1, bG2, bG3, bT, bB, bL, bl])


def gamma_h(y, loops=3):
    """Higgs-field anomalous dimension, d ln(phi)/d ln mu (Landau gauge), Buttazzo et al. App. B."""
    G1, G2, G3, T, B, L, lam = y
    g = (9 / 4 * G2 + 9 / 20 * G1 - 3 * T - 3 * B - L) / K16
    if loops >= 2:
        g += (T * (-3 / 2 * B - 17 / 8 * G1 - 45 / 8 * G2 - 20 * G3 + 27 / 4 * T)
              - L * (15 / 8 * G1 + 15 / 8 * G2 - 9 / 4 * L)
              + B * (-5 / 8 * G1 - 45 / 8 * G2 - 20 * G3 + 27 / 4 * B)
              - 1293 / 800 * G1 ** 2 - 27 / 80 * G2 * G1 + 271 / 32 * G2 ** 2 - 6 * lam ** 2) / K16 ** 2
    if loops >= 3:
        g += (-9 * G1 * lam ** 2 - 45 * G2 * lam ** 2 + 1.07 * G1 ** 2 * lam + 3.57 * G2 * G1 * lam
              + 8.92 * G2 ** 2 * lam + 14.99 * G1 ** 2 * T + 14.13 * G1 * T ** 2 - 13.21 * G2 * G1 * T
              - 8.73 * G3 * G1 * T + 40.11 * G2 * T ** 2 + 79.05 * G3 * T ** 2 + 23.40 * G2 ** 2 * T
              - 178.48 * G3 ** 2 * T - 7.57 * G2 * G3 * T - 5.26 * G1 ** 3 + 1.93 * G2 * G1 ** 2
              + 4.19 * G3 * G1 ** 2 + 1.81 * G2 ** 2 * G1 - 158.51 * G2 ** 3 + 28.57 * G2 ** 2 * G3
              + 36 * lam ** 3 + 135 / 2 * lam ** 2 * T - 45.00 * lam * T ** 2 - 60.13 * T ** 3) / K16 ** 3
    return g


def to_state(c):
    return np.array([c['g1'] ** 2, c['g2'] ** 2, c['g3'] ** 2, c['yt'] ** 2, c['yb'] ** 2, c['ytau'] ** 2,
                     c['lam'], 0.0])


def rhs(t, Y, loops):
    """t = ln(mu/M_t). Y = (G1..L, lam, Gamma). d/d ln mu = 2 d/d ln mu^2; dGamma/d ln mu = gamma."""
    y = Y[:7]
    return np.concatenate([2.0 * betas(y, loops), [gamma_h(y, loops)]])


def run(c0, Mt, mu_max=1.0e20, loops=3, n=600, rtol=1e-10, atol=1e-13):
    """Integrate from mu = M_t to mu_max. Returns dict with arrays over ln(mu)."""
    t1 = np.log(mu_max / Mt)
    ts = np.linspace(0.0, t1, n)
    sol = solve_ivp(rhs, (0.0, t1), to_state(c0), args=(loops,), method='DOP853', t_eval=ts,
                    rtol=rtol, atol=atol, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return dict(t=sol.t, mu=Mt * np.exp(sol.t), Y=sol.y, sol=sol, Mt=Mt, loops=loops)


def couplings_at(res, mu):
    Y = res['sol'].sol(np.log(mu / res['Mt']))
    return dict(g1=np.sqrt(Y[0]), g2=np.sqrt(Y[1]), g3=np.sqrt(Y[2]), yt=np.sqrt(Y[3]), lam=Y[6], Gamma=Y[7])


def beta_lambda_lnmu(Y, loops=3):
    return 2.0 * betas(Y[:7], loops)[6]


# ---------------------------------------------------------------- effective potential (Landau gauge, h >> v)
def _xi(z):
    z = mp.mpc(z)
    s = mp.sqrt(z * z - 4 * z)
    w = (z - s) / (2 * z)
    return complex(s * (2 * mp.log(w) ** 2 - mp.log(z) ** 2 - 4 * mp.polylog(2, w) + mp.pi ** 2 / 3))


def lambda_eff(Y, two_loop=True):
    """lambda_eff(h) at mu = h (Buttazzo et al. App. B), without the overall e^{4 Gamma} (positive) factor."""
    G1, G2, G3, T, B, L, lam, Gam = Y
    gY2 = 3.0 / 5.0 * G1
    g22 = G2
    yt2 = T
    rW = np.log(g22 / 4) + 2 * Gam
    rZ = np.log((g22 + gY2) / 4) + 2 * Gam
    rt = np.log(yt2 / 2) + 2 * Gam
    l1 = (3 * g22 ** 2 / 8 * (rW - 5 / 6) + 3 / 16 * (g22 + gY2) ** 2 * (rZ - 5 / 6)
          - 3 * yt2 ** 2 * (rt - 3 / 2)
          + (3 * lam ** 2 * (4 * np.log(abs(lam)) - 6 + 3 * np.log(3) + 8 * Gam) if lam != 0 else 0.0)) / K16
    if not two_loop:
        return lam + l1
    g32 = G3
    rtW = (rt - rW) * (np.log(abs(yt2 / 2 - g22 / 4)) + 2 * Gam)   # real part if yt^2/2 < g2^2/4
    pi2 = np.pi ** 2
    Li = float(mp.re(mp.polylog(2, g22 / (2 * yt2))))
    l2 = (8 * g32 * yt2 ** 2 * (3 * rt ** 2 - 8 * rt + 9)
          + 0.5 * yt2 ** 3 * (-6 * rt * rW - 3 * rt ** 2 + 48 * rt - 6 * rtW - 69 - pi2)
          + 3 * yt2 * g22 ** 2 / 16 * (8 * rW + 4 * rZ - 3 * rt ** 2 - 6 * rt * rZ - 12 * rt + 12 * rtW + 15 + 2 * pi2)
          + yt2 * gY2 ** 2 / 48 * (27 * rt ** 2 - 54 * rt * rZ - 68 * rt - 28 * rZ + 189)
          + yt2 * g22 * gY2 / 8 * (9 * rt ** 2 - 18 * rt * rZ + 4 * rt + 44 * rZ - 57)
          + g22 ** 3 / 192 * (36 * rt * rZ + 54 * rt ** 2 - 414 * rW * rZ + 69 * rW ** 2 + 1264 * rW + 156 * rZ ** 2
                              + 632 * rZ - 144 * rtW - 2067 + 90 * pi2)
          + g22 ** 2 * gY2 / 192 * (12 * rt * rZ - 6 * rt ** 2 - 6 * rW * (53 * rZ + 50) + 213 * rW ** 2
                                    + 4 * rZ * (57 * rZ - 91) + 817 + 46 * pi2)
          + g22 * gY2 ** 2 / 576 * (132 * rt * rZ - 66 * rt ** 2 + 306 * rW * rZ - 153 * rW ** 2 - 36 * rW
                                    + 924 * rZ ** 2 - 4080 * rZ + 4359 + 218 * pi2)
          + gY2 ** 3 / 576 * (6 * rZ * (34 * rt + 3 * rW - 470) - 102 * rt ** 2 - 9 * rW ** 2 + 708 * rZ ** 2
                              + 2883 + 206 * pi2)
          + yt2 ** 2 / 6 * (4 * gY2 * (3 * rt ** 2 - 8 * rt + 9) - 9 * g22 * (rt - rW + 1))
          + 3 / 4 * (g22 ** 3 - 3 * g22 ** 2 * yt2 + 4 * yt2 ** 3) * Li
          + yt2 / 48 * _xi((g22 + gY2) / (2 * yt2)).real
          * (9 * g22 ** 2 - 6 * g22 * gY2 + 17 * gY2 ** 2 + 2 * yt2 * (7 * gY2 - 73 * g22 + 64 * g22 ** 2 / (gY2 + g22)))
          + g22 / 64 * _xi((g22 + gY2) / g22).real
          * (18 * g22 * gY2 + gY2 ** 2 - 51 * g22 ** 2 - 48 * g22 ** 3 / (gY2 + g22))) / K16 ** 2
    return lam + l1 + l2
