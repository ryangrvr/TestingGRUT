#!/usr/bin/env python3
"""m1_holder_skeptic (Claude, adversarial reviewer; own code). Small extra checks:
  (a) the 94f660d sextic holder point re-solved with scipy (own V, own guess; not the track's holder_point.py) and its Hessian;
  (b) effective suppression scales: Option 1A  Lambda_eff = M_Y/xi;  sextic  Lambda_eff = Lambda/xi6^{1/5}; the tree-level
      Y-exchange sextic kappa^2 u/M_Y^4 already present in the renormalisable model;
  (c) summary of the slow N_X <-> N_Sigma exchange from skeptic_evolve.json (partner identification, periods, amplitudes).
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, root

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
H_GEV = 134284273491.30815
LX, LY, LS = 1e-3, 1.0, 1.0
KAP = 0.0034586847211042523
M, HV, MU = 4028528204739.244 / H_GEV, 5371370939652.326 / H_GEV, 961402368973.181 / H_GEV
XIX = -5.0 / 6.0
LOCK = (10j - MU) / 12
LAMBDA6_GEV, XI6 = 1e15, 4e8
out = {}

# (a) sextic holder of 94f660d: V += 2 Re[i c6 Sigma* X^6], c6 = xi6 R/Lambda^5 at R = 12, lock = 0
lam6 = LAMBDA6_GEV / H_GEV
c6 = XI6 * 12.0 / lam6 ** 5


def V6(r):
    X = (r[0] + 1j * r[1]) / np.sqrt(2); Y = (r[2] + 1j * r[3]) / np.sqrt(2); S = (r[4] + 1j * r[5]) / np.sqrt(2)
    R = 12.0; u = MU
    return float(LX * abs(X) ** 4 + LY * abs(Y) ** 4 + M * M * abs(Y) ** 2 + LS * (abs(S) ** 2 - HV * HV / 2) ** 2
                 + 2 * np.real(KAP * np.conj(Y) * X ** 3 + u * np.conj(S) * Y * Y) + XIX * R * abs(X) ** 2
                 + 2 * np.real(1j * c6 * np.conj(S) * X ** 6))


def grad6(r, h=1e-5):
    g = np.zeros(6)
    for i in range(6):
        e = np.zeros(6); e[i] = h
        g[i] = (-V6(r + 2 * e) + 8 * V6(r + e) - 8 * V6(r - e) + V6(r - 2 * e)) / (12 * h)
    return g


def hess6(r, h=2e-3):
    Hm = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            ei = np.zeros(6); ei[i] = h; ej = np.zeros(6); ej[j] = h
            f = lambda a, b: V6(r + a * ei + b * ej)
            if i == j:
                d = (-f(2, 0) + 16 * f(1, 0) - 30 * f(0, 0) + 16 * f(-1, 0) - f(-2, 0)) / (12 * h * h)
            else:
                d1 = (f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1)) / (4 * h * h)
                d2 = (f(2, 2) - f(2, -2) - f(-2, 2) + f(-2, -2)) / (16 * h * h)
                d = (4 * d1 - d2) / 3
            Hm[i, j] = Hm[j, i] = d
    return Hm


# own guess: X real ~ sqrt(5000), Y ~ -kap x^3/M^2, Sigma phase from the two competing holders (mu wants Y^2 Sigma* real negative;
# sextic wants i Sigma* X^6 real negative -> Sigma ~ -i s * e^{i small})
x0 = np.sqrt(5000.0); y0 = -KAP * x0 ** 3 / M ** 2
r0 = np.sqrt(2) * np.array([x0, 0.0, y0, 0.0, 0.0, -HV / np.sqrt(2)])
r0 = r0 + np.sqrt(2) * np.array([0, 0, 0, 0.3, 2.0, 0])            # kick off the symmetric axis
m = minimize(V6, r0, method='BFGS', options=dict(gtol=1e-10, maxiter=5000))
sol = root(grad6, m.x, method='lm', options=dict(xtol=1e-14, ftol=1e-14, maxiter=20000))
r6 = sol.x
X = (r6[0] + 1j * r6[1]) / np.sqrt(2); Y = (r6[2] + 1j * r6[3]) / np.sqrt(2); S = (r6[4] + 1j * r6[5]) / np.sqrt(2)
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi
ev6 = np.linalg.eigvalsh(hess6(r6))
out['sextic_94f660d'] = dict(c6_over_H5=c6, amplitudes_own=[abs(X), abs(Y), abs(S)], delta_own=float(wrap(np.angle(S) - 6 * np.angle(X))),
                              delta_y_own=float(wrap(np.angle(Y) - 3 * np.angle(X))), V_own=V6(r6), grad_norm=float(np.linalg.norm(grad6(r6))),
                              eigenvalues_fd=ev6.tolist(), lightest_physical_own=float(np.sort(ev6)[1]),
                              codex=dict(amps=[86.657, 3.823, 28.338], delta=-2.09993, lightest=10.669), track=dict(amps=[86.6572, 3.8228, 28.3376], delta=-2.099927, lightest=10.66925),
                              method='scipy BFGS + LM on the finite-difference gradient of the own V; Hessian 4th-order FD (expect 1e-6 accuracy)')

# (b) effective suppression scales
xi = abs(LOCK) * M
out['effective_scales'] = dict(
    option1A_Lambda_eff_over_H=float(M / xi), option1A_Lambda_eff_GeV=float(M / xi * H_GEV),
    sextic_Lambda_eff_over_H=float(lam6 / XI6 ** 0.2), sextic_Lambda_eff_GeV=float(LAMBDA6_GEV / XI6 ** 0.2),
    field_amplitudes_over_H=[87.19556, 5.75425, 28.38718], sqrt_R_over_H=float(np.sqrt(12.0)),
    ratios_to_option1A_scale=dict(Y_over_Lambda_eff=float(5.75425 * xi / M), Sigma_over_Lambda_eff=float(28.38718 * xi / M), X_over_Lambda_eff=float(87.19556 * xi / M),
                                  R_over_Lambda_eff2=float(12 * (xi / M) ** 2)),
    y_exchange_sextic_coefficient_kappa2_u_over_M4_at_R0=float(KAP ** 2 * MU / M ** 4), sextic_holder_c6=c6,
    ratio_c6_holder_to_Y_exchange=float(c6 / (KAP ** 2 * MU / M ** 4)),
    comment='[E] Option 1A with xi = 30.75 is a dimension-5 operator with effective scale M_Y/xi = 0.976 H_i, below every field amplitude and below sqrt(R); '
            'the abandoned sextic had Lambda/xi6^(1/5) = 142 H_i. The renormalisable model already contains a Y-exchange sextic kappa^2 u Sigma* X^6/M_Y^4 '
            'of half the 94f660d c6, with phase arg u(R): this is the channel of the late X <-> Sigma charge exchange.')

# (c) slow exchange summary
ev = json.load(open(RES / 'skeptic_evolve.json'))
for tag in ('original', 'holder'):
    r = ev[tag]['running_80_means']
    t0 = np.array([q['t0'] for q in r]); nx = np.array([q['NX'] for q in r]); ns = np.array([q['NS'] for q in r]); ny = np.array([q['NY'] for q in r])
    # extrema of nx for a crude period
    ext = [float(t0[i] + 40) for i in range(1, len(nx) - 1) if (nx[i] - nx[i - 1]) * (nx[i + 1] - nx[i]) < 0]
    out[f'slow_exchange_{tag}'] = dict(window_centres_of_extrema=ext, NX_pp_rel=ev[tag]['slow_modulation']['NX_peak_to_peak_rel'],
                                        corr_NX_6NS=ev[tag]['slow_modulation']['corr_NX_6NS'], corr_NX_3NY=ev[tag]['slow_modulation']['corr_NX_3NY'],
                                        max_abs_dev_NX_plus_6NS_minus_mean=float(np.max(np.abs((nx + 6 * ns) - (nx + 6 * ns).mean()))),
                                        max_abs_dev_3NY=float(np.max(np.abs(3 * ny - (3 * ny).mean()))),
                                        NS_mean=float(ns.mean()), NX_mean=float(nx.mean()),
                                        late_NY_fast_std=ev[tag]['late_instantaneous']['NY_std'], late_absX_max=ev[tag]['late_instantaneous']['absX_max'])
    # quartic-oscillator orbital period estimate at the late amplitude: T = 4*1.311/(x_m sqrt(2 lx))
    xm = ev[tag]['late_instantaneous']['absX_max']
    out[f'slow_exchange_{tag}']['X_orbit_period_estimate_Hi_t'] = float(4 * 1.31103 / (xm * np.sqrt(2 * LX)))
json.dump(out, open(RES / 'skeptic_extra.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
