#!/usr/bin/env python3
"""m1_ledger (Claude, own code): the included cold X^6 charge-loss normalisation on the OLD and the M1 profile.
Item 4 of the brief.  Independent re-derivation and quadrature of the vacuum-Born pair channel; the coherent forced
channel (6e-5 of the total) is taken from Codex's numbers [Q] and not recomputed.

Operator in the light theory: (A X^6 + h.c.), A = kappa^2 mu (h/sqrt2) / (M_Y^4 (1 - r^2)) GeV^-2 (dimension-6 remnant
of the dimension-7 kappa^2 mu X^6 Sigma^*/M_Y^4 after Sigma -> h/sqrt2).  With X = (v/sqrt2) F(r) e^{-i Omega t} + chi,
X^6 -> C(6,4) X0^4 chi^2 = 15 X0^4 chi^2 creates two anti-X from four condensate quanta (core loses 4, -2 escapes).
Amplitude  <Xbar k1, Xbar k2| = 2 * 15 * A (v/sqrt2)^4 Ft4(k1+k2) * 2 pi delta(E1+E2-4 Omega) / sqrt(2E1 2E2),
Ft4(q) = int d^3x F^4 e^{-iq.x}.  Two-body phase space at fixed total momentum q and energy W = 4 Omega gives
(pi/2) beta(s), s = W^2 - q^2, beta = sqrt(1 - 4 m^2/s); identical particles 1/2.  Rate per object
  Gamma = (900 / (32 pi^3)) |A|^2 (v/sqrt2)^8 m_X^-3 * I,   I = int_0^{qmax} q^2 |Ft4(q)|^2 beta dq  (m_X units).
This reproduces Codex's prefactor 225/(8 pi^3) exactly; the quadrature here is a different substitution
(q = qmax (1 - u^2), Gauss-Legendre in u) on a different radial grid, so the integrals are an independent check.
Output: ../results/charge_loss.json
"""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/m1_ledger/data')
OUT = HERE.parent / 'results'

HBAR = 6.582119569e-25      # GeV s
MRED = 2.435e18             # GeV (reduced Planck mass used by Codex's c6 normalisation)
V = 0.924089637375362
MX = 0.1 * V
Q = 200000
AGE = 4.3549488e17          # s, 13.8 Gyr as in Codex

PROFILES = {
    'old': dict(csv=SP / 'profile_old.csv', Omega=0.8611951308217225, sha='7a737d670a8da512f3a333d22176c28dcac89bd490e858b513ad5adf9c6a8b15',
                src='d71cc01 grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv; Omega [Q FIXED_YIELD_RESULTS refinement]'),
    'M1': dict(csv=SP / 'data' / 'PROFILE_M1.csv' if False else SP / 'PROFILE_M1.csv', Omega=0.8817204919524214, sha='c08176e9be8628f99ef7c2f863e5c2d574c575c286fb82f0befe402b16fe29c7',
               src='d4747cb grut_matched_candidate_20261027/PROFILE_M1.csv; Omega [Q DEFORMATION_RESULTS primary]'),
}
# Codex's numbers [Q]
CODEX = dict(old=dict(coherent_c1=3.91958726382564e-50, born_c1=1.0399701293581432e-46, included_c1=1.0403620880845258e-46,
                      A=2.10886e-31, Omega_codex=0.861195130748691),
             M1=dict(coherent_c1=6.099550997404576e-51, born_c1=9.471343360338484e-47, included_c1=None,
                     dimensionless_integral=2868466.5109355855, events_c1=4.735671680169242e-42, A=2.079605512546138e-31,
                     included_rate_times_age=6.271629359413795e-17, relative_A_change=0.9861279158192362, uniform_scale=1.0070090495735164))
CODEX['M1']['included_c1'] = CODEX['M1']['coherent_c1'] + CODEX['M1']['born_c1']


def ft4(r, F4, q):
    """4 pi int r^2 F^4 j0(q r) dr on a uniform grid with Simpson weights (own implementation)."""
    n = len(r); h = r[1] - r[0]
    w = np.ones(n); w[1:-1:2] = 4; w[2:-1:2] = 2; w *= h / 3
    base = w * r * r * F4
    out = np.empty(len(q))
    for i in range(0, len(q), 64):
        qq = q[i:i + 64, None]
        out[i:i + 64] = 4 * np.pi * np.sum(base[None, :] * np.sinc(qq * r[None, :] / np.pi), axis=1)
    return out


def born(prof, dr=0.01, nu=240):
    a = np.loadtxt(prof['csv'], delimiter=',', skiprows=1)
    r0, F0 = a[:, 0], a[:, 1]
    spl = CubicSpline(r0, F0)
    rmax = 80.0
    n = int(round(rmax / dr)); n += (n % 2)             # even number of intervals for Simpson
    r = np.linspace(0, rmax, n + 1)
    F = spl(r); F[r > r0[-1]] = 0.0
    F4 = F ** 4
    W = 4 * prof['Omega']; qmax = np.sqrt(W * W - 4)
    u, wu = np.polynomial.legendre.leggauss(nu)
    u = (u + 1) / 2; wu = wu / 2                        # u in [0, 1]
    q = qmax * (1 - u * u); dq_du = 2 * qmax * u
    ft = ft4(r, F4, q)
    beta = np.sqrt(np.maximum(1 - 4 / (W * W - q * q), 0.0))
    I = float(np.sum(wu * q * q * ft * ft * beta * dq_du))
    pref = 900 / (32 * np.pi ** 3) * (V / np.sqrt(2)) ** 8 / (MRED ** 4 * MX ** 3)     # per object, A -> c6/M_red^2 with c6 = 1
    events = pref * I / HBAR
    frac = 4 * events / Q
    V8 = 4 * np.pi * np.trapezoid(r * r * F ** 8, r)
    lda_events = 900 / (16 * np.pi) * (V / np.sqrt(2)) ** 8 / (MRED ** 4 * MX ** 3) * V8 * np.sqrt(1 - 4 / W ** 2) / HBAR
    return dict(W_over_mX=W, qmax=qmax, I=I, events_per_second_c1=events, born_fractional_core_rate_c1=frac,
                lda_fractional_core_rate_c1=4 * lda_events / Q, F0=float(F0[0]), int_r2_F8=float(V8 / (4 * np.pi)),
                int_r2_F4_sq_check=float(ft[0] ** 2), grid=dict(dr=dr, rmax=rmax, gauss_nodes=nu))


res = {}
for key, prof in PROFILES.items():
    fine = born(prof, 0.01, 240)
    coarse = born(prof, 0.02, 120)
    fine['quadrature_relative_change'] = fine['I'] / coarse['I'] - 1
    fine['source'] = prof['src']; fine['sha256'] = prof['sha']
    res[key] = fine
    print(key, 'W', fine['W_over_mX'], 'I', fine['I'], 'born frac', fine['born_fractional_core_rate_c1'], 'lda', fine['lda_fractional_core_rate_c1'],
          'conv', fine['quadrature_relative_change'])

# ---- the chain -----------------------------------------------------------------------------------------------------
lam = CODEX['M1']['uniform_scale']
A_old = CODEX['old']['A']
A_new_pred = A_old / lam ** 2                       # [D] A propto kappa^2 mu h / M_Y^4 -> lambda^-2 under the coordinated scale change
c6_old = A_old * MRED ** 2
c6_new = CODEX['M1']['A'] * MRED ** 2
inc_old_ours = CODEX['old']['coherent_c1'] + res['old']['born_fractional_core_rate_c1']
inc_new_ours = CODEX['M1']['coherent_c1'] + res['M1']['born_fractional_core_rate_c1']
chain = dict(
    A_old_GeV_minus2=A_old, A_new_codex=CODEX['M1']['A'], A_new_predicted_lambda_minus2=A_new_pred,
    A_ratio_predicted=1 / lam ** 2, A_ratio_codex=CODEX['M1']['relative_A_change'],
    c6_old=c6_old, c6_new=c6_new,
    born_profile_ratio_ours=res['M1']['born_fractional_core_rate_c1'] / res['old']['born_fractional_core_rate_c1'],
    born_profile_ratio_codex=CODEX['M1']['born_c1'] / CODEX['old']['born_c1'],
    born_vs_codex=dict(old=res['old']['born_fractional_core_rate_c1'] / CODEX['old']['born_c1'] - 1,
                       M1=res['M1']['born_fractional_core_rate_c1'] / CODEX['M1']['born_c1'] - 1,
                       I_M1=res['M1']['I'] / CODEX['M1']['dimensionless_integral'] - 1),
    coherent_share_of_total=dict(old=CODEX['old']['coherent_c1'] / CODEX['old']['included_c1'], M1=CODEX['M1']['coherent_c1'] / CODEX['M1']['included_c1']),
    included_c1=dict(old_ours=inc_old_ours, old_codex=CODEX['old']['included_c1'], M1_ours=inc_new_ours, M1_codex=CODEX['M1']['included_c1']),
    Gamma_t0_old_ours=inc_old_ours * c6_old ** 2 * AGE, Gamma_t0_old_review_quote=7.08e-17,
    Gamma_t0_new_ours=inc_new_ours * c6_new ** 2 * AGE, Gamma_t0_new_codex=CODEX['M1']['included_rate_times_age'],
    Gamma_t0_new_from_old_times_factors=inc_old_ours * c6_old ** 2 * AGE * (inc_new_ours / inc_old_ours) * (1 / lam ** 2) ** 2,
    ratio_new_over_old_ours=(inc_new_ours * c6_new ** 2) / (inc_old_ours * c6_old ** 2),
    ratio_new_over_old_codex=CODEX['M1']['included_rate_times_age'] / 7.08e-17,
    decomposition=dict(profile_factor=inc_new_ours / inc_old_ours, coupling_factor=(1 / lam ** 2) ** 2,
                       note='profile factor = Born(M1)/Born(old) at fixed c6 (W = 4 Omega rises from 3.4448 to 3.5269 and the '
                            'F^4 source changes); coupling factor = (A_new/A_old)^2 = lambda^-4 from the coordinated heavy-scale change'),
    neutral_response_remark='the modified neutral (s) response enters only the coherent forced channel, which is 6e-5 of the total; '
                            'it cannot move Gamma t0 at the quoted precision',
)
print(json.dumps(chain, indent=1, default=float))
OUT.mkdir(exist_ok=True)
(OUT / 'charge_loss.json').write_text(json.dumps(dict(method=__doc__, profiles=res, codex_quoted=CODEX, chain=chain), indent=1, default=float))
