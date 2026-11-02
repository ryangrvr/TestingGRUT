#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code, written fresh, 2 November 2026): an independent relaxed-Hartree solver for the
M1r release, to test dq_entropy item 2 (S/Q(R), the W-universality, W_crit and R_crit).

Same physical model as dq_entropy/code/hartree_dq.py (copied lineage fx_chain skeptic -> wt_critic -> dq_entropy ->
dq_heat): conformal H_i units; X = (chi1 + i chi2)/sqrt 2 with r^2 = chi1^2 + chi2^2; V = lX r^4/4 + g^2 phi^2 r^2/2 +
lP phi^4/4; a circular condensate (r_c^2, frequency mu) plus Hartree quasiparticles X+, X- (mass M, chemical potential
+-mu) and phi (mass Mp, mu = 0); fixed comoving charge N and conformal energy (1 + R) E_X.

What is different (independent numerics):
  * Bose integrals by exponentially-scaled Bessel series (scipy.special.kve), not adaptive quadrature;
  * entropy from the quasiparticle identity s = (eps + P - mu n)/T per species, not from occupation numbers;
  * nested solve: (D, P) gap equations by fixed-point iteration at given (T, rc2), inside a 2-unknown hybrid root in
    (ln T, ln rc2) for (Q, E).  Not the 4-unknown joint root of hartree_dq.py.
  * an analytic massless-ideal-gas law for comparison:
      S/Q = (4/3) (3 pi^2/30)^(1/4) (3/4)^(3/4) lX^(1/4) (W - 1)^(3/4) = 0.19046 (W - 1)^(3/4),  W = (E/E_c)(1 + R).
Usage: sk_hartree_bessel.py OUT.json
Labels: [C] computed, [D] derived.
"""
import json, math, sys, time
import numpy as np
from scipy.special import kve
from scipy.optimize import brentq, root

LX, G2, LP = 1e-3, 1e-2, 3e-2
CC = 2 * math.sqrt(2 * LX) * (3 * LX) ** -0.75      # N_c = CC E^{3/4}
# M1r frozen release (dq_entropy freeze.json, base_1e5) and the H_i-frame clock (dq_entropy convention)
NX_M1R, EX_M1R = 8047.840302004927, 188737.70931829128
NX_PL, EX_PL = 10603.761966396556, 233801.82224588047
H_FREEX = 1.293e11
V_GEV = 0.924089637375362
MX = math.sqrt(G2) * V_GEV
JMAX = 200000
JALL = np.arange(1, JMAX + 1, dtype=float)
_WARM = {}


def bose_series(M, mu, T):
    """n, P, eps, D-tadpole (int f/omega) for one bosonic species with omega = sqrt(k^2 + M^2), chemical potential mu
    (|mu| < M), by the Bessel series; kve(nu, z) = K_nu(z) e^z, so each term carries exp(-j (M - mu)/T)."""
    x = M / T
    nj = int(min(JMAX, max(50, 45.0 * T / (M - abs(mu)) + 10)))     # truncation where exp(-j (M - |mu|)/T) < e^-45
    J = JALL[:nj]
    z = J * x
    w = np.exp(-J * (M - mu) / T)                   # e^{j mu/T} e^{-j x} (the e^{-jx} is put back from kve)
    k1, k2 = kve(1, z), kve(2, z)
    n = (M * M * T / (2 * math.pi ** 2)) * np.sum(w * k2 / J)
    P = (M * M * T * T / (2 * math.pi ** 2)) * np.sum(w * k2 / J ** 2)
    eps = (M * M * T / (2 * math.pi ** 2)) * np.sum(w / J * (3 * T * k2 / J + M * k1))
    Dt = (M * T / (2 * math.pi ** 2)) * np.sum(w * k1 / J)
    tail = w[-1]
    return n, P, eps, Dt, tail


def gaps(T, rc2, D0=None, P0=None):
    """solve D = D_X(T; M, mu), P = D_phi(T; Mp) at fixed (T, rc2) by damped fixed-point iteration."""
    D = D0 if D0 is not None else _WARM.get('D', T * T / 6.0)
    P = P0 if P0 is not None else _WARM.get('P', T * T / 12.0)
    for it in range(400):
        mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P)
        M = math.sqrt(2 * LX * (rc2 + D) + G2 * P)
        Mp = math.sqrt(G2 * (rc2 + D) + 3 * LP * P)
        _, _, _, Dp, _ = bose_series(M, mu, T)
        _, _, _, Dm, _ = bose_series(M, -mu, T)
        _, _, _, Pp, _ = bose_series(Mp, 0.0, T)
        Dn, Pn = Dp + Dm, Pp
        if abs(Dn - D) <= 1e-15 * D and abs(Pn - P) <= 1e-15 * P:
            D, P = Dn, Pn
            break
        D, P = Dn, Pn
    _WARM['D'], _WARM['P'] = D, P
    return D, P, it


def state(T, rc2, D=None, P=None):
    D, P, it = gaps(T, rc2, D, P)
    mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P)
    M2 = 2 * LX * (rc2 + D) + G2 * P
    Mp2 = G2 * (rc2 + D) + 3 * LP * P
    M, Mp = math.sqrt(M2), math.sqrt(Mp2)
    npl, Ppl, epl, _, t1 = bose_series(M, mu, T)
    nmi, Pmi, emi, _, t2 = bose_series(M, -mu, T)
    nph, Pph, eph, _, t3 = bose_series(Mp, 0.0, T)
    # entropy per species at fixed quasiparticle mass: s = (eps + P - mu n)/T
    S = ((epl + Ppl - mu * npl) + (emi + Pmi + mu * nmi) + (eph + Pph)) / T
    Q = mu * rc2 + npl - nmi
    d = D / 2
    pot = LX / 4 * (rc2 ** 2 + 8 * rc2 * d + 8 * d * d) + G2 / 2 * P * (rc2 + D) + 3 * LP / 4 * P * P
    E = 0.5 * mu * mu * rc2 + pot + (epl + emi - 0.5 * M2 * D) + (eph - 0.5 * Mp2 * P)
    return dict(T=T, rc2=rc2, D=D, P=P, mu=mu, M=M, Mp=Mp, n_plus=npl, n_minus=nmi, n_phi=nph, Q=Q, E=E, S=S,
                cond_frac=mu * rc2 / Q, it=it, series_tail=max(t1, t2, t3))


def solve(N, Etot, guess=None):
    """(T, rc2) from (Q = N, E = Etot) by a 2-unknown hybrid root in logs, with the (D, P) gaps nested by fixed point
    (differs from the 4-unknown joint root of hartree_dq.py); start from the massless-gas estimate."""
    if guess is None:
        Tid = (Etot / (3 * math.pi ** 2 / 30)) ** 0.25
        rc2 = (0.92 * N / math.sqrt(LX)) ** (2 / 3)
        guess = [math.log(0.97 * Tid), math.log(rc2)]
    def F(z):
        st = state(math.exp(z[0]), math.exp(z[1]))
        return [st['Q'] / N - 1, st['E'] / Etot - 1]
    r = root(F, guess, method='hybr', options=dict(xtol=1e-14))
    st = state(math.exp(r.x[0]), math.exp(r.x[1]))
    st['root_resid'] = [float(v) for v in F(r.x)]
    st['root_ok'] = bool(r.success)
    return st


def onset(st, N, H=H_FREEX):
    kap = LP * V_GEV ** 2 / H ** 2
    a = math.sqrt((G2 * (st['rc2'] + st['D']) + 3 * LP * st['P']) / kap)
    mc = MX * a / H
    return dict(a_on=a, T_h_over_mX=st['T'] / mc, nQ_mX3=N / mc ** 3)


def ideal(W):
    return (4 / 3) * (3 * math.pi ** 2 / 30) ** 0.25 * 0.75 ** 0.75 * LX ** 0.25 * (W - 1) ** 0.75


def e_over_ec(N, E):
    return E / (N / CC) ** (4 / 3)


def main():
    t0 = time.time()
    out = dict(method=__doc__, release=dict(N_X=NX_M1R, E_X=EX_M1R, E_over_Ec=e_over_ec(NX_M1R, EX_M1R)), points=[])
    w0 = e_over_ec(NX_M1R, EX_M1R)
    # dq_entropy's S/Q(R) table (Hartree) for comparison
    ref = {0.0: 1.4193, 1.0: 2.4498, 2.0: 3.3489, 3.0: 4.1731, 4.0: 4.9460, 5.0: 5.6805, 8.0: 7.7215}
    Rs = [0.0, 1.0, 2.0, 2.1442, 3.0, 3.6484, 4.0, 4.00517, 5.0, 8.0]
    for R in Rs:
        st = solve(NX_M1R, (1 + R) * EX_M1R)
        on = onset(st, NX_M1R)
        W = w0 * (1 + R)
        row = dict(R=R, W=W, S_per_Q=st['S'] / NX_M1R, Xbar_per_Q=st['n_minus'] / NX_M1R, cond_frac=st['cond_frac'],
                   T_conf=st['T'], gap_iters=st['it'], series_tail=st['series_tail'], Q_resid=st['Q'] / NX_M1R - 1,
                   E_resid=st['E'] / ((1 + R) * EX_M1R) - 1, ideal_massless=ideal(W),
                   closed_form_dq=0.190255 * (W - 0.97494) ** 0.75013 - 0.00428, dq_table=ref.get(R), **on)
        out['points'].append(row)
        print(json.dumps({k: (round(v, 7) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    # W-universality with an independent release normalisation (portals zero at M_P release, same W)
    wpl = e_over_ec(NX_PL, EX_PL)
    inv = []
    for R in (0.0, 2.1442, 4.0):
        W = w0 * (1 + R); Rp = W / wpl - 1
        a = solve(NX_M1R, (1 + R) * EX_M1R); b = solve(NX_PL, (1 + Rp) * EX_PL)
        inv.append(dict(W=W, S_M1r=a['S'] / NX_M1R, S_planck_release=b['S'] / NX_PL, rel=(b['S'] / NX_PL) / (a['S'] / NX_M1R) - 1))
    out['universality'] = inv
    print('universality', inv, flush=True)
    # W_crit at the survival thresholds (fx_chain skeptic sigma_min) and R_crit for this release
    sig = dict(bulk=3.4716191256166837, Q_2e5=2.9283778774213727, Q_3e7=3.3839721242149268, binding=2.5530624453282833)
    Wc = {}
    for k, s in sig.items():
        Rc = brentq(lambda R: solve(NX_M1R, (1 + R) * EX_M1R)['S'] / NX_M1R - s, 0.5, 4.0, xtol=1e-9)
        Wc[k] = dict(sigma=s, R_crit=Rc, W_crit=w0 * (1 + Rc),
                     W_ideal=1 + (s / ((4 / 3) * (3 * math.pi ** 2 / 30) ** 0.25 * 0.75 ** 0.75 * LX ** 0.25)) ** (4 / 3))
        print(k, Wc[k], flush=True)
    out['thresholds'] = Wc
    out['seconds'] = time.time() - t0
    json.dump(out, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
