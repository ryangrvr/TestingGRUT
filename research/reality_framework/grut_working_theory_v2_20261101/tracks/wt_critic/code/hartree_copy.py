"""fx_chain_skeptic (Claude, own code, written fresh): independent Hartree equilibrium of the relaxed pre-onset plasma.

Different numerics from fx_chain/heat_onset.py (which copied pt_relax/skeptic): adaptive scipy quad in x = k/T with
explicit Bose tails, a joint 4-unknown root solve (T, rc2, D, P) with scipy.optimize.root, and the entropy computed
from the occupation numbers.  Optionally the tachyonic term -kappa a^2 is kept inside the phi Hartree mass
self-consistently (fx_chain drops it in the state and only uses it in the onset criterion).
Conventions as fx_chain: conformal H_i units, r^2 = chi1^2 + chi2^2, V = lX r^4/4 + g^2 phi^2 r^2/2 + lP phi^4/4,
fixed comoving N_X and energy (1 + R) E_X.  Labels: [C] computed.
"""
import json, math, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import root, brentq

LX, G2, LP = 1e-3, 1e-2, 3e-2
TR = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks'
SRC = json.load(open(TR + '/kin_source/results/source_interface.json'))
FP = SRC['meta']['fixed_point']
H = FP['H_GeV']; V = FP['v_GeV']; MX = math.sqrt(G2) * V
KAPPA = LP * V * V / H ** 2
NX, EX = abs(FP['N_X']), FP['E_X']
A_REF = SRC['release_population']['onset']['a_ref_from_charge']


def bint(fun, T):
    """int_0^inf dk k^2 fun(k)/(2 pi^2), integrand given in k; split for accuracy."""
    tot = 0.0
    edges = [0, 0.05 * T, 0.5 * T, 2 * T, 8 * T, 30 * T, 120 * T]
    for lo, hi in zip(edges[:-1], edges[1:]):
        v, _ = quad(lambda k: k * k * fun(k), lo, hi, limit=200, epsabs=0, epsrel=1e-12)
        tot += v
    return tot / (2 * math.pi ** 2)


def nb(x):
    return 0.0 if x > 700 else 1.0 / math.expm1(x)


def sb(f):
    return 0.0 if f <= 0 else (1 + f) * math.log1p(f) - f * math.log(f)


def thermo(T, rc2, D, P, ka2=0.0, phi=True):
    mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P)
    M2 = 2 * LX * (rc2 + D) + G2 * P
    Mp2 = G2 * (rc2 + D) + 3 * LP * P - ka2
    w = lambda k: math.sqrt(k * k + M2)
    Dn = bint(lambda k: (nb((w(k) - mu) / T) + nb((w(k) + mu) / T)) / w(k), T)
    if phi:
        W = lambda k: math.sqrt(k * k + max(Mp2, 0.0))
        Pn = bint(lambda k: nb(W(k) / T) / W(k) if W(k) > 0 else 0.0, T)
    else:
        Pn = 0.0
    return mu, M2, Mp2, Dn, Pn


def full(T, rc2, D, P, ka2=0.0, phi=True):
    mu, M2, Mp2, Dn, Pn = thermo(T, rc2, D, P, ka2, phi)
    w = lambda k: math.sqrt(k * k + M2)
    npl = bint(lambda k: nb((w(k) - mu) / T), T); nmi = bint(lambda k: nb((w(k) + mu) / T), T)
    eX = bint(lambda k: w(k) * (nb((w(k) - mu) / T) + nb((w(k) + mu) / T)), T)
    sX = bint(lambda k: sb(nb((w(k) - mu) / T)) + sb(nb((w(k) + mu) / T)), T)
    d = D / 2
    pot = LX / 4 * (rc2 ** 2 + 8 * rc2 * d + 8 * d * d) + G2 / 2 * P * (rc2 + D) + 3 * LP / 4 * P * P - 0.5 * ka2 * P
    E = 0.5 * mu * mu * rc2 + pot + eX - 0.5 * M2 * D
    sP = 0.0; Nphi = 0.0
    if phi:
        Mp = math.sqrt(max(Mp2, 0.0)); W = lambda k: math.sqrt(k * k + Mp * Mp)
        eP = bint(lambda k: W(k) * nb(W(k) / T), T); sP = bint(lambda k: sb(nb(W(k) / T)), T)
        Nphi = bint(lambda k: nb(W(k) / T), T)
        E += eP - 0.5 * (Mp2 + ka2) * P + 0.5 * ka2 * P   # mode energy minus the double-counted mass term
    Q = mu * rc2 + npl - nmi
    return dict(T=T, rc2=rc2, D=D, P=P, mu=mu, M2=M2, Mp2=Mp2, n_plus=npl, n_minus=nmi, n_phi=Nphi, Q=Q, E=E,
                S=sX + sP, cond_frac=mu * rc2 / Q)


def solve(R, phi=True, ka2=0.0, guess=None):
    rho = EX * (1 + R)
    def F(z):
        lT, lrc, lD, lP = z
        T, rc2, D, P = map(math.exp, z)
        mu, M2, Mp2, Dn, Pn = thermo(T, rc2, D, P if phi else 0.0, ka2, phi)
        st = full(T, rc2, D, P if phi else 0.0, ka2, phi)
        r = [st['Q'] / NX - 1, st['E'] / rho - 1, Dn / D - 1, (Pn / P - 1) if phi else lP]
        return r
    g = guess or [math.log(21.9 * (1 + R) ** 0.25 / 5.97 ** 0.25), math.log(1.5e3), math.log(70.0), math.log(40.0)]
    s = root(F, g, method='hybr', options=dict(xtol=1e-13))
    T, rc2, D, P = map(math.exp, s.x)
    st = full(T, rc2, D, P if phi else 0.0, ka2, phi)
    st['resid'] = [float(v) for v in F(s.x)]
    st['success'] = bool(s.success)
    return st, list(s.x)


def onset_a(st):
    return math.sqrt((G2 * (st['rc2'] + st['D']) + 3 * LP * st['P']) / KAPPA)


def main():
    out = dict(rows=[])
    for R in (0.0, 0.9127346451110923, 1.1444541138521052, 4.966535664240221, 5.703279002838951, 41.707076999111514):
        for phi in (True, False):
            st, g = solve(R, phi)
            a_on = onset_a(st); mc = MX * a_on / H
            row = dict(R=R, phi=phi, T_conf=st['T'], Xbar_per_net=st['n_minus'] / NX, cond_frac=st['cond_frac'],
                       S_per_net=st['S'] / NX, a_on_fxchain_criterion=a_on, a_on_over_aref=a_on / A_REF,
                       T_over_mX=st['T'] / mc, nQ_over_mX3=NX / mc ** 3, resid=st['resid'])
            # self-consistent tachyonic phi mass: find a where Mphi^2 -> 0 with -kappa a^2 inside the Hartree P
            if phi:
                def gap(a):
                    s2, _ = solve(R, True, KAPPA * a * a, g)
                    return s2['Mp2'] / (KAPPA * a * a)
                try:
                    a_sc = brentq(gap, a_on * 0.999, a_on * 1.2, xtol=1e-6 * a_on)
                    s2, _ = solve(R, True, KAPPA * a_sc * a_sc, g); mcs = MX * a_sc / H
                    row.update(a_on_selfconsistent=a_sc, a_sc_over_a_on=a_sc / a_on, T_over_mX_sc=s2['T'] / mcs,
                               nQ_over_mX3_sc=NX / mcs ** 3, Xbar_per_net_sc=s2['n_minus'] / NX, V0_per_charge_sc=95.0 / (NX / mcs ** 3),
                               E_change_note='energy held at (1+R)E_X; the explicit a-dependence changes it by ~0.2%')
                except ValueError as e:
                    row['a_on_selfconsistent'] = None; row['sc_error'] = repr(e)
            out['rows'].append(row)
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items() if k != 'E_change_note'}), flush=True)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
