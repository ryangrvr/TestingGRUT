"""Independent CMB check of the HF benchmarks [C, given the track's potential].

Method (differs from the track's inflation.py):
- the potential is re-implemented from the exported parameters (sk_common), derivatives by Richardson differences;
- background in e-folds with LSODA; Mukhanov variable v_k = z zeta_k (not zeta) with z''/z built from the exact
  background including V''; Bunch-Davies at k/aH = 400, read at k/aH = 1e-3; 7 wavenumbers, quadratic fit in ln k;
- tensor modes integrated too (exact r, not 2H^2/pi^2 mP^2 / P_zeta);
- N_* recomputed from K with my own bookkeeping, and the n_s sensitivity to K (dn_s/dN_*) for the Pauli-blocking K range.
Usage: python3 -B sk_cmb.py -> ../results/sk_cmb.json
"""
import math, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import sk_common as C

mP = C.mP


SEND = 1e-4      # stop the single-field background at sigma = sigma_c (1 + SEND); the remaining e-folds are added
                 # in slow roll without the bias term (the integrated-out bias -c^2/2M^2 is invalid as M^2 -> 0)


def setup(P):
    lo = P['sc'] * (1 + 0.5 * SEND)

    def V(s):
        return P['V0'] + C.dV1(P, max(s, lo))

    def V1(s):
        return C.dV1_d(P, max(s, lo), 1)

    def V2(s):
        return C.dV1_d(P, max(s, lo), 2)
    return V, V1, V2


def N_sr(P, V, V1, s):
    grid = np.linspace(P['sc'] * (1 + SEND), s, 3001)
    f = np.array([V(z) / (mP ** 2 * V1(z)) for z in grid])
    return float(np.trapezoid(f, grid))


def N_tail(P):
    """slow-roll e-folds from sigma_c(1+SEND) to sigma_c without the bias term (loop and tree only)."""
    grid = np.linspace(P['sc'] * (1 + 1e-9), P['sc'] * (1 + SEND), 2001)
    f = np.array([(P['V0'] + C.dV1(P, z, bias=False)) / (mP ** 2 * C.dV1_d(P, z, 1, bias=False)) for z in grid])
    return float(np.trapezoid(f, grid))


def background(P, Nstar, extra=8.0):
    V, V1, V2 = setup(P)
    Nt = Nstar + extra
    hi = P['sc'] + 0.002 * mP
    while N_sr(P, V, V1, hi) < Nt:
        hi = P['sc'] + 1.5 * (hi - P['sc'])
    s0 = brentq(lambda z: N_sr(P, V, V1, z) - Nt, P['sc'] * (1 + 1e-6), hi, xtol=1e-6)
    P0 = -mP ** 2 * V1(s0) / V(s0)

    def rhs(N, y):
        s, p = y
        eps = 0.5 * p * p / mP ** 2
        return [p, -(3 - eps) * p - (3 - eps) * mP ** 2 * V1(s) / V(s)]
    ev = lambda N, y: y[0] - P['sc'] * (1 + SEND)
    ev.terminal = True
    b = solve_ivp(rhs, (0, Nt + 20), [s0, P0], method='LSODA', rtol=1e-11, atol=[1e-6, 1e-12], dense_output=True,
                  events=ev)
    Nc = b.t_events[0][0] + N_tail(P)
    return b, Nc, (V, V1, V2)


def mode_power(P, b, Vs, k, Nstart, Nend, tensor=False):
    V, V1, V2 = Vs

    def bgq(N):
        s, p = b.sol(N)
        eps = 0.5 * p * p / mP ** 2
        v, v1, v2 = V(s), V1(s), V2(s)
        H = math.sqrt(v / (mP ** 2 * (3 - eps)))
        pN = -(3 - eps) * p - (3 - eps) * mP ** 2 * v1 / v
        epsN = p * pN / mP ** 2
        pNN = (epsN * p - (3 - eps) * pN + epsN * mP ** 2 * v1 / v
               - (3 - eps) * mP ** 2 * (v2 / v - v1 * v1 / v ** 2) * p)
        f = pN / p
        Z = (1 + f) ** 2 + pNN / p - f * f + (1 - eps) * (1 + f)
        Zt = 2 - eps
        return eps, H, Z, Zt, p

    def rhs(N, y):
        eps, H, Z, Zt, p = bgq(N)
        kaH2 = (k / (math.exp(N) * H)) ** 2
        ZZ = Zt if tensor else Z
        vr, vi, dr, di = y
        return [dr, di, -(1 - eps) * dr - (kaH2 - ZZ) * vr, -(1 - eps) * di - (kaH2 - ZZ) * vi]
    eps, H, Z, Zt, p = bgq(Nstart)
    kaH = k / (math.exp(Nstart) * H)
    v0 = 1 / math.sqrt(2 * k)
    y0 = [v0, 0.0, 0.0, -kaH * v0]          # v = e^{-ik tau}/sqrt(2k): dv/dN = -i (k/aH) v
    m = solve_ivp(rhs, (Nstart, Nend), y0, method='DOP853', rtol=1e-10, atol=1e-40)
    vf = m.y[0][-1] + 1j * m.y[1][-1]
    eps, H, Z, Zt, p = bgq(Nend)
    a = math.exp(Nend)
    if tensor:
        return 2 * k ** 3 / (2 * math.pi ** 2) * 4 * abs(vf) ** 2 / (a * a * mP * mP)
    z = a * p
    return k ** 3 * abs(vf) ** 2 / (2 * math.pi ** 2 * z * z)


def cmb(P, Nstar, dks=(-2, -1, -0.5, 0, 0.5, 1, 2)):
    b, Nc, Vs = background(P, Nstar)
    Np = Nc - Nstar
    s, p = b.sol(Np)
    eps = 0.5 * p * p / mP ** 2
    H = math.sqrt(Vs[0](s) / (mP ** 2 * (3 - eps)))
    kp = math.exp(Np) * H
    rows = []
    for dk in dks:
        k = kp * math.exp(dk)
        Ns, Ne = Np + dk - math.log(400.0), Np + dk + math.log(1e3)
        rows.append((dk, mode_power(P, b, Vs, k, Ns, Ne)))
    Pt = mode_power(P, b, Vs, kp, Np - math.log(400.0), Np + math.log(1e3), tensor=True)
    lk = np.array([r[0] for r in rows])
    lP = np.log([r[1] for r in rows])
    c = np.polyfit(lk, lP, 3)
    As = rows[3][1]
    return dict(As=As, ns=1 + c[2], alpha_s=2 * c[1], r=Pt / As, H_star=H, eps_star=eps, sigma_star_MP=s / mP,
                N_total_from_start=Nc, Nstar=Nstar)


def nstar(Hstar_GeV, Hi_GeV, K, g=106.75, gs0=3.909):
    kstar = 0.05 * 1.97327e-14 / 3.08568e24
    T0 = 2.7255 * 8.617333e-14
    rho1 = 3 * K * K * C.MP_GEV ** 2 * Hi_GeV ** 2
    Ta = (30 * rho1 / (math.pi ** 2 * g)) ** 0.25
    a0 = (g / gs0) ** (1 / 3.) * Ta / T0
    return math.log(Hstar_GeV / (kstar * a0))


if __name__ == '__main__':
    out = {}
    for nm in (sys.argv[1:] or ['HFa', 'HFb']):
        d, P = C.load_bench(nm)
        K = d['reheating']['K']
        Hi = d['H_i_GeV']
        Hs = d['cmb']['ms']['H_star'] * C.H_REF
        Nst_track = d['cmb']['Nstar']
        Nst_mine = nstar(Hs, Hi, K, gs0=3.909)
        Nst_mine_3931 = nstar(Hs, Hi, K, gs0=3.931)
        r0 = cmb(P, Nst_track)
        # sensitivity to K: N_* shifts by -(1/2) ln(K'/K)
        sens = {}
        for Kp in (3.26, 5.0, 8.0):
            Nn = nstar(Hs, Hi, Kp, gs0=3.931)
            rr = cmb(P, Nn)
            sens['K=%.2f' % Kp] = dict(Nstar=Nn, ns=rr['ns'], As=rr['As'], dns=rr['ns'] - r0['ns'])
        out[nm] = dict(track=dict(ns=d['cmb']['ms']['ns'], As=d['cmb']['ms']['As'], r=d['cmb']['ms']['r'],
                                  alpha_s=d['cmb']['ms']['alpha_s'], Nstar=Nst_track, K=K),
                       mine=r0, Nstar_mine_gs0_3909=Nst_mine, Nstar_mine_gs0_3931=Nst_mine_3931,
                       pulls=dict(planck=(r0['ns'] - 0.9649) / 0.0042, act=(r0['ns'] - 0.9743) / 0.0034),
                       r_over_BK=r0['r'] / 0.036, K_sensitivity=sens)
        print(nm, out[nm])
    C.dump('sk_cmb.json', dict(label='C (skeptic, independent MS on the track potential)', note=__doc__, **out))
