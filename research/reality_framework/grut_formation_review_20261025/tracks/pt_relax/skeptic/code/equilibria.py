"""pt_relax skeptic (Claude, own code): independent end points and onset.  Imports nothing from Codex or pt_relax.

1. Classical (Rayleigh-Jeans) self-consistent Hartree equilibrium ON the FD4 lattice (N = 32, 40, 48, 64; L = 12)
   at the measured conformal (nQ, rho), circular condensate, X + phi equipartition.  Tests pt_relax's harmonic
   estimate 0.41-0.61 (fixed zero-T mass and mu, interaction energy dropped, thermal charge not removed).
2. Quantum Hartree equilibrium in the continuum (vacuum parts dropped), different numerics from pt_relax:
   fixed Gauss-Legendre grid in x = k/T, nested 1-D root finds (T from energy, r_c^2 from charge).
3. Onset a where kappa a^2 = g^2 <r^2> + 3 lambda_phi <phi^2> (pt_relax drops the phi self term).
4. Massive-phase equilibria after the onset at fixed (nQ, rho + V0): X, Xbar of mass m, charge reservoir at
   mu_c (= omega = 0.861195 m, the case-L cap, or = m, a free BEC), optional radial scalar s (m_s = sqrt(6) m).
   Pilot (nQ, rho) and the heavy-spot release (N_X, E_X (1 + R)).

Conventions (pt_relax / Codex): psi = chi1 + i chi2 with canonical real components, so r^2 = chi1^2 + chi2^2,
X mass from g^2 phi^2 r^2/2 is g phi, phi mass from it is g r.  Conformal H_i units.
"""
import json, math, os
import numpy as np
from scipy.optimize import brentq

LX, G2, LP = 1e-3, 1e-2, 3e-2
AM, C = 8.869218738806, 1.9857962750634526
H = 133349619398.33566; V = 0.924089637375362; MX = math.sqrt(G2) * V
KAPPA = LP * V * V / H ** 2
NQ, RHO = 8112.665383318339, 105481.24746711778
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')

# ------------------------------------------------------------------ 1. classical lattice equilibrium
def fd4_k2(N, L=12.0):
    i = np.fft.fftfreq(N, 1.0 / N); s2 = np.sin(np.pi * i / N) ** 2
    f = 4 * s2 * (1 + s2 / 3) / (L / N) ** 2
    k2 = (f[:, None, None] + f[None, :, None] + f[None, None, :]).ravel()[1:]   # drop k = 0 (condensate)
    u, c = np.unique(np.round(k2, 9), return_counts=True)                        # degenerate eigenvalues
    return u, c.astype(float)


def rj_state(T, rc2, k2c, V3, it=400):
    """classical Hartree at temperature T (energy per oscillator), condensate rc2. Returns dict."""
    k2, cw = k2c
    D, P = 0.0, 0.0
    for _ in range(it):
        mu2 = LX * (rc2 + 2 * D) + G2 * P; mu = math.sqrt(mu2)
        MX2 = 2 * LX * (rc2 + D) + G2 * P
        w = np.sqrt(k2 + MX2); Wp = np.sqrt(k2 + G2 * (rc2 + D) + 3 * LP * P)
        npl = T / (w - mu); nmi = T / (w + mu); nph = T / Wp
        Dn = np.sum(cw * (npl + nmi) / w) / V3; Pn = np.sum(cw * nph / Wp) / V3
        if abs(Dn - D) < 1e-12 * (1 + D) and abs(Pn - P) < 1e-12 * (1 + P): D, P = Dn, Pn; break
        D, P = 0.5 * D + 0.5 * Dn, 0.5 * P + 0.5 * Pn
    mu2 = LX * (rc2 + 2 * D) + G2 * P; mu = math.sqrt(mu2); MX2 = 2 * LX * (rc2 + D) + G2 * P
    Mp2 = G2 * (rc2 + D) + 3 * LP * P
    w = np.sqrt(k2 + MX2); Wp = np.sqrt(k2 + Mp2)
    npl = T / (w - mu); nmi = T / (w + mu); nph = T / Wp
    d = D / 2
    pot = LX / 4 * (rc2 ** 2 + 8 * rc2 * d + 8 * d * d) + G2 / 2 * P * (rc2 + D) + 3 * LP / 4 * P * P
    E = 0.5 * mu2 * rc2 + pot + np.sum(cw * w * (npl + nmi)) / V3 - 0.5 * MX2 * D + np.sum(cw * Wp * nph) / V3 - 0.5 * Mp2 * P
    Q = mu * rc2 + np.sum(cw * (npl - nmi)) / V3
    return dict(T_cl=T, rc2=rc2, mu=mu, MX2=MX2, Mphi2=Mp2, D=D, P=P, charge=Q, energy=E,
                Xbar_action=np.sum(cw * nmi) / V3, Nphi=np.sum(cw * nph) / V3,
                phi_energy_frac=(np.sum(cw * Wp * nph) / V3 - 0.5 * Mp2 * P + 3 * LP / 4 * P * P + G2 / 2 * P * (rc2 + D)) / E)


def rj_solve(N, nQ=NQ, rho=RHO, L=12.0):
    k2 = fd4_k2(N, L); V3 = L ** 3
    def rc2_for(T):
        f = lambda x: rj_state(T, math.exp(x), k2, V3)['charge'] / nQ - 1
        return math.exp(brentq(f, math.log(1e-3), math.log(2e5), xtol=1e-13))
    g = lambda lt: rj_state(math.exp(lt), rc2_for(math.exp(lt)), k2, V3)['energy'] / rho - 1
    sc = (32.0 / N) ** 3
    lt = brentq(g, math.log(50.0 * sc), math.log(1900.0 * sc), xtol=1e-12)
    st = rj_state(math.exp(lt), rc2_for(math.exp(lt)), k2, V3)
    st.update(N=N, L=L, Xbar_per_net=st['Xbar_action'] / nQ, Nphi_per_net=st['Nphi'] / nQ,
              condensate_charge_frac=st['mu'] * st['rc2'] / nQ)
    return st

# ------------------------------------------------------------------ 2. quantum Hartree, continuum
XG, WG = np.polynomial.legendre.leggauss(400)
def _grid(T, Mmax):
    kmax = 80 * T + 20 * Mmax
    k = 0.5 * kmax * (XG + 1); wk = 0.5 * kmax * WG
    return k, wk

def bose(x):
    return np.where(x > 700, 0.0, 1.0 / np.expm1(np.minimum(x, 700)))

def q_state(T, rc2, phi=True, it=500):
    D, P = 0.0, 0.0
    for _ in range(it):
        mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P); M2 = 2 * LX * (rc2 + D) + G2 * P
        Mp2 = G2 * (rc2 + D) + 3 * LP * P
        k, wk = _grid(T, math.sqrt(max(M2, Mp2)))
        w = np.sqrt(k * k + M2); fp = bose((w - mu) / T); fm = bose((w + mu) / T)
        Dn = np.sum(wk * k * k * (fp + fm) / w) / (2 * math.pi ** 2)
        if phi:
            W = np.sqrt(k * k + Mp2); fph = bose(W / T); Pn = np.sum(wk * k * k * fph / W) / (2 * math.pi ** 2)
        else:
            Pn = 0.0
        if abs(Dn - D) < 1e-13 * (1 + D) and abs(Pn - P) < 1e-13 * (1 + P): D, P = Dn, Pn; break
        D, P = 0.5 * (D + Dn), 0.5 * (P + Pn)
    mu = math.sqrt(LX * (rc2 + 2 * D) + G2 * P); M2 = 2 * LX * (rc2 + D) + G2 * P; Mp2 = G2 * (rc2 + D) + 3 * LP * P
    k, wk = _grid(T, math.sqrt(max(M2, Mp2))); c = wk * k * k / (2 * math.pi ** 2)
    w = np.sqrt(k * k + M2); fp = bose((w - mu) / T); fm = bose((w + mu) / T)
    npl, nmi = np.sum(c * fp), np.sum(c * fm)
    d = D / 2
    pot = LX / 4 * (rc2 ** 2 + 8 * rc2 * d + 8 * d * d) + G2 / 2 * P * (rc2 + D) + 3 * LP / 4 * P * P
    E = 0.5 * mu * mu * rc2 + pot + np.sum(c * w * (fp + fm)) - 0.5 * M2 * D
    Ephi = 0.0
    if phi:
        W = np.sqrt(k * k + Mp2); Ephi = np.sum(c * W * bose(W / T)) - 0.5 * Mp2 * P
        E += Ephi
    return dict(T=T, rc2=rc2, mu=mu, M_X=math.sqrt(M2), M_phi=math.sqrt(Mp2) if phi else None, D=D, P=P,
                n_plus=npl, n_minus=nmi, charge=mu * rc2 + npl - nmi, energy=E, Ephi_frac=Ephi / E)

def q_solve(nQ, rho, phi=True):
    def rc2_for(T):
        f = lambda x: q_state(T, math.exp(x), phi)['charge'] / nQ - 1
        return math.exp(brentq(f, math.log(1e-2), math.log(1e6), xtol=1e-14))
    g = lambda lt: q_state(math.exp(lt), rc2_for(math.exp(lt)), phi)['energy'] / rho - 1
    lt = brentq(g, math.log(5.0), math.log(60.0), xtol=1e-13)
    st = q_state(math.exp(lt), rc2_for(math.exp(lt)), phi)
    st.update(Xbar_per_net=st['n_minus'] / nQ, r2_total=st['rc2'] + st['D'], phi_included=phi)
    return st

# ------------------------------------------------------------------ 4. massive phase after the onset (free species)
def n_e(T, m, mu):
    k, wk = _grid(T, m); w = np.sqrt(k * k + m * m); f = bose((w - mu) / T); c = wk * k * k / (2 * math.pi ** 2)
    return np.sum(c * f), np.sum(c * w * f)

def massive_eq(nQ, rho, m, mu_c, with_s):
    """charge reservoir at energy mu_c per unit charge; gas X (mu = mu_c), Xbar (-mu_c); optional s (sqrt6 m, mu 0)."""
    ms = math.sqrt(6.0) * m
    def parts(T):
        npl, epl = n_e(T, m, mu_c * (1 - 1e-12)); nmi, emi = n_e(T, m, -mu_c)
        ns, es = n_e(T, ms, 0.0) if with_s else (0.0, 0.0)
        nc = nQ - (npl - nmi)
        return nc, mu_c * nc + epl + emi + es, nmi, ns, es
    g = lambda lt: parts(math.exp(lt))[1] / rho - 1
    lt = brentq(g, math.log(1e-3 * m), math.log(1e4 * m), xtol=1e-14)
    T = math.exp(lt); nc, E, nmi, ns, es = parts(T)
    return dict(T=T, T_over_m=T / m, condensed_frac=nc / nQ, Xbar_per_net=nmi / nQ, s_per_net=ns / nQ,
                s_energy_frac=es / rho, valid=bool(nc >= 0))


def main():
    res = {}
    # 1. classical lattice end points
    res['rj_lattice'] = []
    for N in [32, 40, 48, 64]:
        st = rj_solve(N); res['rj_lattice'].append(st)
        print('RJ N%d T_cl=%.1f Xbar_action/net=%.4f phi_frac=%.3f Nphi/net=%.3f cond_charge=%.3f MX2=%.2f Mphi2=%.1f mu=%.3f'
              % (N, st['T_cl'], st['Xbar_per_net'], st['phi_energy_frac'], st['Nphi_per_net'], st['condensate_charge_frac'],
                 st['MX2'], st['Mphi2'], st['mu']))
    # 2. quantum Hartree (pilot)
    qa = q_solve(NQ, RHO, True); qb = q_solve(NQ, RHO, False)
    res['quantum_hartree'] = dict(phi_equilibrium=qa, phi_absent=qb)
    for nm, st in [('phi eq', qa), ('phi absent', qb)]:
        print('QH', nm, 'T=%.4f Xbar/net=%.5f rc2=%.2f mu=%.4f M_X=%.4f Ephi=%.4f D=%.3f P=%.3f' %
              (st['T'], st['Xbar_per_net'], st['rc2'], st['mu'], st['M_X'], st['Ephi_frac'], st['D'], st['P']))
    # 3. onset with and without the phi self term
    on = []
    for nm, r2, P in [('brief |chi| = 60', 3600.0, 0.0), ('brief |chi| = 110', 12100.0, 0.0),
                      ('quantum eq, g^2 r^2 only', qa['rc2'] + qa['D'], 0.0),
                      ('quantum eq, + 3 lambda_phi <phi^2>', qa['rc2'] + qa['D'], qa['P'])]:
        a = math.sqrt((G2 * r2 + 3 * LP * P) / KAPPA); s = (a - AM) / C
        on.append(dict(case=nm, a_onset=a, s_onset=s, Gamma_crit=1 / s, m_conf=MX * a / H))
        print('ONSET', nm, 'a=%.4e s=%.4e Gcrit=%.3e mconf=%.3f' % (a, s, 1 / s, MX * a / H))
    res['onset'] = on
    # 4. massive phase at the onset (conformal units: m_conf = m_X a/H; V0_conf = lambda v^4/4 (a/H)^4)
    rows = []
    a_eq = on[2]['a_onset']
    NXs, EXs = 2788.1767940228406, 38386.25063336385
    cases = [('pilot', NQ, RHO, a_eq, 0.0), ('pilot', NQ, RHO, 7.04671e12, 0.0)]
    for R in [0.0, 1.14445, 4.0538, 40.56262, 48.49912]:
        cases.append(('heavy spot', NXs, EXs * (1 + R), 7.04671e12, R))
    for nm, n, rho, a, R in cases:
        m = MX * a / H; V0 = LP * V ** 4 / 4 * (a / H) ** 4
        for with_V0 in [False, True]:
            for mu_c in [0.861195 * m, m]:
                for with_s in [False, True]:
                    st = massive_eq(n, rho + (V0 if with_V0 else 0.0), m, mu_c, with_s)
                    st.update(case=nm, R=R, a=a, m_conf=m, V0_conf=V0, V0_included=with_V0, mu_c_over_m=mu_c / m, with_s=with_s)
                    rows.append(st)
                    print('MASSIVE %-10s R=%7.3f a=%.3e V0=%d mu_c/m=%.3f s=%d T/m=%.3f cond=%.3f Xbar/net=%.4f s/net=%.3f valid=%s' %
                          (nm, R, a, with_V0, mu_c / m, with_s, st['T_over_m'], st['condensed_frac'], st['Xbar_per_net'], st['s_per_net'], st['valid']))
    res['massive_after_onset'] = rows
    json.dump(res, open(os.path.join(OUT, 'equilibria.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
