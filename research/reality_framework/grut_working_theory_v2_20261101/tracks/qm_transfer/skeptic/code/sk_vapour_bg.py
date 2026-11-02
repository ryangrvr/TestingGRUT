"""qm_transfer SKEPTIC: independent re-implementation of the vapour kinetic background (own code, skeptic pass).
Different numerics from qt_hidden.py: own cosmology (rho_c from G), 160-bin comoving-q grid spanning 1e-2..60 q_th
(vs 24 bins), first-order exponential integrator with T_o re-solved every step (vs BDF), non-relativistic X.
Same physics as the track (so that the numbers can be compared): two-phase equilibrium adiabat N_v(5/2 + D/T) = S/Q
with n_sat = (T/2pi)^{3/2} e^{-D/T}; from Gamma_a(vbar)/H = 60 the per-charge distribution N(q) obeys
dN/dlna = (Gamma_a(q)/H)(N_eq(q;T_o) - N), Gamma_a = zeta K n_Q C_o^{2/3} v, K = pi (3/4 pi n_B)^{2/3}, with T_o from
sum Gamma_a (D + E_k)(N - N_eq) = 0 ('thermo'), or exchange switched off at Gamma_a(actual)/H = 1 ('sudden').
Outputs sigma_p = sqrt(<q^2>/3) today, f_v, a_f, the vapour/T_o ratio at Gamma/H = 2, and the naive estimate
sigma_naive = a_f,eq sqrt(T_f,eq) (freeze-out at Gamma(vbar_eq) = H, no lag, no re-emission).  [C]
usage: sk_vapour_bg.py OUT.json"""
import json, math, sys
import numpy as np
from scipy.optimize import brentq

# --- constants and cosmology (Planck 2018 base values as in the track; rho_c from G)
MX = 0.0924089637                         # GeV
MPL = 1.220890e19                         # GeV
h, omb, omc, Neff, T0K = 0.6736, 0.02237, 0.1200, 3.046, 2.7255
H0 = 100 * h / 2.99792458e5 * 299792.458 / 3.0856775814913673e19 * 6.582119569e-25    # km/s/Mpc -> 1/s -> GeV
RHOC = 3 * H0**2 * MPL**2 / (8 * math.pi)                                               # GeV^4
Og = 2.4728e-5 / h**2 * (T0K / 2.7255)**4
Or = Og * (1 + Neff * 7 / 8 * (4 / 11)**(4 / 3))
Om = (omb + omc) / h**2
OL = 1 - Om - Or
NB = 449.31336788
K = math.pi * (3 / (4 * math.pi * NB))**(2 / 3)


def H_mX(a):
    return H0 * math.sqrt(Or / a**4 + Om / a**3 + OL) / MX


def Omega_fit(Q):
    return 0.775865 + (2 / 3) * 8.1429 * Q**(-1 / 3) + (1 / 3) * 125.79 * Q**(-2 / 3)


def model(SQ, zeta, s=1.0, mode='thermo', nq=160, dl=0.002, qhi_fac=1800.0, g_stop=1e-7):
    Q = (s / zeta)**3
    D = 1 - Omega_fit(Q)
    rho_dm = omc / h**2 * RHOC

    def eq(a, NQ0):
        nQ = NQ0 / a**3
        f = lambda lT: 1.5 * (lT - math.log(2 * math.pi)) - D / math.exp(lT) - math.log(nQ) - math.log(SQ / (2.5 + D / math.exp(lT)))
        lT = brentq(f, math.log(1e-5), math.log(0.3))
        T = math.exp(lT); return T, SQ / (2.5 + D / T)

    def Gbar_eq(a, NQ0):
        T, Nv = eq(a, NQ0)
        return zeta * K * NQ0 / a**3 * (1 - Nv)**(2 / 3) * math.sqrt(8 * T / math.pi) / H_mX(a), T, Nv

    # n_Q: rho_dm = n_Q m_X [(1 - N_v)(1 - D) + N_v]; iterate on the late vapour fraction
    Nv_late = 0.03
    for it in range(3):
        NQ0 = rho_dm / (MX**4 * ((1 - Nv_late) * (1 - D) + Nv_late))
        la_f = brentq(lambda la: math.log(Gbar_eq(math.exp(la), NQ0)[0]), math.log(1e-11), math.log(1e-2))
        a_f_eq = math.exp(la_f); _, T_f_eq, Nv_f_eq = Gbar_eq(a_f_eq, NQ0)
        Nv_late = Nv_f_eq
    la0 = brentq(lambda la: math.log(Gbar_eq(math.exp(la), NQ0)[0] / 60.0), la_f - 8.0, la_f)
    a = math.exp(la0); T, Nv = eq(a, NQ0)
    qth = a * math.sqrt(T)
    q = np.exp(np.linspace(math.log(1e-2 * qth), math.log(qhi_fac * a_f_eq * math.sqrt(T_f_eq)), nq))
    w = np.gradient(q)                                      # simple trapezoid-like weights on the log grid
    pref = 4 * math.pi * q * q / ((2 * math.pi)**3 * NQ0)

    def Neq(a, T):
        return pref * np.exp(-(D + q * q / (2 * a * a)) / T)
    N = Neq(a, T); N *= Nv / float(w @ N)
    To = T; frozen = False; a_f_kin = None; lag = None; hist = []
    la = la0
    while la < 0.0:
        a = math.exp(la + 0.5 * dl)
        if not frozen:
            Co = 1 - float(w @ N)
            v = q / a
            G = zeta * K * NQ0 / a**3 * Co**(2 / 3) * v / H_mX(a)            # Gamma_a / H per bin
            Ek = q * q / (2 * a * a)
            B = float(w @ (G * (D + Ek) * N))
            lpre = np.log(pref * w * G * (D + Ek) + 1e-300)

            def F(lT):
                z = lpre - (D + Ek) / math.exp(lT); zm = z.max()
                return zm + math.log(np.exp(z - zm).sum()) - math.log(B)
            To = math.exp(brentq(F, math.log(1e-5), math.log(0.5)))
            Ne = Neq(a, To)
            N = Ne + (N - Ne) * np.exp(-G * dl)
            Nv = float(w @ N)
            Gbar = float(w @ (G * N)) / Nv
            Tv = float(w @ (N * q * q)) / (3 * Nv * a * a)
            if lag is None and Gbar < 2.0:
                lag = dict(a=a, Tv_over_To=Tv / To, Gbar=Gbar)
            if a_f_kin is None and Gbar < 1.0:
                a_f_kin = a; T_f_kin_o = To; Tv_f = Tv
                if mode == 'sudden':
                    frozen = True
            if mode == 'thermo' and Gbar < g_stop:
                frozen = True
        la += dl
    Nv = float(w @ N)
    sig = math.sqrt(float(w @ (N * q * q)) / Nv / 3)
    fv_mass = Nv / (Nv + (1 - Nv) * (1 - D))
    qmean = float(w @ (N * q)) / Nv
    cdf = np.cumsum(w * N) / float(w @ N)
    q50, q90, q99 = (float(np.interp(p, cdf, q)) for p in (0.5, 0.9, 0.99))
    return dict(SQ=SQ, zeta=zeta, Q=Q, Delta=D, mode=mode, nq=nq, dl=dl, qhi_fac=qhi_fac, g_stop=g_stop,
                q_mean=qmean, q50=q50, q90=q90, q99=q99, NQ0=NQ0, a_f_eq=a_f_eq, T_f_eq=T_f_eq, Nv_f_eq=Nv_f_eq,
                a_f_kin=a_f_kin, T_o_at_af=T_f_kin_o, Tv_at_af=Tv_f, lag_at_Gamma2=lag, f_v_charge=Nv, f_v_mass=fv_mass,
                sigma_p=sig, sigma_naive=a_f_eq * math.sqrt(T_f_eq), C=sig / (zeta * T_f_eq), C_naive=a_f_eq * math.sqrt(T_f_eq) / (zeta * T_f_eq))


if __name__ == '__main__':
    out = []
    def show(r, tag=''):
        print('%-14s SQ=%g zeta=%.0e %-6s a_f eq %.3e kin %.3e T_f eq %.4f Tv/To(G/H=2)=%.3f f_v=%.4f (mass %.4f) sigma_p=%.3e C=%.4f C_naive=%.4f  <q>=%.3e q50=%.3e q90=%.3e q99=%.3e' % (
            tag, r['SQ'], r['zeta'], r['mode'], r['a_f_eq'], r['a_f_kin'], r['T_f_eq'], r['lag_at_Gamma2']['Tv_over_To'], r['f_v_charge'],
            r['f_v_mass'], r['sigma_p'], r['C'], r['C_naive'], r['q_mean'], r['q50'], r['q90'], r['q99']), flush=True)
    for SQ, z in ((1.235, 1e-4), (1.235, 1e-3), (1.235, 1e-5), (3.4, 1e-4)):
        for mode in ('thermo', 'sudden'):
            r = model(SQ, z, mode=mode); r['variant'] = 'base'; out.append(r); show(r, 'base')
    # grid / duration dependence of the thermo rms (the track's grid ends at 25 q_f,eq; its thermo stops at Gbar = 1e-6)
    for kw, tag in ((dict(qhi_fac=25.0, g_stop=1e-6), 'track-like'), (dict(qhi_fac=100.0), 'qhi=100'), (dict(qhi_fac=400.0), 'qhi=400'),
                    (dict(qhi_fac=1800.0, g_stop=1e-9), 'stop=1e-9'), (dict(qhi_fac=8000.0, g_stop=1e-9, nq=260), 'qhi=8000'),
                    (dict(nq=320, dl=0.001), 'nq320,dl/2')):
        r = model(1.235, 1e-4, mode='thermo', **kw); r['variant'] = tag; out.append(r); show(r, tag)
    r = model(1.235, 1e-4, mode='sudden', nq=320, dl=0.001); r['variant'] = 'nq320,dl/2'; out.append(r); show(r, 'nq320,dl/2')
    json.dump(dict(scope=__doc__.strip().splitlines()[0], H0_GeV=H0, rho_c_GeV4=RHOC, runs=out), open(sys.argv[1], 'w'), indent=1)
