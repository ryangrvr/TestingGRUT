#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code; issue #13 inflaton-sector pass).  Independent of the track's code paths
(no import of inf_common / calibrate / ch_common): own background in e-folds, own Mukhanov-Sasaki scalar and tensor
modes, own e-fold bookkeeping for N_* from the exact late radiation constant K = a^2 H, own decay fluid, own
oscillation-averaged x0 = rho a^3.

Checks
  T  (TB1 alpha = 0.3, TB2 alpha = 1, S1 slice = one e-fold before eps_H = 1, Gamma_I = Gamma_cal from the track):
     x0, m_I/H_i, eps_H(S1), C = K Gamma^{1/3}/x0^{2/3}, K, T_RH (sudden and from K), N_*, H_i = H(S1), H_*, n_s, r, alpha_s.
  K  (HT6t, T1 with Gamma_I = 0.5 H_i): K = a^2 H late, a = 1 at S1, with the Gamma phi' friction from S1 (engine
     convention) and from eps_H = 1 (inf_potential's Python convention).
Units: M_P (reduced) = 1, V0 = 1 for the shape, H_i := H(S1).  Constants: k_* = 0.05/Mpc, T_0 = 2.7255 K, g = 106.75,
g_s0 = 3.931, A_s = 2.100e-9.
Output: ../results/sk_cmb_fluid.json.  Labels: [C] own code on the stated model choices.
"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MP = 2.435e18                       # GeV, reduced Planck mass
AS = 2.100e-9
G_STAR, GS0 = 106.75, 3.931
T0 = 2.7255 * 8.617333262e-14       # GeV
MPC_INV = 1.973269804e-14 / 3.0856775814913673e24   # GeV
KSTAR = 0.05 * MPC_INV              # GeV
OUT = {}


# ------------------------------------------------------------------------------------------------ potentials
class TModel:
    def __init__(self, alpha):
        self.M = math.sqrt(6 * alpha); self.alpha = alpha
    def V(self, f):  return math.tanh(f / self.M) ** 2
    def dV(self, f):
        y = f / self.M; c = math.cosh(y); return 2 * math.tanh(y) / (c * c) / self.M
    def m2(self):    return 2 / self.M ** 2
    def phi_start(self, N=80.0):        # slow-roll N = (M^2/8) cosh(2 phi/M)
        return 0.5 * self.M * math.acosh(8 * N / self.M ** 2)
    pmin = 0.0


class Hilltop:
    def __init__(self, vev, p, beta):
        self.v, self.p, self.b = vev, p, beta
        self.xtop = -(beta / p) ** (1.0 / (p - 1))
        self.ftop = 1 - beta * self.xtop - self.xtop ** p
        self.xmin = brentq(lambda x: 1 - beta * x - x ** p, 0.5, 1.0)
        self.pmin = self.xmin * vev
    def V(self, f):
        x = f / self.v; return (1 - self.b * x - x ** self.p) ** 2 / self.ftop ** 2
    def dV(self, f):
        x = f / self.v; fx = 1 - self.b * x - x ** self.p
        return 2 * fx * (-self.b - self.p * x ** (self.p - 1)) / self.v / self.ftop ** 2
    def m2(self):
        x = self.xmin; fp = -self.b - self.p * x ** (self.p - 1)
        return 2 * fp * fp / self.v ** 2 / self.ftop ** 2


# ------------------------------------------------------------------------------------------------ background in N
def bg_rhs(N, y, m):
    f, fN = y
    eps = 0.5 * fN * fN
    return [fN, -(3 - eps) * fN - (3 - eps) * m.dV(f) / m.V(f)]


def background(m, f0):
    """slow-roll start at f0; integrate in N to eps_H = 1; dense output."""
    fN0 = -m.dV(f0) / m.V(f0)                       # slow-roll attractor (M_P = 1)
    ev = lambda N, y, m: 0.5 * y[1] ** 2 - 1.0
    ev.terminal = True; ev.direction = 1
    s = solve_ivp(bg_rhs, (0, 400), [f0, fN0], args=(m,), method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True,
                  events=ev)
    Nend = s.t_events[0][0]
    return s, Nend


def H_of(m, f, fN):
    eps = 0.5 * fN * fN
    return math.sqrt(m.V(f) / (3 - eps)), eps


# ------------------------------------------------------------------------------------------------ MS modes
def mode_power(m, sol, Nend, k, aS1_N):
    """P_zeta and P_t for comoving k (units: a = 1 at N = aS1_N, H in shape units).  BD at k/aH = 300."""
    # find start N: k/(aH) = 300
    def kaH(N):
        f, fN = sol.sol(N); H, _ = H_of(m, f, fN); return k / (math.exp(N - aS1_N) * H)
    Na = brentq(lambda N: kaH(N) - 300.0, sol.t[0] + 1e-6, Nend)
    Nb = brentq(lambda N: kaH(N) - 1e-3, Na, Nend + 0.0 - 1e-9) if kaH(Nend - 1e-9) < 1e-3 else Nend - 1e-6
    f, fN = sol.sol(Na); H, eps = H_of(m, f, fN); a = math.exp(Na - aS1_N)
    epsN_over_eps = -6 - 2 * m.dV(f) * (3 - eps) / (fN * m.V(f)) + 2 * eps
    z = a * math.sqrt(2 * eps)
    v = 1 / math.sqrt(2 * k)
    zeta = v / z
    zetaN = zeta * (-1j * k / (a * H) - (1 + 0.5 * epsN_over_eps))
    h = 2 * v / a
    hN = h * (-1j * k / (a * H) - 1)

    def rhs(N, y):
        f, fN = y[0], y[1]
        eps = 0.5 * fN * fN
        Vv = m.V(f); dV = m.dV(f)
        H = math.sqrt(Vv / (3 - eps)); a = math.exp(N - aS1_N)
        q = (k / (a * H)) ** 2
        cz = -3 + eps - 2 * dV * (3 - eps) / (fN * Vv)
        ch = 3 - eps
        zr, zi, zNr, zNi, hr, hi, hNr, hNi = y[2:]
        return [fN, -(3 - eps) * fN - (3 - eps) * dV / Vv,
                zNr, zNi, -cz * zNr - q * zr, -cz * zNi - q * zi,
                hNr, hNi, -ch * hNr - q * hr, -ch * hNi - q * hi]
    y0 = [f, fN, zeta.real, zeta.imag, zetaN.real, zetaN.imag, h.real, h.imag, hN.real, hN.imag]
    s = solve_ivp(rhs, (Na, Nb), y0, method='DOP853', rtol=1e-10, atol=1e-30)
    y = s.y[:, -1]
    Pz = k ** 3 * (y[2] ** 2 + y[3] ** 2) / (2 * math.pi ** 2)
    Pt = 2 * k ** 3 * (y[6] ** 2 + y[7] ** 2) / (2 * math.pi ** 2)
    return Pz, Pt


def spectrum(m, sol, Nend, Nstar):
    """n_s, r, alpha_s at the mode exiting Nstar e-folds before eps_H = 1; P_zeta in shape units (V0 = 1)."""
    aS1_N = Nend - 1.0
    Np = Nend - Nstar
    f, fN = sol.sol(Np); H, _ = H_of(m, f, fN)
    kst = math.exp(Np - aS1_N) * H
    lk = np.array([-0.4, -0.2, 0.0, 0.2, 0.4])
    P = [mode_power(m, sol, Nend, kst * math.exp(x), aS1_N) for x in lk]
    lPz = np.log([p[0] for p in P]); lPt = np.log([p[1] for p in P])
    cz = np.polyfit(lk, lPz, 2); ct = np.polyfit(lk, lPt, 2)
    Pz = math.exp(cz[2]); Pt = math.exp(ct[2])
    return dict(Pz_shape=Pz, r=Pt / Pz, ns=1 + cz[1], alpha_s=2 * cz[0], nt=ct[1], h_star=H, kstar_shape=kst)


# ------------------------------------------------------------------------------------------------ post-inflation
def oscillate(m, sol, Nend, gam_from=None, gam=0.0, t_end=None, ncyc=60):
    """Field from S1 in cosmic time (H_i = H(S1) = 1 units, a = 1 at S1) with optional Gamma phi' friction switched on at
    'S1' or 'end'; radiation fluid sourced by Gamma phidot^2.  Returns x0 (cycle-averaged rho_I a^3, Gamma = 0) or K."""
    NS1 = Nend - 1.0
    f, fN = sol.sol(NS1); Hs, eps = H_of(m, f, fN)
    # rescale: time in 1/H_S1, rho in M_P^2 H_S1^2: V -> V/Hs^2
    sc = 1.0 / Hs ** 2
    fd = fN * Hs / Hs                               # dphi/dt in units H_S1 M_P: phidot = H fN -> /H_S1
    fd = fN * 1.0                                   # H(S1) = 1 after rescaling: phidot = H fN = fN
    mI = math.sqrt(m.m2() * sc)
    tend_eps = None

    def rhs(t, y, g_on):
        f, p, lna, rr = y
        rho = 0.5 * p * p + sc * m.V(f) + rr
        H = math.sqrt(rho / 3)
        g = gam if g_on else 0.0
        return [p, -(3 * H + g) * p - sc * m.dV(f), H, g * p * p - 4 * H * rr]
    y0 = [f, fd, 0.0, 0.0]
    if gam_from == 'end' or gam == 0.0:
        ev = lambda t, y, g: 0.5 * y[1] ** 2 / ((0.5 * y[1] ** 2 + sc * m.V(y[0]) + y[3]) / 3) - 1.0
        ev.terminal = True; ev.direction = 1
        s1 = solve_ivp(rhs, (0, 50), y0, args=(False,), method='DOP853', rtol=1e-11, atol=1e-14, events=ev)
        tend_eps = s1.t_events[0][0]; y0 = list(s1.y[:, -1]); t0 = tend_eps
    else:
        t0 = 0.0
    if gam == 0.0:
        T = 2 * math.pi / mI
        tE = t0 + ncyc * T * 4
        s2 = solve_ivp(rhs, (t0, tE), y0, args=(False,), method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True)
        tt = np.linspace(tE - 8 * T, tE, 4001)
        Y = s2.sol(tt)
        x = (0.5 * Y[1] ** 2 + sc * np.array([m.V(v) for v in Y[0]])) * np.exp(3 * Y[2])
        tt2 = np.linspace(tE - 40 * T, tE - 32 * T, 4001); Y2 = s2.sol(tt2)
        x2 = (0.5 * Y2[1] ** 2 + sc * np.array([m.V(v) for v in Y2[0]])) * np.exp(3 * Y2[2])
        return dict(x0=float(np.mean(x)), x0_earlier=float(np.mean(x2)), m_I_over_Hi=mI, eps_S1=eps, t_end_eps1=tend_eps,
                    a_end=math.exp(s1.y[2, -1]))
    # with decay: integrate until rho_I/rho_r < 1e-12
    tE = t0 + 40.0 / gam
    s2 = solve_ivp(rhs, (t0, tE), y0, args=(True,), method='DOP853', rtol=1e-11, atol=1e-16)
    f, p, lna, rr = s2.y[:, -1]
    rhoI = 0.5 * p * p + sc * m.V(f)
    K = math.sqrt(rr * math.exp(4 * lna) / 3)
    return dict(K=K, rhoI_over_rhor=rhoI / rr, m_I_over_Hi=mI, t_end_eps1=tend_eps, gam_from=gam_from)


def fluid_K(x0, gam):
    """matter (rho a^3 -> x0 at early times) decaying into radiation; K = a^2 H late (a = 1 at S1, H_i = 1)."""
    a0 = (1e-3 * gam) ** (-2.0 / 3.0) * 1e-2     # start deep in matter domination, H >> Gamma
    a0 = min(a0, 1e2)
    def rhs(N, y):
        rm, rr = math.exp(y[0]), math.exp(y[1])
        H = math.sqrt((rm + rr) / 3)
        return [-3 - gam / H, -4 + gam * rm / (H * rr)]
    rm0 = x0 / a0 ** 3
    H0 = math.sqrt(rm0 / 3)
    rr0 = rm0 * 0.4 * gam / H0                      # quasi-static radiation seed rho_r = (2/5) Gamma rho_m / H
    N0 = math.log(a0)
    Nf = N0 + 2 / 3 * math.log(H0 / gam) + 20
    s = solve_ivp(rhs, (N0, Nf), [math.log(rm0), math.log(rr0)], method='LSODA', rtol=1e-11, atol=1e-13)
    rr = math.exp(s.y[1, -1]); N = s.t[-1]
    K = math.sqrt(rr * math.exp(4 * N) / 3)
    return K


def nstar_loop(m, sol, Nend, K_of_none, Kdim):
    """iterate N_* with V0 from A_s (exact MS) and the exact late radiation constant K (dimensionless)."""
    Ns = 52.5
    for it in range(8):
        sp = spectrum(m, sol, Nend, Ns)
        V0 = AS / sp['Pz_shape']                         # M_P^4
        f, fN = sol.sol(Nend - 1.0); hS1, _ = H_of(m, f, fN)
        Hi = math.sqrt(V0) * hS1 * MP
        Hst = math.sqrt(V0) * sp['h_star'] * MP
        rr_a4 = 3 * Kdim ** 2 * MP ** 2 * Hi ** 2       # GeV^4 (a = 1 at S1)
        Ta = (30 * rr_a4 / (math.pi ** 2 * G_STAR)) ** 0.25
        a0 = Ta * (G_STAR / GS0) ** (1 / 3) / T0      # a_today in S1 units
        Nn = math.log(math.e * Hst / (KSTAR * a0))
        if abs(Nn - Ns) < 1e-5:
            Ns = Nn; break
        Ns = Nn
    sp = spectrum(m, sol, Nend, Ns)
    V0 = AS / sp['Pz_shape']
    f, fN = sol.sol(Nend - 1.0); hS1, epsS1 = H_of(m, f, fN)
    Hi = math.sqrt(V0) * hS1 * MP; Hst = math.sqrt(V0) * sp['h_star'] * MP
    return dict(Nstar=Ns, H_i_GeV=Hi, H_star_GeV=Hst, Hstar_over_Hi=Hst / Hi, ns=sp['ns'], r=sp['r'], alpha_s=sp['alpha_s'],
                r_from_Hstar=2 * Hst ** 2 / (math.pi ** 2 * AS * MP ** 2), V0=V0, eps_H_S1=epsS1, iterations=it + 1)


def run_T():
    track = {'TB1': dict(alpha=0.3, gam=7.512248636886893e-09, x0=21.8506048158377, mI=2.207748554153159, eps=0.13469877444547865,
                         K=2000.4722069213312, Nstar=52.34459311737356, Hi=7466455428038.676, Hst=8989265703919.133,
                         ns=0.961850295628048, r=0.001315031850875881, TRH=199789092385.92355, NMD=12.220747969015472),
             'TB2': dict(alpha=1.0, gam=2.1126537975537015e-09, x0=17.062791831474147, mI=1.541287456486854, eps=0.2271371355723044,
                         K=2589.224491651785, Nstar=52.60718086004074, Hi=10604555120818.209, Hst=16112365581655.734,
                         ns=0.9621531821789163, r=0.004224373813520002, TRH=126268040613.06866, NMD=13.000587492618338)}
    for key, tr in track.items():
        m = TModel(tr['alpha'])
        sol, Nend = background(m, m.phi_start(80.0))
        osc = oscillate(m, sol, Nend)
        K = fluid_K(osc['x0'], tr['gam'])
        C = K * tr['gam'] ** (1 / 3) / osc['x0'] ** (2 / 3)
        ns = nstar_loop(m, sol, Nend, None, K)
        Hi = ns['H_i_GeV']
        gG = tr['gam'] * Hi
        TRH_sudden = (90 / (math.pi ** 2 * G_STAR)) ** 0.25 * math.sqrt(gG * MP)
        # T a from K and the 'H = Gamma' epoch in a matter era (a = 1 at S1)
        rr_a4 = 3 * K ** 2 * MP ** 2 * Hi ** 2
        Ta = (30 * rr_a4 / (math.pi ** 2 * G_STAR)) ** 0.25
        a_HG = (math.sqrt(osc['x0'] / 3) / tr['gam']) ** (2 / 3)
        T_at_HG = Ta / a_HG
        N_MD = math.log(a_HG / osc['a_end'])
        y_maj = math.sqrt(16 * math.pi * tr['gam'] / osc['m_I_over_Hi'])
        row = dict(alpha=tr['alpha'], Gamma_over_Hi=tr['gam'], x0=osc['x0'], x0_earlier_cycles=osc['x0_earlier'], m_I_over_Hi=osc['m_I_over_Hi'],
                   eps_H_S1_osc=osc['eps_S1'], a_end_over_aS1=osc['a_end'], C_K=C, K=K, T_RH_sudden_GeV=TRH_sudden, T_at_H_eq_Gamma_GeV=T_at_HG,
                   N_MD_end_to_H_eq_Gamma=N_MD, y_majorana=y_maj, **ns, track=tr)
        row['rel'] = dict(x0=osc['x0'] / tr['x0'] - 1, K=K / tr['K'] - 1, Nstar=ns['Nstar'] - tr['Nstar'], Hi=Hi / tr['Hi'] - 1,
                          Hstar=ns['H_star_GeV'] / tr['Hst'] - 1, ns=ns['ns'] - tr['ns'], r=ns['r'] / tr['r'] - 1,
                          TRH=TRH_sudden / tr['TRH'] - 1, mI=osc['m_I_over_Hi'] / tr['mI'] - 1)
        OUT[key] = row
        print(key, json.dumps({k: (v if not isinstance(v, dict) else None) for k, v in row.items()}, default=float))
        print('   rel:', row['rel'], flush=True)


def run_K():
    benches = {'HT6t': (Hilltop(0.3469421539393091, 6, 1.5529238776548362e-05), 0.03, 12.340784016982916, 1.5163631380626654),
               'T1': (TModel(6.0967783372331265e-05), None, 15.200283662768815, 1.8538702665196651)}
    for key, (m, f0, Ktr, Kend_ip) in benches.items():
        if f0 is None:
            f0 = m.phi_start(70.0)
        sol, Nend = background(m, f0)
        f, fN = sol.sol(Nend - 1.0)
        rS1 = oscillate(m, sol, Nend, gam_from='S1', gam=0.5)
        rE = oscillate(m, sol, Nend, gam_from='end', gam=0.5)
        row = dict(phi_S1=f, K_S1_switch=rS1['K'], K_end_switch=rE['K'], K_end_switch_endnorm=rE['K'] / math.e ** 2,
                   ratio=rS1['K'] / rE['K'], m_I_over_Hi=rS1['m_I_over_Hi'], track_engine_K=Ktr, inf_potential_K_endnorm=Kend_ip,
                   rel_engine=rS1['K'] / Ktr - 1, rel_ip=(rE['K'] / math.e ** 2) / Kend_ip - 1, resid=[rS1['rhoI_over_rhor'], rE['rhoI_over_rhor']])
        OUT['K_' + key] = row
        print(key, row, flush=True)


if __name__ == '__main__':
    what = sys.argv[1:] or ['T', 'K']
    if 'T' in what:
        run_T()
    if 'K' in what:
        run_K()
    import pathlib
    p = pathlib.Path(__file__).resolve().parent.parent / 'results' / 'sk_cmb_fluid.json'
    old = json.load(open(p)) if p.exists() else {}
    old.update(OUT)
    json.dump(old, open(p, 'w'), indent=1, default=float)
