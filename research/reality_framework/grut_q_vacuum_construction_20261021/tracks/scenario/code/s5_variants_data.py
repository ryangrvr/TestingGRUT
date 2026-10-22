"""Part 5: variants and retrospective comparisons (all data already exposed: every comparison is retrospective).

(a) QV1 (P_q): no late vacuum term. H0 follows from the CMB physical densities and T0; distances, age, q0;
    Pantheon+ SN (shape, M marginalised), DESI DR2 BAO (shape with H0 r_d free, and absolute with r_d from omega_b,
    omega_c), Planck 100 theta_*, SH0ES H0, globular-cluster age, SH0ES+Pantheon+ q0.
(b) QV1-cal: mu = mu0[eps_inf] - rho_Lambda/q0 (calibrated offset; the calibration fixes the residual).
(c) Rejected readings: initial-equilibrium mu (AdS turnaround/crunch after the EW step); fixed-number mu with a
    heavy-threshold constant c (history moves one-for-one) versus P_q (c re-absorbed into mu0[eps_inf + c]: the
    history does not move; only G is renormalised by c/M^4).
(d) Full-equation confirmation of (c) at a compressed hierarchy (m/H_c = 100, KV shape).
Data files (downloaded to the scratchpad, never to the repository):
  DESI DR2 BAO: CobayaSampler/bao_data desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt (DESI collaboration
  release accompanying arXiv:2503.14738).
  Pantheon+: PantheonPlusSH0ES/DataRelease Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov (Scolnic et al. 2022,
  Brout et al. 2022, arXiv:2202.04077).
Run: python3 -I s5_variants_data.py [outdir] [data_dir] [eos_dir]
"""
import json
import math
import os
import sys
import time
import urllib.request

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qv1_core import KV_SHAPE, QV1, ToyMatter, make_rhs, quasi_static_toy, reduced_toy  # noqa: E402
from s4_realistic_history import (GYR, H100, HBAR, K0, M, M_4MEV, M_NAT, OMEGA_B, OMEGA_C, RHO_CRIT100,  # noqa: E402
                                  RHO_LAMBDA, T0, Background)
from sm_eos import DEFAULT_DIR as EOS_DIR, SMEoS  # noqa: E402

mp.mp.dps = 300
C_KMS = 299792.458
DH100 = C_KMS / 100.0             # Mpc, c/H100
DATA_DIR = ('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/'
            'scenario_dl')
URLS = {'desi/desi_gaussian_bao_ALL_GCcomb_mean.txt':
        'https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt',
        'desi/desi_gaussian_bao_ALL_GCcomb_cov.txt':
        'https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt',
        'pplus/Pantheon+SH0ES.dat':
        'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/'
        'Pantheon%2BSH0ES.dat',
        'pplus/cov.cov':
        'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/'
        'Pantheon%2BSH0ES_STAT%2BSYS.cov'}
# exposed data quoted at abstract level [Q]
SH0ES_H0 = (73.04, 1.04)          # Riess et al. 2022, arXiv:2112.04510 abstract
SH0ES_Q0 = (-0.51, 0.024)         # same abstract: "H0 = 73.30 +- 1.04 with q0 = -0.51 +- 0.024"
PLANCK_THETA = (1.0411, 0.0003)   # Planck 2018 VI abstract, 100 theta_*
GC_AGE = (13.5, 0.16, 0.14, 0.5)  # Valcin et al. 2020, arXiv:2007.06594 abstract: t_U = 13.5 +0.16 -0.14 (stat) +- 0.5 (sys) Gyr
DELTA_EW = 125.25 ** 2 * 246.22 ** 2 / 8.0   # GeV^4, tree level m_H^2 v^2/8 (equations track section 7) [R inputs]
DELTA_QCD = (0.145 ** 4, 0.235 ** 4)         # bag-constant range (equations track section 7) [R]


def fetch(d=DATA_DIR):
    for rel, url in URLS.items():
        p = os.path.join(d, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if not os.path.exists(p):
            urllib.request.urlretrieve(url, p)


# ------------------------------------------------------------------------------------------------ late cosmology
class Cosmo:
    """Flat FRW with omega_r (from T0 and the Saikawa-Shirai asymptotic g_rho, i.e. N_eff = 3.044 equivalent),
    omega_m (dust) and a constant vacuum term omega_L; h fixed by flatness: h^2 = omega_r + omega_m + omega_L."""

    def __init__(self, omega_m, omega_L, g_rho_late=3.3830836):
        self.omega_r = math.pi ** 2 / 30 * g_rho_late * T0 ** 4 / RHO_CRIT100
        self.omega_g = math.pi ** 2 / 15 * T0 ** 4 / RHO_CRIT100
        self.omega_m, self.omega_L = omega_m, omega_L
        self.h = math.sqrt(self.omega_r + omega_m + omega_L)

    def E100(self, z):
        zp = 1 + z
        return np.sqrt(self.omega_r * zp ** 4 + self.omega_m * zp ** 3 + self.omega_L)

    def DM(self, z):
        z = np.atleast_1d(z)
        zz = np.linspace(0, z.max(), 200001)
        f = 1 / self.E100(zz)
        cum = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(zz))])
        return DH100 * np.interp(z, zz, cum)

    def DH(self, z):
        return DH100 / self.E100(z)

    def age_Gyr(self):
        f = lambda lna: 1.0 / (self.E100(math.exp(-lna) - 1))
        val, _ = quad(f, -60, 0, limit=500)
        return val / (H100 / HBAR) / GYR

    def q0(self):
        # q = -a addot/adot^2 = (omega_r*2 + omega_m - 2 omega_L)/(2 h^2)
        return (2 * self.omega_r + self.omega_m - 2 * self.omega_L) / (2 * self.h ** 2)

    def z_star(self, ob):
        om = self.omega_m
        g1 = 0.0783 * ob ** -0.238 / (1 + 39.5 * ob ** 0.763)
        g2 = 0.560 / (1 + 21.1 * ob ** 1.81)
        return 1048 * (1 + 0.00124 * ob ** -0.738) * (1 + g1 * om ** g2)

    def rs(self, z, ob):
        R0 = 0.75 * ob / self.omega_g
        f = lambda lnzp: (math.exp(lnzp) / math.sqrt(3 * (1 + R0 / math.exp(lnzp)))) / self.E100(math.exp(lnzp) - 1)
        val, _ = quad(f, math.log(1 + z), math.log(1 + z) + 25, limit=500)
        return DH100 * val

    def theta_star_100(self, ob):
        zs = self.z_star(ob)
        zz = np.linspace(0, zs, 400001)
        f = 1 / self.E100(zz)
        DMs = DH100 * np.trapezoid(f, zz)
        return 100 * self.rs(zs, ob) / DMs, zs


def rd_aubourg(ob, ocb, onu=0.0):
    """Aubourg et al. 2015 (arXiv:1411.1074) eq. (16) fitting formula [R], standard N_eff: 0.02% accuracy."""
    return 55.154 * math.exp(-72.3 * (onu + 0.0006) ** 2) / (ocb ** 0.25351 * ob ** 0.12807)


# ------------------------------------------------------------------------------------------------ data
def load_desi(d):
    rows = [ln.split() for ln in open(os.path.join(d, 'desi/desi_gaussian_bao_ALL_GCcomb_mean.txt'))
            if not ln.startswith('#') and ln.strip()]
    z = np.array([float(r[0]) for r in rows])
    v = np.array([float(r[1]) for r in rows])
    kind = [r[2] for r in rows]
    cov = np.loadtxt(os.path.join(d, 'desi/desi_gaussian_bao_ALL_GCcomb_cov.txt'))
    return z, v, kind, cov


def bao_model(cos, z, kind):
    """Distances in units of c/H0 (i.e. D/(DH100/h))."""
    DM = cos.DM(z) * cos.h / DH100
    DH = cos.DH(z) * cos.h / DH100
    out = []
    for i, k in enumerate(kind):
        if k == 'DM_over_rs':
            out.append(DM[i])
        elif k == 'DH_over_rs':
            out.append(DH[i])
        elif k == 'DV_over_rs':
            out.append((z[i] * DM[i] ** 2 * DH[i]) ** (1 / 3))
        else:
            raise ValueError(k)
    return np.array(out)


def chi2_bao(cos, data, rd=None):
    z, v, kind, cov = data
    Ci = np.linalg.inv(cov)
    f = bao_model(cos, z, kind)
    if rd is None:          # amplitude A = c/(H0 r_d) free
        A = (f @ Ci @ v) / (f @ Ci @ f)
        r = v - A * f
        return float(r @ Ci @ r), float(A)
    A = DH100 / cos.h / rd
    r = v - A * f
    return float(r @ Ci @ r), float(A)


def load_pplus(d):
    tab = np.genfromtxt(os.path.join(d, 'pplus/Pantheon+SH0ES.dat'), names=True, dtype=None, encoding='utf-8')
    with open(os.path.join(d, 'pplus/cov.cov')) as fh:
        n = int(fh.readline())
        cov = np.array(fh.read().split(), dtype=float).reshape(n, n)
    sel = tab['zHD'] > 0.01
    Ci = np.linalg.inv(cov[np.ix_(sel, sel)])
    return dict(zHD=tab['zHD'][sel], zHEL=tab['zHEL'][sel], mb=tab['m_b_corr'][sel], Ci=Ci, n=int(sel.sum()))


def chi2_sn(cos, sn):
    DM = cos.DM(sn['zHD'])
    mu = 5 * np.log10((1 + sn['zHEL']) * DM) + 25
    r = sn['mb'] - mu
    Ci = sn['Ci']
    a = r @ Ci @ r
    b = np.sum(Ci @ r)
    c = np.sum(Ci)
    return float(a - b * b / c), float(b / c)


def flat_lcdm_shape(Om):
    """Flat LambdaCDM shape with Omega_m (radiation negligible at z < 3): use h = 0.7 nominal (shape only)."""
    h = 0.7
    om = Om * h * h
    return Cosmo(om, h * h - om - 0.0, g_rho_late=0.0)


# ------------------------------------------------------------------------------------------------ (c) realistic
def turnaround(bg, const):
    """Flat FRW with plasma + dust + constant 'const' (< 0): N, T and t at H = 0 (radiation-era start at 1 TeV)."""
    rho = bg.rho + const
    i = np.where(rho <= 0)[0]
    if len(i) == 0:
        return None
    i = i[0]
    Nt = float(np.interp(0.0, -rho[i - 1:i + 1], bg.N[i - 1:i + 1]))
    Tt = float(np.interp(Nt, bg.N[i - 1:i + 1], bg.T[i - 1:i + 1]))
    rhoN = lambda N: np.interp(N, bg.N, rho)
    # t = int dN/H; substitute N = Nt - u^2 to remove the 1/sqrt singularity
    H = lambda N: math.sqrt(max(rhoN(N), 0.0) / (3 * K0))
    N0 = bg.N[0]
    t0 = 1.0 / (2 * math.sqrt(bg.rho[0] / (3 * K0)))
    U = math.sqrt(Nt - N0)
    g = lambda u: 2 * u / H(Nt - u * u) if u > 0 else 2.0 / math.sqrt(-np.gradient(rho, bg.N)[i] / (3 * K0))
    val, err = quad(g, 0, U, limit=400)
    t_turn = t0 + val
    return dict(N_turn=Nt, T_turn_GeV=Tt, t_turn_s=t_turn * HBAR, t_crunch_s=2 * t_turn * HBAR,
                note='crunch time = 2 t_turn for a time-symmetric (reversible, entropy-conserving) recollapse [A]',
                analytic_RD_t_turn_s=(math.pi / 4) * math.sqrt(3 * K0 / abs(const)) * HBAR)


def kv_lambda(chi_c):
    """P_q re-absorption with the KV quartic: (q0'/q0)^2 = (1 + sqrt(1 + 8 chi c))/2 (exact)."""
    mp.mp.dps = 300
    return mp.sqrt((1 + mp.sqrt(1 + 8 * mp.mpf(chi_c))) / 2)


def fixed_mu_static_root(c_over_M4, m2_over_M2=mp.mpf(2) / 3):
    """Fixed numerical mu (= mu0 of eps_inf), eps -> eps + c (type i): static dS/AdS root in Y-variables (KV shape,
    K0 = M^2 = 1): m^2 F(Y) = (4c)/3 ... from S_const = 4c, and rho_V,tot = rho_V(Y) + c; H^2 = rho/(3(1+Y))."""
    mp.mp.dps = 300
    c = mp.mpf(c_over_M4)
    m2 = mp.mpf(m2_over_M2)
    b3, b4 = mp.mpf(4) / 3, mp.mpf(1) / 3
    F = lambda Y: Y + mp.mpf(1.5) * b3 * Y ** 2 + (b3 / 2 + 2 * b4) * Y ** 3 + b4 * Y ** 4
    dF = lambda Y: 1 + 3 * b3 * Y + 3 * (b3 / 2 + 2 * b4) * Y ** 2 + 4 * b4 * Y ** 3
    Y = 4 * c / (3 * m2)
    for _ in range(200):                  # Newton with a relative stopping rule (Y can be ~1e-120)
        dY = (m2 * F(Y) - 4 * c / 3) / (m2 * dF(Y))
        Y -= dY
        if abs(dY) <= abs(Y) * mp.mpf(10) ** (-280):
            break
    rhoV = mp.mpf(3) / 4 * m2 * Y ** 2 * (1 + b3 * Y + b4 * Y ** 2) + c
    H2 = rhoV / (3 * (1 + Y))
    return Y, rhoV, H2


# ------------------------------------------------------------------------------------------------ (d) compressed
def compressed_run(r, c_mu=0.0, c_eps=0.0, shape=None, w=0.3, fD=0.05, fdust=0.02, Ni=-1.0, Nend=2.5,
                   K0scale=1.0, rtol=1e-11):
    """Toy (radiation + dust + step) at m/H_c = r with the q-offsets: c_mu (type ii, mu shifted) and c_eps
    (type i, eps + c at fixed mu). Returns N, t, H, rho_V,tot = 3H^2(1+Y) + 3HV - (rho_m - c_eps) - ... ."""
    sh = dict(KV_SHAPE) if shape is None else shape
    q = QV1(r, c=c_mu, **sh)
    mat = ToyMatter(3 * fD, 0.0, w, extra_const=c_eps)
    rd_i = 3 * fdust * math.exp(-3 * Ni)
    rr_i = 3 * (1 - fdust - fD) * math.exp(-4 * Ni)
    z = (rr_i, rd_i)
    H = math.sqrt(mat.rho(z, Ni) / 3)
    for _ in range(4):
        Yq, Vq, _, _ = quasi_static_toy(q, mat, Ni, rr_i, rd_i, H)
        H = q.H_from_constraint(Yq, Vq, mat.rho(z, Ni))
    rhs = make_rhs(q, mat)

    def ev_turn(t, y):
        return y[2]
    ev_turn.terminal = False

    def ev_crunch(t, y):
        return y[2] + 3.0
    ev_crunch.terminal = True

    def ev_end(t, y):
        return y[3] - Nend
    ev_end.terminal = True
    sol = solve_ivp(rhs, (0, 1e6), [Yq, Vq, H, Ni, 0.0, rr_i, rd_i], method='DOP853', rtol=rtol, atol=1e-20,
                    events=[ev_turn, ev_crunch, ev_end], dense_output=True)
    Y, V, Hs, N, J, rr, rd = sol.y
    rho_m_noconst = rr + rd + np.array([mat.dV(n) for n in N])
    rhoV_tot = 3 * Hs ** 2 * (1 + Y) + 3 * Hs * V - rho_m_noconst
    return dict(sol=sol, q=q, mat=mat, t=sol.t, N=N, H=Hs, Y=Y, rhoV_tot=rhoV_tot,
                t_turn=[float(x) for x in sol.t_events[0]], status=int(sol.status))


def compressed_demo():
    r = 100.0
    out = {}
    Delta = 3 * 0.05
    chi = 2.0 / (3 * r ** 2)            # KV shape in units K0 = M^2 = 1 with m = r
    # (1) P_q
    a = compressed_run(r)
    late = a['N'] > 2.0
    out['P_q'] = dict(rhoV_tot_late_over_3H2=float(np.max(np.abs(a['rhoV_tot'][late] / (3 * a['H'][late] ** 2)))),
                      note='late rho_V,tot/(3H^2) -> 0: Minkowski attractor (matter dilutes)', N_end=float(a['N'][-1]))
    # (2) initial-equilibrium mu: mu = mu0[eps_inf + Delta] (type ii, exact KV): c = -(mu - mu0) q0
    lam = float(kv_lambda(mp.mpf(chi * Delta)))
    dmu_q0 = (lam * (-1 + (2.0 / 3.0) * lam ** 2) + 1.0 / 3.0) / chi
    b = compressed_run(r, c_mu=-dmu_q0, Nend=6.0)
    Ys, rhoVs, H2s = None, None, None
    # static AdS root for type ii: m^2 F(Y) = (c/3)(2 + 2Y), rho_V = (3/4)m^2Y^2(1+b3Y+b4Y^2) + c(1+Y)
    qb = b['q']
    cc = -dmu_q0
    Yr = brentq(lambda Y: qb.m2 * qb.F(Y) - (cc / 3) * (2 + 2 * Y), -0.5, 0.5)
    rhoV_root = qb.rhoV(Yr)
    i_t = int(np.argmin(np.abs(b['H']))) if b['t_turn'] else None
    red = reduced_toy(ToyMatter(Delta, 0.0, 0.3, extra_const=0.0), qb, (3 * 0.93 * math.exp(4.0), 3 * 0.02 *
                                                                        math.exp(3.0)), -1.0, 0.5, order=0)
    out['initial_equilibrium'] = dict(c_mu=cc, minus_Delta=-Delta, t_turn_full=b['t_turn'][:1],
                                      rhoV_at_turnaround=float(b['rhoV_tot'][i_t]) if i_t is not None else None,
                                      static_AdS_root_rhoV=float(rhoV_root),
                                      series_minus_Delta_plus=float(-Delta + (4.0 / 3.0) * Delta ** 2 / r ** 2),
                                      status=b['status'], H_min=float(b['H'].min()))
    # (3) fixed mu + heavy threshold c_eps = +0.2 Delta (type i): de Sitter
    ch = 0.2 * Delta
    d = compressed_run(r, c_eps=ch, Nend=6.0)
    Yd, rVd, H2d = fixed_mu_static_root(ch, r ** 2)
    late = d['N'] > d['N'][-1] - 1
    out['fixed_mu_plus_c'] = dict(c=ch, rhoV_tot_late=float(np.mean(d['rhoV_tot'][late])),
                                  static_root_rhoV=float(rVd), static_root_H2=float(H2d),
                                  H2_late=float(np.mean(d['H'][late] ** 2)),
                                  one_for_one_check=float(np.mean(d['rhoV_tot'][late]) / ch),
                                  series=float(1 + (4.0 / 3.0) * ch / r ** 2))
    # (4) fixed mu + c_eps = -0.2 Delta: AdS
    e = compressed_run(r, c_eps=-ch, Nend=6.0)
    out['fixed_mu_minus_c'] = dict(c=-ch, t_turn=e['t_turn'][:1], status=e['status'], H_min=float(e['H'].min()))
    # (5) P_q with c re-absorbed: eps + c, mu = mu0[eps + c] -> new q0' = lambda q0; in units of the NEW K0' the
    #     equations are those of (1) with m'^2 = lambda(2 lambda^2 - 1) m^2 (K0 units), b3', b4' (KV shape)
    lam = float(kv_lambda(mp.mpf(chi * ch)))
    m_new = r * math.sqrt(lam * (2 * lam ** 2 - 1) / lam)   # m'^2 in units of K0' (divide by K0'/K0 = lambda)
    shape_new = dict(b3=(4.0 / 3.0) * lam ** 2 / (2 * lam ** 2 - 1), b4=(1.0 / 3.0) * lam ** 2 / (2 * lam ** 2 - 1),
                     beta=1.0)
    # matter densities in K0' units are rho/lambda: rescale fD, fdust consistently by keeping physical rho fixed
    f = compressed_run(m_new, shape=shape_new)
    # compare H(N) physical: 3 K0' H^2 = ... ; H_phys(N) = H_f(N) * sqrt(K0'/K0)... in our toy units both runs use
    # K0 = 1 for their own sector, so physical H with the same physical matter is H_f / sqrt(lambda)
    def binned(run):
        tt = np.linspace(run['t'][0], run['t'][-1], int((run['t'][-1] - run['t'][0]) * r / (2 * math.pi) * 32))
        yy = run['sol'].sol(tt)
        edges = np.arange(-0.8, 2.41, 0.05)
        idx = np.digitize(yy[3], edges)
        return np.array([yy[2][idx == k].mean() for k in range(1, len(edges))])
    Ha, Hf = binned(a), binned(f)
    out['P_q_c_reabsorbed'] = dict(c=ch, chi_c=chi * ch, lambda_minus_1=lam - 1, delta_G_over_G=-(lam - 1),
                                   m_over_sqrtK0_ratio=m_new / r,
                                   max_rel_diff_binned_H_own_G_units=float(np.max(np.abs(Hf / Ha - 1))),
                                   physical_H_shift_same_matter=float(1 / math.sqrt(lam) - 1),
                                   late_rhoV_tot_over_3H2=float(np.max(np.abs(f['rhoV_tot'][f['N'] > 2.0] /
                                                                              (3 * f['H'][f['N'] > 2.0] ** 2)))),
                                   note='same dimensionless history up to the G renormalisation K0 -> lambda K0 '
                                        'and O(chi c H^2/m^2) shape terms; Minkowski attractor kept')
    return out


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                                                  'results')
    data_dir = sys.argv[2] if len(sys.argv) > 2 else DATA_DIR
    eos_dir = sys.argv[3] if len(sys.argv) > 3 else EOS_DIR
    T00 = time.time()
    fetch(data_dir)
    res = {}
    om = OMEGA_B + OMEGA_C
    # ---------------------------------------------------------------------------------- models
    qv1 = Cosmo(om, 0.0)
    omega_L_cal = RHO_LAMBDA / RHO_CRIT100
    cal = Cosmo(om, omega_L_cal)
    planck = Cosmo(0.02237 + 0.1200, 0.6736 ** 2 - (0.02237 + 0.1200) - Cosmo(0, 0).omega_r)  # code check [R]
    desi = load_desi(data_dir)
    sn = load_pplus(data_dir)
    rd = rd_aubourg(OMEGA_B, om)
    models = {}
    for name, cos in (('QV1', qv1), ('QV1-cal', cal), ('Planck2018-LCDM code check', planck)):
        th, zs = cos.theta_star_100(OMEGA_B if name != 'Planck2018-LCDM code check' else 0.02237)
        c2s, Msn = chi2_sn(cos, sn)
        c2b_shape, A = chi2_bao(cos, desi)
        rd_use = rd if name != 'Planck2018-LCDM code check' else rd_aubourg(0.02237, 0.1200 + 0.02237)
        c2b_abs, _ = chi2_bao(cos, desi, rd=rd_use)
        models[name] = dict(h=cos.h, H0=100 * cos.h, omega_m=cos.omega_m, omega_r=cos.omega_r, omega_L=cos.omega_L,
                            Omega_m=cos.omega_m / cos.h ** 2, age_Gyr=cos.age_Gyr(), q0=cos.q0(),
                            theta_star_100=th, z_star=zs, r_d_Mpc=rd_use,
                            chi2_SN_shape=c2s, n_SN=sn['n'], chi2_BAO_shape=c2b_shape, chi2_BAO_absolute=c2b_abs,
                            n_BAO=len(desi[0]), d_L_Mpc={str(z): float((1 + z) * cos.DM(z)[0]) for z in (0.1, 0.5,
                                                                                                         1.0, 2.0)})
    # flat LambdaCDM best fits (shape)
    fsn = minimize_scalar(lambda Om: chi2_sn(flat_lcdm_shape(Om), sn)[0], bounds=(0.05, 0.99), method='bounded',
                          options=dict(xatol=1e-6))
    fbao = minimize_scalar(lambda Om: chi2_bao(flat_lcdm_shape(Om), desi)[0], bounds=(0.05, 0.99), method='bounded',
                           options=dict(xatol=1e-6))

    def sig1(fun, x0, c0):
        lo = brentq(lambda x: fun(x) - c0 - 1, 0.05, x0)
        hi = brentq(lambda x: fun(x) - c0 - 1, x0, 0.99)
        return (hi - lo) / 2
    sn_fun = lambda Om: chi2_sn(flat_lcdm_shape(Om), sn)[0]
    bao_fun = lambda Om: chi2_bao(flat_lcdm_shape(Om), desi)[0]
    fits = dict(SN_flatLCDM=dict(Omega_m=fsn.x, sigma=sig1(sn_fun, fsn.x, fsn.fun), chi2_min=fsn.fun,
                                 published='Brout et al. 2022 abstract: Omega_M = 0.334 +- 0.018 (SN only) [Q]'),
                BAO_flatLCDM=dict(Omega_m=fbao.x, sigma=sig1(bao_fun, fbao.x, fbao.fun), chi2_min=fbao.fun,
                                  published='DESI DR2: Omega_m = 0.2975 +- 0.0086 (BAO only; ledger row 9) [Q]'))
    cmp = {}
    th_check = models['Planck2018-LCDM code check']['theta_star_100']
    for name in models:
        # remove the code's systematic (fitting formula for z_*, massless neutrinos): scale by 1.0411/theta(check)
        models[name]['theta_star_100_calibrated'] = models[name]['theta_star_100'] * PLANCK_THETA[0] / th_check
    for name in ('QV1', 'QV1-cal'):
        mdl = models[name]
        cmp[name] = dict(
            dchi2_SN_vs_bestfit=mdl['chi2_SN_shape'] - fsn.fun,
            dchi2_BAO_shape_vs_bestfit=mdl['chi2_BAO_shape'] - fbao.fun,
            BAO_absolute_chi2=mdl['chi2_BAO_absolute'],
            H0_vs_SH0ES_sigma=(mdl['H0'] - SH0ES_H0[0]) / SH0ES_H0[1],
            q0_vs_SH0ES_sigma=(mdl['q0'] - SH0ES_Q0[0]) / SH0ES_Q0[1],
            theta_vs_Planck_sigma_raw=(mdl['theta_star_100'] - PLANCK_THETA[0]) / PLANCK_THETA[1],
            theta_vs_Planck_sigma_calibrated=(mdl['theta_star_100_calibrated'] - PLANCK_THETA[0]) / PLANCK_THETA[1],
            age_minus_GC_Gyr=mdl['age_Gyr'] - GC_AGE[0])
    # QV1 forced to the SH0ES H0 (EdS with Omega_m = 1): age
    h_sh = SH0ES_H0[0] / 100
    eds_sh = Cosmo(h_sh ** 2, 0.0)
    cmp['QV1_with_SH0ES_H0'] = dict(age_Gyr=eds_sh.age_Gyr(), omega_m_required=h_sh ** 2,
                                    age_deficit_sigma=(GC_AGE[0] - eds_sh.age_Gyr()) / math.hypot(GC_AGE[2], GC_AGE[3]),
                                    omega_m_vs_Planck=(h_sh ** 2) / om)
    res['late_universe'] = dict(models=models, flat_LCDM_fits=fits, comparisons=cmp,
                                data=dict(DESI='DESI DR2 BAO (13 numbers, full covariance), arXiv:2503.14738',
                                          SN='Pantheon+ (zHD > 0.01, STAT+SYS covariance, M marginalised), '
                                             'arXiv:2202.04077', H0='SH0ES 73.04 +- 1.04, arXiv:2112.04510',
                                          q0='SH0ES+Pantheon+ q0 = -0.51 +- 0.024 (same abstract)',
                                          theta='Planck 2018 100 theta_* = 1.0411 +- 0.0003 (abstract)',
                                          age='globular clusters t_U = 13.5 +0.16/-0.14 (stat) +- 0.5 (sys) Gyr, '
                                              'arXiv:2007.06594 (lower bound on the age)'),
                                status='retrospective: all data exposed')
    # ---------------------------------------------------------------------------------- (c) realistic
    eos = SMEoS(eos_dir)
    bg = Background(eos, T_hi=1e3)
    m_h = 1e4
    c_heavy = -1.5 * m_h ** 4 / (64 * math.pi ** 2)
    rej = {}
    rej['initial_equilibrium_EW'] = turnaround(bg, -DELTA_EW)
    rej['initial_equilibrium_EW']['Delta_EW_GeV4'] = DELTA_EW
    rej['initial_equilibrium_EW']['rho_V_inf_over_rhoLambda'] = -DELTA_EW / RHO_LAMBDA
    rej['initial_equilibrium_EW_plus_QCD'] = turnaround(bg, -(DELTA_EW + DELTA_QCD[1]))
    # fixed mu + heavy threshold (10 TeV real scalar, one loop at mu_RG = m_h) [R: MS-bar CW]
    rej['fixed_mu_heavy_threshold_negative'] = turnaround(bg, c_heavy)
    rej['fixed_mu_heavy_threshold_negative']['c_GeV4'] = c_heavy
    T_dom = brentq(lambda lT: eos.rho(math.exp(lT)) - abs(c_heavy), math.log(1.0), math.log(1e4))
    rej['fixed_mu_heavy_threshold_positive'] = dict(c_GeV4=-c_heavy, T_vacuum_domination_GeV=math.exp(T_dom),
                                                    H_dS_GeV=math.sqrt(-c_heavy / (3 * K0)),
                                                    note='de Sitter domination below T ~ 0.9 TeV: the electroweak '
                                                         'crossover never reaches a radiation-dominated plasma')
    # small c = +- rho_Lambda in the late universe: one-for-one
    small = {}
    for sgn in (+1, -1):
        if sgn > 0:
            cos = Cosmo(om, omega_L_cal)
            small['+rho_Lambda'] = dict(H0=100 * cos.h, age_Gyr=cos.age_Gyr(), q0=cos.q0(),
                                        note='fixed-number mu with c = +rho_Lambda is QV1-cal: a free Lambda')
        else:
            # flat matter + negative Lambda: turnaround and crunch (exact): a^3 ~ sin^2(1.5 sqrt(|L|/3) t)
            HL = math.sqrt(RHO_LAMBDA / (3 * K0))
            t_crunch = 2 * math.pi / (3 * HL)
            a_turn = (om / omega_L_cal) ** (1.0 / 3.0)
            small['-rho_Lambda'] = dict(a_turn_over_a0=a_turn, T_turn_K=2.7255 / a_turn,
                                        t_turn_Gyr=t_crunch / 2 * HBAR / GYR, t_crunch_Gyr=t_crunch * HBAR / GYR,
                                        note='matter + (-rho_Lambda) with the CMB omega_m: turnaround when rho_m = '
                                             'rho_Lambda, before the photon temperature reaches 2.7255 K (T0 is never '
                                             'reached); crunch at t = 2 pi/(3 H_L), radiation neglected [D]')
    rej['fixed_mu_small_c'] = small
    # one-for-one with fixed mu vs P_q re-absorption, natural chi = 1/M^4 (exact KV formulas)
    oneforone = []
    for c in (RHO_LAMBDA, -RHO_LAMBDA, DELTA_EW, c_heavy, -c_heavy):
        Yf, rVf, H2f = fixed_mu_static_root(mp.mpf(c) / mp.mpf(M) ** 4)
        lam = kv_lambda(mp.mpf(c) / mp.mpf(M) ** 4)
        oneforone.append(dict(c_GeV4=c, fixed_mu_rhoV_inf_over_c_minus_1=float(rVf / (mp.mpf(c) / mp.mpf(M) ** 4) - 1),
                              fixed_mu_rhoV_inf_GeV4=float(rVf * mp.mpf(M) ** 4),
                              P_q_rhoV_inf=0.0, P_q_delta_K_over_K=float(lam - 1),
                              P_q_delta_m2_over_K0_relative=float(lam * (2 * lam ** 2 - 1) / lam - 1),
                              P_q_delta_m2_relative=float(lam * (2 * lam ** 2 - 1) - 1)))
    rej['one_for_one_vs_P_q'] = oneforone
    rej['existence_P_q_KV_shape'] = 'c >= -M^4/8 (1 + 8 chi c >= 0) on the stable branch; natural chi = 1/M^4'
    res['rejected_readings'] = rej
    # ---------------------------------------------------------------------------------- (d) compressed demo
    res['compressed_full_equations'] = compressed_demo()
    res['meta'] = dict(wall_s=time.time() - T00)
    with open(os.path.join(outdir, 'variants_comparison.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=float)[:12000])


if __name__ == '__main__':
    main()
