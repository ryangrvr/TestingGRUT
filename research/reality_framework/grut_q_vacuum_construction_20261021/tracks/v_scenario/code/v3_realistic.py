"""v_scenario / v3: independent recomputation of the scenario's realistic-scale and late-universe numbers.

Inputs (each with provenance):
  M = 2.435e18 GeV (brief), K0 = M^2; natural KV scales -> m_q = sqrt(2/3) M (v1 check C)
  rho_Lambda = 2.52e-47 GeV^4 (brief)
  omega_b = 0.0224, omega_c = 0.120 (Planck 2018 VI abstract, as used by the scenario; retrospective)
  T0 = 2.7255 K (Fixsen 2009)
  Saikawa-Shirai tables (arXiv:1803.01038 ancillary standardmodel2018.dat: g_rho, g_s; author page eos_w.dat: w),
      downloaded independently by this verifier to the scratchpad v_scenario_dl/ss (byte-identical to the
      scenario's copies)
  Delta_EW = m_H^2 v^2/8 (125.25, 246.22 GeV), heavy 10 TeV real scalar c = -(3/2) m^4/(64 pi^2)
  Planck 2018 100 theta_* = 1.0411 +- 0.0003; SH0ES H0 = 73.04 +- 1.04, q0 = -0.51 +- 0.024 (as quoted by the
      scenario; abstract level)
Method differences from the scenario: own SS loader and interpolation in log T; R1 corrections evaluated with
mpmath (no double-precision floor); turnaround times integrated from T = 1e10 GeV (not from 1 TeV with
t = 1/(2 H_GR)); theta_* by the ratio D_M^LCDM(z*)/D_M^QV1(z*) at fixed r_s (no z* or r_s fitting formula).
Output: ../results/v3_realistic.json
"""
import json
import math
import os

import mpmath as mp
import numpy as np
from scipy.integrate import quad

mp.mp.dps = 260
DL = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/v_scenario_dl/ss'
M = 2.435e18
K0 = M * M
m_nat = math.sqrt(2.0 / 3.0) * M
HBAR = 6.582119569e-25          # GeV s
KB = 8.617333262e-14            # GeV/K
MPC_KM = 3.0856775814913673e19
H100 = 100.0 / MPC_KM / 1.0 * HBAR   # GeV
YR = 3.15576e7
GYR = 1e9 * YR
RHO_C100 = 3 * K0 * H100 ** 2
T0 = 2.7255 * KB
OB, OC = 0.0224, 0.120
RHO_L = 2.52e-47


def load(fn):
    rows = [l.split() for l in open(fn) if l.strip() and not l.startswith('#')]
    return np.array([[float(x) for x in r] for r in rows])


G = load(os.path.join(DL, 'x/anc/standardmodel2018.dat'))
W = load(os.path.join(DL, 'eos_w.dat'))
lTg, grho, gs = np.log(G[:, 0]), G[:, 1], G[:, 3]
lTw, wv, werr = np.log(W[:, 0]), W[:, 1], W[:, 2]


def g_rho(T):
    return np.interp(np.log(T), lTg, grho)


def g_s(T):
    return np.interp(np.log(T), lTg, gs)


def S_over_rho(T):
    """plasma trace / plasma density; SS w above 10 MeV; below, e+- ideal-gas trace at the photon temperature
    (own Fermi-Dirac integral; muons negligible) divided by the SS rho."""
    if T >= 0.999 * W[0, 0]:
        return 1 - 3 * np.interp(math.log(T), lTw, wv)
    trace = 0.0
    for mf in (0.51099895e-3, 0.1056583755):      # e+-, mu+- (g = 4 each), Fermi-Dirac, photon temperature
        f = lambda p: p * p / (math.sqrt(p * p + mf * mf) * (math.exp(min(math.sqrt(p * p + mf * mf) / T, 700)) + 1))
        trace += 4 / (2 * math.pi ** 2) * mf * mf * quad(f, 0, 60 * T + 30 * mf, limit=200)[0]
    return trace / (math.pi ** 2 / 30 * g_rho(T) * T ** 4)


s0 = 2 * math.pi ** 2 / 45 * g_s(T0) * T0 ** 3
rho_d0 = (OB + OC) * RHO_C100


def state(T):
    rho_pl = math.pi ** 2 / 30 * g_rho(T) * T ** 4
    s = 2 * math.pi ** 2 / 45 * g_s(T) * T ** 3
    rho_d = rho_d0 * s / s0
    N = -math.log(T / T0) - math.log(g_s(T) / g_s(T0)) / 3
    return rho_pl, rho_d, N


def S_tot(T):
    rp, rd, N = state(T)
    return S_over_rho(T) * rp + rd, rp + rd, N


def q_numbers(T, m=m_nat, h=1e-3):
    S, rho, N = S_tot(T)
    Sp, _, Np = S_tot(T * math.exp(-h))
    Sm, _, Nm = S_tot(T * math.exp(h))
    SN = (Sp - Sm) / (Np - Nm)
    S, rho, SN = mp.mpf(S), mp.mpf(rho), mp.mpf(SN)
    k0m2 = mp.mpf(K0) * mp.mpf(m) ** 2
    H2_GR = rho / (3 * K0)
    H2 = (rho + S * S / (12 * k0m2)) / (3 * K0 * (1 + (S + SN) / (3 * k0m2)))
    H = mp.sqrt(H2_GR)
    return dict(T_GeV=T, S_over_rho=float(S / rho), m_over_H=float(m / H),
                H_over_HGR_minus_1=mp.nstr(mp.sqrt(H2 / H2_GR) - 1, 4),
                rhoV_GeV4=mp.nstr(S * S / (12 * k0m2), 4),
                deltaG_over_G=mp.nstr(-S / (3 * k0m2), 4),
                Gdot_over_G_per_yr=mp.nstr(-H * SN / (3 * k0m2) / HBAR * YR, 4),
                t_s_radiation=float(1 / (2 * H) * HBAR))


out = {}
# ------------------------------------------------------------------ realistic epochs (natural m_q)
ep = {}
for lab, T in [('1 TeV', 1e3), ('EW crossover 160 GeV', 160.0), ('100 GeV', 100.0), ('QCD 155 MeV', 0.155),
               ('10 MeV', 0.010), ('9.99 MeV (ideal-gas branch)', 0.00999), ('e+e- 0.2 MeV', 2e-4),
               ('BBN 70 keV', 7e-5)]:
    ep[lab] = q_numbers(T)
# today: S = rho_d0 (radiation traceless), S_N = -3 rho_d0
rho_r0 = math.pi ** 2 / 30 * g_rho(T0) * T0 ** 4
S0, rho0 = mp.mpf(rho_d0), mp.mpf(rho_d0 + rho_r0)
k0m2 = mp.mpf(K0) * mp.mpf(m_nat) ** 2
H2GR = rho0 / (3 * K0)
H2 = (rho0 + S0 ** 2 / (12 * k0m2)) / (3 * K0 * (1 + (S0 - 3 * S0) / (3 * k0m2)))
H0QV1 = mp.sqrt(H2GR)
ep['today'] = dict(rhoV_GeV4=mp.nstr(S0 ** 2 / (12 * k0m2), 4), deltaG_over_G=mp.nstr(-S0 / (3 * k0m2), 4),
                   H_over_HGR_minus_1=mp.nstr(mp.sqrt(H2 / H2GR) - 1, 4),
                   Gdot_over_G_per_yr=mp.nstr(-H0QV1 * (-3 * S0) / (3 * k0m2) / HBAR * YR, 4),
                   m_over_H=float(m_nat / H0QV1))
ep['max_rhoV_after_1TeV_over_rhoLambda'] = float(max(float(mp.mpf(ep[k]['rhoV_GeV4'])) for k in ep
                                                      if 'rhoV_GeV4' in ep[k]) / RHO_L)
out['natural_mq_epochs'] = ep

# ------------------------------------------------------------------ late universe QV1, QV1-cal
omega_r = rho_r0 / RHO_C100


def cosmo(omega_L, omega_m=OB + OC):
    h2 = omega_r + omega_m + omega_L
    h = math.sqrt(h2)
    E = lambda z: math.sqrt(omega_r * (1 + z) ** 4 + omega_m * (1 + z) ** 3 + omega_L)
    age = quad(lambda x: 1 / E(math.expm1(x)), 0, 40, limit=400)[0]     # x = ln(1+z): dt = dz/((1+z)H)
    age_Gyr = age / H100 * HBAR / GYR
    q0 = (omega_m / 2 + omega_r - omega_L) / h2
    return dict(h=h, H0=100 * h, age_Gyr=age_Gyr, q0=q0, Omega_m=omega_m / h2, E=E)


qv1 = cosmo(0.0)
qcal = cosmo(RHO_L / RHO_C100)
lcdm = cosmo(0.6736 ** 2 - (OB + OC) - omega_r)     # Planck best fit h = 0.6736 (abstract), flat


def DM(c, z):
    return quad(lambda zz: 1 / c['E'](zz), 0, z, limit=400, epsrel=1e-11)[0]   # units c/H100


zstar = 1089.92   # Planck 2018 VI table 2 [R]; the ratio below varies by < 1e-4 over z* +- 10
theta_qv1 = 1.0411 * DM(lcdm, zstar) / DM(qv1, zstar)
theta_qcal = 1.0411 * DM(lcdm, zstar) / DM(qcal, zstar)
out['QV1_late'] = dict(H0=qv1['H0'], age_Gyr=qv1['age_Gyr'], q0=qv1['q0'], Omega_m=qv1['Omega_m'],
                       omega_r=omega_r,
                       H0_tension_sigma=(qv1['H0'] - 73.04) / 1.04, q0_tension_sigma=(qv1['q0'] + 0.51) / 0.024,
                       theta100_ratio_method=theta_qv1, theta_tension_sigma=(theta_qv1 - 1.0411) / 0.0003,
                       H0_forced_7304=dict(omega_m=0.7304 ** 2 - omega_r,
                                           age_Gyr=cosmo(0.0, 0.7304 ** 2 - omega_r)['age_Gyr']),
                       dL_Mpc={str(z): (1 + z) * DM(qv1, z) * 2997.92458 for z in (0.5, 1.0, 2.0)},
                       dL_Mpc_QV1cal={str(z): (1 + z) * DM(qcal, z) * 2997.92458 for z in (0.5, 1.0, 2.0)})
ag = out['QV1_late']['H0_forced_7304']['age_Gyr']
out['QV1_late']['H0_forced_7304']['age_tension_sigma'] = (ag - 13.5) / math.sqrt(0.15 ** 2 + 0.5 ** 2)
out['QV1cal_late'] = dict(H0=qcal['H0'], age_Gyr=qcal['age_Gyr'], q0=qcal['q0'], Omega_m=qcal['Omega_m'],
                          theta100_ratio_method=theta_qcal, H0_tension_sigma=(qcal['H0'] - 73.04) / 1.04,
                          c_over_M4=RHO_L / M ** 4, second_order_rel=(4.0 / 3.0) * RHO_L / (K0 * m_nat ** 2))

# ------------------------------------------------------------------ rejected readings: turnaround and crunch
DEW = 125.25 ** 2 * 246.22 ** 2 / 8
CH = -1.5 * 1e4 ** 4 / (64 * math.pi ** 2)


def turnaround(const, T_start=1e10):
    """flat FRW: 3K H^2 = rho_pl(T) + rho_d(T) + const (const < 0); t from T_start (t = 1/(2H) there) by
    t = int dN/H, N(T) from entropy conservation; returns T, t at H = 0 (crunch = 2 t_turn, reversible)."""
    lT = np.linspace(math.log(T_start), math.log(1e-3), 400001)
    rho = np.empty_like(lT)
    Nn = np.empty_like(lT)
    for i, l in enumerate(lT):
        rp, rd, N = state(math.exp(l))
        rho[i], Nn[i] = rp + rd + const, N
    i = int(np.argmax(rho <= 0))
    # turnaround N by linear interpolation
    Nt = Nn[i - 1] + (Nn[i] - Nn[i - 1]) * rho[i - 1] / (rho[i - 1] - rho[i])
    Tt = math.exp(lT[i - 1] + (lT[i] - lT[i - 1]) * rho[i - 1] / (rho[i - 1] - rho[i]))
    t0 = 1 / (2 * math.sqrt(rho[0] / (3 * K0)))
    # integrate dN/H on the grid up to i-1, then the square-root end analytically:
    H = np.sqrt(rho[:i] / (3 * K0))
    tgrid = t0 + np.concatenate([[0], np.cumsum(0.5 * (1 / H[1:] + 1 / H[:-1]) * np.diff(Nn[:i]))])
    # last interval: rho ~ k (Nt - N): int dN/sqrt(k (Nt-N)/(3K)) = 2 sqrt(3K (Nt-N_{i-1})/k)
    k = rho[i - 1] / (Nt - Nn[i - 1])
    t_turn = tgrid[-1] + 2 * math.sqrt(3 * K0 * (Nt - Nn[i - 1]) / k)
    return dict(T_turn_GeV=Tt, t_turn_s=t_turn * HBAR, t_crunch_s=2 * t_turn * HBAR,
                analytic_RD_from_a0_s=(math.pi / 4) * math.sqrt(3 * K0 / abs(const)) * HBAR,
                g_rho_at_turn=float(g_rho(Tt)))


rej = dict(EW=turnaround(-DEW), heavy10TeV=turnaround(CH))
# the scenario's convention: start at 1 TeV with t = 1/(2 H_GR) ignoring the constant
rp, rd, _ = state(1e3)
rej['heavy10TeV']['t_at_1TeV_exact_s'] = None
rej['heavy10TeV']['scenario_start_t_1TeV_s'] = 1 / (2 * math.sqrt((rp + rd) / (3 * K0))) * HBAR
rej['heavy10TeV']['true_t_at_1TeV_s'] = None
# true t at 1 TeV in the c < 0 universe (same grid integration)
lT = np.linspace(math.log(1e10), math.log(1e3), 200001)
rr = np.array([sum(state(math.exp(l))[:2]) + CH for l in lT])
NN = np.array([state(math.exp(l))[2] for l in lT])
Hh = np.sqrt(rr / (3 * K0))
rej['heavy10TeV']['true_t_at_1TeV_s'] = float((1 / (2 * Hh[0]) + np.sum(0.5 * (1 / Hh[1:] + 1 / Hh[:-1]) * np.diff(NN))) * HBAR)
rej['plus_heavy_10TeV'] = dict(H_dS_GeV=math.sqrt(-CH / (3 * K0)),
                               T_equal_GeV=None)
# T where rho_pl = |c|
from scipy.optimize import brentq
rej['plus_heavy_10TeV']['T_equal_GeV'] = brentq(lambda T: state(T)[0] - abs(CH), 100, 1e4)
# -rho_Lambda: matter + negative Lambda (radiation neglected), exact
HL = math.sqrt(RHO_L / (3 * K0))
rej['minus_rhoLambda'] = dict(a_turn=(rho_d0 / RHO_L) ** (1 / 3), T_turn_K=2.7255 / (rho_d0 / RHO_L) ** (1 / 3),
                              t_turn_Gyr=math.pi / (3 * HL) * HBAR / GYR, t_crunch_Gyr=2 * math.pi / (3 * HL) * HBAR / GYR)
rej['Delta_EW_GeV4'] = DEW
rej['Delta_EW_over_rhoL'] = DEW / RHO_L
rej['c_heavy_GeV4'] = CH
out['rejected'] = rej

# ------------------------------------------------------------------ light-q bookkeeping numbers
m4 = 4e-12
Gam = lambda m: 4 * m ** 3 / (192 * math.pi * K0)
Hrad = lambda T: math.sqrt(math.pi ** 2 / 30 * g_rho(T) * T ** 4 / (3 * K0))
out['light_q'] = dict(m_over_H_1TeV=m4 / Hrad(1e3), m_over_H_160GeV=m4 / Hrad(160.0),
                      Gamma_over_H_1MeV=Gam(m4) / Hrad(1e-3),
                      m_for_Gamma_eq_H_4MeV=brentq(lambda lm: math.log(Gam(math.exp(lm)) / Hrad(4e-3)), 0, 40),
                      Gamma_nat_over_m=Gam(m_nat) / m_nat)
out['light_q']['m_for_Gamma_eq_H_4MeV'] = math.exp(out['light_q']['m_for_Gamma_eq_H_4MeV'])
out['inputs'] = dict(M=M, m_nat=m_nat, RHO_C100=RHO_C100, T0_GeV=T0, rho_d0=rho_d0, rho_r0=rho_r0, zstar=zstar,
                     g_rho_T0=float(g_rho(T0)), g_s_T0=float(g_s(T0)))


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if k != 'E'}
    return o


here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v3_realistic.json'), 'w') as f:
    json.dump(clean(out), f, indent=1, default=str)
print(json.dumps(clean(out), indent=1, default=str))
