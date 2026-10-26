#!/usr/bin/env python3
"""kin_source (Claude, own code): fixed inputs, unit conventions and shared helpers for the kinetic source terms.

Inputs are embedded with their provenance.  Nothing here imports or executes Codex code; the Codex JSON files are
read only by check scripts through `git show <commit>:<path>` into the scratchpad.

Conventions (ours = Codex's, stated once):
  * a = 1 at the end of inflation / instant radiation release; H(a) = H_i / a^2 while g* is constant.
  * Comoving numbers N are a^3 n / H_i^3 (dimensionless); N_X is the comoving Noether charge of the released X
    condensate in these units; one bosonic degree per species (X and anti-X are separate degrees).
  * Comoving momentum q = a p, quoted as q/H_i; the conformal energy of a relativistic quantum is q.
  * E_X = a^4 rho_X / H_i^4 is the conformal energy of the released X condensate (conserved in the quartic era);
    R = (conformal energy injected) / E_X.
  * Phase-space convention n = int dp p^2 f / (2 pi^2) per degree, so comoving f(q) = 2 pi^2 (dN/dln q) / q^3.
  * Per net charge = divided by |N_X|; per visible entropy = multiplied by Y_Q/|N_X|; physical density at a is
    N H_i^3 / a^3.
"""
import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRACKS = TRACK.parent
REPO = TRACKS.parents[2]
SCR = Path(os.environ.get('KIN_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                          'scratchpad/review/kin_source'))

# ---------------------------------------------------------------- fixed point [Q, Codex 225028b INPUTS.json / BACKGROUND_RESULTS.json]
COMMIT_PARENT = '225028b'
COMMIT_KINETIC = '9877518'
H_GEV = 133349619398.33566          # H_i
H_GEV_HWGATES = 133350496019.93945  # our hw_gates spot normalisation (differs by 6.6e-6)
MY_GEV = 4000488581950.07
HH_GEV = 5333984775933.427          # h (Sigma vev scale, h/H_i = 40)
MU_GEV = 954710753970.2344
KAPPA = 0.0034586847211042523
LAM_X, LAM_Y, LAM_S = 1e-3, 1.0, 1.0
GX, LAM_PHI, V_PHI = 0.1, 0.03, 0.924089637375362
MX = GX * V_PHI                      # 0.09240896373753621 GeV
MS = math.sqrt(2 * LAM_PHI) * V_PHI  # 0.22635 GeV (radial phi quantum s)
Q_REF = 200000
YQ = 4.8294009076988175e-9
DM_PER_S = 4.4e-10                   # GeV of calibrated dark matter per unit visible entropy
F_ALLOW = 0.01
OMEGA_OVER_M = 0.8611951308217225    # dE/dQ of the reference object (mu_cap)
E_OVER_QM = 0.9233161353726057       # E/Q of the reference object
MRED = 2.435e18
GSTAR_RH = 106.75

# release (primary window a = 30-40) [Q 225028b BACKGROUND_RESULTS.json 'late']
N_X = -2788.1767940228406
E_X = 38386.25063336385
EPS = -0.14571744282991028
N_YH = 0.01725303422930909
N_YL = 0.3344649635088846
N_SIG = 3.158315013358995e-09

# our hw_gates release of the same point (own hwgrel, window a = 30-60) [Q own track]
HW = dict(N_X=-2788.161, E_X=38398.22, EPS=-0.14568, N_YH=0.01725, N_YL=0.33446, N_SIG=3.161e-9,
          N_A_P1=2.4691556381365065, N_A_P0=2.4684882597394293, N_A_P1_e1e6=0.7712575811986991)

# Codex finite-band vector numbers (k/k* = 0.2-5) [Q 225028b VECTOR_RESULTS.json]
CODEX_BAND = {1e-7: 2.0650977334710046, 1e-6: 0.6526449887956494}

# coarse SM g* table (same nodes as our ls_threshold track; reuse stated)
_TS = np.log([1e-6, 5e-4, 1e-3, .1, .15, .3, 1., 5., 100., 200., 1e19])
_G = [3.36, 3.36, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]
_GS = [3.91, 3.91, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]


def gstar(T):
    return float(np.interp(math.log(T), _TS, _G))


def gs(T):
    return float(np.interp(math.log(T), _TS, _GS))


def entropy(T):
    return 2 * math.pi**2 / 45 * gs(T) * T**3


def T_RH(H=H_GEV):
    return (90 * MRED**2 * H**2 / (math.pi**2 * GSTAR_RH))**0.25


def s_RH(H=H_GEV):
    return 2 * math.pi**2 / 45 * GSTAR_RH * T_RH(H)**3


def T_vis_of_a(a, H=H_GEV):
    """Visible temperature at scale factor a with entropy conservation: a T g_s^(1/3) = const (iterated)."""
    T = T_RH(H) / a
    for _ in range(60):
        Tn = T_RH(H) / a * (GSTAR_RH / gs(T))**(1 / 3)
        if abs(Tn / T - 1) < 1e-14:
            break
        T = Tn
    return T


def hubble_of_T(T):
    return math.sqrt(math.pi**2 * gstar(T) / 90) * T * T / MRED


# ---------------------------------------------------------------- light-sector onset quantities [D]
XC = math.sqrt(LAM_PHI * V_PHI**2 / (2 * GX**2))     # amplitude at the phi transition
NC = 2 * math.sqrt(2 * LAM_X) * XC**3                   # circular-orbit charge density at that amplitude
RHO_C = 3 * LAM_X * XC**4                                # its energy density
V0 = LAM_PHI * V_PHI**4 / 4                              # phi vacuum energy released at the transition
N_Q = abs(EPS) * NC                                      # case-L density coordinate (net charge density)
A_REF = H_GEV * (abs(N_X) / N_Q)**(1 / 3)                # scale factor of that coordinate
A_REF_E = H_GEV * (E_X / RHO_C)**0.25


def vector_params(e, H=H_GEV, h=HH_GEV):
    """m_A = 6 e h (charge-6 Higgs), Gamma(A -> X Xbar) = e^2 m_A/(48 pi) beta^3, p* and the decay scales [D]."""
    mA = 6 * e * h
    beta3 = (1 - 4 * MX**2 / mA**2)**1.5
    Gam = e**2 * mA / (48 * math.pi) * beta3
    gam = Gam / H
    pstar = math.sqrt(mA**2 / 4 - MX**2)
    a_unit = 1 / math.sqrt(gam)
    return dict(e=e, mA_GeV=mA, Gamma_A_GeV=Gam, gamma=gam, p_star_GeV=pstar, a_unit=a_unit,
                a_mean=a_unit * math.sqrt(math.pi / 2), q0_over_H=pstar * a_unit / H, mA_over_H=mA / H)


def y_masses(H=H_GEV):
    """Real Majorana-like eigenmodes of Y in the Sigma background: m^2 = M_Y^2 +- sqrt2 mu h [D, sp_yledger]."""
    mL = math.sqrt(MY_GEV**2 - math.sqrt(2) * MU_GEV * HH_GEV)
    mH = math.sqrt(MY_GEV**2 + math.sqrt(2) * MU_GEV * HH_GEV)
    GamL = 3 * KAPPA**2 * mL / (256 * math.pi**3)        # Y_L -> XXX or XbarXbarXbar (half each)
    # Goldstone-equivalence estimate of Y_H -> Y_L + A_L with coupling sqrt2 mu G Y1 Y2 [E, sp_yledger method]
    pA = (mH**2 - mL**2) / (2 * mH)
    GamHLA = 2 * MU_GEV**2 * pA / (8 * math.pi * mH**2)
    msig = math.sqrt(2 * LAM_S) * HH_GEV
    return dict(mYL_GeV=mL, mYH_GeV=mH, Gamma_L_GeV=GamL, gamma_L=GamL / H, a_unit_L=1 / math.sqrt(GamL / H),
                a_mean_L=math.sqrt(math.pi * H / (2 * GamL)), Gamma_H_to_LA_GeV_estimate=GamHLA,
                E_A_from_YH_GeV=pA, msigma_GeV=msig)


def rayleigh_pdf_x(x):
    return x * np.exp(-x * x / 2)


def rayleigh_cdf_x(x):
    return -np.expm1(-x * x / 2)


def save(path, obj):
    def conv(o):
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o).__name__)
    Path(path).write_text(json.dumps(obj, indent=1, default=conv, allow_nan=False) + '\n')


def git_show_json(commit, path):
    """Read a Codex JSON file at a stated commit through git (read-only) and cache it in the scratchpad."""
    out = SCR / 'gitshow' / commit / path.replace('/', '__')
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        data = subprocess.check_output(['git', 'show', f'{commit}:{path}'], cwd=REPO)
        out.write_bytes(data)
    return json.loads(out.read_text())


def translation_table():
    """Multiplicative factors from 'comoving per H_i^3' to the other conventions used in this review."""
    return dict(
        per_net_charge=1 / abs(N_X),
        per_visible_entropy=YQ / abs(N_X),
        physical_density_GeV3_at_a='N * H_i^3 / a^3',
        physical_density_GeV3_at_a_ref=H_GEV**3 / A_REF**3,
        conformal_energy_to_R='E / E_X',
        conformal_energy_to_physical_GeV4_at_a='E * H_i^4 / a^4',
        conformal_energy_to_physical_GeV4_at_a_ref=H_GEV**4 / A_REF**4,
        q_over_H_to_p_over_mX_at_a_ref=H_GEV / (A_REF * MX),
        q_over_H_to_p_GeV_at_a='q * H_i / a',
        f_rest_per_pair_per_net_charge=2 * MX * YQ / DM_PER_S,
        rest_mass_of_net_charge_in_DM_units=MX * YQ / DM_PER_S,
    )
