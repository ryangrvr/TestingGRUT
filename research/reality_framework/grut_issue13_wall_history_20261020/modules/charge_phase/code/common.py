"""Shared parameters, unit conversions and helpers for the charge_phase module.

Unchanged Candidate 1 potential (GRUT_CANDIDATE_1.md section 2):
  U = (lam_phi/4)(phi^2 - v^2)^2 + gX^2 phi^2 |X|^2 + lam_X |X|^4
      + 1/2 (mD0^2 + hD^2 phi^2) D^2 + (lam_D/4) D^4.

Dimensionless (Codex) variables used in parts 1, 2, 4, 6:
  S = phi/v,  x = F^2 = 2|X|^2/v^2,  a = gX^2/lam_phi,  b = lam_X/lam_phi,
  length/time unit  ell = 1/(sqrt(lam_phi) v)          (= 1/m_X at the benchmark, a = 1)
  energy density unit  lam_phi v^4,  charge density unit  sqrt(lam_phi) v^3,
  frequency unit  sqrt(lam_phi) v.
  Vbar(S,x) = (S^2-1)^2/4 + a S^2 x/2 + b x^2/4.
For a homogeneous rotating condensate X = |X| e^{i w t}:  nbar = Omega x,  Omega = w/(sqrt(lam_phi) v).

Thermal part (part 3) works in units of v (T/v, n/v^3, f/v^4).
"""
import json
import os
import numpy as np

# ---- benchmark (BRIEF.md) ----
V_GEV = 1000.0
G_X = 0.1
LAM_PHI = 0.01
LAM_X = 0.001
H_D = 0.01
LAM_D = 0.001
M_D0_GEV = G_X * V_GEV

A = G_X**2 / LAM_PHI          # = 1
B = LAM_X / LAM_PHI           # = 0.1
UNIT_FREQ_OVER_V = np.sqrt(LAM_PHI)     # frequency/length^-1 unit in units of v (=0.1)
UNIT_N_OVER_V3 = np.sqrt(LAM_PHI)       # charge-density unit in v^3
UNIT_RHO_OVER_V4 = LAM_PHI              # energy-density unit in v^4
M_X_OVER_UNIT = np.sqrt(A)              # m_X = gX v = sqrt(a) * sqrt(lam_phi) v
MSTAR_OVER_MX = min(1.0, (LAM_PHI * LAM_X) ** 0.25 / G_X)   # 0.5623413252
GEV_TO_FM = 0.1973269804          # hbar c in GeV fm

# archived soliton (grut_candidate_1_20261018/code/results_notes.md, results.json)
ARCHIVED_SOLITON = {
    "Omega": 0.8000000019763,
    "S0": 0.0490799788706,
    "F0": 1.979050307129842,
    "Qbar": 399.4384624223,
    "Ebar": 363.9229985517,
    "E_over_mXQ": 0.9110865196725,
    "R90bar": 5.10343221347,
    "Q_phys": 399.4384624223 / LAM_PHI,
}

# Codex production preparation (479f111 protocol.json / REPORT.md section 5)
CODEX_PRODUCTION = {
    "F2": 0.25, "S": np.sqrt(0.75), "Omega": np.sqrt(0.775),
    "Qbar_per_volume": 0.2200852108, "box_side": 64.0, "seed_rms": 0.005,
    "grids_h": [1.0, 2.0 / 3.0, 0.5],
    "first_negative_time": [73.35, 72.45, 72.35],
    "minS_every_step": [-0.2975431625086389, -0.23768193687699504, -0.26177603005144995],
    "endpoint_negative_F2_q0": [1.0720348358154297, 1.2872283458709717, 1.7876859903335571],
    "endpoint_negative_F2_quartiles": [[3.152033805847168, 3.651457190513611, 4.003987789154053],
                                       [3.0402438640594482, 3.4697275161743164, 3.9329371452331543],
                                       [2.9853484630584717, 3.4360060691833496, 3.9856009483337402]],
    "endpoint_negative_cells": [1228, 549, 1291],
    "endpoint_retained_spacing": [1.0, 4.0 / 3.0, 1.0],
    "min_meanF2_on_negative_cells_every_step": [2.116223047225572, 1.9054749052302398, 1.8919151284606477],
    "visible_T_GeV_codex": 1.7563e5,
}

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.normpath(os.path.join(HERE, "..", "results"))


def Vbar(S, x, a=A, b=B):
    return (S**2 - 1.0) ** 2 / 4.0 + a * S**2 * x / 2.0 + b * x**2 / 4.0


def save_json(name, obj):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, name)

    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(str(type(o)))

    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=conv, sort_keys=False)
    return path
