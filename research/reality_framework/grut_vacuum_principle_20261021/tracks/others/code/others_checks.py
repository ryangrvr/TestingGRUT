"""Checks for track `others` (grut_vacuum_principle_20261021): remaining classes of
vacuum-energy proposals scored against Codex's four-question gate.

Every number produced here is [C] (computed) unless the JSON record labels it as an input.
Nothing here is a new theory, an empirical validation or a confirmatory test.  All
comparisons with Planck/DESI numbers are RETROSPECTIVE: those data are already exposed.

Conventions: signature (-,+,+,+); reduced Planck mass M = (8 pi G)^(-1/2) = 2.435e18 GeV;
hbar = c = 1.  Target quantity: total renormalized rho_vac / M^4.

Sections
  A  target value rho_DE0/M^4 (inputs h, Omega_DE stated)
  B  well-tempered self-tuning (Appleby-Linder 1805.00470): symbolic degeneracy check of the
     simplest example; de Sitter rate h = kappa is a Lagrangian parameter, independent of rho_L
  C  ADDG non-local average-curvature equation (hep-th/0209227 eq. 8, 12): reproduction of
     the quoted M-bar values; sensitivity of the residual to a heavy threshold
  D  DHK degravitation (hep-th/0703027 sec. XI): filter scale L required for a residual M^2/L^2
  E  holographic / IR-cutoff (CKN hep-th/9803132; Li hep-th/0403127; Hsu hep-th/0403052):
     CKN UV cutoff; known thresholds versus the holographic density; Hsu w = 0 identity;
     Li w(z) reproduction; Omega_de0 as a free integration constant; DESI DR2 BAO-only fits
     (retrospective) of LCDM (validation against DESI Table 5), HDE c = 1, HDE free c and the
     Hubble-cutoff case; CPL projection of HDE
  F  Euclidean / wormhole measures (Baum 1983, Hawking 1984, Coleman 1988, Duff 1989):
     on-shell Euclidean de Sitter action, weight ratios at the observed value
  G  anthropic / landscape (Weinberg 1987; Bousso-Polchinski hep-th/0004134): symbolic
     derivation of the 500/729 spherical-collapse bound; EdS window estimate; BP charge scale
     reproduction (J = 100)
  H  global-constraint proposals: Barrow-Shaw (1007.3086, 1105.3105) reproduction of
     Omega_k0 = -0.0055 and update with current inputs, compared with DESI DR2+CMB Omega_K
     (retrospective); everpresent-Lambda order of magnitude
"""
import json
import math
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid, quad
from scipy.optimize import minimize_scalar, brentq

T0 = time.time()
CHECKS = []
OUT = {}


def check(name, cond, detail=""):
    CHECKS.append({"name": name, "pass": bool(cond), "detail": detail})


# ------------------------------------------------------------------ A. target and units
M_GEV = 2.435e18            # reduced Planck mass [GeV] (input, convention)
H_LITTLE = 0.674            # input (Planck-2018-like), not predicted
OMEGA_DE = 0.685            # input, not predicted
H100_GEV = 2.1332e-42       # 100 km/s/Mpc in GeV (hbar = 1)
H0_GEV = H100_GEV * H_LITTLE
RHO_DE0 = 3 * OMEGA_DE * M_GEV**2 * H0_GEV**2          # GeV^4
TARGET = RHO_DE0 / M_GEV**4
OUT["A_target"] = {
    "inputs": {"M_GeV": M_GEV, "h": H_LITTLE, "Omega_DE": OMEGA_DE},
    "rho_DE0_GeV4": RHO_DE0,
    "rho_DE0_quarter_meV": RHO_DE0**0.25 * 1e12,
    "rho_DE0_over_M4": TARGET,
    "(H0/M)^2": (H0_GEV / M_GEV) ** 2,
    "label": "retrospective inputs (Planck 2018-like); matches sequestering track 7.16e-121",
}
check("A1_target_matches_sibling_track", abs(TARGET / 7.16e-121 - 1) < 0.01, f"{TARGET:.4e}")

# ------------------------------------------------------------------ B. well-tempered example
# Appleby-Linder 1805.00470 eqs. (16),(17) with A = const, B = (c1/9k^2)(ln y + 1/y), y = psi'
c1, k, y, ypp, mu, rhoL, psi = sp.symbols("c1 kappa y ypp mu rhoL psi", real=True)
B = c1 / (9 * k**2) * (sp.log(y) + 1 / y)
By, Byy = sp.diff(B, y), sp.diff(B, y, 2)
Ay = 0
Ayy = 0
h, hp = k, 0          # on-shell de Sitter ansatz h = kappa, h' = 0
eq16 = -Ay * y + By * y * (ypp - 3 * h * y)            # times f(y); 2 mu h' f = ... with h' = 0
eq17 = 3 * h * Ay + Ayy * ypp - c1 + 3 * By * (hp * y + h * ypp + 3 * h**2 * y) + 3 * h * Byy * y * ypp
sol16 = sp.solve(sp.Eq(eq16, 0), ypp)
eq17_on = sp.simplify(eq17.subs(ypp, 3 * k * y))
F = -3 * k / (1 - y)
eq28_rhs = sp.simplify(-(3 * k + F) * F / sp.diff(F, y))
OUT["B_well_tempered"] = {
    "eq16_on_shell_solution_psi''": [str(s) for s in sol16],
    "eq17_after_substituting_psi''=3 kappa psi'": str(eq17_on),
    "paper_eq28_rhs_for_F=-3k/(1-psi')": str(eq28_rhs),
    "statement": "On h = kappa both scalar and scale-factor equations reduce to psi'' = 3 kappa psi' "
                 "for any rho_Lambda; rho_Lambda enters only the Hamiltonian constraint, which fixes psi. "
                 "H_dS = kappa*M is a Lagrangian parameter (supplied).",
}
check("B1_eq16_gives_psi''=3k psi'", any(sp.simplify(s - 3 * k * y) == 0 for s in sol16))
check("B2_eq17_identically_satisfied_on_shell", eq17_on == 0, str(eq17_on))
check("B3_matches_paper_eq28", sp.simplify(eq28_rhs - 3 * k * y) == 0, str(eq28_rhs))
# Hamiltonian constraint (paper eq. 34): 3 mu k^2 = mu ML2 - c1 psi + (c1/(3 k)) (psi'-1)*... solve for psi:
ML2 = sp.symbols("ML2", positive=True)
ham = sp.Eq(3 * mu * k**2, mu * ML2 - c1 * psi + c1 / (3 * k**2) * k * (y - 1))
psi_sol = sp.solve(ham, psi)[0]
dpsi_dML2 = sp.diff(psi_sol, ML2)
OUT["B_well_tempered"]["psi_from_Hamiltonian_constraint"] = str(psi_sol)
OUT["B_well_tempered"]["d psi / d(M_Lambda^2)"] = str(dpsi_dML2)
check("B4_vacuum_shift_absorbed_by_psi_not_H", sp.simplify(dpsi_dML2 - mu / c1) == 0, str(dpsi_dML2))

# ------------------------------------------------------------------ C. ADDG average curvature
# eq. (8): M^2 G_mn - (1/4) Mbar^2 g_mn Rbar = T_mn ; vacuum E: rho_eff = E M^2/(M^2+Mbar^2)
def addg_Mbar(E):
    return M_GEV * math.sqrt(E / RHO_DE0 - 1.0)


Mbar_TeV = addg_Mbar(1e3**4)
Mbar_Pl = addg_Mbar(M_GEV**4)
OUT["C_ADDG"] = {
    "Mbar_for_E=(1TeV)^4_GeV": Mbar_TeV,
    "Mbar_for_E=M^4_GeV": Mbar_Pl,
    "paper_quotes": "M-bar ~ 1e48 GeV (TeV^4) and ~1e80 GeV (M_Pl^4) [Q, hep-th/0209227 p.6]",
}
check("C1_reproduces_ADDG_1e48", 47 < math.log10(Mbar_TeV) < 49, f"{Mbar_TeV:.3e}")
check("C2_reproduces_ADDG_order_1e78-80", 77 < math.log10(Mbar_Pl) < 81, f"{Mbar_Pl:.3e}")
# heavy-threshold sensitivity: one-loop m^4/(64 pi^2) ln(m^2/mu^2) per real scalar dof, mu = m/e
def cw(m, ndof=1, lnfac=1.0):
    return ndof * m**4 / (64 * math.pi**2) * lnfac


E0 = 1e3**4                     # supplied 'natural' vacuum energy, TeV^4, fitted via Mbar
rho_eff0 = E0 * M_GEV**2 / (M_GEV**2 + Mbar_TeV**2)
rows = []
for m in (1e3, 1e4, 1e16):
    dE = cw(m)
    rho_eff = (E0 + dE) * M_GEV**2 / (M_GEV**2 + Mbar_TeV**2)
    rows.append({"threshold_mass_GeV": m, "delta_E_GeV4": dE, "rho_eff_over_rho_obs": rho_eff / RHO_DE0})
OUT["C_ADDG"]["threshold_sensitivity_fixed_Mbar"] = rows
OUT["C_ADDG"]["rho_eff_before_over_obs"] = rho_eff0 / RHO_DE0
check("C3_residual_scales_with_E", rows[-1]["rho_eff_over_rho_obs"] > 1e40,
      "adding a 1e16 GeV threshold multiplies the residual; Mbar would need retuning")

# ------------------------------------------------------------------ D. DHK filter scale
L_needed_invGeV = M_GEV / math.sqrt(RHO_DE0)        # rho_res ~ M^2/L^2
hubble_len_invGeV = 1 / H0_GEV
OUT["D_DHK"] = {
    "residual_form": "rho_res ~ M^2 L^-2 (DHK sec. XI: 'leave a small effective cosmological term of order L^-2')",
    "L_needed_GeV^-1": L_needed_invGeV,
    "L_needed_over_1/H0": L_needed_invGeV / hubble_len_invGeV,
    "note": "L is a supplied filter scale; order-one coefficient not computed in DHK",
}
check("D1_L_is_Hubble_scale", 0.1 < L_needed_invGeV / hubble_len_invGeV < 10)

# ------------------------------------------------------------------ E. holographic / IR cutoff
# CKN eq. (2) saturated: L^3 Lam^4 ~ L M_P^2  (CKN use M_P; we show both reduced and non-reduced)
MP_NONRED = M_GEV * math.sqrt(8 * math.pi)
L_H = 1 / H0_GEV
lam_red = (M_GEV**2 / L_H**2) ** 0.25
lam_nonred = (MP_NONRED**2 / L_H**2) ** 0.25
me, mmu, mtop = 0.51099895e-3, 0.1056584, 172.5
rho_e = cw(me, ndof=4)            # Dirac fermion magnitude m^4/(16 pi^2) up to log; order estimate
OUT["E_CKN"] = {
    "UV_cutoff_for_L=1/H0_reduced_eV": lam_red * 1e9,
    "UV_cutoff_for_L=1/H0_nonreduced_eV": lam_nonred * 1e9,
    "CKN_quote": "Lambda ~ 1e-2.5 eV for L ~ horizon [Q]",
    "electron_loop_magnitude_over_rho_obs": rho_e / RHO_DE0,
    "top_loop_magnitude_over_rho_obs": cw(mtop, ndof=12) / RHO_DE0,
    "note": "Known thresholds are below/above this UV cutoff; a holographic IR cutoff does not "
            "subtract the m^4 log terms of particles that exist; CKN state it does not predict the value.",
}
check("E1_CKN_cutoff_meV", 1e-3 < lam_red * 1e9 < 1e-1, f"{lam_red*1e9:.3e} eV")
check("E2_electron_threshold_exceeds_target", rho_e / RHO_DE0 > 1e29, f"{rho_e/RHO_DE0:.2e}")

# Hsu: L = 1/H -> rho_de = 3c^2 M^2 H^2, rho_m = 3(1-c^2) M^2 H^2 -> rho_de/rho_m const -> w = 0
cc, Hs = sp.symbols("c H", positive=True)
ratio = sp.simplify((3 * cc**2 * Hs**2) / (3 * (1 - cc**2) * Hs**2))
check("E3_Hsu_ratio_constant_so_w0", sp.diff(ratio, Hs) == 0, str(ratio))


def w_li(Om, c):
    return -1.0 / 3.0 - 2.0 * np.sqrt(Om) / (3.0 * c)


def li_slope(Om, c):          # Li eq. (28): dw/dz at z = 0
    return np.sqrt(Om) * (1 - Om) * (1 + 2 * np.sqrt(Om) / c) / (6 * c)


OUT["E_Li"] = {
    "w0(Omega=0.73,c=1)": float(w_li(0.73, 1.0)),
    "dw/dz(Omega=0.73,c=1)": float(li_slope(0.73, 1.0)),
    "paper": "w = -0.903 + 0.104 z [Q, hep-th/0403127 eq. 29]",
}
check("E4_reproduces_Li_eq29", abs(w_li(0.73, 1) + 0.903) < 1e-3 and abs(li_slope(0.73, 1) - 0.104) < 1e-3)


def hde_E(z, Om0, c):
    """E(z) = H/H0 for Li HDE + dust, flat. dOmega/dx = Omega(1-Omega)(1+2 sqrt(Omega)/c)."""
    zmax = float(np.max(z))
    xmin = -math.log(1 + zmax) - 0.01
    f = lambda x, O: O * (1 - O) * (1 + 2 * np.sqrt(np.clip(O, 0, 1)) / c)
    s = solve_ivp(f, (0.0, xmin), [1 - Om0], rtol=1e-10, atol=1e-12, dense_output=True)
    x = -np.log1p(z)
    Ode = s.sol(x)[0]
    return np.sqrt(Om0 * (1 + z) ** 3 / (1 - Ode)), Ode


# integration constant: same c = 1, different Omega_de0, all solve the same ODE
consts = []
for Od0 in (0.5, 0.685, 0.9):
    _, Ode = hde_E(np.array([0.0, 1.0, 3.0]), 1 - Od0, 1.0)
    consts.append({"Omega_de0": Od0, "Omega_de(z=1)": float(Ode[1]), "Omega_de(z=3)": float(Ode[2])})
OUT["E_Li"]["integration_constant_family_c=1"] = consts
check("E5_Omega_de0_is_free_initial_condition", len({round(r['Omega_de(z=1)'], 6) for r in consts}) == 3)

# --- DESI DR2 BAO (arXiv:2503.14738 Table 4; [Q] values, retrospective use)
BAO = [  # tracer, z, DV/rd (BGS only) or (DM/rd, sM, DH/rd, sH, r)
    ("BGS", 0.295, ("DV", 7.942, 0.075)),
    ("LRG1", 0.510, (13.588, 0.167, 21.863, 0.425, -0.459)),
    ("LRG2", 0.706, (17.351, 0.177, 19.455, 0.330, -0.404)),
    ("LRG3+ELG1", 0.934, (21.576, 0.152, 17.641, 0.193, -0.416)),
    ("ELG2", 1.321, (27.601, 0.318, 14.176, 0.221, -0.434)),
    ("QSO", 1.484, (30.512, 0.760, 12.817, 0.516, -0.500)),
    ("Lya", 2.330, (38.988, 0.531, 8.632, 0.101, -0.431)),
]
data, cov_blocks, kinds = [], [], []
for name, z, v in BAO:
    if v[0] == "DV":
        data.append(v[1]); cov_blocks.append(np.array([[v[2] ** 2]])); kinds.append(("DV", z))
    else:
        dm, sm, dh, sh, r = v
        data += [dm, dh]
        cov_blocks.append(np.array([[sm**2, r * sm * sh], [r * sm * sh, sh**2]]))
        kinds += [("DM", z), ("DH", z)]
data = np.array(data)
cov = np.zeros((len(data), len(data)))
i0 = 0
for b in cov_blocks:
    n = b.shape[0]; cov[i0:i0 + n, i0:i0 + n] = b; i0 += n
icov = np.linalg.inv(cov)
zs = np.array(sorted({z for _, z in kinds}))


def model_vector(Efun):
    """Distances in units of c/H0 (beta = c/(H0 rd) multiplies all)."""
    zgrid = np.linspace(0, 2.4, 4801)
    E = Efun(zgrid)
    DMg = cumulative_trapezoid(1 / E, zgrid, initial=0)
    out = []
    for kind, z in kinds:
        DM = np.interp(z, zgrid, DMg); DH = 1 / np.interp(z, zgrid, E)
        if kind == "DM": out.append(DM)
        elif kind == "DH": out.append(DH)
        else: out.append((z * DM**2 * DH) ** (1 / 3))
    return np.array(out)


def chi2_profile(Efun):
    m = model_vector(Efun)
    beta = (m @ icov @ data) / (m @ icov @ m)
    r = data - beta * m
    return float(r @ icov @ r), float(beta)


OR_H2 = 4.18e-5       # photons + 3 massless nu, Omega_r h^2 (input)
def E_lcdm(Om, rad=True, h=0.68):
    Orad = OR_H2 / h**2 if rad else 0.0
    return lambda z: np.sqrt(Orad * (1 + z) ** 4 + Om * (1 + z) ** 3 + 1 - Om - Orad)


def fit_1d(fun, lo, hi):
    r = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
    return r.x, r.fun


lcdm_x, lcdm_chi = fit_1d(lambda Om: chi2_profile(E_lcdm(Om))[0], 0.2, 0.4)
# Delta chi2 = 1 interval
lo = brentq(lambda Om: chi2_profile(E_lcdm(Om))[0] - lcdm_chi - 1, 0.2, lcdm_x)
hi = brentq(lambda Om: chi2_profile(E_lcdm(Om))[0] - lcdm_chi - 1, lcdm_x, 0.4)
lcdm_norad_x, lcdm_norad_chi = fit_1d(lambda Om: chi2_profile(E_lcdm(Om, rad=False))[0], 0.2, 0.4)
hde1_x, hde1_chi = fit_1d(lambda Om: chi2_profile(lambda z: hde_E(z, Om, 1.0)[0])[0], 0.15, 0.45)
# free c: profile over grid
c_grid = np.round(np.arange(0.40, 2.01, 0.05), 3)
prof = []
for cval in c_grid:
    xo, ch = fit_1d(lambda Om: chi2_profile(lambda z: hde_E(z, Om, cval)[0])[0], 0.15, 0.45)
    prof.append({"c": float(cval), "Om_best": float(xo), "chi2": float(ch)})
best = min(prof, key=lambda r: r["chi2"])
eds_chi = chi2_profile(lambda z: (1 + z) ** 1.5)[0]   # Hubble-cutoff HDE: total ~ a^-3 exactly
# CPL projection of HDE c=1 at its BAO best fit, and at Omega_m = 0.3
def cpl_project(Om0, c):
    a = np.linspace(1 / 3.33, 1, 400); z = 1 / a - 1
    _, Ode = hde_E(z, Om0, c)
    w = w_li(Ode, c)
    A = np.vstack([np.ones_like(a), 1 - a]).T
    (w0, wa), *_ = np.linalg.lstsq(A, w, rcond=None)
    return float(w0), float(wa), float(w[-1]), float(w[0])
cpl_rows = [{"Om0": om, "c": cv, "w0_fit": cpl_project(om, cv)[0], "wa_fit": cpl_project(om, cv)[1],
             "w(z=0)": cpl_project(om, cv)[2], "w(z=2.33)": cpl_project(om, cv)[3]}
            for om, cv in ((0.30, 1.0), (float(hde1_x), 1.0), (best["Om_best"], best["c"]), (0.30, 0.6))]
OUT["E_DESI_BAO_retrospective"] = {
    "data": "DESI DR2 Table 4 (arXiv:2503.14738), 13 numbers, per-bin DM-DH correlations; "
            "H0*rd profiled analytically; flat; BAO only",
    "LCDM_best_Om": lcdm_x, "LCDM_Om_1sigma": [lo, hi], "LCDM_chi2": lcdm_chi,
    "LCDM_no_radiation_best_Om": lcdm_norad_x, "LCDM_no_radiation_chi2": lcdm_norad_chi,
    "DESI_Table5_LCDM_DESI_only": "Omega_m = 0.2975 +/- 0.0086 [Q]",
    "HDE_c1_best_Om": hde1_x, "HDE_c1_chi2": hde1_chi,
    "HDE_free_c_best": best, "HDE_free_c_profile": prof,
    "Hubble_cutoff_EdS_shape_chi2": eds_chi,
    "dof_note": "13 data; LCDM and HDE(c=1) fit 2 params (Om, H0 rd); HDE free c fits 3",
    "CPL_projection_of_HDE": cpl_rows,
    "DESI_DR2_w0wa_quadrant": "DESI+CMB+{Pantheon+,Union3,DESY5}: w0 = -0.838, -0.667, -0.752; wa = -0.62, -1.09, -0.86 [Q, Table 5]",
}
check("E6_LCDM_BAO_fit_reproduces_DESI_Om", abs(lcdm_x - 0.2975) < 0.004 and abs((hi - lo) / 2 - 0.0086) < 0.002,
      f"Om = {lcdm_x:.4f} [{lo:.4f},{hi:.4f}]")
check("E7_Hubble_cutoff_excluded_by_BAO_shape", eds_chi - lcdm_chi > 100, f"dchi2 = {eds_chi-lcdm_chi:.1f}")
check("E8_HDE_c1_CPL_wa_positive", cpl_rows[0]["wa_fit"] > 0, f"wa = {cpl_rows[0]['wa_fit']:.3f}")

# ------------------------------------------------------------------ F. Euclidean measures
S_dS = 24 * math.pi**2 / TARGET            # -I_E = 24 pi^2 M^4/rho  (reduced Planck units)
OUT["F_Euclidean"] = {
    "minus_I_E_at_rho_obs": S_dS,
    "ln[P(rho_obs/2)/P(rho_obs)]_Hawking_weight": S_dS,
    "Coleman_weight": "exp(exp(+24 pi^2 M^4/rho)) -> peaked at rho -> 0+ more strongly",
    "sign_flipped_weight_(tunnelling/Duff-sign)": "exp(-24 pi^2 M^4/rho) -> favours largest rho",
    "formula_check": "I_E = -(1/16 pi G) int (R - 2 Lambda) = -3 pi/(G Lambda) = -24 pi^2 M^4/rho",
}
# symbolic: S^4 radius r, Lambda = 3/r^2, R = 4 Lambda, vol = 8 pi^2 r^4/3
G_, Lam_ = sp.symbols("G Lambda", positive=True)
r2 = 3 / Lam_
IE = -(1 / (16 * sp.pi * G_)) * (4 * Lam_ - 2 * Lam_) * sp.Rational(8, 3) * sp.pi**2 * r2**2
check("F1_IE_equals_-3pi/(G Lambda)", sp.simplify(IE + 3 * sp.pi / (G_ * Lam_)) == 0, str(sp.simplify(IE)))
check("F2_matches_obstruction_track_lnd_3.3e122", abs(S_dS / 3.3e122 - 1) < 0.02, f"{S_dS:.3e}")

# ------------------------------------------------------------------ G. anthropic / landscape
# spherical shell: (dr/dt)^2 = 2GM/r + (Lam/3) r^2 - K; turnaround iff Lam/3 <= K^3/(27 (GM)^2)
GM, Kc, rr, Hh, dd, Lm = sp.symbols("GM K r H delta Lam", positive=True)
V = -2 * GM / rr - Lm / 3 * rr**2
rstar = sp.solve(sp.diff(V, rr), rr)[0]
Vmax = sp.simplify(V.subs(rr, rstar))
Lam_crit = sp.solve(sp.Eq(Kc, -Vmax), Lm)[0]
# growing-mode shell: K = (5/3) delta H^2 r^2, GM = H^2 r^3/2 ; H^2 = (8 pi G/3) rho_m
Lam_crit_sub = sp.simplify(Lam_crit.subs({Kc: sp.Rational(5, 3) * dd * Hh**2 * rr**2, GM: Hh**2 * rr**3 / 2}))
# rho_V,max / rho_m = (Lam/3) / H^2  (since Lam/3 = 8 pi G rho_V/3 and H^2 = 8 pi G rho_m/3)
coef = sp.simplify(Lam_crit_sub / 3 / Hh**2 / dd**3)
OUT["G_anthropic"] = {"rho_V_max_over_rho_m_delta^3": str(coef)}
check("G1_Weinberg_500_over_729", sp.simplify(coef - sp.Rational(500, 729)) == 0, str(coef))
zc = 4.0
ratio_max = 500 / 729 * (1.686 * (1 + zc)) ** 3                 # rho_V,max / rho_m0 (EdS extrapolation)
OMEGA_M = 0.315
OUT["G_anthropic"].update({
    "z_c": zc, "rho_V_max_over_rho_m0_EdS": ratio_max,
    "rho_V_max_over_rho_obs": ratio_max * OMEGA_M / OMEGA_DE,
    "BP_window_quote": "-1e-120 < lambda < 1e-118 M^4, upper bound ~100x observed [Q, hep-th/0004134 eq. 3.13]",
    "flat_prior_P(rho<=rho_obs | 0<rho<rho_max)": OMEGA_DE / (ratio_max * OMEGA_M),
    "label": "[E] order of magnitude; MSW 1998 compute a collapsed-fraction-weighted distribution",
})
check("G2_window_order_100", 30 < ratio_max * OMEGA_M / OMEGA_DE < 1000, f"{ratio_max*OMEGA_M/OMEGA_DE:.0f}")
# BP shell estimate eq. (2.26) with D = 1, lambda_bare = 1, J = 100: q such that dlambda_min = 1e-120
def log10_dlam(J, log10q, lam_bare=1.0, D=1.0):
    log10_omega = math.log10(2) + (J / 2) * math.log10(math.pi) - math.lgamma(J / 2) / math.log(10)
    return math.log10(D) + J * log10q - log10_omega - (J / 2 - 1) * math.log10(2 * lam_bare)
lq = brentq(lambda t: log10_dlam(100, t) + 120, -5, 0)
# BP closing argument: at nucleosynthesis the present Lambda is dynamically irrelevant.
gstar_bbn = 10.75
bbn_rows = []
for T_MeV in (1.0, 0.1):
    rho_r = math.pi**2 / 30 * gstar_bbn * (T_MeV * 1e-3) ** 4
    bbn_rows.append({"T_MeV": T_MeV, "rho_Lambda_over_rho_rad": RHO_DE0 / rho_r})
OUT["G_BBN_fraction"] = {"rows": bbn_rows, "BP_quote": "'only a fraction 1e-30' [Q, hep-th/0004134 sec. 4] (loose; computed value is smaller)"}
check("G4_BBN_fraction_below_1e-30", bbn_rows[0]["rho_Lambda_over_rho_rad"] < 1e-30, f"{bbn_rows[0]['rho_Lambda_over_rho_rad']:.2e}")
OUT["G_BP"] = {"J": 100, "log10_q_needed_for_1e-120": lq, "BP_quote": "q ~ 1e-1.6 (J = 100, lambda_bare ~ 1, D small) [Q]"}
check("G3_reproduces_BP_charge_scale", abs(lq + 1.6) < 0.3, f"log10 q = {lq:.2f}")

# ------------------------------------------------------------------ H. global constraints
def barrow_shaw(OL, Ob, h, zeta=0.5, Tcmb=2.725, Neff=3.046, n=200001):
    """Omega_k0 from Barrow-Shaw eq. (3) (1007.3086), flat background to leading order.
    Units H0 = 1, a0 = 1, conformal time tau."""
    Og = 2.4728e-5 * (Tcmb / 2.725) ** 4 / h**2
    Orad = Og * (1 + 0.22711 * Neff)
    Om = 1 - OL - Orad
    a = np.logspace(-9, 0, n)
    Hf = np.sqrt(Orad / a**4 + Om / a**3 + OL)
    dtau = 1 / (a**2 * Hf)
    tau = cumulative_trapezoid(dtau, a, initial=0) + a[0] / np.sqrt(Orad)
    I2 = cumulative_trapezoid(1 / (a**2 * Hf**3), a, initial=0) + (a[0] ** 3 / 3) / Orad**1.5
    A = a * Hf / 6 * I2
    t0 = tau[-1]
    num = np.trapezoid((t0 - tau) ** 3 * a**4 * zeta * 3 * Ob * a**-3 * A * dtau, a)
    den = np.trapezoid(a**2 * (t0 - tau) ** 2 * (4 * (t0 - tau) + 6 * tau) * A * dtau, a)
    return -num / den


bs_orig = barrow_shaw(0.73, 0.0423, 0.71)
bs_planck = barrow_shaw(0.6847, 0.0493, 0.6736)
bs_desi = barrow_shaw(1 - 0.3027, 0.02218 / 0.6817**2, 0.6817)
OK_DESI, sOK_DESI = 0.0023, 0.0011          # DESI DR2 + CMB, LCDM + Omega_K [Q, Table 5]
OUT["H_Barrow_Shaw"] = {
    "paper_prediction": "Omega_k0 = -0.0055 (2 zeta_b), zeta_b = 1/2 chiral bag [Q, 1007.3086]",
    "reproduced_with_paper_inputs": bs_orig,
    "updated_Planck2018_like_inputs(OL=0.6847,Ob=0.0493,h=0.6736)": bs_planck,
    "updated_DESI+CMB_LCDM_inputs(Om=0.3027,h=0.6817)": bs_desi,
    "DESI_DR2+CMB_Omega_K": [OK_DESI, sOK_DESI],
    "tension_sigma_paper_value": (OK_DESI - (-0.0055)) / sOK_DESI,
    "tension_sigma_updated": (OK_DESI - bs_desi) / sOK_DESI,
    "sign_only_tension_sigma(zeta_b->0+)": OK_DESI / sOK_DESI,
    "label": "retrospective; prediction published 2010-2011 before Planck 2013 and DESI 2025",
}
check("H1_reproduces_BS_-0.0055", abs(bs_orig + 0.0055) < 0.0002, f"{bs_orig:.5f}")
check("H2_updated_prediction_negative", bs_desi < 0 and bs_planck < 0)
# everpresent Lambda: Lambda ~ +/- V^(-1/2) in Planck units (l_P^4 = (8 pi G)^2 ...) order only
V4 = (1 / H0_GEV) ** 4                      # crude 4-volume ~ H0^-4 [GeV^-4]
Lam_ep = 1 / math.sqrt(V4 / (1 / MP_NONRED) ** 4) * MP_NONRED**2   # Lambda ~ l_P^-2 (V/l_P^4)^-1/2 [GeV^2]
rho_ep = Lam_ep * M_GEV**2
OUT["H_everpresent"] = {
    "rho_estimate_over_rho_obs": rho_ep / RHO_DE0,
    "note": "[E] Lambda ~ V^(-1/2) with V ~ H0^-4 gives rho ~ M_P^2 H0^2 up to O(1-10) factors; mean of Lambda assumed zero; sign fluctuates",
}
check("H3_everpresent_order_H0^2", 0.01 < rho_ep / RHO_DE0 < 1e3, f"{rho_ep/RHO_DE0:.2f}")

OUT["checks"] = CHECKS
OUT["n_checks"] = len(CHECKS)
OUT["all_checks_pass"] = all(c["pass"] for c in CHECKS)
OUT["versions"] = {"python": sys.version.split()[0], "numpy": np.__version__, "sympy": sp.__version__,
                   "scipy": __import__("scipy").__version__}
OUT["elapsed_s"] = time.time() - T0
OUT["scope"] = ("Reproductions of published formulas and numbers, symbolic identities, order-of-magnitude "
                "estimates and a BAO-only retrospective fit. No new theory, no empirical validation, "
                "no confirmatory test.")
if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "others_checks.json"
    with open(path, "w") as fh:
        json.dump(OUT, fh, indent=1, default=float)
    for c in CHECKS:
        print(("PASS " if c["pass"] else "FAIL ") + c["name"] + "  " + c["detail"])
    print("all pass:", OUT["all_checks_pass"], f"({len(CHECKS)} checks)")
