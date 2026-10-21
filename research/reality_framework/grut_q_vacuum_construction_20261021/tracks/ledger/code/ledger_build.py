#!/usr/bin/env python3
"""Track `ledger`: numbers for the GRUT mapping (NOTES.md sections 3-4) and the machine-readable
exposed-data ledger (results/ledger.json).

Run:  python3 -I code/ledger_build.py results/
Writes results/ledger_numbers.json and results/ledger.json. Deterministic, standard library + numpy.

Every observational number below is a quoted datum (already exposed, retrospective). The source
identifier and the read level achieved are recorded with it. Nothing here is fitted.
"""
import json
import math
import sys
import platform

import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/"
if not OUT.endswith("/"):
    OUT += "/"

checks = []


def check(name, ok, detail=""):
    checks.append({"name": name, "pass": bool(ok), "detail": detail})


# ---------------------------------------------------------------- constants (BRIEF conventions)
M = 2.435e18                    # reduced Planck mass, GeV (BRIEF)
K = M**2                        # GeV^2
M4 = M**4
rho_L = 2.52e-47                # GeV^4, rho_Lambda,obs (BRIEF; Planck 2018 inputs)
hbar_GeV_s = 6.582119569e-25    # GeV s
hbarc_eV_m = 1.973269804e-7     # eV m
Mpc_m = 3.0856775814913673e22
yr_s = 3.15576e7
R_sun_m = 6.957e8               # nominal solar radius (IAU 2015 B3) [R]

H0_kms = 67.4                   # Planck 2018 VI abstract (exposed)
Om = 0.315                      # Planck 2018 VI abstract (exposed)
OL = 1 - Om
H0_s = H0_kms * 1e3 / Mpc_m
H0_GeV = H0_s * hbar_GeV_s
H0_eV = H0_GeV * 1e9
H0_yr = H0_s * yr_s

num = {"constants": {"M_GeV": M, "rho_Lambda_obs_GeV4": rho_L, "rho_Lambda_over_M4": rho_L / M4,
                     "H0_GeV": H0_GeV, "H0_eV": H0_eV, "H0_per_yr": H0_yr}}
check("rho_Lambda^(1/4) = 2.24 meV", abs(rho_L**0.25 * 1e12 - 2.2405) < 0.002, f"{rho_L**0.25*1e12:.4f} meV")

# natural q-sector scales of the `equations` track section 7: q0^2 eps'' = M^4, q0 K' = M^2
mq2_nat = (2.0 / 3.0) * M**2     # m_q^2 = (2/3) K eps''/K'^2 at natural scales
mq_nat = math.sqrt(mq2_nat)
num["constants"]["m_q_natural_GeV"] = mq_nat

# ---------------------------------------------------------------- Candidate 1 inputs (repository, quoted)
v = 1000.0                       # GeV, benchmark v = 1 TeV (ledger_s08 13c; wall-history thermal module)
c1loop = -1.22652632365e-6 * v**4          # GeV^4, quantum REPORT.md section 6
lam_phi = 0.01                              # at mu_RG = 0.1 v (thermal_potential NOTES section 1)
C_tree = 0.0025 * v**4                      # polynomial constant = lam_phi v^4 / 4
check("C = lam_phi v^4/4 (tree broken minimum at V = 0)", abs(C_tree - lam_phi * v**4 / 4) < 1e-6 * C_tree)
Tc_lo, Tc_hi = 1.557140 * v, 1.564147 * v   # thermal_potential adopted range (AE, Parwani, DR schemes)
Tc_mid = 0.5 * (Tc_lo + Tc_hi)

# SM vacuum-energy drops (equations track section 7 inputs, recalled values)
mH, vEW = 125.25, 246.22
D_EW = mH**2 * vEW**2 / 8.0
D_QCD = (4.4e-4, 3.0e-3)                    # bag constant range [R], GeV^4

c1 = {}
c1["rho_vac_1loop_GeV4"] = c1loop
c1["rho_vac_1loop_over_rhoL"] = c1loop / rho_L
c1["rho_vac_1loop_over_M4"] = c1loop / M4
# Under P_q the constant is absorbed by re-selecting mu: d mu0/dc = 1/q0, d q0/dc = 1/(q0 eps'') (equations 5.1)
# natural scales: delta q0/q0 = c chi0 with chi0^-1 = q0^2 eps'' = M^4; delta mu0/mu0 = c/(mu0 q0) = c/eps(q0) ~ c/M^4
c1["delta_q0_over_q0_natural"] = c1loop / M4
c1["delta_mu0_over_mu0_natural"] = c1loop / M4
c1["delta_K_over_K_natural"] = c1loop / M4          # (q0 K'/K) (delta q0/q0) with q0 K' ~ K
c1["Z2_DeltaV_tree_GeV4"] = C_tree
c1["Z2_DeltaV_tree_over_rhoL"] = C_tree / rho_L
c1["Z2_DeltaV_over_Delta_EW"] = C_tree / D_EW
c1["Delta_EW_GeV4"] = D_EW
c1["Z2_1loop_rel_to_tree"] = abs(c1loop) / C_tree
c1["Tc_GeV"] = [Tc_lo, Tc_hi]

# plateau fraction and Hubble rate at the hidden transition, for hidden/SM temperature ratio xi
g_SM, g_hid = 106.75, 4.0      # hidden bosons phi, q1, q2, D (thermal_potential NOTES section 1)
plateau = []
for xi in (1.0, 0.1, 0.01):
    rows = []
    for Tc in (Tc_lo, Tc_hi):
        T_SM = Tc / xi
        rho_r = (math.pi**2 / 30.0) * (g_SM * T_SM**4 + g_hid * Tc**4)
        H = math.sqrt((rho_r + C_tree) / (3 * K))
        rows.append({"Tc_GeV": Tc, "T_SM_GeV": T_SM, "rho_r_GeV4": rho_r,
                     "DeltaV_over_rho_r": C_tree / rho_r, "H_GeV": H, "H4_over_rhoL": H**4 / rho_L,
                     "m_q_nat_over_H": mq_nat / H})
    plateau.append({"xi": xi, "rows": rows})
c1["plateau_at_Tc"] = plateau
check("Z2 plateau fraction at T_c below 2e-5 for xi = 1",
      all(r["DeltaV_over_rho_r"] < 2e-5 for r in plateau[0]["rows"]))

# second-order transient (equations 5.3): (4/3) Delta^2 / (K m_q^2), not a constant remnant
sec = (4.0 / 3.0) * C_tree**2 / (K * mq2_nat)
c1["Z2_second_order_transient_GeV4"] = sec
c1["Z2_second_order_transient_over_rhoL"] = sec / rho_L
secEW = (4.0 / 3.0) * D_EW**2 / (K * mq2_nat)
c1["EW_second_order_transient_GeV4_check"] = secEW
check("EW second-order transient reproduces equations track 8.0e-58 GeV^4", abs(secEW / 8.0e-58 - 1) < 0.02,
      f"{secEW:.3e}")

# initial-equilibrium reading (mu fixed before the transitions): remnant = -(sum of later drops)
for lab, dq in (("lo", D_QCD[0]), ("hi", D_QCD[1])):
    rem = -(C_tree + D_EW + dq)
    c1[f"initial_reading_remnant_GeV4_{lab}"] = rem
    c1[f"initial_reading_remnant_over_rhoL_{lab}"] = rem / rho_L
premise = -2.72e-47      # rho_vac,tot >= -(2.28 meV)^4 (ledger_s08 13b, stress corner)
c1["H0_exclusion_premise_GeV4"] = premise
check("premise = -(2.28 meV)^4", abs((-premise) ** 0.25 * 1e12 - 2.28) < 0.01)
c1["initial_reading_violates_premise_by"] = c1["initial_reading_remnant_GeV4_lo"] / premise

# final reading: today rho_V = 0 + matter-induced rho_m^2/(12 K m_q^2) (equations 6.4)
rho_m0 = 3 * K * H0_GeV**2 * Om
c1["final_reading_rhoV_today_matter_induced_GeV4"] = rho_m0**2 / (12 * K * mq2_nat)
c1["final_reading_satisfies_premise"] = True
num["candidate1"] = c1

# ---------------------------------------------------------------- bounds on the q-mode from exposed data
b = {}
# Cassini: gamma - 1 = -(2/3) e^{-m b}/(1 + e^{-m b}/3) for a Yukawa of strength 1/3 (equations 4.3)
sig_gamma = 2.3e-5      # Bertotti et al. 2003, gamma - 1 = (2.1 +- 2.3)e-5
for lab, bimp in (("b=R_sun", R_sun_m), ("b=1.6R_sun", 1.6 * R_sun_m)):
    # solve (2/3) e/(1+e/3) = sig  ->  e = sig/(2/3 - sig/3)
    e = sig_gamma / (2.0 / 3.0 - sig_gamma / 3.0)
    mr = -math.log(e)
    b[f"cassini_m_min_eV_{lab}"] = mr / bimp * hbarc_eV_m
b["cassini_m_r_required"] = -math.log(sig_gamma / (2.0 / 3.0 - sig_gamma / 3.0))
# Eot-Wash: |alpha| = 1 excluded for lambda >= 38.6 um (95%); alpha = 1/3 point not read (estimate)
b["eotwash_m_min_eV_alpha1_point"] = hbarc_eV_m / 38.6e-6
# LLR Gdot/G: heavy branch Gdot/G = 3 Om H0^3/m^2 (quasi-static, m >> H)
Gdot_max = 2 * 9.6e-15   # per yr, twice the quoted LLR uncertainty (Biskupek et al. 2021)
b["llr_Gdot_max_per_yr_used"] = Gdot_max
b["llr_m_min_over_H0"] = math.sqrt(3 * Om * H0_yr / Gdot_max)
b["llr_m_min_eV"] = b["llr_m_min_over_H0"] * H0_eV
b["heavy_branch_Gdot_over_G_today_natural_per_yr"] = 3 * Om * H0_yr * (H0_GeV / mq_nat) ** 2
# BBN: delta K/K ~ kappa H^2/m^2 with trace fraction kappa ~ 0.1 near e+e- annihilation [E]
g_bbn = 10.75
T_bbn = 1e-3   # GeV
rho_bbn = (math.pi**2 / 30) * g_bbn * T_bbn**4
H_bbn = math.sqrt(rho_bbn / (3 * K))
b["H_BBN_eV"] = H_bbn * 1e9
b["bbn_m_min_eV_kappa0.1_dG0.06"] = H_bbn * 1e9 * math.sqrt(0.1 / 0.06)
# EDE-equivalent at BBN from Yeh et al. 2022: N_nu < 3.180 (2 sigma), SM 3 -> Delta N < 0.18
dN = 3.180 - 3.0
b["bbn_DeltaN_max"] = dN
b["bbn_max_extra_fraction"] = 1.75 * dN / (g_bbn + 1.75 * dN)
b["rhoL_over_rho_r_at_1MeV"] = rho_L / rho_bbn
# QCD bag plateau in GR (context, same in GR and in q-theory with final mu): B / rho_r at T = 155 MeV, g* = 61.75
rho_qcd = (math.pi**2 / 30) * 61.75 * 0.155**4
b["qcd_bag_over_rho_r"] = [D_QCD[0] / rho_qcd, D_QCD[1] / rho_qcd]
# heavy-branch deviation of the vacuum from w = -1 today: rho_m^2/(12 K m^2) / rho_Lambda
def dev_w(m_eV):
    m_GeV = m_eV * 1e-9
    return (rho_m0**2 / (12 * K * m_GeV**2)) / (3 * K * H0_GeV**2 * OL)
b["heavy_branch_rhoV_shift_over_rhoL_natural"] = dev_w(mq_nat * 1e9)
b["heavy_branch_rhoV_shift_over_rhoL_at_cassini_bound"] = dev_w(b["cassini_m_min_eV_b=R_sun"])
# Planck-mass q-oscillation lifetime, scalaron-type decay Gamma ~ m^3/(192 pi M^2) per scalar [R/E]
Gam = mq_nat**3 / (192 * math.pi * M**2)
b["planck_mass_q_decay_rate_GeV_estimate"] = Gam
b["planck_mass_q_lifetime_s_estimate"] = hbar_GeV_s / Gam
num["bounds"] = b
check("Cassini bound ~3e-15 eV (equations track)", 2e-15 < b["cassini_m_min_eV_b=R_sun"] < 4e-15,
      f"{b['cassini_m_min_eV_b=R_sun']:.3e}")
check("Eot-Wash alpha=1 point ~5 meV", 4e-3 < b["eotwash_m_min_eV_alpha1_point"] < 6e-3)
check("natural m_q passes all bounds",
      mq_nat * 1e9 > max(b["cassini_m_min_eV_b=R_sun"], b["eotwash_m_min_eV_alpha1_point"], b["llr_m_min_eV"]))

# ---------------------------------------------------------------- the ledger (machine-readable)
L = []


def entry(**kw):
    kw.setdefault("exposed", True)
    kw.setdefault("status_of_comparison", "retrospective")
    L.append(kw)


entry(id="Gdot_LLR", observable="dG/dt / G today (lunar laser ranging)",
      construction_quantity="K(q) drift: Gdot/G = -Kdot/K; heavy branch 3 Om H0^3/m_q^2",
      value="(-5.0 +- 9.6)e-15 /yr", numbers={"central": -5.0e-15, "sigma": 9.6e-15, "unit": "1/yr"},
      cl="as quoted in the abstract (convention not stated there)",
      source="Biskupek, Mueller, Torre 2021, arXiv:2012.12032 (Universe 7, 34 [R])", read_level="abstract",
      also=[{"value": "(4 +- 9)e-13 /yr", "source": "Williams, Turyshev, Boggs, PRL 93, 261101 (2004), gr-qc/0411113",
             "read_level": "abstract"},
            {"value": "(7.1 +- 7.6)e-14 /yr", "source": "Hofmann, Mueller, CQG 35, 035015 (2018), DOI 10.1088/1361-6382/aa8f7a",
             "read_level": "recalled [R] (publisher page blocked)"}],
      implication=f"heavy branch: Gdot/G ~ {b['heavy_branch_Gdot_over_G_today_natural_per_yr']:.1e}/yr at natural m_q; "
                  f"bound only m_q > {b['llr_m_min_over_H0']:.0f} H0 (weaker than Cassini)")
entry(id="Gdot_pulsar", observable="dG/dt / G (binary pulsar J1713+0747)",
      construction_quantity="as Gdot_LLR, for strongly self-gravitating bodies",
      value="(-0.1 +- 0.9)e-12 /yr", numbers={"central": -0.1e-12, "sigma": 0.9e-12, "unit": "1/yr"},
      cl="68% (stated as +-)", source="Zhu et al., MNRAS 482, 3249 (2019), arXiv:1802.09206, DOI 10.1093/mnras/sty2905",
      read_level="abstract", implication="weaker than LLR; no additional constraint")
entry(id="Gdot_helio", observable="dG/dt / G (helioseismology)",
      construction_quantity="as Gdot_LLR, integrated over the solar age",
      value="(1.25 +- 0.30)e-13 /yr (claimed 4 sigma nonzero)", numbers={"central": 1.25e-13, "sigma": 0.30e-13, "unit": "1/yr"},
      cl="Bayesian, as quoted", source="Bonanno, Froehlich, arXiv:1707.01866 (no journal reference in arXiv metadata)",
      read_level="abstract",
      implication="a single-analysis claim of nonzero Gdot, in mild tension with LLR at its central value; "
                  "the heavy branch predicts Gdot/G ~ 0 and cannot produce 1e-13/yr without m_q ~ few H0 (excluded by Cassini)")
entry(id="G_BBN", observable="G at BBN relative to today",
      construction_quantity="K(q_BBN)/K(q_0); heavy branch delta K/K ~ kappa H_BBN^2/m_q^2",
      value="G_BBN/G_0 = 0.99 (+0.06 -0.05)", numbers={"central": 0.99, "plus": 0.06, "minus": 0.05},
      cl="2 sigma", source="Alvey, Sabti, Escudero, Fairbairn, EPJC 80, 148 (2020), arXiv:1910.10730, DOI 10.1140/epjc/s10052-020-7727-y",
      read_level="abstract",
      implication=f"m_q > ~{b['bbn_m_min_eV_kappa0.1_dG0.06']:.0e} eV [E]; trivially met by natural m_q")
entry(id="G_rec", observable="G at recombination (CMB + BAO + BBN), relative to laboratory G_N",
      construction_quantity="K(q_rec)^-1 vs local G (G_lab = (4/3) G_* if q light, G_* if heavy)",
      value="lambda_G = sqrt(G/G_N) = 1.0056 +- 0.0091, i.e. G/G_N = 1.011 +- 0.018",
      numbers={"lambda_G": 1.0056, "sigma": 0.0091, "G_over_GN": 1.0056**2, "G_over_GN_sigma": 2 * 1.0056 * 0.0091},
      cl="68%", source="Lamine et al., A&A 697, A109 (2025), arXiv:2407.15553 (Table 1, P20 TTTEEE + DESI DR1 BAO + BBN)",
      read_level="section (Table 1 and abstract)",
      also=[{"value": "Delta = -0.022 +- 0.023 (68%); |G_rad - G_now|/G < 3% (95%) in scalar-tensor models",
             "source": "Ballardini, Finelli, Sapone, JCAP 06 (2022) 004, arXiv:2111.09168", "read_level": "abstract"}],
      implication="a light q would give G_cos/G_lab = 3/4 (excluded at > 10 sigma); heavy branch: no difference")
entry(id="Yukawa_EotWash", observable="Yukawa deviation from 1/r^2 (torsion balance)",
      construction_quantity="q-mode Yukawa: strength alpha = 1/3, range 1/m_q (equations 4.3)",
      value="|alpha| = 1 excluded for lambda >= 38.6 um", numbers={"lambda_max_um_alpha1": 38.6},
      cl="95%", source="Lee, Adelberger, Cook, Fleischer, Heckel, PRL 124, 101101 (2020), arXiv:2002.11761",
      read_level="abstract",
      implication=f"m_q > ~{b['eotwash_m_min_eV_alpha1_point']*1e3:.1f} meV at alpha = 1; for alpha = 1/3 the bound is "
                  "somewhat weaker (~4 meV) [E: alpha = 1/3 crossing not read from the exclusion curve]")
entry(id="PPN_gamma_Cassini", observable="PPN gamma (Shapiro delay, Cassini conjunction)",
      construction_quantity="gamma(r) = (1 - e^{-m r}/3)/(1 + e^{-m r}/3); gamma = 1/2 for m r << 1",
      value="gamma - 1 = (2.1 +- 2.3)e-5", numbers={"central": 2.1e-5, "sigma": 2.3e-5},
      cl="1 sigma", source="Bertotti, Iess, Tortora, Nature 425, 374 (2003), DOI 10.1038/nature01997",
      read_level="abstract (publisher page)",
      also=[{"value": "gamma - 1 = (1.7 +- 1.6)e-4", "source": "Biskupek et al. 2021, arXiv:2012.12032 (LLR)", "read_level": "abstract"}],
      implication=f"m_q > {b['cassini_m_min_eV_b=R_sun']:.1e} eV (impact parameter R_sun) or "
                  f"{b['cassini_m_min_eV_b=1.6R_sun']:.1e} eV (1.6 R_sun) unless screened; excludes every ultralight q")
entry(id="DESI_w0wa", observable="dark-energy equation of state, w0-wa (CPL)",
      construction_quantity="w of rho_V(t); heavy branch w = -1 up to O(H^2/m_q^2)",
      value="DESI+CMB: w0 = -0.42 +- 0.21, wa = -1.75 +- 0.58 (3.1 sigma); +Pantheon+: -0.838 +- 0.055, -0.62 (+0.22 -0.19) (2.8 sigma); "
            "+Union3: -0.667 +- 0.088, -1.09 (+0.31 -0.27) (3.8 sigma); +DESY5: -0.752 +- 0.057, -0.86 (+0.23 -0.20) (4.2 sigma)",
      numbers={"DESI+CMB": {"w0": -0.42, "w0_sig": 0.21, "wa": -1.75, "wa_sig": 0.58, "pref_sigma": 3.1},
               "DESI+CMB+PantheonPlus": {"w0": -0.838, "w0_sig": 0.055, "wa": -0.62, "wa_plus": 0.22, "wa_minus": 0.19, "pref_sigma": 2.8},
               "DESI+CMB+Union3": {"w0": -0.667, "w0_sig": 0.088, "wa": -1.09, "wa_plus": 0.31, "wa_minus": 0.27, "pref_sigma": 3.8},
               "DESI+CMB+DESY5": {"w0": -0.752, "w0_sig": 0.057, "wa": -0.86, "wa_plus": 0.23, "wa_minus": 0.20, "pref_sigma": 4.2}},
      cl="68% (posterior means); preference from Delta chi^2_MAP (Table 6)",
      source="DESI Collaboration, DESI DR2 Results II, PRD 112, 083515 (2025), arXiv:2503.14738v3, Tables 5-6",
      read_level="section (Tables 5 and 6, HTML)",
      implication=f"heavy branch sits at the LambdaCDM point (rho_V shift today {b['heavy_branch_rhoV_shift_over_rhoL_at_cassini_bound']:.0e} "
                  "rho_Lambda even at the Cassini bound on m_q); it cannot produce the preferred w0 > -1, wa < 0 without screening")
entry(id="LCDM_Omega_m", observable="matter density Omega_m in flat LambdaCDM",
      construction_quantity="background matter density with a constant rho_V",
      value="DESI: 0.2975 +- 0.0086; DESI+CMB: 0.3027 +- 0.0036 (H0 = 68.17 +- 0.28); Planck 2018: 0.315 +- 0.007",
      numbers={"DESI": [0.2975, 0.0086], "DESI+CMB": [0.3027, 0.0036], "Planck2018": [0.315, 0.007]},
      cl="68%", source="arXiv:2503.14738v3 Table 5; Planck 2018 VI, A&A 641, A6 (2020), arXiv:1807.06209 abstract",
      read_level="section (DESI Table 5); abstract (Planck)",
      implication="calibration input for any rho_Lambda used by the construction; not a test")
entry(id="Omega_K", observable="spatial curvature Omega_K",
      construction_quantity="curvature of the FRW solution (the construction assumes k = 0)",
      value="DESI+CMB (LambdaCDM+Omega_K): 0.0023 +- 0.0011; w0wa+Omega_K with SNe: 0.0011/0.0006/0.0008 +- 0.0013",
      numbers={"DESI+CMB": [0.0023, 0.0011], "w0wa_PantheonPlus": [0.0011, 0.0013], "w0wa_Union3": [0.0006, 0.0013],
               "w0wa_DESY5": [0.0008, 0.0013]},
      cl="68%", source="arXiv:2503.14738v3 Table 5", read_level="section (Table 5, HTML)",
      implication="no curvature statement follows from P_q; the Candidate 1 exclusion premise of flat FRW is consistent")
entry(id="EDE_recombination", observable="early dark energy fraction near recombination (axion-like EDE, n = 3)",
      construction_quantity="pre-transition vacuum plateaus Delta + (4/3)Delta^2/(K m_q^2) (equations 5.3)",
      value="f_EDE < 0.088 (ACT), < 0.12 (P-ACT; P-ACT-LB), Planck alone < 0.087 (Hill et al. 2020)",
      numbers={"ACT": 0.088, "P-ACT": 0.12, "P-ACT-LB": 0.12, "Planck_Hill2020": 0.087, "log10_zc_prior": [3.0, 4.3]},
      cl="95%", source="ACT DR6, Calabrese et al., arXiv:2503.14454v2, eq. (14), (16); Hill et al., PRD 102, 043507 (2020), arXiv:2003.07355 (number via ACT DR6 text)",
      read_level="section (ACT DR6 EDE section, HTML); abstract (Hill et al.)",
      also=[{"value": "constant early DE: Omega_e < 0.0071 (Planck TT+lensing+BSH), < 0.0036 (Planck TT,TE,EE+lowP+BSH)",
             "cl": "95%", "source": "Planck 2015 XIV, A&A 594, A14 (2016), arXiv:1502.01590, conclusions (TeX)",
             "read_level": "section (conclusions)"}],
      implication="no SM or Candidate 1 plateau survives below T ~ 150 MeV; rho_Lambda/rho_r ~ 1e-35 at 1 MeV; limits met trivially, not a test")
entry(id="EDE_BBN", observable="expansion rate at BBN (extra energy density, N_nu)",
      construction_quantity="any vacuum or q-oscillation energy present at T ~ 1 MeV",
      value="N_nu = 2.898 +- 0.141 (BBN+CMB), N_nu < 3.180 (2 sigma)",
      numbers={"N_nu": 2.898, "sigma": 0.141, "N_nu_2sigma_max": 3.180,
               "equivalent_max_fraction_at_1MeV": b["bbn_max_extra_fraction"]},
      cl="68% / 2 sigma", source="Yeh, Shelton, Olive, Fields, JCAP 10 (2022) 046, arXiv:2207.13133",
      read_level="abstract",
      also=[{"value": "N_eff = 2.86 +- 0.13 (ACT DR6), 2.89 +- 0.11 with BBN", "source": "arXiv:2503.14454 abstract",
             "read_level": "abstract"},
            {"value": "N_eff = 2.99 +- 0.17 (Planck+BAO)", "source": "arXiv:1807.06209 abstract", "read_level": "abstract"}],
      implication=f"an extra component at 1 MeV must be < ~{b['bbn_max_extra_fraction']:.3f} of the total [E: treated as extra radiation]; "
                  "the construction has none")
entry(id="Omega_c", observable="cold dark matter density",
      construction_quantity="dust-like q-oscillation energy (equations 5.3; KV 0806.2805 V.E; 1612.02326)",
      value="Omega_c h^2 = 0.120 +- 0.001", numbers={"central": 0.120, "sigma": 0.001},
      cl="68%", source="Planck 2018 VI, A&A 641, A6 (2020), arXiv:1807.06209", read_level="abstract",
      implication="transition-deposited q-oscillation energy is suppressed by (x/sinh x)^2 with x ~ m_q/H ~ 1e30 (zero); "
                  "a primordial amplitude would be an input, calibrated if matched to 0.120; a Planck-mass oscillation "
                  f"decays in ~{b['planck_mass_q_lifetime_s_estimate']:.0e} s if any non-conformal matter couples [E]")
entry(id="ultralight_scalar_DM", observable="oscillating (ultralight) scalar dark-matter component",
      construction_quantity="q-oscillation with m_q in the ultralight range",
      value="m_a > 2e-20 eV for all DM (Lyman-alpha); percent-level fraction limits for 1e-33 <= m_a <= 1e-24 eV",
      numbers={"m_min_eV_all_DM": 2e-20, "fraction_limit_mass_range_eV": [1e-33, 1e-24]},
      cl="95%", source="Rogers, Peiris, PRL 126, 071302 (2021), arXiv:2007.12705; Hlozek, Marsh, Grin, MNRAS 476, 3063 (2018) [R journal], arXiv:1708.05681",
      read_level="abstract",
      implication="moot for an unscreened q-mode: Cassini already requires m_q > ~3e-15 eV; above that the "
                  "oscillation is CDM-like on observed scales and only Omega_c h^2 bounds its amount")
entry(id="spectral_distortion", observable="CMB spectral distortions (energy injection into photons)",
      construction_quantity="vacuum-energy decay into photons (only in dissipative q-theory variants)",
      value="|y| < 1.5e-5, |mu| < 9e-5", numbers={"y_max": 1.5e-5, "mu_max": 9e-5},
      cl="95%", source="Fixsen et al., ApJ 473, 576 (1996), astro-ph/9605054", read_level="abstract",
      also=[{"value": "T(z) = T0 (1+z)^(1-beta), beta_max ~ 3e-3 for vacuum decay into CMB photons",
             "source": "Opher, Pelinson, MNRAS 362, 167 (2005), astro-ph/0409451", "read_level": "abstract"}],
      implication="the non-dissipative construction (m_q >> H) has no late decay; KV's dissipative variants "
                  "(0905.1919 eq 4.1; 1601.00601) produce gravitons/neutrinos, bounded by N_eff (EDE_BBN entry)")

ledger = {"track": "ledger", "date": "2026-10-21",
          "note": "All data already exposed; every comparison is retrospective. Read levels as achieved in this track.",
          "entries": L}

num["checks"] = checks
num["n_checks"] = len(checks)
num["n_pass"] = sum(c["pass"] for c in checks)
num["versions"] = {"python": platform.python_version(), "numpy": np.__version__}

with open(OUT + "ledger_numbers.json", "w") as f:
    json.dump(num, f, indent=1)
with open(OUT + "ledger.json", "w") as f:
    json.dump(ledger, f, indent=1)
print(f"checks {num['n_pass']}/{num['n_checks']} pass; ledger entries {len(L)}")
for c in checks:
    print(("PASS " if c["pass"] else "FAIL ") + c["name"], c["detail"])
