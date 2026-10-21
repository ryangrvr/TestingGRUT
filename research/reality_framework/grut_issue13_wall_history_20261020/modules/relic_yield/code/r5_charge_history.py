"""Items 3-5 (support): hidden charge in a thermal hidden plasma, number needed for thermal restoration,
charge supplied at T_Q into a cold hidden sector, the floor inequality with a symmetric relic, and the
FLS-solitogenesis prior-art table. Reads results/boltzmann.json. Writes results/charge_history.json.
"""
import json
import math
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from common import (GX, LAM_X, M_B_HIGH, M_P_PROTON, M_X, MASER_PRODUCT_PRIMARY, MASER_PRODUCT_STRESS, MSTAR,
                    N_C_PLACEHOLDER, N_S1, OMEGA_B, OMEGA_B_LOW, OMEGA_DM_KNOWN, OMEGA_PER_YIELD, RESULTS, S_SM_COEF,
                    TC_LEADING, TC_RESUMMED, V, Y_B, ZETA3, dump)

BZ = json.load(open(os.path.join(RESULTS, "boltzmann.json")))
out = dict(status="derived/computed; placeholders and outside-action inputs labelled")
PI_X = GX ** 2 / 12 + LAM_X / 3       # symmetric-phase Debye coefficient per real X mode (thermal_potential, derived)
MX_T = math.sqrt(PI_X)                # M_X(T)/T at phi = 0


# ---------------------------------------------------------------- (4) charge in a thermal plasma
def nQ_over_T3(mu_over_T, M_over_T):
    """Exact BE charge density of one complex scalar of mass M: n_Q/T^3 (|mu| <= M)."""
    def f(p):
        E = math.sqrt(p * p + M_over_T ** 2)
        a = E - mu_over_T
        b = E + mu_over_T
        return p * p * (1.0 / math.expm1(a) - 1.0 / math.expm1(b)) if a > 0 else 0.0
    val, _ = quad(f, 0.0, 60.0, limit=400, points=[1e-6, 0.1, 1.0])
    return val / (2 * math.pi ** 2)


fmax = nQ_over_T3(MX_T * (1 - 1e-12), MX_T)
fsmall = nQ_over_T3(1e-3 * MX_T, MX_T) / (1e-3 * MX_T)
out["thermal_charge_capacity"] = dict(
    M_X_over_T_symmetric=MX_T, nQ_max_over_T3=fmax, leading_estimate_M_over_3T=MX_T / 3,
    dnQ_dmu_over_T2_small_mu=fsmall, massless_value_1_over_3=1 / 3,
    note="a complex scalar holds at most n_Q,max = n_Q(mu = M) in the gas; excess charge Bose-condenses")

rows4 = []
coefQ = Y_B * S_SM_COEF          # n_Q = r * coefQ * T_SM^3
for r in (1.0, 0.639344, 0.0897, 1e-2, 1e-3):
    xi_bec = (r * coefQ / fmax) ** (1 / 3)
    for xi in (BZ["xi_max"]["all_species"], 1e-3, 1e-2, 0.1, 1.0):
        need = r * coefQ / xi ** 3          # n_Q/T_h^3
        if need < fmax:
            mu = brentq(lambda m: nQ_over_T3(m, MX_T) - need, 0.0, MX_T * (1 - 1e-12))
            cond = 0.0
        else:
            mu = MX_T
            cond = 1 - fmax / need
        rows4.append(dict(r=r, xi=xi, nQ_over_Th3=need, mu_over_T=mu, condensed_fraction_at_Tc=cond, xi_BEC=xi_bec))
out["charge_in_thermal_plasma"] = dict(
    rows=rows4, xi_BEC_formula="xi_BEC = (r Y_B (2 pi^2/45) g*s / f_max)^(1/3) = %.4g r^(1/3)" % ((coefQ / fmax) ** (1 / 3)),
    note="evaluated at T_h ~ T_c in the restored phase (M = sqrt(Pi_X) T). For xi < xi_BEC the excess charge forms a "
         "homogeneous condensate whose phase structure is the charge_phase module's (H1/H2)")

# charge holds phi = 0 after thermal restoration ends?  n_Q(T_SM = T_c/xi) vs n_c (placeholder)
Tn = lambda r, n_thr: (n_thr * V ** 3 / (r * coefQ)) ** (1 / 3)   # GeV
rows_hold = []
for r in (1.0, 0.639344, 0.0897, 1e-2, 1e-3):
    for Tc in (TC_LEADING, TC_RESUMMED[0]):
        xi_hold = Tc * V / Tn(r, N_C_PLACEHOLDER)
        rows_hold.append(dict(r=r, Tc_over_v=Tc, T_n_GeV=Tn(r, N_C_PLACEHOLDER), xi_below_which_charge_keeps_phi_0=xi_hold,
                              nQ_over_nc_at_restoration_end_xi_max=r * coefQ * (Tc * V / BZ["xi_max"]["all_species"]) ** 3 / (N_C_PLACEHOLDER * V ** 3)))
out["charge_after_thermal_restoration"] = dict(
    rows=rows_hold,
    note="if T_c/xi > T_n(r), i.e. xi < T_c/T_n(r), the charge density at the end of thermal restoration exceeds n_c "
         "(cold-sector placeholder) and phi = 0 persists until T_SM = T_n(r); the exit is then the charge route of "
         "charge_phase (spinodal at n_c, p(+) = p(-))")

# solitogenesis scale estimate (charge_phase KZ estimate quoted)
ell = 1.0 / 100.0
for r, xiKZ in ((1.0, 850.0), (0.01, 300.0)):
    out.setdefault("solitogenesis_estimate", []).append(dict(
        r=r, xi_KZ_ell=xiKZ, charge_per_correlation_volume=N_C_PLACEHOLDER * V ** 3 * (xiKZ * ell) ** 3,
        thin_wall_binding_threshold_Q=1.3e4,
        note="estimate: Q per KZ volume at the spinodal exit (charge_phase xi_KZ, quoted) >> binding threshold, so "
             "charge-restored regions can leave bound Q-balls; with exact Z2 the surrounding vacuum signs are random"))

# ---------------------------------------------------------------- (3 support) number needed for thermal restoration


def Li(s, z, kmax=20000):
    if z >= 1.0:
        return {2: math.pi ** 2 / 6, 3: ZETA3}[s]
    k = np.arange(1, kmax + 1, dtype=float)
    return float(np.sum(z ** k / k ** s))


zs = [1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]
Fz = [Li(2, z) / Li(3, z) ** (2 / 3) for z in zs]
run = BZ["xi_max"]["run_at_xi_max"]
Ytot = sum(run["Y"])
mbar = run["Omega_h2"] / (Ytot * OMEGA_PER_YIELD)
YNmax = OMEGA_DM_KNOWN / (mbar * OMEGA_PER_YIELD)
rows3 = []
for label, Tc in (("leading", TC_LEADING), ("resummed_low", TC_RESUMMED[0]), ("resummed_high", TC_RESUMMED[1])):
    n_th = 4 * ZETA3 / math.pi ** 2 * (Tc * V) ** 3
    Tmin_number = (n_th / (S_SM_COEF * YNmax)) ** (1 / 3)
    Tmin_cond = (N_C_PLACEHOLDER * V ** 3 / (S_SM_COEF * YNmax)) ** (1 / 3)
    rows3.append(dict(Tc_label=label, Tc_over_v=Tc, n_thermal_at_Tc_v3=n_th / V ** 3, T_SM_min_thermal_number_GeV=Tmin_number,
                      T_SM_min_condensate_placeholder_GeV=Tmin_cond))
out["restoration_number_bound"] = dict(
    F_of_z=dict(z=zs, F=Fz, monotonic_increasing=bool(all(np.diff(Fz) > 0)),
                meaning="at fixed number density, Pi_phi ~ n^(2/3) Li2(z)/Li3(z)^(2/3) is largest at z = 1 (mu = 0): "
                        "a bath with mu_N < 0 needs MORE quanta per T_c^3 to restore"),
    mean_mass_at_xi_max_GeV=mbar, Y_N_max=YNmax, rows=rows3,
    note="T_SM,min = minimum visible temperature at which a hidden population with relic yield <= Y_N_max can restore "
         "Z2 thermally (number conserved after restoration because 3->2 is inactive for xi <~ 3e-3, boltzmann.json). "
         "The condensate column uses the cold charged-condensate density n_c as a PLACEHOLDER for a coherent X "
         "population and is an estimate only")

# ---------------------------------------------------------------- (2/5) floor inequality with a symmetric relic
coef_sym = (M_X / MSTAR) * (M_B_HIGH / OMEGA_B_LOW)     # GeV per unit omega_sym
omega_asym_per_r = MSTAR * Y_B * OMEGA_PER_YIELD
out["floor_inequality"] = dict(
    statement="omega_h >= (m*|r|/m_b) omega_b + omega_sym (late-time dilute population: bound charge carries >= m* per "
              "unit, every other hidden quantum is a free NR particle of mass >= m_X); rigorous weaker form "
              "omega_h >= max[(m*|r|/m_b) omega_b, omega_sym]. Hence H_0 >= 100 sqrt(omega_b + omega_h,min + omega_sym) "
              "and the maser-envelope necessary condition becomes m_X|r| + (m_X/m*)(m_b/omega_b) omega_sym < 63.93 GeV "
              "(primary), 74.64 GeV (stress), with the envelope's conservative omega_b = 0.02125, m_b = 0.94738 GeV",
    coefficient_GeV_per_unit_omega_sym=coef_sym,
    r_max_primary_no_sym=MASER_PRODUCT_PRIMARY / M_X, r_max_primary_sym_012=(MASER_PRODUCT_PRIMARY - coef_sym * OMEGA_DM_KNOWN) / M_X,
    r_max_stress_no_sym=MASER_PRODUCT_STRESS / M_X, r_max_stress_sym_012=(MASER_PRODUCT_STRESS - coef_sym * OMEGA_DM_KNOWN) / M_X,
    H0_floor_r1_no_sym=100 * math.sqrt(OMEGA_B_LOW * (1 + MSTAR / M_B_HIGH)),
    H0_floor_r1_sym_012=100 * math.sqrt(OMEGA_B_LOW * (1 + MSTAR / M_B_HIGH) + OMEGA_DM_KNOWN),
    omega_asym_min_per_unit_r=omega_asym_per_r,
    retrospective_DM_bound_r_max=OMEGA_DM_KNOWN / omega_asym_per_r,
    label="retrospective: omega_DM = 0.12 and the maser product bounds are known data (Codex 2ee7f82); not a prediction")

# ---------------------------------------------------------------- (5) charge supplied at T_Q into a cold hidden sector
rmax_p = MASER_PRODUCT_PRIMARY / M_X
rows5 = []
for label, r in (("r=1 benchmark", 1.0), ("maser primary r_max", rmax_p), ("maser stress r_max", MASER_PRODUCT_STRESS / M_X),
                 ("maser primary with omega_sym=0.12", (MASER_PRODUCT_PRIMARY - coef_sym * OMEGA_DM_KNOWN) / M_X),
                 ("retrospective DM (omega_h<=0.12)", OMEGA_DM_KNOWN / omega_asym_per_r), ("r=0.01", 0.01)):
    rows5.append(dict(label=label, r=r, T_Q_max_nc_placeholder_GeV=Tn(r, N_C_PLACEHOLDER), T_Q_max_ns1_GeV=Tn(r, N_S1)))
out["charge_supply_cold_sector"] = dict(
    condition="no charge-driven restoration iff n_Q(T_Q) = r Y_B s_SM(T_Q) < n_c, i.e. T_Q < T_n(r) = "
              "[n_c/(r Y_B (2 pi^2/45) g*s)]^(1/3) = %.4g r^(-1/3) (v/TeV) GeV" % Tn(1.0, N_C_PLACEHOLDER),
    n_c_used_v3=N_C_PLACEHOLDER, n_c_label="PLACEHOLDER (charge_phase T_h = 0 lower spinodal; its NOTES confirm 0.0316228 v^3)",
    rows=rows5,
    v_scaling="T_n ~ v r^(-1/3); with r_max = 63.93 GeV/(gX v): T_n(r_max) = %.4g GeV (v/TeV)^(4/3)" % Tn(rmax_p, N_C_PLACEHOLDER),
    note="valid for a homogeneously supplied charge in a cold (T_h << T_c) hidden sector with a pre-selected phi sign; "
         "supply below T_n is the charge_phase H2 case (fragmentation into Q-balls inside the sign-memory radius)")

# ---------------------------------------------------------------- prior art (read levels in NOTES.md)
out["prior_art"] = [
    dict(ref="Frieman, Gelmini, Gleiser, Kolb, PRL 60 (1988) 2101; SLAC-PUB-4538/FERMILAB-PUB-88-013-A",
         read="full preprint (OCR text of the Fermilab scan)", Z2_handling="(b)+(a)",
         detail="FLS potential with h|phi|^2(sigma - sigma0)^2 coupling and a sigma^3 term (a = lam2/lam1 = 0.15): Z2 "
                "explicitly broken, vacua non-degenerate; states that the exact-Z2 limit (h = lam2 = 0) gives stable "
                "domain walls that dominate unless sigma0 <~ 10 MeV; solitogenesis relies on the bias (pressure pinches "
                "false-vacuum regions); percolating case analysed only with bias"),
    dict(ref="Griest, Kolb, Massarotti, PRD 40 (1989) 3529; FERMILAB-PUB-89-110-A",
         read="preprint introduction and model section (OCR of scan; equations partly illegible)", Z2_handling="(b)",
         detail="same FGGK potential with lam2/lam1 = 0.15 (non-degenerate true vacuum); wall (surface) energy of the "
                "soliton ignored; statistical charge fluctuations as the charge source; no exact-Z2 wall analysis"),
    dict(ref="Frieman, Olinto, Gleiser, Alcock, PRD 40 (1989) 3241; FERMILAB-PUB-89-128-A",
         read="sections I-III and figure captions (OCR text) + abstract", Z2_handling="(b)",
         detail="generic NTS evolution (evaporation/accretion/freeze-out); formation discussed with schematic "
                "potentials having a unique true vacuum (metastable or unstable NTS phase at sigma = 0); 'pinch-off in "
                "the asymmetric vacuum model'; exact-Z2 domain walls not analysed"),
    dict(ref="Macpherson, Campbell, PLB 347 (1995) 205, arXiv:hep-ph/9408387",
         read="abstract + introduction and formation section", Z2_handling="(a)+(c)",
         detail="explicitly 'approximate' discrete symmetry: bias term A(phi) with vacuum splitting Lambda; notes exact "
                "Z2 walls would dominate (Zel'dovich-Kobzarev-Okun); analyses the biased (percolating) wall network "
                "following Gelmini-Gleiser-Kolb and its fragmentation into Fermi balls"),
    dict(ref="Krylov, Levin, Rubakov, PRD 87 (2013) 083528, arXiv:1301.0354",
         read="abstract + section 2 (model) incl. footnote 1", Z2_handling="(a) assumed, not analysed",
         detail="tree potential lam(phi^2 - v^2)^2 + h^2 phi^2 |chi|^2 (same structure as Candidate 1); footnote: the "
                "phi -> -phi symmetry makes the model non-viable because of the domain-wall problem; assume explicit "
                "breaking or a charged phi; 'irrelevant for what follows'"),
    dict(ref="Ponton, Bai, Jain, JHEP 09 (2019) 011, arXiv:1906.10739",
         read="abstract + section list + model section (searched)", Z2_handling="(b)",
         detail="confining field is the SM Higgs doublet (continuous gauged vacuum manifold, no discrete-vacuum walls); "
                "first-order EW transition with a singlet extension; no Z2 wall issue arises"),
    dict(ref="Bai, Long, Lu, PRD 99 (2019) 055047, arXiv:1810.04360",
         read="abstract + introduction (searched)", Z2_handling="(b)",
         detail="dark-QCD first-order confinement/chiral transition; no real Z2 scalar; walls mentioned only for "
                "other authors' axion quark nuggets"),
    dict(ref="Hong, Jung, Xie, PRD 102 (2020) 075028, arXiv:2008.04430",
         read="abstract + sections II and V (toy model)", Z2_handling="(b)",
         detail="generic Yukawa-trapped fermions; toy model U = (mu^2 + cT^2)phi^2/2 + mu3 phi^3/3 + lam phi^4/4 with a "
                "tree-level phi^3 (no Z2); a Z2 appears only as an R-parity-like stabiliser of the fermion, not of phi"),
]
out["prior_art_summary"] = (
    "None of the eight works analyses a persistent wall network of an EXACT Z2 in an FLS/Q-ball formation history. "
    "FGGK/GKM/FOGA, Ponton-Bai-Jain, Bai-Long-Lu and Hong-Jung-Xie use potentials without a degenerate Z2 (b); "
    "Macpherson-Campbell and Krylov-Levin-Rubakov keep the Z2 structure but break it explicitly (a), the former "
    "analysing the biased network (c). FGGK and Macpherson-Campbell state that the exact limit leads to the "
    "domain-wall problem. Candidate 1 forbids every such bias (BRIEF rule 5), so no prior-art mechanism removes its walls.")
dump("charge_history.json", out)
print(json.dumps({k: out[k] for k in ("thermal_charge_capacity", "floor_inequality")}, indent=1))
for r in rows3:
    print(r)
for r in rows5:
    print(r)
print(out["restoration_number_bound"]["F_of_z"])
