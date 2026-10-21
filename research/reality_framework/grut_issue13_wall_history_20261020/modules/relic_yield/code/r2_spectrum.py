"""Item 2 (part): spectrum, exact stability and gravitational decay estimates in the +v vacuum.
Writes results/spectrum_stability.json.

All rates are ESTIMATES unless stated. Inputs marked 'recalled' were not fetched in this session.
"""
import math

from common import (GEV_TO_INV_S, GX, HD, LAM_PHI, M_D, M_H, M_X, MD0, MP, MU0, T_UNIVERSE_S, V, dump)

out = dict(status="derived (masses, thresholds, symmetry stability) + estimate (gravitational decay rates)")
out["masses_GeV"] = dict(m_X=M_X, m_h=M_H, m_D=M_D, formulas="m_X = gX v, m_h = sqrt(2 lam_phi) v, m_D = sqrt(mD0^2 + hD^2 v^2)")
out["thresholds_GeV"] = dict(two_m_X=2 * M_X, two_m_D=2 * M_D, h_to_XXbar_open=M_H > 2 * M_X, h_to_DD_open=M_H > 2 * M_D)
out["broken_phase_vertices"] = dict(
    hhh="lam_phi v h^3 (Feynman 6 lam_phi v = %.4g GeV)" % (6 * LAM_PHI * V),
    hhhh="lam_phi/4 h^4 (6 lam_phi = %.4g)" % (6 * LAM_PHI),
    hXX="2 gX^2 v h |X|^2 (vertex %.4g GeV)" % (2 * GX ** 2 * V),
    hhXX="gX^2 h^2 |X|^2 (vertex %.4g)" % (2 * GX ** 2),
    hDD="hD^2 v h D^2 (vertex %.4g GeV)" % (2 * HD ** 2 * V),
    hhDD="hD^2/2 h^2 D^2 (vertex %.4g)" % (2 * HD ** 2),
    note="symmetric phase (phi = 0) has no cubic vertex: number changes only by 2->4 there (thermal_potential)")
out["stability"] = [
    dict(species="X, Xbar", mass=M_X, verdict="exactly stable in the action", reason="lightest U(1)_X-charged state; U(1)_X exact"),
    dict(species="D", mass=M_D, verdict="exactly stable in the action", reason="lightest Z2_D-odd state; Z2_D exact"),
    dict(species="h (phi radial mode)", mass=M_H, verdict="no non-gravitational decay; gravitational lifetime >> age",
         reason="h -> X Xbar needs m_h > 200 GeV, h -> D D needs m_h > 201.0 GeV (both closed); Higgs and other "
                "portals are zero at mu_0 and stay zero under non-gravitational running (hidden_heating)"),
]

# ---------------------------------------------------------------- gravitational decay estimates of h
G_naive = M_H ** 5 / MP ** 4
# (i) loop-induced h Riem^2 from heavy-field one-loop (Schwinger-DeWitt a_2 Riem^2/180 per real scalar):
#     L = (c/v) h R_{abcd}R^{abcd}, c = sum_i n_i v d ln m_i^2/dh /(32 pi^2 180); Gamma(h -> gg) = c^2 m_h^7/(4 pi v^2 M_P^4)
contrib = dict(X=2 * 2.0, D=2 * HD ** 2 * V ** 2 / M_D ** 2, h_self=3.0)
c_loop = sum(contrib.values()) / (32 * math.pi ** 2 * 180)
G_gg = c_loop ** 2 * M_H ** 7 / (4 * math.pi * V ** 2 * MP ** 4)
# (ii) nonminimal xi R phi^2 -> -xi v h R -> (via R = -T/M_P^2) h couples to the trace T with 1/f = xi v / M_P^2:
#      a dilaton-like coupling; Gamma ~ (xi v v_H / M_P^2)^2 Gamma_SM-Higgs(m_h) x A (A <= 10 anomaly allowance)
V_H = 246.22
GAMMA_SM_H_141 = 8.6e-3   # GeV, SM-Higgs-like total width at 141 GeV (recalled from LHC HXSWG tables, not fetched)
beta_xi_phi = -1.0575e-4  # per e-fold (hidden_heating module, quoted)
xi_run_mh = beta_xi_phi * math.log(M_H / MU0)
rows_xi = []
for xi in (1e-5, 1e-4, 1e-3):
    kap = xi * V * V_H / MP ** 2
    G = kap ** 2 * GAMMA_SM_H_141
    rows_xi.append(dict(xi_eff=xi, kappa=kap, Gamma_GeV=G, Gamma_upper_GeV=10 * G,
                        tau_s=1 / (G * GEV_TO_INV_S), tau_lower_s=1 / (10 * G * GEV_TO_INV_S)))


def tau(G):
    return 1.0 / (G * GEV_TO_INV_S)


out["h_decay"] = dict(
    naive_dimensional=dict(Gamma_GeV=G_naive, tau_s=tau(G_naive),
                           note="m_h^5/M_P^4 with no small factors: generous upper estimate of any gravitational 2-body rate"),
    loop_h_to_graviton_pair=dict(contributions_v_dlnm2_dh=contrib, c=c_loop, Gamma_GeV=G_gg, tau_s=tau(G_gg),
                                 band_tau_s=[tau(10 * G_gg), tau(0.1 * G_gg)],
                                 note="heavy-loop low-energy theorem; m_h^2/(4 m_X^2) = 0.5 so O(1)-O(10) form-factor "
                                      "uncertainty; tree-level h -> gg vanishes for minimal coupling (no h-linear vertex; "
                                      "h R term vanishes on shell)"),
    xi_trace_coupling=dict(xi_phi_running_at_m_h=xi_run_mh, Gamma_SM_Higgs_141GeV_recalled=GAMMA_SM_H_141, rows=rows_xi,
                           note="xi_phi = 0 at mu_0 by matching; |xi_eff| ~ 1e-5..1e-3 band covers running + finite parts"),
    age_of_universe_s=T_UNIVERSE_S,
    verdict="every estimate gives tau_h > 1e36 s >> 4.4e17 s: h is a cosmologically stable relic species of the action")

# ---------------------------------------------------------------- outside-action comparisons (NOT in the action)
G_dim5 = M_H ** 3 / (64 * math.pi * MP ** 2)              # (phi/(4 M_P)) F F, Z2_phi-odd: forbidden by exact Z2_phi
G_dim6 = (2 * V) ** 2 * M_H ** 3 / (64 * math.pi * MP ** 4)  # phi^2 F F/(4 M_P^2), zero at mu_0 by matching
dlam = 3.4e-34                                          # graviton-induced Higgs-portal estimate (hidden_heating, quoted)
M_HIGGS = 125.25
theta = dlam * V * V_H / abs(M_H ** 2 - M_HIGGS ** 2)
G_mix = theta ** 2 * GAMMA_SM_H_141
out["outside_action_comparison"] = dict(
    dim5_Z2odd_phiFF_over_MP=dict(Gamma_GeV=G_dim5, tau_s=tau(G_dim5),
                                  note="NOT in the action (Z2_phi exact by postulate P4); such an operator would also bias the walls"),
    dim6_phi2FF_over_MP2=dict(Gamma_GeV=G_dim6, tau_s=tau(G_dim6), note="zero at mu_0 by matching; shown for scale"),
    graviton_induced_higgs_portal=dict(delta_lambda=dlam, mixing=theta, Gamma_GeV=G_mix, tau_s=tau(G_mix)))
dump("spectrum_stability.json", out)
print("tau naive %.3g s, loop gg %.3g s, xi(1e-4) %.3g s, dim5(outside) %.3g s, dim6 %.3g s, portal %.3g s" % (
    tau(G_naive), tau(G_gg), rows_xi[1]["tau_s"], tau(G_dim5), tau(G_dim6), tau(G_mix)))
