"""Part 6: consistency of Codex's production replays (135d535 / 479f111, read with git show; not rerun)
with the fixed-|X| and fixed-n phase structure derived here.  Also a trivial numerical Z2 evenness check.
Run after t0_phase.py and z2_memory.py:  python3 codex_consistency.py
"""
import json
import os
import numpy as np
from common import A, B, CODEX_PRODUCTION as CP, ARCHIVED_SOLITON, RESULTS_DIR, LAM_PHI, V_GEV, GEV_TO_FM, save_json
import thermal_phase as tp

a, b = A, B
t0 = json.load(open(os.path.join(RESULTS_DIR, "t0_phase.json")))
zm = json.load(open(os.path.join(RESULTS_DIR, "z2_memory.json")))
HERE = os.path.dirname(os.path.abspath(__file__))
prof = np.loadtxt(os.path.normpath(os.path.join(HERE, "../../../../grut_candidate_1_20261018/code/results_profile.csv")),
                  delimiter=",", skiprows=1)
rp, Fp, Sp = prof[:, 0], prof[:, 1], prof[:, 3]


def S_minima_fixed_x(x):
    """minima of Vbar(S, x) in S at fixed x (fixed |X|)."""
    if a * x < 1:
        return [np.sqrt(1 - a * x), -np.sqrt(1 - a * x)]
    return [0.0]


# 1. fixed-|X| threshold
xs = np.linspace(0, 3, 3001)
nmin = np.array([len(S_minima_fixed_x(x)) for x in xs])
x_threshold = xs[np.argmax(nmin == 1)]
endpoint_rows = []
for k, name in enumerate(["coarse h=1", "medium h=2/3", "fine h=1/2"]):
    q0 = CP["endpoint_negative_F2_q0"][k]
    qs = CP["endpoint_negative_F2_quartiles"][k]
    endpoint_rows.append(dict(grid=name, min_F2_negative_cells=q0, F2_minus_threshold=q0 - 1 / a,
                              S_mass2_at_min_F2=a * q0 - 1, F2_quartiles=qs,
                              S_mass2_at_median=a * qs[1] - 1,
                              median_over_x_star=qs[1] / t0["binodal_x_star"],
                              every_step_min_meanF2_negative=CP["min_meanF2_on_negative_cells_every_step"][k],
                              minS_every_step=CP["minS_every_step"][k],
                              equivalent_single_ball_radius_of_negative_cells=float(
                                  (3 * CP["endpoint_negative_cells"][k] * CP["endpoint_retained_spacing"][k] ** 3
                                   / (4 * np.pi)) ** (1 / 3))))

# 2. archived soliton geometry
r_F2_1 = float(rp[np.argmax(Fp**2 < 1)])
r_S_half = float(rp[np.argmax(Sp > 0.5)])

# 3. R_* comparisons (from z2_memory.json)
sudden = {(s["x_after"], s["eps"]): s["R_first_flip"] for s in zm["sudden_quench_scans"]}
prev_R = {}
for s in zm["sudden_quench_scans"]:
    rows = sorted(s["rows"], key=lambda r: r["R"])
    last_inherit = max([r["R"] for r in rows if not r["flipped"] and (s["R_first_flip"] is None or r["R"] < s["R_first_flip"])],
                       default=None)
    prev_R[(s["x_after"], s["eps"])] = last_inherit
arch = zm["archived_soliton_charge_removed"]

# 4. Z2 evenness of the leading-order thermal free energy (trivial by construction)
phis = np.linspace(0, 1.2, 25)
dev = 0.0
for T in [0.0, 0.5, 1.2]:
    for p in phis:
        dev = max(dev, abs(tp.UT(p, T) - tp.UT(-p, T)) + abs((tp.g2 * p * p) - (tp.g2 * (-p) * (-p))))

cpl = t0["codex_production_linear"]
out = dict(
    status="derived (fixed-|X| threshold), quoted (Codex numbers via git show 135d535/479f111), computed (toy R_*)",
    fixed_X_threshold_F2=float(x_threshold), fixed_X_threshold_exact="F^2 = 1/a = lam_phi/gX^2 = 1 (|X|^2 = lam_phi v^2/(2 gX^2))",
    statement="At fixed |X|, Vbar(S,x) = (S^2-1)^2/4 + a x S^2/2 + const has minima S=+-sqrt(1-a x) for a x<1 and the "
              "single minimum S=0 for a x>=1. F^2>1 is exactly the region where phi=0 is the unique fixed-|X| minimum.",
    codex_endpoint_negative_cells=endpoint_rows,
    production_background=dict(x=CP["F2"], S=CP["S"], Omega=CP["Omega"], n_over_nc=cpl["n_over_nc"],
                               branch="phi!=0 (branch 1), fragmentation-unstable",
                               cs2=cpl["cs2_longwave"], gamma_max=cpl["gamma_max"],
                               k_at_gamma_max=cpl["k_at_gamma_max"], wavelength=cpl["wavelength_at_gamma_max"],
                               t_nl_estimate=cpl["t_nonlinear_estimate_ln_1_over_eps"],
                               codex_first_negative_time=CP["first_negative_time"],
                               box_side=CP["box_side"], unstable_harmonics_per_axis=cpl["n_unstable_box_harmonics_per_axis"]),
    archived_soliton=dict(r_where_F2_eq_1=r_F2_1, r_where_S_eq_half=r_S_half, R90=ARCHIVED_SOLITON["R90bar"],
                          F2_center=ARCHIVED_SOLITON["F0"] ** 2, S_center=ARCHIVED_SOLITON["S0"],
                          charge_removed_toy=arch),
    R_star=dict(dirichlet_full_charge_loss=np.pi, dirichlet_full_charge_loss_fm=np.pi / (np.sqrt(LAM_PHI) * V_GEV) * GEV_TO_FM,
                sudden_coherent_seed_first_flip_R={"x_after=%g eps=%g" % k: v for k, v in sudden.items()},
                sudden_coherent_seed_last_inherit_R={"x_after=%g eps=%g" % k: v for k, v in prev_R.items()},
                ramp_phase_randomised=[dict(tauQ=s["tauQ"], R=s["R"], frac_negative=s["fraction_negative_center_decision"])
                                       for s in zm["ramp_quench_scans"]]),
    patch_vs_Rstar="Archived/Codex-scale charged patches have phi~0 radius ~3.3 ell (F^2>1 out to 3.38). This is "
                   "comparable to the bare Dirichlet pi/m_tach = 3.14 ell (full charge loss), and below the toy "
                   "sign-memory radius (first flip at R=8 ell for a coherent 1e-2 seed, 14 ell for 1e-4). The archived "
                   "profile with its charge removed instantly inherits +v for coherent negative seeds up to 0.3.",
    verdict="consistent: Codex's negative cells sit where F^2>1, i.e. where phi=0 is the unique fixed-|X| minimum and "
            "the local phi mass^2 a F^2 - 1 > 0; their F^2 quartiles (3.0-4.0) match Q-matter (x*=3.16) and the archived "
            "core (3.92). The shallow negative phi (>-0.30) are oscillations about phi~0 inside charge-restored "
            "fragments, not opposite-vacuum domains. These patches are within the sign-memory radius.",
    z2_evenness_max_deviation=dev,
)
save_json("codex_consistency.json", out)
if __name__ == "__main__":
    print(json.dumps(out, indent=1, default=float)[:6000])
