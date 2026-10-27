"""Implication checks, not a full RG fit or cosmological abundance solver."""
import math
import numpy as np
from common import ROOT, save
import json


def vacuum():
    delta = .01
    q0 = math.sqrt(2.)
    mu = q0

    def potential(s):
        return (s*s - 1)**2/4 + delta*(2 - 3*s + s**3)/4

    phases = []
    for name, s in (("false_present_candidate", -1.), ("matched_true_reference", 1.)):
        energy = q0*q0/2 + 1 + potential(s)
        phases.append({"name": name, "s": s, "matter_potential": potential(s),
                       "matter_derivative": (s*s-1)*(s+3*delta/4),
                       "matter_second_derivative": 3*s*s-1+1.5*delta*s,
                       "gravitating_density_epsilon_minus_mu_q": energy-mu*q0})
    save("VACUUM_COUNTEREXAMPLE.json", {
        "scope": "fixed-mu ground-reference equilibrium matching, not a computed future transition",
        "epsilon": "q^2/2+1+V(s)", "q0": q0, "mu": mu, "delta": delta,
        "epsilon_q_minus_mu": q0-mu, "epsilon_qq": 1., "phases": phases,
        "present_vacuum_equals_final_vacuum_is_additional": True,
        "actual_future_exit_history": None, "observed_dark_energy_prediction": None})


def portal():
    """Fixed light masses and constant-portal leading log; no running improvement."""
    mh, v, gamma_sm, inv_max = 125.20, 246.22, .0041, .107
    mu0, cutoff = 173., 2.435e18
    inp = json.loads((ROOT.parent / "grut_heavy_vector_budget_20261025/INPUTS.json").read_text())
    light = inp["light_sector"]
    mass = math.sqrt(2*light["lambdaPhi"]) * light["vPhi_GeV"]
    phase = math.sqrt(1 - 4*mass*mass/mh**2)
    allowed_width = inv_max*gamma_sm/(1-inv_max)
    log_span = math.log(cutoff/mu0)
    rows = []
    for count in (1, 2):
        width_coefficient = count*v*v*phase/(32*math.pi*mh)
        lp = math.sqrt(allowed_width/width_coefficient)
        beta = count*lp*lp/(32*math.pi**2)
        targets = []
        for desired in (1e-5, 1e-4, 1e-3, .01):
            coupling = math.sqrt(desired*32*math.pi**2/(count*log_span))
            width = width_coefficient*coupling**2
            targets.append({"illustrative_shift": desired, "constant_portal": coupling,
                            "invisible_branching_fraction": width/(gamma_sm+width)})
        rows.append({"real_scalar_count": count, "portal_max_in_this_width_proxy": lp,
                     "width_coefficient_GeV": width_coefficient,
                     "delta_beta_lambda_H_at_max": beta,
                     "frozen_leading_log_shift_at_max": beta*log_span,
                     "illustrative_targets_not_mass_fit": targets})
    # Independent Hessian construction for V=lambda_H*(h.h)^2/4+lp*(h.h)*S^2/4.
    # At S=0: 4 Higgs eigenvalues and N scalar eigenvalues, no mixed blocks.
    rng = np.random.default_rng(20261027)
    max_error = 0.
    for _ in range(20):
        h = rng.normal(size=4)
        lh, lp = .13, .27
        r2 = h @ h
        for count in (1, 2):
            hessian = np.zeros((4+count, 4+count))
            hessian[:4, :4] = lh*(r2*np.eye(4)+2*np.outer(h, h))
            hessian[4:, 4:] = lp*r2*np.eye(count)/2
            actual = np.trace(hessian @ hessian)
            expected = (12*lh**2 + count*lp**2/4)*r2**2
            max_error = max(max_error, abs(actual/expected-1))
    save("PORTAL_JOINT_GATE.json", {
        "scope": "light unmixed scalar width proxy and frozen one-loop leading logarithm only",
        "portal_convention": "V=(lambda_p/2)*(H_dagger_H)*sum_a(S_a^2)",
        "mH_GeV": mh, "vEW_GeV": v, "SM_width_GeV": gamma_sm,
        "light_scalar_mass_GeV": mass, "phase_space_factor": phase,
        "invisible_BR_benchmark": inv_max,
        "BR_source": "ATLAS arXiv:2301.10731v2, published 2023; not a latest global fit",
        "reference_scale_GeV": mu0, "cutoff_GeV": cutoff, "log_span": log_span,
        "allowed_extra_width_GeV": allowed_width, "rows": rows,
        "hessian_trace_max_relative_error": max_error,
        "full_RG_matching_result": None, "thermal_source_erasure_result": None,
        "heavy_threshold_stabilization_excluded": False})


def tcc():
    """One explicit TCC + constant-H + immediate-reheating history, not finite-state lore."""
    inp = json.loads((ROOT.parent / "grut_heavy_vector_budget_20261025/INPUTS.json").read_text())
    h = inp["normalized_parameters"]["H_GeV"]
    mp = 2.435e18
    h0 = 67.4 / 3.0856775814913673e19 * 6.582119569e-25
    t0 = 2.7255 * 8.617333262145e-14
    gr, gsr, gs0 = 106.75, 106.75, 43/11
    a = (90/(math.pi**2*gr))**.25
    tr = a*math.sqrt(mp*h)
    ae = t0/tr*(gs0/gsr)**(1/3)
    n_needed = math.log(ae*h/h0)
    n_max = math.log(mp/h)
    h_max = mp*(a*h0/t0*(gsr/gs0)**(1/3))**(2/3)
    # Verify equality independently by solving the difference in log H by bisection.
    def gap(logh):
        hv = math.exp(logh)
        end = t0/(a*math.sqrt(mp*hv))*(gs0/gsr)**(1/3)
        return math.log(end*hv/h0)-math.log(mp/hv)
    lo, hi = math.log(1e-10), math.log(1e12)
    for _ in range(100):
        mid = (lo+hi)/2
        if gap(mid) > 0:
            hi = mid
        else:
            lo = mid
    root = math.exp((lo+hi)/2)
    save("TCC_JOINT_GATE.json", {
        "assumptions": ["TCC adopted as an additional conjecture", "constant-H single-stage inflation",
                        "present Hubble scale inside initial Hubble radius", "instantaneous thermal reheating",
                        "standard subsequent entropy history"],
        "finite_Hilbert_space_implies_TCC": False,
        "M_reduced_GeV": mp, "H0_GeV": h0, "T0_GeV": t0,
        "gstar_reheat": gr, "gstar_entropy_reheat": gsr, "gstar_entropy_today": gs0,
        "heavy_point_H_GeV": h, "reheating_temperature_GeV": tr,
        "minimum_horizon_efolds": n_needed, "maximum_TCC_efolds": n_max,
        "efold_gap": n_needed-n_max, "H_upper_GeV": h_max,
        "heavy_H_over_upper": h/h_max, "bisection_formula_relative_error": abs(root/h_max-1),
        "all_inflation_histories_excluded": False})


def uv():
    ts = (0., 1., 5., 20.)
    rows = [{"B": b, "flow": [{"t": t, "g": 1., "y": 2+b*math.exp(-t)} for t in ts]}
            for b in (-.5, 0., .5)]
    save("UV_COUNTEREXAMPLE.json", {
        "scope": "two-coupling mathematical flow, not a beta function for SM plus gravity",
        "beta_g": "g-1", "beta_y": "-(y-2)", "fixed_point": [1., 2.],
        "general_solution": "g=1+A*exp(t), y=2+B*exp(-t)",
        "UV_completeness_requires": "A=0; B remains free",
        "flows": rows, "GRUT_fixed_point_or_matching_calculated": False})


if __name__ == "__main__":
    vacuum()
    portal()
    tcc()
    uv()
    print("Wrote four separate implication diagnostics; no mass likelihood or new cosmological point.")
