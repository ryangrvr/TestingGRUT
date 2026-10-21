"""Audit of the Lorentz-weighting statements in Codex a69bd86 DERIVATION.md section 2 and of the energy-loss identity
Q = [1 - Ht(1 + 3<v^2>)] rho_w / t, plus a synthetic test of what network_sim's velocity estimator weights.

Read-only use of Claude's earlier packet (network_sim wallsim.measure; cosmo_observable core.history) by import.

Usage:  python3 -I lorentz_bounds.py <claude_packet_dir> <out_json>
"""
import json
import math
import sys

import numpy as np
from scipy.optimize import linprog


def R_E_bounds(s):
    """Bounds on R = E/(sigma A_area) given energy-weighted s = <v^2>_E."""
    return 1.0 / math.sqrt(1.0 - s), 1.0 / (1.0 - s)


def R_E_upper_capped(s, vmax):
    """Sharp upper bound when every speed is <= vmax: two-point {0, vmax} in energy measure."""
    x = vmax ** 2
    if s > x:
        return None
    return 1.0 / (1.0 - (s / x) * (1.0 - math.sqrt(1.0 - x)))


def lp_extremes_energy(s, grid):
    """Min and max of <sqrt(1-x)>_E over energy-measure distributions on grid with <x>_E = s (LP)."""
    f = np.sqrt(1.0 - grid)
    A_eq = np.vstack([np.ones_like(grid), grid])
    b_eq = np.array([1.0, s])
    lo = linprog(f, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    hi = linprog(-f, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    def support(res):
        idx = np.nonzero(res.x > 1e-9)[0]
        return [[float(math.sqrt(grid[i])), float(res.x[i])] for i in idx]
    # max R = 1/min<sqrt(1-x)>, min R = 1/max<sqrt(1-x)>
    return {"R_max": float(1.0 / lo.fun), "R_max_support_[v,energy_weight]": support(lo),
            "R_min": float(1.0 / (-hi.fun)), "R_min_support_[v,energy_weight]": support(hi)}


def lp_extremes_area(s, grid):
    """Min and max of <gamma>_A over area-measure distributions on grid with <x>_A = s."""
    g = 1.0 / np.sqrt(1.0 - grid)
    A_eq = np.vstack([np.ones_like(grid), grid])
    b_eq = np.array([1.0, s])
    lo = linprog(g, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    hi = linprog(-g, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    return {"R_min": float(lo.fun), "R_max_on_grid": float(-hi.fun), "vmax_grid": float(math.sqrt(grid[-1]))}


def Q_coeff(v2, Ht):
    """Q t / rho_w for a scaling network rho_w = A_E sigma / t with wall stress w = <v^2>_E - 2/3."""
    return 1.0 - Ht * (1.0 + 3.0 * v2)


# ---------------------------------------------------------------------------------------------------------------
# Synthetic test of network_sim's estimator: separated boosted phi^4 kinks with chosen speeds on a 1D lattice
# ---------------------------------------------------------------------------------------------------------------
DELTA = 2.0  # network_sim wall half-width (tanh(z/delta))


def kink_segment(v, dx, offset, half_len, sign):
    g = 1.0 / math.sqrt(1.0 - v * v)
    n = int(round(2 * half_len / dx))
    x = (np.arange(n) - n / 2 + offset) * dx
    u = g * x / DELTA
    phi = sign * np.tanh(u)
    pi = -sign * v * g / DELTA / np.cosh(u) ** 2      # d phi / d tau at the instant the kink is at x = 0
    return phi, pi


def estimator_on_segments(speeds, dx, offsets, half_len=40.0, use_wallsim=None):
    """Sum pi^2 / sum (central grad)^2 over |phi| < 1/2, exactly as wallsim.measure (core mask, central gradient)."""
    num = den = 0.0
    phis, pis = [], []
    for i, v in enumerate(speeds):
        sgn = 1.0 if i % 2 == 0 else -1.0
        phi, pi = kink_segment(v, dx, offsets[i], half_len, sgn)
        phis.append(phi)
        pis.append(pi)
    phi = np.concatenate(phis)
    pi = np.concatenate(pis)
    # periodic central gradient as in wallsim.measure
    c = (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * dx)
    core = np.abs(phi) < 0.5
    num = float(np.sum(pi[core] ** 2))
    den = float(np.sum(c[core] ** 2))
    out = {"own": num / den}
    if use_wallsim is not None:
        m = use_wallsim(phi.astype(np.float32), pi.astype(np.float32), dx, 10.0, 1.0, 0.0, clusters=False)
        out["wallsim_measure"] = m["v2_wall"]
    return out


def synthetic_tests(measure):
    rng = np.random.default_rng(20261021)
    cases = {
        "one_speed_0.456": [0.456] * 8,
        "mix_rest_and_0.9_equal_area": [0.0, 0.9] * 4,
        "mix_rest_and_0.99_equal_area": [0.0, 0.99] * 4,
        "uniform_speeds_0_to_0.8": list(np.linspace(0.0, 0.8, 8)),
    }
    out = {}
    for name, sp in cases.items():
        sp = np.array(sp)
        g = 1 / np.sqrt(1 - sp ** 2)
        s_A = float(np.mean(sp ** 2))
        s_E = float(np.sum(g * sp ** 2) / np.sum(g))
        R_true = float(np.mean(g))           # E/(sigma A_area) for equal-area walls
        res = {"s_area_weighted": s_A, "s_energy_weighted": s_E, "R_true_E_over_sigmaA": R_true,
               "R_bounds_from_s_E": R_E_bounds(s_E), "R_one_speed_from_s_A": 1 / math.sqrt(1 - s_A)}
        for dx in (1.0, 0.5, 0.25):
            vals, vals_ws = [], []
            for _ in range(16):   # average over sub-cell kink positions
                offs = rng.uniform(0, 1, size=len(sp))
                r = estimator_on_segments(sp, dx, offs, use_wallsim=measure)
                vals.append(r["own"])
                if "wallsim_measure" in r:
                    vals_ws.append(r["wallsim_measure"])
            res[f"estimator_dx{dx}"] = {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
                                        "rel_to_s_E": float(np.mean(vals) / s_E - 1) if s_E > 0 else None,
                                        "rel_to_s_A": float(np.mean(vals) / s_A - 1) if s_A > 0 else None,
                                        "wallsim_measure_mean": float(np.mean(vals_ws)) if vals_ws else None}
        out[name] = res
    return out


def main(packet, out_path):
    sys.path.insert(0, packet + "/modules/network_sim/code")
    try:
        from wallsim import measure
    except Exception as e:  # pragma: no cover
        print("wallsim import failed:", e)
        measure = None

    v_rms = 0.456
    s = v_rms ** 2
    lo, hi = R_E_bounds(s)
    res = {"inputs": {"v_rms_3D_PRS_256": v_rms, "seed_std": 0.006, "s_E": s}}
    res["energy_weighted"] = {
        "R_lower_(1-s)^-1/2": lo, "R_upper_(1-s)^-1": hi,
        "lower_attained_by": "single speed sqrt(s) (Jensen equality, sqrt concave)",
        "upper_approached_by": "energy fraction s at v -> 1 (zero area), energy fraction 1-s at rest; "
                               "sqrt(1-x) >= 1-x with equality only at x in {0,1}; not attained for v < 1",
        "upper_capped": {str(vm): R_E_upper_capped(s, vm) for vm in (0.8, 0.9, 0.95, 0.99, 0.999)},
        "A_increase_percent": [100 * (lo - 1), 100 * (hi - 1)],
        "Tpar_factor_R^1/2": [math.sqrt(lo), math.sqrt(hi)],
        "t_par_factor_R^-1": [1 / lo, 1 / hi],
    }
    grid = np.linspace(0.0, 0.999999, 20001) ** 1  # x = v^2 grid
    res["energy_weighted"]["LP_check"] = lp_extremes_energy(s, grid)
    # uncertainty and alternative 3D values
    alt = {}
    for label, v in {"0.456-0.006": 0.450, "0.456": 0.456, "0.456+0.006": 0.462, "384^3 (1 seed) 0.468": 0.468,
                     "verifier 0.457": 0.457, "verifier 0.476": 0.476}.items():
        a, b = R_E_bounds(v * v)
        alt[label] = {"s_E": v * v, "R_lo": a, "R_hi": b, "Q_coeff_radiation": Q_coeff(v * v, 0.5)}
    res["sensitivity"] = alt
    # area-weighted: Codex's counterexample and LP
    vf = 0.99
    p = s / vf ** 2
    R_mix = (1 - p) + p / math.sqrt(1 - vf ** 2)
    gf = 1 / math.sqrt(1 - vf ** 2)
    sE_of_mix = p * gf * vf ** 2 / ((1 - p) + p * gf)
    res["area_weighted"] = {
        "lower_(1-s_A)^-1/2": lo,
        "codex_mix_v0.99_area_fraction": p, "codex_mix_R": R_mix,
        "energy_weighted_s_of_that_mix": sE_of_mix,
        "R_of_mix_as_vfast_to_1": {str(v): (1 - s / v ** 2) + (s / v ** 2) / math.sqrt(1 - v ** 2)
                                   for v in (0.99, 0.999, 0.9999, 0.999999)},
        "LP_grid": lp_extremes_area(s, grid),
    }
    # Q identity
    res["Q_identity"] = {
        "formula": "Q = [1 - Ht(1 + 3<v^2>_E)] rho_w/t from d rho_w/dt = -H(1+3<v^2>_E) rho_w - Q, rho_w = A_E sigma/t",
        "radiation_Ht": 0.5,
        "Q_coeff_s=0.456^2": Q_coeff(s, 0.5),
        "Q_coeff_s=0.456_literal": Q_coeff(0.456, 0.5),
        "sign_change_at_<v^2>": 1.0 / 3.0,
        "Q_coeff_range_for_s_0.450^2..0.476^2": [Q_coeff(0.476 ** 2, 0.5), Q_coeff(0.450 ** 2, 0.5)],
        "free_network_n_=1+3s": 1 + 3 * s,
    }
    # T_par shift with cosmo_observable's own history integrator
    sys.path.insert(0, packet + "/modules/cosmo_observable/code")
    try:
        import core
        tp = {}
        for A in (0.7, 0.8, 0.9):
            row = {}
            for lab, fac in (("area", 1.0), ("E_lo", lo), ("E_hi", hi)):
                for br in (False, True):
                    h = core.history(core.SIGMA_BENCH, A=A * fac, backreaction=br)
                    row[f"{lab}_br{int(br)}_keV"] = h["T_switch_GeV"] * 1e6
            tp[str(A)] = row
        res["T_par_cosmo_observable"] = tp
    except Exception as e:  # pragma: no cover
        res["T_par_cosmo_observable"] = f"failed: {e}"
    res["synthetic_estimator_tests"] = synthetic_tests(measure)
    json.dump(res, open(out_path, "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
