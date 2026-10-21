"""Adversarial re-derivation of the estimators track: Lorentz-weighting bounds, area-weighted counterexample,
thin-wall stress / energy-loss identity Q, kernel calibrations.  Written independently of audit/estimators/code.

Methods differ from the audited track on purpose:
  * bounds checked by Monte-Carlo over random discrete energy-measure distributions (not LP), plus a brute-force
    two-point search for the capped sharp bound (not the chord formula);
  * Q identity derived symbolically (sympy) from the boosted thin-wall stress tensor, not quoted;
  * kernel fractions by quadrature.

Usage: python3 -I verify_bounds_Q.py <out_json>
"""
import json
import math
import sys

import numpy as np
import sympy as sp
from scipy import integrate


def R_from_energy_dist(x, wE):
    """R = E/(sigma A_area) for energy weights wE (sum 1) at squared speeds x: A_area = sum wE sqrt(1-x)."""
    return 1.0 / np.sum(wE * np.sqrt(1.0 - x), axis=-1)


def mc_bounds(n_samples=400000, k=6, seed=7):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.0, 1.0 - 1e-12, size=(n_samples, k)) ** rng.uniform(0.2, 5.0, size=(n_samples, 1))
    w = rng.dirichlet(np.ones(k) * 0.3, size=n_samples)
    s = np.sum(w * x, axis=1)
    R = R_from_energy_dist(x, w)
    lo = (1 - s) ** -0.5
    hi = 1 / (1 - s)
    return {
        "n": n_samples,
        "min_R_over_lower": float(np.min(R / lo)),
        "max_R_over_upper": float(np.max(R / hi)),
        "violations_lower": int(np.sum(R < lo * (1 - 1e-12))),
        "violations_upper": int(np.sum(R > hi * (1 + 1e-12))),
    }


def two_point_sup(s, xmax, n=4001):
    """Brute force over energy-measure two-point distributions {x1 <= s <= x2 <= xmax}."""
    best = (0.0, None)
    x1 = np.linspace(0.0, s, n)
    x2 = np.linspace(s, xmax, n)
    X1, X2 = np.meshgrid(x1, x2, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        p2 = np.where(X2 > X1, (s - X1) / (X2 - X1), 0.0)
    val = (1 - p2) * np.sqrt(1 - X1) + p2 * np.sqrt(1 - X2)
    i = np.unravel_index(np.argmin(val), val.shape)
    return float(1.0 / val[i]), float(math.sqrt(X1[i])), float(math.sqrt(X2[i])), float(p2[i])


def area_mix(s, vf):
    p = s / vf ** 2
    g = 1 / math.sqrt(1 - vf ** 2)
    R = (1 - p) + p * g
    sE = p * g * vf ** 2 / R
    return p, R, sE


def stress_and_Q():
    v, sig, H, t, s = sp.symbols("v sigma H t s", positive=True)
    g = 1 / sp.sqrt(1 - v ** 2)
    # rest-frame surface stress of a wall in the x-y plane: diag(sigma, -sigma, -sigma, 0) for (t, x, y, z)
    T_rest = sp.diag(sig, -sig, -sig, 0)
    # boost along the normal z with speed v (Lambda^mu_nu acting on contravariant indices; metric (+,-,-,-))
    L = sp.Matrix([[g, 0, 0, g * v], [0, 1, 0, 0], [0, 0, 1, 0], [g * v, 0, 0, g]])
    T_lab = sp.simplify(L * T_rest * L.T)
    # surface density per unit lab-frame area: the delta function delta(z') = delta(z - v t) / gamma
    T_area = sp.simplify(T_lab / g)
    rho = sp.simplify(T_area[0, 0])
    p_iso = sp.simplify((T_area[1, 1] + T_area[2, 2] + T_area[3, 3]) / 3)
    w = sp.simplify(p_iso / rho)
    # population: rho = sum gamma sigma dA, p = sum gamma sigma (v^2-2/3) dA -> w = <v^2>_E - 2/3
    # continuity d rho/dt = -3 H (1 + w) rho - Q with rho = A sigma / t
    A = sp.symbols("A", positive=True)
    rho_w = A * sig / t
    Q = sp.simplify(-sp.diff(rho_w, t) - 3 * H * (1 + (s - sp.Rational(2, 3))) * rho_w)
    coeff = sp.simplify(Q / (rho_w / t))
    mu = sp.symbols("mu", real=True)
    rho_mu = A * sig * t ** ((mu - 1) / 2)
    Qmu = sp.simplify((-sp.diff(rho_mu, t) - 3 * H * (1 + (s - sp.Rational(2, 3))) * rho_mu) / (rho_mu / t))
    return {
        "T00_per_area": str(rho), "Tzz_per_area": str(sp.simplify(T_area[3, 3])),
        "Txx_per_area": str(sp.simplify(T_area[1, 1])), "w_single_wall": str(w),
        "Q_coeff_general": str(sp.expand(coeff)),
        "Q_coeff_radiation": str(sp.expand(coeff.subs(H, 1 / (2 * t)))),
        "Q_coeff_radiation_offscaling_mu": str(sp.expand(Qmu.subs(H, 1 / (2 * t)))),
    }, (lambda sv, Ht=0.5: 1 - Ht * (1 + 3 * sv)), (lambda sv, m: (-m - 3 * sv) / 2)


def kernel_fractions():
    e = lambda S: 1 - S * S
    tot = integrate.quad(e, -1, 1)[0]
    sharp = integrate.quad(e, -0.5, 0.5)[0] / tot
    smooth = integrate.quad(lambda S: e(S) * e(S), -1, 1)[0] / tot
    norm = integrate.quad(lambda S: 0.75 * e(S), -1, 1)[0]
    # cross-check rest energy measure from the kink profile itself: S = tanh(z/delta); energy density 2 U = S'^2
    d = math.sqrt(2.0)
    edens = lambda z: (1 / d / math.cosh(z / d) ** 2) ** 2
    zc = d * math.atanh(0.5)
    sharp_z = integrate.quad(edens, -zc, zc)[0] / integrate.quad(edens, -40, 40)[0]
    return {"kernel_norm": norm, "sharp_fraction": sharp, "sharp_fraction_from_z_profile": sharp_z,
            "smooth_weighted_fraction": smooth, "11/16": 11 / 16, "4/5": 0.8}


def main(out):
    s = 0.456 ** 2
    res = {"s": s}
    res["bounds_closed_form"] = {"lower": (1 - s) ** -0.5, "upper": 1 / (1 - s)}
    res["monte_carlo_bounds"] = mc_bounds()
    caps = {}
    for vm in (0.8, 0.9, 0.95, 0.99, 0.999, 0.999999):
        R, v1, v2, p2 = two_point_sup(s, vm ** 2)
        chord = 1.0 / (1.0 - (s / vm ** 2) * (1.0 - math.sqrt(1.0 - vm ** 2)))
        caps[str(vm)] = {"brute_force_R": R, "support_v": [v1, v2], "fast_energy_fraction": p2, "chord_formula": chord}
    res["capped_upper"] = caps
    mix = {}
    for vf in (0.99, 0.999, 0.9999, 0.999999):
        p, R, sE = area_mix(s, vf)
        mix[str(vf)] = {"fast_area_fraction": p, "R": R, "s_E_of_mix": sE}
    res["area_weighted_mix"] = mix
    sym, Qc, Qmu = stress_and_Q()
    res["stress_Q_symbolic"] = sym
    res["Q_numbers"] = {
        "s_E=0.456^2": Qc(s), "s_E=0.456_literal": Qc(0.456), "sign_change_s": 1 / 3,
        "rms_0.450..0.476": [Qc(0.476 ** 2), Qc(0.450 ** 2)],
        "mu_-0.896": Qmu(s, -0.896), "mu_-0.975": Qmu(s, -0.975),
        # downward-shifted s_E if the dx = delta/2 estimator reads high by 10 % or 20 %
        "s_E/1.10": Qc(s / 1.10), "s_E/1.20": Qc(s / 1.20),
    }
    shifted = {}
    for lab, f in (("x1.00", 1.0), ("/1.05", 1.05), ("/1.10", 1.10), ("/1.20", 1.20)):
        ss = s / f
        shifted[lab] = {"s_E": ss, "R_lo": (1 - ss) ** -0.5, "R_hi": 1 / (1 - ss), "Q_rad": Qc(ss)}
    res["if_estimator_reads_high"] = shifted
    res["kernels"] = kernel_fractions()
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
