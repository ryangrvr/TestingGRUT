"""Item 1 (thermal function): tabulate the exact bosonic J_B by quadrature, check against
(a) the convergent high-T series derived in NOTES.md (denominator (l+2)!), (b) the Bessel sum,
(c) the truncated high-T expansion through the cubic and log terms (fit of the non-analytic coefficients).
Writes results/jb_checks.json and results/jb_table.npz (cache used by transition.py)."""
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy import special

from common import (JBSeries, JBTable, jb_quad, jb_nonanalytic, _jb_bessel, LN_AB, PI, JB0, to_builtin)

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(exist_ok=True)


def main():
    t0 = time.time()
    ser = JBSeries()
    tab = JBTable()                      # 5001 quadratures on [-25, 25]
    np.savez(OUT / "jb_table.npz", grid=tab.grid, jq=tab.jq)
    t_tab = time.time() - t0
    res = {"table": {"y2_min": tab.grid[0], "y2_max": tab.grid[-1], "step": 0.01, "n_points": len(tab.grid),
                     "seconds": t_tab,
                     "method": "scipy.integrate.quad of the defining integral; real part ln|1-exp(-i theta)| for x<|y| when y2<0"}}

    # (a) quadrature vs convergent series at every table node with |y2|<=10
    m = np.abs(tab.grid) <= 10
    d_ser = tab.jq[m] - ser(tab.grid[m])
    res["quad_vs_series_max_abs_on_nodes_|y2|<=10"] = float(np.max(np.abs(d_ser)))
    # quadrature vs Bessel sum for 10 < y2 <= 25
    mb = tab.grid > 10
    d_bes = tab.jq[mb] - _jb_bessel(tab.grid[mb])
    res["quad_vs_bessel_max_abs_on_nodes_10<y2<=25"] = float(np.max(np.abs(d_bes)))

    # spline interpolant vs series off-grid (midpoints) and derivative check
    mid = 0.5 * (tab.grid[1:] + tab.grid[:-1])
    mid = mid[np.abs(mid) <= 10]
    res["table_spline_vs_series_max_abs_midpoints"] = float(np.max(np.abs(tab(mid) - ser(mid))))
    res["table_spline_vs_series_max_abs_derivative_midpoints"] = float(np.max(np.abs(tab.d(mid) - ser.d(mid))))
    # independent direct quadrature at random off-grid points
    rng = np.random.default_rng(20261020)
    pts = np.concatenate([rng.uniform(-24.9, 24.9, 40), 10 ** rng.uniform(-6, 0, 20)])
    dq = np.array([jb_quad(x) for x in pts])
    res["random_offgrid_quad_vs_table_max_abs"] = float(np.max(np.abs(dq - tab(pts))))
    res["random_offgrid_quad_vs_series_or_bessel_max_abs"] = float(np.max(np.abs(dq - ser(pts))))

    # (c) truncated high-T expansion: error table and coefficient fit
    rows = []
    for y2 in [1e-4, 1e-3, 0.01, 0.04, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, -0.01, -0.1, -0.5, -1.0]:
        ex = jb_quad(y2)
        tr4 = float(jb_nonanalytic(y2))
        no_cubic = tr4 + PI / 6 * max(y2, 0) ** 1.5
        y6 = special.zeta(3) * y2**3 / (384 * PI**2)
        rows.append({"y2": y2, "J_exact_quad": ex, "highT_through_log": tr4, "err_highT_through_log": ex - tr4,
                     "err_if_cubic_dropped": ex - no_cubic, "next_term_zeta3_y6_over_384pi2": y6})
    res["highT_truncation_table"] = rows
    y = np.linspace(0.02, 0.35, 60)
    y2 = y**2
    rhs = np.array([jb_quad(t) for t in y2]) - (JB0 + PI**2 / 12 * y2)
    A = np.stack([y**3, y**4 * np.log(y2), y**4, y**6], axis=1)
    coef, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    res["fit_nonanalytic_coefficients"] = {
        "basis": "J - (-pi^4/45 + pi^2 y^2/12) = a y^3 + b y^4 ln y^2 + c y^4 + e y^6, y in [0.02,0.35]",
        "a_fit": coef[0], "a_expected_-pi/6": -PI / 6,
        "b_fit": coef[1], "b_expected_-1/32": -1 / 32,
        "c_fit": coef[2], "c_expected_ln(a_b)/32": LN_AB / 32,
        "e_fit": coef[3], "e_expected_zeta3/(384pi^2)": special.zeta(3) / (384 * PI**2),
        "ln_a_b": LN_AB,
    }
    # commonly quoted alternative with (l+1)! in the series denominator: its y^6 coefficient
    res["series_denominator_note"] = {
        "derived_y6_coefficient_(l+2)!": special.zeta(3) / (384 * PI**2),
        "alternative_(l+1)!_y6_coefficient": special.zeta(3) / (128 * PI**2),
        "quad_minus_alternative_at_y2=0.1": jb_quad(0.1) - (float(jb_nonanalytic(0.1)) + special.zeta(3) * 0.1**3 / (128 * PI**2)),
        "quad_minus_derived_at_y2=0.1": jb_quad(0.1) - float(ser(0.1)[0]),
    }
    res["checks_pass"] = {
        "series_matches_quadrature_1e-11": res["quad_vs_series_max_abs_on_nodes_|y2|<=10"] < 1e-11,
        "bessel_matches_quadrature_1e-11": res["quad_vs_bessel_max_abs_on_nodes_10<y2<=25"] < 1e-11,
        "spline_matches_series_1e-11": res["table_spline_vs_series_max_abs_midpoints"] < 1e-11,
        "cubic_coefficient_within_1e-4": abs(coef[0] + PI / 6) < 1e-4,
        "log_coefficient_within_1e-3": abs(coef[1] + 1 / 32) < 1e-3,
    }
    res["seconds_total"] = time.time() - t0
    (OUT / "jb_checks.json").write_text(json.dumps(to_builtin(res), indent=1))
    print("J_B checks:", json.dumps(to_builtin(res["checks_pass"])), "fit", coef, "time %.1fs" % res["seconds_total"])


if __name__ == "__main__":
    main()
