"""Stationary classical profiles, tangent BVPs and reduced stability checks."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson, solve_bvp
from scipy.interpolate import CubicSpline, PPoly
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import spsolve

ARCHIVE = Path(__file__).resolve().parent.parent
REFERENCE = ARCHIVE / "grut_candidate_1_20261018/code/results_profile.npz"
CHARGE_FACTOR = 4 * np.pi / .1**2
Q_REF = 39943.84624223011


@dataclass
class Profile:
    label: str
    omega: float
    spline: PPoly
    statistics: dict


def fields(r, y, omega):
    f, fp, s, sp = y[:4]
    return fp, (s*s + .1*f*f - omega*omega)*f, sp, (s*s - 1 + f*f)*s


def reference_profile():
    with np.load(REFERENCE) as raw:
        return Profile("archived", float(raw["Omega"]), CubicSpline(raw["rho"], raw["y"], axis=1), {})


def stationary(seed, radius, tolerance, settings, *, charge=None, omega=None, label):
    r = np.linspace(0, radius, int(radius * settings["initial_nodes_per_unit_radius"]) + 1)
    scale = (charge / Q_REF)**(1/3) if charge is not None else 1.
    eval_r = r / scale
    guess = seed.spline(np.minimum(eval_r, seed.spline.x[-1]))[:4]
    guess[[1, 3]] /= scale
    guess[:, eval_r >= seed.spline.x[-1]] = np.array([[0.], [0.], [1.], [0.]])
    if charge is not None:
        guess = np.vstack([guess, cumulative_trapezoid(r*r*guess[0]**2, r, initial=0)])

        def ode(r, y, p):
            return np.vstack([*fields(r, y, p[0]), r*r*y[0]**2])

        def bc(a, b, p):
            return np.array([a[1], a[3], a[4], b[0], b[2]-1, CHARGE_FACTOR*p[0]*b[4]-charge])

        sol = solve_bvp(ode, bc, r, guess, p=[.8 if charge > 30000 else .87],
                        S=np.diag([0, -2, 0, -2, 0]), tol=tolerance, bc_tol=tolerance,
                        max_nodes=settings["maximum_nodes"])
        omega = float(sol.p[0])
        bc_error = np.max(abs(bc(sol.y[:, 0], sol.y[:, -1], sol.p)))
    else:
        def ode(r, y):
            return np.array(fields(r, y, omega))

        def bc(a, b):
            return np.array([a[1], a[3], b[0], b[2]-1])

        sol = solve_bvp(ode, bc, r, guess, S=np.diag([0, -2, 0, -2]),
                        tol=tolerance, bc_tol=tolerance, max_nodes=settings["maximum_nodes"])
        bc_error = np.max(abs(bc(sol.y[:, 0], sol.y[:, -1])))

    r = np.linspace(0, radius, settings["charge_quadrature_intervals"] + 1)
    f, fp, s, sp = sol.sol(r)[:4]
    potential = .25*(s*s-1)**2 + .5*s*s*f*f + .025*f**4
    grad = .5*(fp*fp + sp*sp)
    kinetic = .5*omega*omega*f*f
    energy, q, scaling = (
        CHARGE_FACTOR * simpson(r*r*v, x=r)
        for v in (potential + grad + kinetic, omega*f*f, grad + 3*(potential - kinetic))
    )
    active = f > 1e-7
    stats = {
        "target_charge": charge, "charge": q, "omega": omega, "energy": energy, "radius": radius,
        "tolerance": tolerance, "success": bool(sol.success), "message": sol.message,
        "nodes": len(sol.x), "max_rms": np.max(sol.rms_residuals), "max_bc": bc_error,
        "relative_virial": scaling/energy, "central_F": f[0], "central_S": s[0],
        "min_F": min(f), "min_S": min(s), "max_active_F_prime": max(fp[active], default=np.inf),
        "min_active_S_prime": min(sp[active], default=-np.inf),
    }
    spline = PPoly.construct_fast(np.ascontiguousarray(sol.sol.c[:, :, :4]), sol.x, axis=sol.sol.axis)
    return Profile(label, omega, spline, stats)


def tangent(profile, settings):
    """Differentiate the stationary equations before solving them."""
    omega = profile.omega
    r = profile.spline.x

    def ode(x, y):
        f, _, s, _ = profile.spline(x)[:4]
        u, up, v, vp = y
        mixing = 2*f*s
        return np.array([up, (s*s+.3*f*f-omega*omega)*u + mixing*v - 2*omega*f,
                         vp, (3*s*s-1+f*f)*v + mixing*u])

    def bc(a, b):
        return np.array([a[1], a[3], b[0], b[2]])

    sol = solve_bvp(ode, bc, r, np.zeros((4, len(r))), S=np.diag([0, -2, 0, -2]),
                    tol=2e-9, bc_tol=2e-9, max_nodes=settings["maximum_nodes"])
    r = np.linspace(0, r[-1], settings["charge_quadrature_intervals"] + 1)
    f = profile.spline(r)[0]
    d_q = CHARGE_FACTOR * simpson(r*r*(f*f+2*omega*f*sol.sol(r)[0]), x=r)
    return {"success": bool(sol.success), "dQ_dOmega": d_q, "max_rms": np.max(sol.rms_residuals),
            "max_bc": np.max(abs(bc(sol.y[:, 0], sol.y[:, -1])))}


def amplitude_inertia(a, b, mixing, step, test):
    """Block LDL by congruence. Count negative pivots of K_amp - test I."""
    off_square = 1/step**4
    ia = ib = ic = 0.
    negative = 0
    min_pivot = np.inf
    for x, y, z in zip(a, b, mixing):
        da, db, dc = x-test-off_square*ia, y-test-off_square*ib, z-off_square*ic
        determinant = da*db - dc*dc
        vals = np.linalg.eigvalsh([[da, dc], [dc, db]])
        negative += int(np.sum(vals < 0))
        min_pivot = min(min_pivot, np.min(abs(vals)))
        ia, ib, ic = db/determinant, da/determinant, -dc/determinant
    return {"negative": negative, "minimum_absolute_pivot": min_pivot}


def stability(profile, step, settings):
    gap = settings["stability_gap"]
    r = np.arange(1, round(40/step)) * step
    f, _, s, _ = profile.spline(r)[:4]
    n = len(r)
    fp = r * f
    fp /= np.linalg.norm(fp)
    mixing = 2*f*s
    result = {"h": step, "rho_box": 40., "sectors": {}}
    for ell in [0, 1, 2]:
        centrifugal = ell*(ell+1)/r**2
        lap = 2/step**2 + centrifugal
        a = lap + s*s+.3*f*f-profile.omega**2
        b = lap + 3*s*s-1+f*f
        phase = lap + s*s+.1*f*f-profile.omega**2
        off = np.full(n-1, -1/step**2)
        phase_eigen = eigh_tridiagonal(phase, off, select="i", select_range=(0, 1))[0]
        tests = [-gap, 0., gap]
        inertia = [amplitude_inertia(a, b, mixing, step, t) for t in tests]
        row = {"test_eigenvalues": tests, "amplitude_inertia": inertia, "phase_lowest": phase_eigen}
        if ell == 0:
            corrections = []
            for test in [0., gap]:
                matrix = bmat([[diags([off, a-test, off], [-1, 0, 1]), diags(mixing)],
                               [diags(mixing), diags([off, b-test, off], [-1, 0, 1])]],
                              format="csc")
                e = np.r_[fp, np.zeros(n)]
                sol = spsolve(matrix, e)
                corrections.append({
                    "test": test,
                    "secular": 1+4*profile.omega**2*np.dot(e, sol),
                    "residual": np.max(abs(matrix @ sol - e)),
                })
            row["rank_one_test"] = corrections
            row["passes"] = (
                all(x["negative"] == 1 for x in inertia[1:]) and
                all(x["secular"] < 0 and x["residual"] < 1e-8 for x in corrections) and
                abs(phase_eigen[0]) < settings["translation_phase_zero_tolerance"] and phase_eigen[1] > gap
            )
        else:
            if ell == 1:
                ztol = settings["translation_phase_zero_tolerance"]
                translation = [amplitude_inertia(a, b, mixing, step, t) for t in [-ztol, ztol]]
                row["translation_zero_window"] = translation
                row["passes"] = (
                    [x["negative"] for x in translation] == [0, 1] and
                    inertia[-1]["negative"] == 1 and phase_eigen[0] > gap
                )
            else:
                row["passes"] = inertia[-1]["negative"] == 0 and phase_eigen[0] > gap
        result["sectors"][str(ell)] = row
    result["passes"] = all(row["passes"] for row in result["sectors"].values())
    return result
