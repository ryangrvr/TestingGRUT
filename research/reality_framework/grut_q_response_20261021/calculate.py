#!/usr/bin/env python3
"""Fixed-mu q/f(R) response calculations; see DERIVATION.md for scope.

This independently implements a local quadratic representative. It neither
imports nor replays Claude's solver. The theory and its known f(R) limits
are not new. Finite masses below are hierarchy controls, not allowed fits.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
SPEC = json.loads((ROOT / "INPUTS.json").read_text())
B = SPEC["background"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def densities(N, om, orad):
    return 3 * om * np.exp(-3 * N), 3 * orad * np.exp(-4 * N)


def gr_growth(grid, om, orad, ol, rtol, atol):
    def rhs(N, y):
        d, w = y
        rd, rr = densities(N, om, orad)
        h2 = (rd + rr) / 3 + ol
        hlog = -(rd + 4 * rr / 3) / (2 * h2)
        return [w, -(2 + hlog) * w + rd * d / (2 * h2)]

    a0 = math.exp(grid[0])
    sol = solve_ivp(rhs, (grid[0], grid[-1]), [a0, a0],
                    t_eval=grid, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


def background_and_growth(m, refined=False, ol=None):
    om, orad = B["omega_m_reference"], B["omega_r_reference"]
    ol = B["omega_lambda_reference"] if ol is None else ol
    lam = 3 * ol  # supplied residual, NOT selected by this calculation
    grid = np.linspace(math.log(B["a_initial"]), math.log(B["a_final"]), B["samples"])
    ks = B["k_over_Href"]
    rd0, rr0 = densities(grid[0], om, orad)
    x0 = (rd0 + 4 * lam) / (3 * m * m)
    # v0=-H0*rd0/m^2, with H0 fixed by the exact constraint.
    h0 = math.sqrt((rd0 + rr0 + lam + .75 * m*m*x0*x0)
                   / (3 * (1 + x0 - rd0 / (m*m))))
    v0 = -h0 * rd0 / (m*m)
    a0 = B["a_initial"]
    y0 = [x0, v0, h0] + [value for _ in ks for value in [a0, a0]]

    def rhs(N, y):
        x, v, h = y[:3]  # x=F-1, v=dF/dt, time in H_ref^-1
        f = 1 + x
        rd, rr = densities(N, om, orad)
        hdot = .5 * m*m*x - 2*h*h
        out = [v / h, (-3*h*v - m*m*x + (rd + 4*lam)/3) / h, hdot / h]
        for j, k in enumerate(ks):
            d, w = y[3+2*j:5+2*j]
            p2 = k*k * math.exp(-2*N)
            # Leading quasi-static f(R) expression: nu^2 = m^2 F.
            s = p2 / (p2 + m*m*f)
            gdyn = (1 + s/3) / f
            out += [w, -(2 + hdot/(h*h))*w + rd*gdyn*d/(2*h*h)]
        return out

    rtol = B["refined_rtol"] if refined else B["rtol"]
    atol = B["refined_atol"] if refined else B["atol"]
    sol = solve_ivp(rhs, (grid[0], grid[-1]), y0, t_eval=grid,
                    method="DOP853", rtol=rtol, atol=atol, max_step=min(.01, .8/m))
    if not sol.success:
        raise RuntimeError(sol.message)
    x, v, h = sol.y[:3]
    f = 1 + x
    rd, rr = densities(grid, om, orad)
    constraint = 3*f*h*h + 3*h*v - rd - rr - lam - .75*m*m*x*x
    denom = 3*f*h*h + np.abs(3*h*v) + rd + rr + abs(lam) + .75*m*m*x*x
    hgr = np.sqrt((rd+rr)/3 + ol)
    growth_gr = gr_growth(grid, om, orad, ol, rtol, atol)
    hdot = .5*m*m*x - 2*h*h
    rows = []
    for j, k in enumerate(ks):
        p2 = k*k * np.exp(-2*grid)
        s = p2 / (p2 + m*m*f)
        gdyn = (1+s/3)/f
        eta = (1-s/3)/(1+s/3)
        sigma = 1/f
        d = sol.y[3+2*j]
        # Time derivatives, scalar homogeneous modes and relativistic terms
        # have NOT been solved in the perturbation calculation.
        rows.append({
            "k_over_Href": k,
            "D_final": float(d[-1]),
            "D_GR_final": float(growth_gr[0,-1]),
            "D_ratio_to_GR_final": float(d[-1]/growth_gr[0,-1]),
            "Gdyn_over_Gref_final": float(gdyn[-1]),
            "eta_final": float(eta[-1]),
            "Sigma_final": float(sigma[-1]),
            "minimum_p_over_H": float(np.min(np.sqrt(p2)/h)),
            "maximum_H2_over_p2_plus_nu2": float(np.max(h*h/(p2+m*m*f))),
            "algebraic_Gdyn_eta_Sigma_max_abs_residual":
                float(np.max(np.abs(gdyn*(1+eta)/2-sigma)))
        })
    summary = {
        "m_over_Href": m, "omega_lambda_reference_supplied": ol,
        "H_final_over_Href": float(h[-1]),
        "F_initial": float(f[0]), "F_final": float(f[-1]),
        "max_relative_Friedmann_residual_sum_abs_denominator": float(np.max(np.abs(constraint)/denom)),
        "max_abs_H_over_H_GR_minus_1": float(np.max(np.abs(h/hgr-1))),
        "max_H_over_m": float(np.max(h/m)),
        "min_F": float(np.min(f)),
        "deceleration_final": float(-1-hdot[-1]/h[-1]**2),
        "propagation_dGW_over_dEM_from_initial_to_final": float(math.sqrt(f[-1]/f[0])),
        "growth": rows,
        "nfev": sol.nfev,
    }
    return summary, grid, sol.y


def exact_de_sitter_control():
    """A different exact solution: radiation and dust absent, F constant."""
    m, lam = 7., 2.1
    h0 = math.sqrt(lam/3)
    x0 = 4*lam/(3*m*m)

    def rhs(t, y):
        x, v, h = y
        return [v, -3*h*v-m*m*x+4*lam/3, .5*m*m*x-2*h*h]

    sol = solve_ivp(rhs, (0., 2.), [x0, 0., h0], method="DOP853", rtol=1e-11, atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(np.max(np.abs(sol.y - np.array([x0,0.,h0])[:,None])))


def eds_control():
    """Compare ODE integration to analytic EdS power laws, not to a replay."""
    rows = []
    for g in [1., 4/3]:
        p = (math.sqrt(1+24*g)-1)/4
        def rhs(N, y):
            return [y[1], -.5*y[1]+1.5*g*y[0]]
        sol = solve_ivp(rhs, (0., math.log(10)), [1., p], method="DOP853", rtol=1e-12, atol=1e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        rows.append({"Gdyn_over_G": g, "analytic_growth_exponent": p,
                     "relative_endpoint_error": float(abs(sol.y[0,-1]/10**p-1))})
    return rows


def tensor_control():
    """Direct wave integration against the independently derived WKB law."""
    cfg = SPEC["tensor_control"]
    start, stop = cfg["conformal_interval"]
    def F(eta):
        return 1+.02/eta**2
    def qlog(eta, modified):
        return 2/eta + (-.02/eta**3/F(eta) if modified else 0.)
    rows = []
    for k in cfg["k"]:
        amps=[]
        for modified in [False, True]:
            def rhs(t, y):
                return [y[1], -2*qlog(t,modified)*y[1]-k*k*y[0]]
            initial = np.array([1+0j, -qlog(start,modified)-1j*k])
            sol = solve_ivp(rhs, (start,stop), initial, method="DOP853", rtol=1e-11, atol=1e-13,
                            max_step=.4/k)
            if not sol.success:
                raise RuntimeError(sol.message)
            amps.append(abs(sol.y[0,-1]))
        expected=math.sqrt(F(stop)/F(start))
        observed=amps[0]/amps[1]
        rows.append({"k":k, "WKB_distance_ratio":expected, "direct_wave_distance_ratio":observed,
                     "relative_error":abs(observed/expected-1)})
    return rows


def scales():
    cfg=SPEC["scale_control"]
    hbar,c,mpc=cfg["hbar_eV_s"],cfg["c_m_s"],cfg["Mpc_m"]
    h0=hbar*cfg["H_reference_km_s_Mpc"]*1000/mpc
    m=cfg["assumed_uniform_unscreened_mass_floor_eV"]
    rows=[]
    for z in cfg["redshifts"]:
        for k in cfg["k_comoving_Mpc_inverse"]:
            p=hbar*c*k*(1+z)/mpc
            s=p*p/(p*p+m*m)
            # Record small departures directly, avoiding 1+1e-54 cancellation.
            rows.append({"z":z,"k_Mpc_inverse":k,"physical_momentum_eV":p,
                         "Gdyn_over_Glens_minus_1":s/3,
                         "one_minus_eta":2*s/(3+s)})
    return {"H_reference_eV":h0,"mass_floor_over_H_reference":m/h0,
            "Compton_length_m":hbar*c/m,"rows":rows}


def main():
    summaries=[]
    refinement=[]
    for m in B["m_over_Href"]:
        original, grid, y=background_and_growth(m)
        fine, _, yf=background_and_growth(m,True)
        # Absolute error for x=F-1; fractional for H and the growth factors.
        refinement.append({"m_over_Href":m,
            "max_abs_F_difference":float(np.max(np.abs(y[0]-yf[0]))),
            "max_relative_H_difference":float(np.max(np.abs(y[2]/yf[2]-1))),
            "max_relative_D_difference":float(max(np.max(np.abs(y[j]/yf[j]-1)) for j in range(3,len(y),2)))})
        summaries.append(fine)
    zero,_,_=background_and_growth(B["m_over_Href"][-1],True,ol=0.)
    tensor=tensor_control()
    eds=eds_control()
    ds=exact_de_sitter_control()
    checks={
        "background_constraint_under_1e-7":all(r["max_relative_Friedmann_residual_sum_abs_denominator"]<1e-7 for r in summaries+[zero]),
        "refinement_H_and_growth_under_2e-7":all(max(r["max_relative_H_difference"],r["max_relative_D_difference"])<2e-7 for r in refinement),
        "GR_background_limit_improves_with_mass":all(summaries[j+1]["max_abs_H_over_H_GR_minus_1"]<summaries[j]["max_abs_H_over_H_GR_minus_1"] for j in range(len(summaries)-1)),
        "positive_tensor_kinetic_coefficient":all(r["min_F"]>0 for r in summaries+[zero]),
        "exact_de_sitter_control_under_1e-9":ds<1e-9,
        "EdS_growth_control_under_1e-10":all(r["relative_endpoint_error"]<1e-10 for r in eds),
        "tensor_fine_WKB_error_under_1e-4":tensor[-1]["relative_error"]<1e-4,
        "zero_residual_local_branch_decelerates_at_endpoint":zero["deceleration_final"]>0,
    }
    checks={key:bool(value) for key,value in checks.items()}
    out={"classification":"Standard-theory construction and author calculation; not independent physical validation.",
         "source_commit":SPEC["source_commit"],
         "source_hashes":{"INPUTS.json":digest(ROOT/"INPUTS.json"),"calculate.py":digest(Path(__file__))},
         "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
         "calibrated_residual_backgrounds":summaries,"numerical_refinement":refinement,
         "zero_residual_control":zero,"exact_de_sitter_max_abs_error":ds,
         "EdS_growth_controls":eds,"tensor_propagation_control":tensor,
         "conditional_unscreened_scale_control":scales(),"numerical_checks":checks,
         "not_computed":["Full relativistic scalar perturbations and free scalar isocurvature modes",
             "Nonlinear screening and laboratory likelihood", "Vacuum selection, transitions and quantum corrections",
             "Primordial spectrum or CMB transfer functions", "Observed lensing or gravitational-wave likelihood"]}
    (ROOT/"RESULTS.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"numerical_checks":checks,"backgrounds":[{"m":s["m_over_Href"],
        "H_final":s["H_final_over_Href"],"F_final":s["F_final"],
        "growth_ratios":[g["D_ratio_to_GR_final"] for g in s["growth"]]} for s in summaries]},indent=2))
    if not all(checks.values()):
        raise SystemExit("A numerical control failed; retain RESULTS.json and diagnose.")


if __name__=="__main__":
    main()
