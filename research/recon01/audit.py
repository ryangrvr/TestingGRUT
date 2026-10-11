"""Audit a published standard-model discriminator, not a new GRUT candidate.

Units in the finite model: every generator coefficient has units s^-1.
Published comparison: arXiv:2603.26075v2, B.2, with g=alpha_G dx^2/2,
a=Gamma_1, b=Gamma_2, c=beta dx^2/4, and hbar restored in SI estimates.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
I = np.eye(2, dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
ZS = (np.kron(Z, I), np.kron(I, Z))
ZZ = ZS[0] @ ZS[1]
SIGNS = np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]])
RHO0 = np.ones((4, 4), dtype=complex) / 4
# Entries of H tensor H are exactly +/-1/2, avoiding rounded square roots.
BX = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
               [1, 1, -1, -1], [1, -1, -1, 1]]) / 2
KERNEL_BASIS = BX[:, [2, 1, 3]]
G_SI = F("6.67430e-11")
HBAR_SI = F("1.054571817e-34")


def pt(x):
    return x.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def dissipator(x, matrix):
    out = np.zeros((4, 4), dtype=complex)
    for i in range(2):
        for j in range(2):
            q = ZS[j] @ ZS[i]
            out += matrix[i, j] * (ZS[i] @ x @ ZS[j] - (q @ x + x @ q)/2)
    return out


def generator(x, a, b, c, g):
    return -1j*g*(ZZ@x-x@ZZ) + dissipator(x, np.array([[a, c], [c, b]]))


def scalar_generator(a, b, c, g):
    result = np.zeros((4, 4), dtype=complex)
    matrix = np.array([[a, c], [c, b]])
    for i in range(4):
        for j in range(4):
            delta = SIGNS[i]-SIGNS[j]
            result[i, j] = (-1j*g*(np.prod(SIGNS[i])-np.prod(SIGNS[j]))
                            - .5*delta @ matrix @ delta)
    return result


def negativity(rho):
    ev = np.linalg.eigvalsh(pt(rho))
    return float(np.maximum(-ev, 0).sum())


def abs_lower(interval):
    lo, hi = interval
    return F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))


def interval_decision(a, b, c, g, squared_rate_error=F(0)):
    """Conservative exact box test. Inputs are rational interval endpoints.

    The caller must separately certify model applicability and the error bound.
    Positive error protects against an unmodelled contribution to the margin.
    """
    if (not all(isinstance(q, (tuple, list)) and len(q) == 2
                and all(isinstance(v, F) for v in q) for q in [a,b,c,g])
            or not isinstance(squared_rate_error, F)):
        return {"status": "INVALID_INPUT"}
    if (any(q[0] > q[1] for q in [a, b, c, g])
            or a[0] < 0 or b[0] < 0 or squared_rate_error < 0):
        return {"status": "INVALID_INPUT"}
    cp_possible_margin = a[1]*b[1] - abs_lower(c)**2
    if cp_possible_margin < 0:
        return {"status": "INCONSISTENT_CP_CALIBRATION"}
    margin = abs_lower(g)**2 + abs_lower(c)**2 - a[1]*b[1] - squared_rate_error
    return {"status": "EXCLUDES_NONENTANGLING_GKSL_CLASS" if margin > 0
            else "WITHIN_OR_UNRESOLVED", "certified_margin_s_minus_2": str(margin)}


def run():
    fixtures = [
        ("coherent_only", F(0), F(0), F(0), F(1)),
        ("below_uncorrelated", F(3,4), F(3,4), F(0), F(1)),
        ("boundary_uncorrelated", F(1), F(1), F(0), F(1)),
        ("above_uncorrelated", F(3,2), F(3,2), F(0), F(1)),
        ("correlation_active", F(1), F(1), F(3,4), F(1)),
        ("boundary_correlated", F(5,4), F(5,4), F(3,4), F(1)),
        ("above_correlated", F(3,2), F(3,2), F(3,4), F(1)),
        ("asymmetric_boundary", F(1,2), F(2), F(0), F(1)),
    ]
    times = [0., .001, .02, .1, .5, 1., 2.]
    maximum = {"matrix_unit_ODE_error": 0., "PT_generator_error": 0.,
               "kernel_projection_error": 0., "expm_solution_error": 0.,
               "slope_estimate_error": 0., "trace_error": 0.}
    min_schur_eigenvalue = 0.
    records = []
    for name, a, b, c, g in fixtures:
        assert a*b >= c*c  # ordinary generator must be CP
        det_cp, det_pt = a*b-c*c, a*b-c*c-g*g
        # Exact rational discriminant identity, an arithmetic control for the proof.
        assert (a+b)**2-((a-b)**2+4*(c*c+g*g)) == 4*det_pt
        af, bf, cf, gf = map(float, (a, b, c, g))
        f = scalar_generator(af, bf, cf, gf)
        Cpt = np.array([[af, -cf-1j*gf], [-cf+1j*gf, bf]])
        columns = []
        for k in range(16):
            E = np.zeros((4, 4), dtype=complex)
            E.flat[k] = 1
            out = generator(E, af, bf, cf, gf)
            maximum["matrix_unit_ODE_error"] = max(maximum["matrix_unit_ODE_error"],
                                                     float(np.max(np.abs(out-f*E))))
            maximum["PT_generator_error"] = max(maximum["PT_generator_error"],
                float(np.max(np.abs(pt(generator(pt(E), af,bf,cf,gf))
                                    - dissipator(E, Cpt)))))
            columns.append(out.flatten())
        L = np.array(columns).T
        projected = KERNEL_BASIS.conj().T @ pt(generator(RHO0,af,bf,cf,gf)) @ KERNEL_BASIS
        expected = np.zeros((3, 3), dtype=complex)
        expected[:2,:2] = Cpt
        maximum["kernel_projection_error"] = max(maximum["kernel_projection_error"],
            float(np.max(np.abs(projected-expected))))
        slope = max(0., .5*(math.sqrt((af-bf)**2+4*(cf*cf+gf*gf))-(af+bf)))
        eps = 1e-6
        slope_observed = negativity(np.exp(f*eps)*RHO0)/eps
        maximum["slope_estimate_error"] = max(maximum["slope_estimate_error"],
                                               abs(slope_observed-slope))
        rows = []
        for t in times:
            multiplier = np.exp(t*f)
            independent = expm(t*L)
            maximum["expm_solution_error"] = max(maximum["expm_solution_error"],
                float(np.max(np.abs(independent-np.diag(multiplier.flatten())))))
            ev = np.linalg.eigvalsh(multiplier)
            min_schur_eigenvalue = min(min_schur_eigenvalue, float(ev.min()))
            rho = multiplier*RHO0
            maximum["trace_error"] = max(maximum["trace_error"], abs(complex(np.trace(rho))-1))
            n = negativity(rho)
            if det_pt >= 0:
                assert n < 2e-12
            rows.append({"t_s":t,"negativity":n})
        # Independent preparations/readouts defining the four decay rates.
        assert -f[0,2].real == 2*af
        assert -f[0,1].real == 2*bf
        assert abs(-f[0,3].real-2*(af+bf+2*cf)) < 1e-15
        assert abs(-f[1,2].real-2*(af+bf-2*cf)) < 1e-15
        assert abs(((-f[0,3].real)-(-f[1,2].real))/8-cf) < 1e-15
        records.append({"name":name,"a_s_minus_1":str(a),"b_s_minus_1":str(b),
                        "c_s_minus_1":str(c),"g_s_minus_1":str(g),
                        "CP_determinant":str(det_cp),"PT_determinant":str(det_pt),
                        "negativity_initial_slope":slope,
                        "slope_finite_step":slope_observed,"trajectory":rows})
    assert maximum["matrix_unit_ODE_error"] == 0
    assert maximum["PT_generator_error"] == 0
    assert maximum["kernel_projection_error"] == 0
    assert maximum["expm_solution_error"] < 5e-13
    assert maximum["slope_estimate_error"] < 5e-6
    assert min_schur_eigenvalue > -2e-12

    # SI dimensional audit: k=2|g| is the unnormalized a+b threshold.
    dx, d = F("1e-4"), F("1e-3")
    unit_rows = []
    for name, m in [("10 femtograms",F("1e-17")),
                    ("10 picograms",F("1e-14")), ("1 nanogram",F("1e-12"))]:
        k = G_SI*m*m*dx*dx/(HBAR_SI*d**3)
        k_exact = G_SI*m*m*dx*dx/(HBAR_SI*d*(d*d-dx*dx))
        unit_rows.append({"mass_label":name,"mass_kg":str(m),
                          "a_plus_b_threshold_s_minus_1":float(k),
                          "normalized_decay_sum_threshold_s_minus_1":float(2*k),
                          "exact_point_geometry_threshold_s_minus_1":float(k_exact),
                          "leading_geometry_relative_error":str(k_exact/k-1),
                          "paper_printed_6_divided_by_threshold":float(F(6)/k)})
    assert unit_rows[0]["a_plus_b_threshold_s_minus_1"] < 1e-9
    assert unit_rows[2]["a_plus_b_threshold_s_minus_1"] > 6

    # Time-local is weaker than CP-divisible: ordinary 3-state static phase memory.
    # H_SE=Z_S tensor diag(-1,0,1), sigma_E=diag(1/4,1/2,1/4).
    memory = []
    for t in [0., math.pi/3, 2*math.pi/3, 5*math.pi/6, math.pi]:
        q = .5+.5*math.cos(2*t)
        unitary = expm(-1j*t*np.kron(Z,np.diag([-1.,0.,1.])))
        sigma_e = np.diag([.25,.5,.25])
        plus = np.ones((2,2))/2
        whole = unitary @ np.kron(plus,sigma_e) @ unitary.conj().T
        reduced = np.einsum('ikjk->ij',whole.reshape(2,3,2,3))
        expected = np.array([[1,q],[q,1]])/2
        assert np.max(np.abs(reduced-expected)) < 1e-14
        memory.append({"t":t,"q":q,"normalized_Choi_min_eigenvalue":min(0.,(1-abs(q))/2),
                       "time_local_dephasing_rate":math.tan(t)})
    ratio = F(3,4)/F(1,4)
    assert ratio == 3
    memory_witness = {"q_at_2pi_over_3":"1/4","q_at_5pi_over_6":"3/4",
        "intermediate_multiplier":str(ratio),"normalized_intermediate_Choi_eigenvalue":"-1",
        "full_maps_CPTP":True,"CP_divisible":False,
        "scope":"A scope counterexample to time-local => GKSL with nonnegative rates; not a classical-gravity counterexample."}

    synthetic_intervals = [
        ("exclusion", ((F('.7'),F('.8')),(F('.7'),F('.8')),(F('-.1'),F('.1')),
                       (F('.99'),F('1.01')),F('.01')), "EXCLUDES_NONENTANGLING_GKSL_CLASS"),
        ("uncertain", ((F('.9'),F('1.1')),(F('.9'),F('1.1')),(F('-.1'),F('.1')),
                       (F('.99'),F('1.01')),F('.01')), "WITHIN_OR_UNRESOLVED"),
        ("CP_inconsistent", ((F('.2'),F('.3')),(F('.2'),F('.3')),(F('.9'),F('1')),
                       (F('.99'),F('1.01')),F(0)), "INCONSISTENT_CP_CALIBRATION"),
        ("invalid", ((F('-1'),F('.3')),(F('.2'),F('.3')),(F('.0'),F('.1')),
                       (F('.99'),F('1.01')),F(0)), "INVALID_INPUT"),
    ]
    decisions = []
    for label, args, expected in synthetic_intervals:
        result = interval_decision(*args)
        assert result["status"] == expected
        decisions.append({"name":label,**result,"data_kind":"synthetic; no experimental data"})

    result = {
        "status":"PUBLISHED_DISCRIMINATOR_AUDIT_PASS_NO_NEW_GRUT_LAW",
        "scope":"Known GKSL model / NULL-PiT conventional comparator; no candidate law search",
        "new_GRUT_candidate_formulations":0,"new_GRUT_candidate_evaluations":0,
        "frozen_slots_consumed_total":1,"remaining_slots":2,
        "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
        "counts":{"published_model_fixtures":len(fixtures),"matrix_unit_generator_cells":16*len(fixtures),
                  "matrix_unit_time_cells":16*len(fixtures)*len(times),"finite_horizon_trajectories":len(fixtures)*len(times)},
        "maximum_errors":maximum,"minimum_Schur_eigenvalue":min_schur_eigenvalue,
        "fixtures":records,"SI_benchmark_audit":unit_rows,"memory_trace":memory,
        "memory_scope_witness":memory_witness,"interval_decisions":decisions,
        "uncertainty_note":"Numerical tolerances are diagnostics, not rigorous floating-point enclosures. Interval arithmetic is exact on supplied rational boxes; certification of physical boxes is external.",
        "Card1_ledger_sha256":hashlib.sha256((ROOT.parent/'card01'/'ATTEMPT_LEDGER.json').read_bytes()).hexdigest(),
    }
    (ROOT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    tracked = ['audit.py','REPORT.md','SOURCES.json','SCOPE.md','reproduce.sh','requirements.txt','results.json']
    hashes = {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in tracked}
    (ROOT/'MANIFEST.json').write_text(json.dumps({"sha256":hashes},indent=2)+'\n')
    print(result['status'])
    print(json.dumps({"maximum_errors":maximum,"fixture_count":len(fixtures),
                      "10_fg_threshold_s_minus_1":unit_rows[0]['a_plus_b_threshold_s_minus_1']},indent=2))


if __name__ == '__main__':
    run()
