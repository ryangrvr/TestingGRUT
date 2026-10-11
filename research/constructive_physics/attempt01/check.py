#!/usr/bin/env python3
"""Exact CP-01 algebra and independent tree-amplitude contractions.

Python 3 standard library only. Writes results.json next to this file.
No measured outcomes, fitted coefficients or private data are inputs.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json


def clean(p):
    return {k: v for k, v in p.items() if v}


def add(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out.get(k, Q(0)) + v
    return clean(out)


def scale(p, c):
    return clean({k: v * c for k, v in p.items()})


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i + j] = out.get(i + j, Q(0)) + a * b
    return clean(out)


def deriv(p):
    return clean({k - 1: k * v for k, v in p.items() if k})


def eval_poly(p, x):
    return sum((v * x ** k for k, v in p.items()), Q(0))


def exact_json(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {str(k): exact_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [exact_json(v) for v in x]
    return x


checks = []


def verify(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


# Polynomials in F; U and R are divided by m^2.
F = {1: Q(1)}
U = {3: Q(3, 2), 2: Q(-3), 1: Q(3, 2)}
R = deriv(U)
f = add(mul(F, R), scale(U, -1))
verify("Jordan curvature constraint",
       R == {2: Q(9, 2), 1: Q(-6), 0: Q(3, 2)})
verify("Legendre polynomial and df/dR=F",
       f == {3: Q(3), 2: Q(-3)} and deriv(f) == mul(F, deriv(R)))
trace = add(mul(F, R), scale(U, -2))
verify("Linearized Jordan trace gives scalar mass m",
       eval_poly(trace, Q(1)) == 0 and eval_poly(deriv(trace), Q(1)) == 3)
verify("Legendre fold F=2/3",
       eval_poly(deriv(R), Q(2, 3)) == 0 and
       eval_poly(R, Q(2, 3)) == Q(-1, 2))

# Einstein potential divided by m^2 M_P^2: Laurent polynomial in F.
V_laurent = {k - 2: v / 2 for k, v in U.items()}
verify("Conformal potential is reciprocal cosh",
       V_laurent == {1: Q(3, 4), 0: Q(-3, 2), -1: Q(3, 4)})

order = 8
cosh1 = {k: Q(1, factorial(k)) for k in range(0, order + 1, 2)}
cosh2 = {k: Q(2 ** k, factorial(k)) for k in range(0, order + 1, 2)}
base = add(cosh1, {0: Q(-1)})
hostile = add(cosh2, scale(cosh1, -4), {0: Q(3)})
a2 = Q(2, 3)
verify("Canonical mass and quartic relation",
       base[2] == Q(1, 2) and factorial(4) * base[4] * a2 == Q(2, 3))
verify("Hostile term leaves vacuum, tadpole and mass unchanged",
       all(hostile.get(k, Q(0)) == 0 for k in range(4)) and hostile[4] == Q(1, 2))
# Exact all-x identity, checked as Laurent polynomials in y=exp(x).
cy = {1: Q(1, 2), -1: Q(1, 2)}
c2y = {2: Q(1, 2), -2: Q(1, 2)}
verify("Hostile term is nonnegative: 2(cosh x-1)^2",
       add(c2y, scale(cy, -4), {0: Q(3)}) ==
       scale(mul(add(cy, {0: Q(-1)}), add(cy, {0: Q(-1)})), 2))

comparators = []
for epsilon in [Q(0), Q(1, 12), Q(1)]:
    coeff = add(base, scale(hostile, epsilon))
    quartic = factorial(4) * coeff[4] * a2
    comparators.append({"epsilon": epsilon, "lambda4_M_P2_over_m2": quartic})
verify("Conventional comparator violates the candidate relation",
       comparators[1]["lambda4_M_P2_over_m2"] == Q(4, 3))

# Scalar one-loop potential pole is proportional to cosh(x)^2.
loop_y = mul(cy, cy)
verify("Loop creates second reciprocal harmonic",
       loop_y == {2: Q(1, 4), 0: Q(1, 2), -2: Q(1, 4)})
loop_series = scale(add(cosh2, {0: Q(1)}), Q(1, 2))
# Allow the original mass, width/normalization and vacuum counterterms.
# At orders x^2 and x^4 the only fitting coefficients are A=-4, B=3.
xsinh = {k: Q(1, factorial(k - 1)) for k in range(2, order + 1, 2)}
original_tangent = add(scale(base, -4), scale(xsinh, 3), {0: Q(1)})
verify("Finite-family counterterms fit x^0,x^2,x^4 but fail at x^6",
       all(loop_series[k] == original_tangent[k] for k in [0, 2, 4]) and
       loop_series[6] - original_tangent[6] == Q(1, 40))

# Independently contract the graviton numerator, using (+---) just for
# momentum bookkeeping. Physical Mandelstam s,t,u have the same definitions
# as -q^2 in the LAW.md (-+++) convention.
eta = [1, -1, -1, -1]


def dot(p, q):
    return sum((Q(eta[i]) * p[i] * q[i] for i in range(4)), Q(0))


def sub(p, q):
    return [a - b for a, b in zip(p, q)]


def plus(p, q):
    return [a + b for a, b in zip(p, q)]


def tau(p, q, m2):
    return [[p[i] * q[j] + q[i] * p[j] -
             (eta[i] * (dot(p, q) - m2) if i == j else 0)
             for j in range(4)] for i in range(4)]


def contract(T, S):
    traceT = sum((eta[i] * T[i][i] for i in range(4)), Q(0))
    traceS = sum((eta[i] * S[i][i] for i in range(4)), Q(0))
    contraction = sum((eta[i] * eta[j] * T[i][j] * S[i][j]
                       for i in range(4) for j in range(4)), Q(0))
    return contraction - traceT * traceS / 2


cells = []
for E, p in [(Q(5, 4), Q(3, 4)), (Q(13, 12), Q(5, 12))]:
    for cost, sint in [(Q(0), Q(1)), (Q(3, 5), Q(4, 5))]:
        p1, p2 = [E, 0, 0, p], [E, 0, 0, -p]
        p3, p4 = [E, p * sint, 0, p * cost], [E, -p * sint, 0, -p * cost]
        s, t, u = dot(plus(p1, p2), plus(p1, p2)), dot(sub(p1, p3), sub(p1, p3)), dot(sub(p1, p4), sub(p1, p4))
        T, S = tau(p1, p3, Q(1)), tau(p2, p4, Q(1))
        N = contract(T, S)
        expected = -(s - 2) * (u - 2) - 2
        # Stress conservation: q_mu T^{mu nu}=0.
        q = sub(p1, p3)
        ward = [sum((eta[i] * q[i] * T[i][j] for i in range(4)), Q(0))
                for j in range(4)]
        cells.append({"s_over_m2": s, "t_over_m2": t, "u_over_m2": u,
                      "numerator_over_m4": N, "ward_residual": ward})
        assert s + t + u == 4 and N == expected and all(v == 0 for v in ward)
verify("Independent conserved-stress contractions match exchange numerator",
       len(cells) == 4)


def gravity_amplitude(s, t, u, m2=Q(1)):
    # M_P^2 times the amplitude, in h=g_eta+2h/M_P normalization.
    return sum(((x - 2 * m2) * (y - 2 * m2) + 2 * m2 ** 2) / z
               for x, y, z in [(s, u, t), (s, t, u), (t, u, s)])


spin2 = gravity_amplitude(Q(8), Q(-2), Q(-2))
candidate = spin2 - a2
comparator = spin2 - 2 * a2
verify("Full finite-angle tree amplitude benchmark",
       spin2 == Q(97, 4) and candidate == Q(283, 12) and
       comparator == Q(275, 12))

source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
result = {
    "type": "CP01_AUTHOR_EXACT_AND_TREE_RESULTS",
    "freeze_commit": "9ca24c0d987971bd466d40dd7f58b2657edcd674",
    "evidence": "Exact rational identities and tree-level conventional Feynman controls; no experimental data or independent review",
    "status": "CP01_EXACT_ALGEBRA_AND_TREE_CONTROLS_PASS",
    "quantum_closure_status": "DECLARED_FINITE_FAMILY_NOT_CLOSED_UNDER_SCALAR_LOOP",
    "check_groups": len(checks),
    "checks": checks,
    "Jordan": {"U_over_m2": U, "R_over_m2": R, "f_over_m2": f,
               "fold_F": Q(2, 3), "fold_R_over_m2": Q(-1, 2)},
    "Einstein": {"V_over_m2_M_P2_as_F_Laurent": V_laurent,
                 "mass_squared_over_m2": "1",
                 "lambda3": "0", "lambda4_M_P2_over_m2": a2,
                 "beta_squared": Q(1, 6), "Yukawa_strength": Q(1, 3)},
    "hostile_comparators": comparators,
    "one_loop": {"cosh_squared_as_exp_x_Laurent": loop_y,
                 "new_operator_after_vacuum_and_mass_subtraction": "cosh(2x)-4cosh(x)+3",
                 "x6_residual_after_all_original_tangent_fit": Q(1, 40),
                 "pole_convention": "d=4-2epsilon; V1_pole=-(V'')^2/(64*pi^2*epsilon); counterterm has opposite sign",
                 "full_graviton_and_matter_loops_computed": False},
    "stress_contraction_cells": cells,
    "scattering": {"kinematics": "s=8m^2,t=u=-2m^2, p=m and theta=pi/2",
                   "amplitudes_in_units_m2_over_M_P2": {"spin2": spin2,
                                                       "CP01": candidate,
                                                       "epsilon_1_over_12": comparator},
                   "cross_section_candidate_over_comparator": (candidate / comparator) ** 2,
                   "cross_section_ratio_float": float((candidate / comparator) ** 2),
                   "scope": "Leading tree, same inclusive/readout conventions; no available scalar apparatus certified"},
    "script_sha256": source_hash
}
path = Path(__file__).with_name("results.json")
path.write_text(json.dumps(exact_json(result), indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": result["status"], "exact_check_groups": len(checks),
                  "quantum_closure_status": result["quantum_closure_status"],
                  "tree_cross_section_ratio": result["scattering"]["cross_section_ratio_float"],
                  "results": str(path)}))
