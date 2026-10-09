"""Issue #3 STANDARD RESULT / REPRODUCTION / COUNTEREXAMPLE controls.

No GRUT candidate, fitted parameter, empirical validation, or novelty claim.
Exact rational arithmetic except explicitly tagged transcendental evaluations.
Run: python development/foundations_lab/run.py --check
"""
from fractions import Fraction as F
from itertools import product
import math


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def add(a, b):
    return [[x+y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def tr(a):
    return sum(a[i][i] for i in range(len(a)))


def kron(a, b):
    return [[x*y for x in r for y in s] for r in a for s in b]


def purity(a):
    return tr(mm(a, a))


def partial_b(a):
    return [[sum(a[2*i+k][2*j+k] for k in range(2)) for j in range(2)]
            for i in range(2)]


def control_qubit():
    # The phase distinguishes states despite identical Z probabilities.
    plus = [[F(1, 2), F(1, 2)], [F(1, 2), F(1, 2)]]
    minus = [[F(1, 2), -F(1, 2)], [-F(1, 2), F(1, 2)]]
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "P_Z_plus": [plus[0][0], plus[1][1]],
            "P_Z_minus": [minus[0][0], minus[1][1]],
            "P_X_plus_plus": tr(mm(plus, plus)),
            "P_X_plus_minus": tr(mm(plus, minus)),
            "QFI_theta": 1, "FS_theta": F(1, 4)}


def control_units():
    # Unit vectors are (power of length, power of time). theta dimensionless.
    return {"classification": "COUNTEREXAMPLE",
            "early_cos_formula": F(1, 1)/F(3, 5),
            "v_formula_units_LT": [-1, 1], "velocity_units_LT": [1, -1],
            "Phi_formula_units_LT": [-2, 1], "potential_units_LT": [2, -2],
            "theta_second_derivative_units_LT": [0, -2],
            "acceleration_units_LT": [1, -2]}


def control_spectrum():
    # D=diag(1,2), H=[[0,D],[D,0]]; exact eigenvalues +/-1,+/-2.
    h = [[F(0), F(0), F(1), F(0)], [F(0), F(0), F(0), F(2)],
         [F(1), F(0), F(0), F(0)], [F(0), F(2), F(0), F(0)]]
    gamma = [[F((1 if i < 2 else -1) if i == j else 0)
              for j in range(4)] for i in range(4)]
    samples = []
    for omega in map(F, [0, 1, 2, 14]):
        # Paper uses exponent -z/2: its finite determinant is a square root.
        paper = (1+omega**2)*(1+omega**2/4)
        samples.append({"omega": omega, "paper_ratio": paper,
                        "ordinary_det_ratio": paper**2})
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "H": h, "Gamma": gamma, "H_squared": mm(h, h),
            "anticommutator": add(mm(h, gamma), mm(gamma, h)),
            "trace_H": tr(h), "eigenvalues": [-2, -1, 1, 2],
            "determinant_samples": samples,
            "scope": "Finite invertible self-adjoint H; no numerical zeta evaluation"}


def control_gibbs():
    # A self-adjoint direct sum +/- n on l2 is not bounded below.
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "beta": 1, "arithmetic": "FLOAT: math.cosh, not a convergence proof",
            "cutoffs": [{"N": n, "Z_N": 2*sum(math.cosh(k) for k in range(1, n+1)),
                         "lower_bound": math.exp(n)} for n in [1, 2, 4, 8]],
            "exact_obstruction": "Z_N >= exp(N) -> infinity"}


def control_encoding():
    u = [[F(1), F(0), F(0)], [F(0), F(1), F(0)]]
    gram = mm(transpose(u), u)
    err = add(gram, scale(-1, eye(3)))
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "U": u, "U_adjoint_U": gram,
            "identity_error_squared": sum(x*x for row in err for x in row),
            "domain_dimension": 3, "codomain_dimension": 2,
            "general_bound": "rank(U*U) <= dim(codomain)"}


def control_metrics():
    p, q = F(1, 4), F(3, 4)
    # Tangent A=[[0,1],[1,0]] at diag(p,q), tr(A)=0.
    # SLD c(p,q)=2/(p+q); BKM c(p,q)=(log p-log q)/(p-q).
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "rho_eigenvalues": [p, q], "SLD_metric": 4,
            "BKM_metric_exact": "4*ln(3)", "BKM_metric_float": 4*math.log(3),
            "arithmetic": "Exact metric formulas; FLOAT log evaluation",
            "scope": "Contractive mixed-state quantum metrics are not unique"}


def control_factorization():
    bell = [[F(1, 2) if i in (0, 3) and j in (0, 3) else F(0)
             for j in range(4)] for i in range(4)]
    cnot = [[F(i == [0, 1, 3, 2][j]) for j in range(4)] for i in range(4)]
    mapped = mm(mm(cnot, bell), transpose(cnot))
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "rho": bell, "unitary": cnot, "mapped_rho": mapped,
            "marginal_before": partial_b(bell), "marginal_after": partial_b(mapped),
            "purity_before": purity(partial_b(bell)),
            "purity_after": purity(partial_b(mapped)),
            "scope": "Alternative tensor factorization without fixed local algebra; CNOT is not local"}


def pr_probability(a, b, x, y):
    return F(1, 2) if (a ^ b) == x*y else F(0)


def control_pr():
    cells = [{"a": a, "b": b, "x": x, "y": y,
              "p": pr_probability(a, b, x, y)} for a, b, x, y in product(range(2), repeat=4)]
    corr = [sum((-1)**(a+b)*pr_probability(a, b, x, y)
                for a, b in product(range(2), repeat=2))
            for x, y in product(range(2), repeat=2)]
    chsh_local = [a0*b0+a0*b1+a1*b0-a1*b1
                  for a0, a1, b0, b1 in product((-1, 1), repeat=4)]
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "cells": cells, "correlators_00_01_10_11": corr,
            "PR_CHSH": corr[0]+corr[1]+corr[2]-corr[3],
            "local_max": max(chsh_local), "quantum_bound": "2*sqrt(2)"}


def control_tomography():
    i, x, z = eye(2), [[F(0), F(1)], [F(1), F(0)]], [[F(1), F(0)], [F(0), -F(1)]]
    # Y tensor Y is real, even though each Y is not a real observable.
    yy = [[F(0), F(0), F(0), -F(1)], [F(0), F(0), F(1), F(0)],
          [F(0), F(1), F(0), F(0)], [-F(1), F(0), F(0), F(0)]]
    plus, minus = scale(F(1, 4), add(eye(4), yy)), scale(F(1, 4), add(eye(4), scale(-1, yy)))
    local = [{"a": a, "b": b,
              "plus": tr(mm(plus, kron(aa, bb))), "minus": tr(mm(minus, kron(aa, bb)))}
             for a, aa in zip("IXZ", [i, x, z]) for b, bb in zip("IXZ", [i, x, z])]
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "rho_plus": plus, "rho_minus": minus, "YY": yy,
            "local_expectations": local, "global_plus": tr(mm(plus, yy)),
            "global_minus": tr(mm(minus, yy)), "real_dimensions": [3, 10, 9],
            "complex_dimensions": [4, 16, 16]}


def control_kinetics():
    # Reproduction of the existing two-state control, NOT a new no-go family.
    u, t = 1.0, 1.0
    pi = [math.exp(-u)/(2*math.cosh(u)), math.exp(u)/(2*math.cosh(u))]
    rows = []
    for g in (0.0, u*u):
        forward, reverse = .5*math.exp(g+u), .5*math.exp(g-u)
        gap = forward+reverse
        rows.append({"g": g, "rates": [forward, reverse], "gap": gap,
                     "P_01_at_t1": pi[1]*(1-math.exp(-gap*t)),
                     "detailed_balance_residual": pi[0]*forward-pi[1]*reverse})
    return {"classification": "REPRODUCTION / STANDARD COUNTEREXAMPLE",
            "u": u, "t": t, "stationary": pi, "models": rows,
            "unforced_rates_both": [F(1, 2), F(1, 2)],
            "arithmetic": "FLOAT evaluation; exact derivation in PROOFS.md",
            "scope": "Fixed labelled classical states, reversible Markov dynamics"}


def anomaly_charges(q):
    h, l = 3*q, -3*q
    u, d, e = q+h, q-h, l-h
    return {"q": q, "u": u, "d": d, "l": l, "e": e, "h": h,
            "SU3_squared_U1": 2*q-u-d, "SU2_squared_U1": 3*q+l,
            "gravity_squared_U1": 6*q-3*u-3*d+2*l-e,
            "U1_cubed": 6*q**3-3*u**3-3*d**3+2*l**3-e**3}


def control_anomalies():
    return {"classification": "STANDARD RESULT / REPRODUCTION",
            "families": [anomaly_charges(q) for q in [F(0), F(1, 6), F(1, 3)]],
            "supplied": ["SU(3)xSU(2)xU(1)", "one SM generation", "no RH neutrino",
                         "one Higgs doublet", "all three Yukawa couplings"],
            "freedom_remaining": "Overall hypercharge normalization and gauge coupling; representation choice"}


def control_prime_noise():
    units = [n for n in range(30) if math.gcd(n, 30) == 1]
    x = [[F(0), F(1)], [F(1), F(0)]]
    plus = [[F(1, 2), F(1, 2)], [F(1, 2), F(1, 2)]]
    zero = [[F(1), F(0)], [F(0), F(0)]]
    def channel(rho):
        return scale(F(1, 2), add(rho, mm(mm(x, rho), x)))
    return {"classification": "STANDARD RESULT / COUNTEREXAMPLE",
            "first_four_primes_product": 2*3*5*7, "units_mod_30": units,
            "nonunit_fraction": F(30-len(units), 30), "composite_unit": 49,
            "gcd_49_30": math.gcd(49, 30), "plus_after_X_dephasing": channel(plus),
            "zero_after_X_dephasing": channel(zero),
            "purity_plus_after": purity(channel(plus)),
            "purity_zero_after": purity(channel(zero))}


def control_velocity():
    # May 2026 definitions: vx=cos(omega*t), vy=sin(omega*t), c=1.
    omega, t = F(1, 2), 1.0
    return {"classification": "COUNTEREXAMPLE",
            "omega": omega, "t": t, "vector_norm_squared_exact": 1,
            "vector_norm_squared_float": math.cos(float(omega)*t)**2+math.sin(float(omega)*t)**2,
            "claimed_v_squared": omega**2,
            "scope": "Literal simultaneous definitions; distinct clock and velocity variables would need a replacement derivation"}


CONTROLS = {"C01_qubit": control_qubit, "C02_units": control_units,
            "C03_spectral": control_spectrum, "C04_gibbs": control_gibbs,
            "C05_encoding": control_encoding, "C06_metrics": control_metrics,
            "C07_factorization": control_factorization, "C08_PR_box": control_pr,
            "C09_local_tomography": control_tomography, "C10_kinetics": control_kinetics,
            "C11_anomalies": control_anomalies, "C12_prime_noise": control_prime_noise,
            "C13_velocity": control_velocity}


def run_controls():
    return {name: fn() for name, fn in CONTROLS.items()}
