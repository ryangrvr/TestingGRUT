#!/usr/bin/env python3
"""marginal_compatibility.py -- GRUT closed-whole working theory, calculation 2.

Nested and overlapping restriction compatibility for a closed whole: nested restrictions must agree
(necessary and automatic), overlapping restrictions must agree on the overlap (necessary), and overlap
agreement is NOT sufficient for a common positive global realization (the quantum marginal problem).
Everything here is standard quantum theory, computed as semidefinite programs (cvxpy) with analytic checks.

Sections
  C1  nested restrictions compose and agree (identity check on a random 4-qubit state)
  C2  chain A-B-C with Bell-pair marginals on AB and BC: rho_B = I/2 from both sides, yet no global state
      (SDP infeasible; analytic proof: a pure marginal factorizes the global state)
  C3  quantitative obstruction: max F_AB + F_BC over all 3-qubit states (SDP); max symmetric singlet fidelity
  C4  Werner-state chain and triangle: largest p for which the marginals have a common global state
  C5  classical contrast: chain marginals always extend (Markov gluing); an odd anticorrelated cycle does not (LP)
  C6  inconsistent overlaps are immediately infeasible (necessary condition)
Writes ../results/marginal_compatibility.json and prints a summary.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import cvxpy as cp
import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
OUT.mkdir(exist_ok=True)
SOLVER = cp.CLARABEL

I2 = np.eye(2, dtype=complex)
phi_plus = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
PHI = np.outer(phi_plus, phi_plus.conj())
I4 = np.eye(4) / 4


def ptrace(rho, keep, n):
    keep = sorted(keep)
    r = rho.reshape([2] * (2 * n))
    cur = n
    for q in sorted([k for k in range(n) if k not in keep], reverse=True):
        r = np.trace(r, axis1=q, axis2=q + cur)
        cur -= 1
    d = 2 ** len(keep)
    return r.reshape(d, d)


def werner(p):
    return p * PHI + (1 - p) * I4


def random_density(n, seed):
    rng = np.random.default_rng(seed)
    G = rng.normal(size=(2 ** n, 2 ** n)) + 1j * rng.normal(size=(2 ** n, 2 ** n))
    rho = G @ G.conj().T
    return rho / np.trace(rho)


def tr_C(rho):  # 3-qubit expression -> AB
    return cp.partial_trace(rho, [2, 2, 2], axis=2)


def tr_A(rho):  # -> BC
    return cp.partial_trace(rho, [2, 2, 2], axis=0)


def tr_B(rho):  # -> AC
    return cp.partial_trace(rho, [2, 2, 2], axis=1)


def feasibility(constraints_fn):
    rho = cp.Variable((8, 8), hermitian=True)
    cons = [rho >> 0, cp.real(cp.trace(rho)) == 1] + constraints_fn(rho)
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver=SOLVER)
    return prob.status, rho.value


results = {}
t0 = time.time()

# ----------------------------------------------------------------------------- C1
rho4 = random_density(4, seed=7)
rAB_direct = ptrace(rho4, [0, 1], 4)
rAB_nested = ptrace(ptrace(rho4, [0, 1, 2], 4), [0, 1], 3)
rB_via_AB = ptrace(rAB_direct, [1], 2)
rB_via_BC = ptrace(ptrace(rho4, [1, 2], 4), [0], 2)
results["C1_nested_restriction"] = {
    "max_abs_Tr_C(Tr_D rho) - Tr_CD rho": float(np.max(np.abs(rAB_direct - rAB_nested))),
    "max_abs_rho_B_from_AB_minus_rho_B_from_BC": float(np.max(np.abs(rB_via_AB - rB_via_BC))),
    "note": "Restriction maps compose; nested and overlapping restrictions of ONE global state agree automatically. "
            "This is the necessary direction and it is a theorem (partial trace is a conditional expectation).",
}

# ----------------------------------------------------------------------------- C2
status, _ = feasibility(lambda r: [tr_C(r) == PHI, tr_A(r) == PHI])
rB_AB = ptrace(PHI, [1], 2)
rB_BC = ptrace(PHI, [0], 2)
results["C2_bell_chain"] = {
    "rho_AB": "Phi+ (pure Bell state)", "rho_BC": "Phi+ (pure Bell state)",
    "overlap_rho_B_from_AB": rB_AB.real.tolist(), "overlap_rho_B_from_BC": rB_BC.real.tolist(),
    "overlap_agreement_max_abs_dev": float(np.max(np.abs(rB_AB - rB_BC))),
    "sdp_status_common_global_state": status,
    "analytic_proof": "If rho_AB is pure then any extension satisfies rho_ABC = rho_AB (x) rho_C (purity => no correlations with C), "
                      "so rho_BC = Tr_A(rho_AB) (x) rho_C = (I/2) (x) rho_C, whose fidelity with Phi+ is 1/2 != 1. "
                      "Overlap agreement is necessary, not sufficient (entanglement monogamy).",
}

# ----------------------------------------------------------------------------- C3
rho = cp.Variable((8, 8), hermitian=True)
F_AB = cp.real(cp.trace(PHI @ tr_C(rho)))
F_BC = cp.real(cp.trace(PHI @ tr_A(rho)))
prob = cp.Problem(cp.Maximize(F_AB + F_BC), [rho >> 0, cp.real(cp.trace(rho)) == 1])
prob.solve(solver=SOLVER)
max_sum = float(prob.value)
rho_opt = rho.value
# symmetric: maximize t s.t. F_AB >= t, F_BC >= t
rho = cp.Variable((8, 8), hermitian=True)
tv = cp.Variable()
prob2 = cp.Problem(cp.Maximize(tv), [rho >> 0, cp.real(cp.trace(rho)) == 1,
                                     cp.real(cp.trace(PHI @ tr_C(rho))) >= tv, cp.real(cp.trace(PHI @ tr_A(rho))) >= tv])
prob2.solve(solver=SOLVER)
# analytic candidate: the W-type / symmetric three-qubit states; also the product of a singlet with a free qubit gives (1, 1/2)
results["C3_quantitative_obstruction"] = {
    "max_F_AB_plus_F_BC_over_all_3qubit_states": max_sum,
    "max_min(F_AB,F_BC)_symmetric": float(prob2.value),
    "F_AB_and_F_BC_at_optimum_of_sum": [float(np.real(np.trace(PHI @ ptrace(rho_opt, [0, 1], 3)))), float(np.real(np.trace(PHI @ ptrace(rho_opt, [1, 2], 3))))],
    "value_if_compatible_would_be": 2.0,
    "note": "Linear objective over the convex set of 3-qubit states (SDP). The optimum is strictly below 2: the two requested "
            "Bell marginals cannot coexist; the gap is the quantitative size of the obstruction.",
}

# ----------------------------------------------------------------------------- C4
def max_werner(constraint_sets):
    rho = cp.Variable((8, 8), hermitian=True)
    pv = cp.Variable()
    cons = [rho >> 0, cp.real(cp.trace(rho)) == 1, pv >= 0, pv <= 1]
    for red in constraint_sets:
        cons.append(red(rho) == pv * PHI + (1 - pv) * I4)
    prob = cp.Problem(cp.Maximize(pv), cons)
    prob.solve(solver=SOLVER)
    return float(prob.value), prob.status


p_chain, st_chain = max_werner([tr_C, tr_A])
p_tri, st_tri = max_werner([tr_C, tr_A, tr_B])
# separability threshold of a Werner state (PPT): p <= 1/3
results["C4_werner_thresholds"] = {
    "chain_AB_BC_max_p": p_chain, "chain_status": st_chain,
    "triangle_AB_BC_CA_max_p": p_tri, "triangle_status": st_tri,
    "werner_entangled_iff_p_gt": 1 / 3,
    "note": "W(p) = p Phi+ + (1-p) I/4 on each listed pair, with identical overlaps by construction. The largest p admitting a "
            "common global state is computed by SDP; above it the pairwise-consistent family has no global realization.",
}

# ----------------------------------------------------------------------------- C5
# classical chain: P(a,b) and P(b,c) with equal P(b) always extend via P(a,b,c) = P(a,b) P(c|b)
rng = np.random.default_rng(3)
Pab = rng.random((2, 2)); Pab /= Pab.sum()
Pb = Pab.sum(axis=0)
Pc_given_b = rng.random((2, 2)); Pc_given_b /= Pc_given_b.sum(axis=1, keepdims=True)  # rows: b
Pbc = Pb[:, None] * Pc_given_b
Pabc = Pab[:, :, None] * Pc_given_b[None, :, :]
results["C5_classical_contrast"] = {
    "chain_extension_max_abs_marginal_error": float(max(np.max(np.abs(Pabc.sum(axis=2) - Pab)), np.max(np.abs(Pabc.sum(axis=0) - Pbc)))),
    "chain_note": "Classically, two marginals sharing one variable ALWAYS have a joint (Markov gluing). The quantum Bell chain fails. "
                  "The obstruction in C2 is specifically quantum (monogamy), not a bookkeeping failure.",
}
# classical odd cycle: perfect anticorrelation on AB, BC, CA has no joint distribution
P = cp.Variable(8, nonneg=True)
idx = lambda a, b, c: 4 * a + 2 * b + c
cons = [cp.sum(P) == 1]
for a in range(2):
    for b in range(2):
        cons.append(cp.sum([P[idx(a, b, c)] for c in range(2)]) == (0.5 if a != b else 0.0))
        cons.append(cp.sum([P[idx(c, a, b)] for c in range(2)]) == (0.5 if a != b else 0.0))  # P(b=a, c=b) pattern -> (B,C)
        cons.append(cp.sum([P[idx(b, c, a)] for c in range(2)]) == (0.5 if a != b else 0.0))  # (C,A)
lp = cp.Problem(cp.Minimize(0), cons)
lp.solve(solver=cp.CLARABEL)
results["C5_classical_contrast"]["odd_cycle_anticorrelated_status"] = lp.status
results["C5_classical_contrast"]["cycle_note"] = ("Even classically, overlap agreement on a CYCLE is not sufficient (three pairwise perfectly "
                                                  "anticorrelated bits). Sufficiency of overlap agreement is a property of acyclic (tree) access "
                                                  "structures in the classical case and fails already for a chain in the quantum case.")

# ----------------------------------------------------------------------------- C6
rho_BC_bad = werner(0.0)  # rho_B = I/2
rho_AB_bad = np.kron(np.diag([0.8, 0.2]).astype(complex), np.diag([0.6, 0.4]).astype(complex))  # rho_B = diag(0.6,0.4)
status_bad, _ = feasibility(lambda r: [tr_C(r) == rho_AB_bad, tr_A(r) == rho_BC_bad])
results["C6_inconsistent_overlap"] = {
    "rho_B_from_AB": np.real(np.diag(ptrace(rho_AB_bad, [1], 2))).tolist(),
    "rho_B_from_BC": np.real(np.diag(ptrace(rho_BC_bad, [0], 2))).tolist(),
    "sdp_status": status_bad,
    "note": "Disagreeing overlaps are infeasible at once: overlap agreement is the necessary condition every sector's effective environment must satisfy.",
}

# ----------------------------------------------------------------------------- provenance
results["provenance"] = {
    "script_sha256": hashlib.sha256((HERE / "marginal_compatibility.py").read_bytes()).hexdigest(),
    "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__, "cvxpy": cp.__version__,
    "solver": "CLARABEL", "platform": platform.platform(), "wall_seconds": round(time.time() - t0, 2),
}


def rnd(o):
    if isinstance(o, float):
        return float(f"{o:.12g}")
    if isinstance(o, dict):
        return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list):
        return [rnd(v) for v in o]
    return o


(OUT / "marginal_compatibility.json").write_text(json.dumps(rnd(results), indent=1))
print("C1 nested: max dev %.1e ; overlap agreement %.1e" % (results["C1_nested_restriction"]["max_abs_Tr_C(Tr_D rho) - Tr_CD rho"], results["C1_nested_restriction"]["max_abs_rho_B_from_AB_minus_rho_B_from_BC"]))
print("C2 Bell chain: overlap dev %.1e ; SDP status: %s" % (results["C2_bell_chain"]["overlap_agreement_max_abs_dev"], status))
print("C3 max F_AB + F_BC = %.6f (would be 2 if compatible); max symmetric min(F_AB,F_BC) = %.6f ; at optimum (%.4f, %.4f)" % (max_sum, float(prob2.value), *results["C3_quantitative_obstruction"]["F_AB_and_F_BC_at_optimum_of_sum"]))
print("C4 Werner chain p* = %.6f (%s); triangle p* = %.6f (%s); entangled iff p > 1/3" % (p_chain, st_chain, p_tri, st_tri))
print("C5 classical chain extension error %.1e ; classical anticorrelated cycle LP: %s" % (results["C5_classical_contrast"]["chain_extension_max_abs_marginal_error"], lp.status))
print("C6 inconsistent overlap SDP: %s" % status_bad)
print("wall %.1fs sha256 %s" % (results["provenance"]["wall_seconds"], results["provenance"]["script_sha256"][:16]))
