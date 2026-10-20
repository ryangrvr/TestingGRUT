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
  C4  isotropic-state chain and triangle: largest p for which the marginals have a common global state
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
rAB_nested_DC = ptrace(ptrace(rho4, [0, 1, 2], 4), [0, 1], 3)   # trace D, then C
rAB_nested_CD = ptrace(ptrace(rho4, [0, 1, 3], 4), [0, 1], 3)   # trace C, then D
rB_via_AB = ptrace(rAB_direct, [1], 2)
rB_via_BC = ptrace(ptrace(rho4, [1, 2], 4), [0], 2)
results["C1_nested_restriction"] = {
    "max_abs_Tr_C(Tr_D rho) - Tr_CD rho": float(np.max(np.abs(rAB_direct - rAB_nested_DC))),
    "max_abs_Tr_D(Tr_C rho) - Tr_CD rho": float(np.max(np.abs(rAB_direct - rAB_nested_CD))),
    "max_abs_rho_B_from_AB_minus_rho_B_from_BC": float(np.max(np.abs(rB_via_AB - rB_via_BC))),
    "note": "Restriction maps compose (Tr_CD = Tr_C o Tr_D = Tr_D o Tr_C) and nested or overlapping restrictions of ONE global state "
            "agree automatically. This is the necessary direction; it is an identity, checked here in both trace orders.",
}

# ----------------------------------------------------------------------------- C2
status, _ = feasibility(lambda r: [tr_C(r) == PHI, tr_A(r) == PHI])


def feasibility_scs(constraints_fn):
    rho = cp.Variable((8, 8), hermitian=True)
    prob = cp.Problem(cp.Minimize(0), [rho >> 0, cp.real(cp.trace(rho)) == 1] + constraints_fn(rho))
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=200000)
    return prob.status


def min_violation(targets):
    """Always-feasible epigraph SDP: smallest t with ||Tr_X rho - target_X||_F <= t for every listed reduction."""
    rho = cp.Variable((8, 8), hermitian=True)
    tv = cp.Variable()
    cons = [rho >> 0, cp.real(cp.trace(rho)) == 1]
    for red, target in targets:
        E = red(rho) - target
        cons.append(cp.norm(cp.hstack([cp.vec(cp.real(E), order="F"), cp.vec(cp.imag(E), order="F")])) <= tv)
    prob = cp.Problem(cp.Minimize(tv), cons)
    prob.solve(solver=SOLVER)
    return float(prob.value)


status_scs = feasibility_scs(lambda r: [tr_C(r) == PHI, tr_A(r) == PHI])
viol_c2 = min_violation([(tr_C, PHI), (tr_A, PHI)])
rB_AB = ptrace(PHI, [1], 2)
rB_BC = ptrace(PHI, [0], 2)
results["C2_bell_chain"] = {
    "sdp_status_second_solver_SCS": status_scs,
    "minimum_achievable_marginal_violation_frobenius": viol_c2,
    "rho_AB": "Phi+ (pure Bell state)", "rho_BC": "Phi+ (pure Bell state)",
    "overlap_rho_B_from_AB": rB_AB.real.tolist(), "overlap_rho_B_from_BC": rB_BC.real.tolist(),
    "overlap_agreement_max_abs_dev": float(np.max(np.abs(rB_AB - rB_BC))),
    "sdp_status_common_global_state": status,
    "analytic_proof": "If rho_AB is pure then any extension satisfies rho_ABC = rho_AB (x) rho_C (purity => no correlations with C), "
                      "so rho_BC = Tr_A(rho_AB) (x) rho_C = (I/2) (x) rho_C, whose fidelity with Phi+ is Tr(rho_C)/4 = 1/4 != 1. "
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
# optimum = top eigenspace of Phi+_AB (x) I + I (x) Phi+_BC (lambda_max = 3/2, F_AB = F_BC = 3/4), e.g. (2|000> + |110> + |011>)/sqrt(6);
# the W state gives (1/6, 1/6); Phi+_AB (x) rho_C gives (1, 1/4)
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
# certificates just above the thresholds: the minimum achievable marginal violation is strictly positive
viol_chain_above = min_violation([(tr_C, werner(2 / 3 + 1e-3)), (tr_A, werner(2 / 3 + 1e-3))])
viol_tri_above = min_violation([(tr_C, werner(1 / 3 + 1e-3)), (tr_A, werner(1 / 3 + 1e-3)), (tr_B, werner(1 / 3 + 1e-3))])
# separability threshold of the isotropic state (PPT): p <= 1/3
results["C4_werner_thresholds"] = {
    "state_family": "isotropic state W(p) = p Phi+ + (1-p) I/4 (U x U* invariant); locally equivalent to the Werner state p Psi- + (1-p) I/4 "
                    "pairwise, and for the chain globally via Y on B; the triangle problems are distinct but both give 1/3",
    "chain_AB_BC_max_p": p_chain, "chain_status": st_chain,
    "chain_min_violation_at_p_plus_1e-3": viol_chain_above,
    "triangle_AB_BC_CA_max_p": p_tri, "triangle_status": st_tri,
    "triangle_min_violation_at_p_plus_1e-3": viol_tri_above,
    "isotropic_entangled_iff_p_gt": 1 / 3,
    "note": "W(p) on each listed pair, with identical overlaps by construction. The largest p admitting a common global state is computed "
            "by SDP (interior-point value within 1e-8 of 2/3 and 1/3); above it the pairwise-consistent family has no global realization, "
            "certified by a strictly positive minimum marginal violation.",
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
# parity contrast: a 4-cycle with four anticorrelated edges DOES extend (1/2 on 0101 + 1/2 on 1010)
P4 = cp.Variable(16, nonneg=True)
idx4 = lambda a, b, c, d: 8 * a + 4 * b + 2 * c + d
cons4 = [cp.sum(P4) == 1]
for x in range(2):
    for y in range(2):
        val = 0.5 if x != y else 0.0
        cons4.append(cp.sum([P4[idx4(x, y, c, d)] for c in range(2) for d in range(2)]) == val)  # (A,B)
        cons4.append(cp.sum([P4[idx4(a, x, y, d)] for a in range(2) for d in range(2)]) == val)  # (B,C)
        cons4.append(cp.sum([P4[idx4(a, b, x, y)] for a in range(2) for b in range(2)]) == val)  # (C,D)
        cons4.append(cp.sum([P4[idx4(y, b, c, x)] for b in range(2) for c in range(2)]) == val)  # (D,A)
lp4 = cp.Problem(cp.Minimize(0), cons4)
lp4.solve(solver=cp.CLARABEL)
results["C5_classical_contrast"]["even_cycle_four_anticorrelated_status"] = lp4.status
results["C5_classical_contrast"]["cycle_note"] = ("Even classically, overlap agreement on a cycle need not suffice: three pairwise perfectly "
                                                  "anticorrelated bits have no joint distribution. The obstruction is an odd number of anticorrelated "
                                                  "edges around the cycle, not cyclicity as such: the 4-cycle with four anticorrelated edges extends "
                                                  "(1/2 on 0101 and 1/2 on 1010). Overlap agreement suffices for EVERY consistent family exactly on "
                                                  "acyclic structures (Vorob'ev 1962); in the quantum case it fails already on a chain.")

# ----------------------------------------------------------------------------- C6
rho_BC_bad = werner(0.0)  # rho_B = I/2
rho_AB_bad = np.kron(np.diag([0.8, 0.2]).astype(complex), np.diag([0.6, 0.4]).astype(complex))  # rho_B = diag(0.6,0.4)
status_bad, _ = feasibility(lambda r: [tr_C(r) == rho_AB_bad, tr_A(r) == rho_BC_bad])
status_bad_scs = feasibility_scs(lambda r: [tr_C(r) == rho_AB_bad, tr_A(r) == rho_BC_bad])
viol_bad = min_violation([(tr_C, rho_AB_bad), (tr_A, rho_BC_bad)])
results["C6_inconsistent_overlap"] = {
    "sdp_status_second_solver_SCS": status_bad_scs,
    "minimum_achievable_marginal_violation_frobenius": viol_bad,
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
print("C2 Bell chain: overlap dev %.1e ; SDP status: %s (SCS: %s); min marginal violation %.6f" % (results["C2_bell_chain"]["overlap_agreement_max_abs_dev"], status, status_scs, viol_c2))
print("C3 max F_AB + F_BC = %.6f (would be 2 if compatible); max symmetric min(F_AB,F_BC) = %.6f ; at optimum (%.4f, %.4f)" % (max_sum, float(prob2.value), *results["C3_quantitative_obstruction"]["F_AB_and_F_BC_at_optimum_of_sum"]))
print("C4 isotropic chain p* = %.6f (%s, violation at p*+1e-3: %.2e); triangle p* = %.6f (%s, violation %.2e); entangled iff p > 1/3" % (p_chain, st_chain, viol_chain_above, p_tri, st_tri, viol_tri_above))
print("C5 classical chain extension error %.1e ; odd anticorrelated 3-cycle LP: %s ; 4-cycle four anticorrelated: %s" % (results["C5_classical_contrast"]["chain_extension_max_abs_marginal_error"], lp.status, lp4.status))
print("C6 inconsistent overlap SDP: %s (SCS: %s); min violation %.6f" % (status_bad, status_bad_scs, viol_bad))
print("wall %.1fs sha256 %s" % (results["provenance"]["wall_seconds"], results["provenance"]["script_sha256"][:16]))
