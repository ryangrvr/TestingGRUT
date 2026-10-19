#!/usr/bin/env python3
"""closed_whole_demo.py -- GRUT closed-whole working theory, calculation 1.

A finite closed whole Omega = A u Abar (u R) with a SUPPLIED Hamiltonian, preparation,
access structure (A = qubit 0) and instruments, carried from preparation through global
unitary evolution to local Born-rule measurements with exact energy and correlation
accounting. Exact diagonalization only (numpy/scipy), deterministic, no fitted inputs.

Sections
  S1  preparation -> evolution -> measurement; global invariants (what closure forces)
  S2  non-monotone local relaxation and information backflow (not forced: monotone relaxation)
  S3  reduced dynamics is DERIVED and state-dependent: CPTP map for product preparations,
      time-convolutionless generator, CP-divisibility failure, generator depends on bath state
  S4  initial correlations: identical marginals, different futures (reduced state alone does
      not fix the future; the inhomogeneous Nakajima-Zwanzig term is nonzero)
  S5  non-unique local law: same H_A, same rho_A(0), different complement -> different dynamics
  S6  records and interventions as included subsystems / unitaries; record persistence depends on H
  S7  time reversal: exact recovery and an entropy-decreasing legal preparation (no arrow forced)
  S8  closed-whole invariants and standard bounds that closure DOES force (checked numerically)
Writes ../results/closed_whole_demo.json and prints a summary.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
OUT.mkdir(exist_ok=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PLUS = np.array([1, 1], dtype=complex) / np.sqrt(2)
MINUS = np.array([1, -1], dtype=complex) / np.sqrt(2)
ZERO = np.array([1, 0], dtype=complex)
ONE = np.array([0, 1], dtype=complex)


# ----------------------------------------------------------------------------- utilities
def kron_all(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def site_op(o, k, n):
    return kron_all([o if j == k else I2 for j in range(n)])


def two_site(o1, k1, o2, k2, n):
    return kron_all([o1 if j == k1 else (o2 if j == k2 else I2) for j in range(n)])


def ptrace(rho, keep, n):
    """Partial trace of an n-qubit operator, keeping qubits in `keep` (ascending order)."""
    keep = sorted(keep)
    r = rho.reshape([2] * (2 * n))
    cur = n
    for q in sorted([k for k in range(n) if k not in keep], reverse=True):
        r = np.trace(r, axis1=q, axis2=q + cur)
        cur -= 1
    d = 2 ** len(keep)
    return r.reshape(d, d)


def herm(a):
    return (a + a.conj().T) / 2


def vn_entropy(rho):
    w = np.linalg.eigvalsh(herm(rho))
    w = w[w > 1e-14]
    return float(max(0.0, -np.sum(w * np.log(w))))


def trace_distance(a, b):
    w = np.linalg.eigvalsh(herm(a - b))
    return float(0.5 * np.sum(np.abs(w)))


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def proj(v):
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def thermal(H, beta):
    E, V = eigh(H)
    w = np.exp(-beta * (E - E.min()))
    w /= w.sum()
    return (V * w) @ V.conj().T


def ground(H):
    E, V = eigh(H)
    return proj(V[:, 0])


class Closed:
    """A closed whole: fixed Hamiltonian, exact unitary evolution by eigendecomposition."""

    def __init__(self, H):
        self.H = herm(H)
        self.E, self.V = eigh(self.H)

    def U(self, t):
        return (self.V * np.exp(-1j * self.E * t)) @ self.V.conj().T

    def evolve(self, rho0, t):
        U = self.U(t)
        return U @ rho0 @ U.conj().T


def chain_H(n, sites, J, Jz, hs):
    """XXZ chain on consecutive `sites` with longitudinal fields hs, embedded in n qubits."""
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for a, b in zip(sites[:-1], sites[1:]):
        H += J * (two_site(X, a, X, b, n) + two_site(Y, a, Y, b, n)) + Jz * two_site(Z, a, Z, b, n)
    for s, h in zip(sites, hs):
        H += h * site_op(Z, s, n)
    return H


def embed_bath_state(rho_bath, n_bath, n_total, bath_sites):
    """Return the n_total-qubit operator rho_bath on bath_sites (assumed the trailing qubits)."""
    assert bath_sites == list(range(n_total - n_bath, n_total))
    return rho_bath


PAULI = [I2, X, Y, Z]


def map_matrix(phi):
    """4x4 real matrix of a qubit map in the Pauli basis: phi(P_j) = sum_i M_ij P_i / ... (Hilbert-Schmidt)."""
    M = np.zeros((4, 4))
    for j, Pj in enumerate(PAULI):
        out = phi(Pj)
        for i, Pi in enumerate(PAULI):
            M[i, j] = np.real(np.trace(Pi.conj().T @ out)) / 2
    return M


def choi_from_matrix(M):
    """Choi matrix (4x4) of the qubit map with Pauli-basis matrix M."""
    C = np.zeros((4, 4), dtype=complex)
    for a in range(2):
        for b in range(2):
            Eab = np.zeros((2, 2), dtype=complex)
            Eab[a, b] = 1
            # expand Eab in Paulis, apply M, resum
            coeffs = np.array([np.trace(P.conj().T @ Eab) / 2 for P in PAULI])
            outc = M @ coeffs
            out = sum(c * P for c, P in zip(outc, PAULI))
            C += np.kron(out, Eab)
    return C


# ----------------------------------------------------------------------------- model datum
DATUM = {
    "N_total": 8,
    "A": [0],
    "Abar": [1, 2, 3, 4, 5, 6, 7],
    "H_A": "omega0/2 * Z_0, omega0 = 1.0",
    "H_Abar": "XXZ chain on qubits 1..7: J (X X + Y Y) + Jz Z Z nearest-neighbour, J=1.0, Jz=0.5, "
              "longitudinal fields h = [0.3,-0.7,0.5,-0.2,0.8,-0.4,0.6] (fixed, non-integrable)",
    "H_AAbar": "g X_0 X_1, g = 0.5",
    "preparation_S1": "rho_A(0) = |+><+| ; rho_Abar(0) = exp(-beta H_Abar)/Z, beta = 1.0 (mixed global state)",
    "instruments": "projective Z and X measurements on A; Born rule p = Tr(rho_A E)",
    "external_time": "parameter t of the global unitary; grid t in [0, 40] step 0.1",
    "units": "hbar = 1, energies in units of omega0",
}
omega0, J, Jz, g, beta = 1.0, 1.0, 0.5, 0.5, 1.0
hs = [0.3, -0.7, 0.5, -0.2, 0.8, -0.4, 0.6]
N = 8
A = [0]
ABAR = list(range(1, 8))
NB = 7

t_grid = np.round(np.arange(0.0, 40.0 + 1e-9, 0.1), 10)

H_A = omega0 / 2 * site_op(Z, 0, N)
H_B = chain_H(N, ABAR, J, Jz, hs)
H_AB = g * two_site(X, 0, X, 1, N)
H = H_A + H_B + H_AB
H_bath_only = chain_H(NB, list(range(NB)), J, Jz, hs)  # same chain on 7 qubits
rho_bath_th = thermal(H_bath_only, beta)
rho_bath_gs = ground(H_bath_only)
rho_bath_inf = np.eye(2 ** NB) / 2 ** NB

closed = Closed(H)
results = {"datum": DATUM}
t0 = time.time()

# ----------------------------------------------------------------------------- S1
rho_A0 = proj(PLUS)
rho0 = np.kron(rho_A0, rho_bath_th)
E_tot0 = np.real(np.trace(H @ rho0))
S1 = {"t": [], "E_tot": [], "E_A": [], "E_Abar": [], "E_int": [], "S_A": [], "S_Abar": [],
      "S_tot": [], "purity_tot": [], "MI_A_Abar": [], "pZ_plus": [], "pX_plus": [], "born_global_minus_local_max": 0.0}
PZp = proj(ZERO)
PXp = proj(PLUS)
rhoA_traj = []
rho_at = {}
for t in t_grid:
    rho = closed.evolve(rho0, t)
    rA = ptrace(rho, A, N)
    rB = ptrace(rho, ABAR, N)
    rhoA_traj.append(rA)
    if t in (0.0, 3.0, 6.0, 20.0, 40.0):
        rho_at[t] = rho
    S1["t"].append(float(t))
    S1["E_tot"].append(float(np.real(np.trace(H @ rho))))
    S1["E_A"].append(float(np.real(np.trace(H_A @ rho))))
    S1["E_Abar"].append(float(np.real(np.trace(H_B @ rho))))
    S1["E_int"].append(float(np.real(np.trace(H_AB @ rho))))
    sA, sB, sT = vn_entropy(rA), vn_entropy(rB), vn_entropy(rho)
    S1["S_A"].append(sA)
    S1["S_Abar"].append(sB)
    S1["S_tot"].append(sT)
    S1["purity_tot"].append(purity(rho))
    S1["MI_A_Abar"].append(sA + sB - sT)
    pz_local = float(np.real(np.trace(rA @ PZp)))
    px_local = float(np.real(np.trace(rA @ PXp)))
    pz_global = float(np.real(np.trace(rho @ site_op(PZp, 0, N))))
    px_global = float(np.real(np.trace(rho @ site_op(PXp, 0, N))))
    S1["pZ_plus"].append(pz_local)
    S1["pX_plus"].append(px_local)
    S1["born_global_minus_local_max"] = max(S1["born_global_minus_local_max"], abs(pz_local - pz_global), abs(px_local - px_global))
S1["invariants"] = {
    "max_abs_dE_tot": float(np.max(np.abs(np.array(S1["E_tot"]) - E_tot0))),
    "max_abs_dS_tot": float(np.max(np.abs(np.array(S1["S_tot"]) - S1["S_tot"][0]))),
    "max_abs_dpurity_tot": float(np.max(np.abs(np.array(S1["purity_tot"]) - S1["purity_tot"][0]))),
    "E_tot": E_tot0,
    "S_tot": S1["S_tot"][0],
    "range_E_A": [float(min(S1["E_A"])), float(max(S1["E_A"]))],
    "range_E_int": [float(min(S1["E_int"])), float(max(S1["E_int"]))],
    "range_S_A": [float(min(S1["S_A"])), float(max(S1["S_A"]))],
    "S_A_max_possible": float(np.log(2)),
}
# exact reduced equation check at t = 3: d rho_A/dt = -i Tr_Abar[H_A + H_AAbar, rho]  (the H_Abar term traces to zero)
tt = 3.0
dt = 1e-4
rA_p = ptrace(closed.evolve(rho0, tt + dt), A, N)
rA_m = ptrace(closed.evolve(rho0, tt - dt), A, N)
lhs = (rA_p - rA_m) / (2 * dt)
rho_t = closed.evolve(rho0, tt)
rhs_full = -1j * ptrace(H @ rho_t - rho_t @ H, A, N)
rhs_no_bath = -1j * ptrace((H_A + H_AB) @ rho_t - rho_t @ (H_A + H_AB), A, N)
bath_term = -1j * ptrace(H_B @ rho_t - rho_t @ H_B, A, N)
S1["reduced_equation_check_t3"] = {
    "max_abs_lhs_minus_rhs_full": float(np.max(np.abs(lhs - rhs_full))),
    "max_abs_lhs_minus_rhs_without_H_Abar": float(np.max(np.abs(lhs - rhs_no_bath))),
    "max_abs_Tr_Abar_commutator_H_Abar": float(np.max(np.abs(bath_term))),
    "norm_drhoA_dt": float(np.linalg.norm(lhs)),
}
results["S1_preparation_evolution_measurement"] = S1

# ----------------------------------------------------------------------------- S2
SA = np.array(S1["S_A"])
dSA = np.diff(SA)
sign_changes = int(np.sum(np.diff(np.sign(dSA[np.abs(dSA) > 1e-12])) != 0))
first_max_idx = int(np.argmax(SA))
S2 = {
    "S_A_t0": float(SA[0]),
    "S_A_max": float(SA.max()), "t_of_S_A_max": float(t_grid[first_max_idx]),
    "S_A_min_after_max": float(SA[first_max_idx:].min()),
    "t_of_min_after_max": float(t_grid[first_max_idx + int(np.argmin(SA[first_max_idx:]))]),
    "number_of_sign_changes_of_dS_A_dt_on_grid": sign_changes,
    "S_A_final_t40": float(SA[-1]),
}
# BLP information backflow between two preparations differing only on A (same bath state)
rho0_minus = np.kron(proj(MINUS), rho_bath_th)
D = []
for t in t_grid:
    rA1 = ptrace(closed.evolve(rho0, t), A, N)
    rA2 = ptrace(closed.evolve(rho0_minus, t), A, N)
    D.append(trace_distance(rA1, rA2))
D = np.array(D)
dD = np.diff(D)
S2["BLP"] = {
    "D_t0": float(D[0]), "D_min": float(D.min()), "t_of_D_min": float(t_grid[int(np.argmin(D))]),
    "D_t40": float(D[-1]),
    "backflow_sum_of_positive_increments": float(np.sum(dD[dD > 0])),
    "largest_single_increase": float(dD.max()),
    "note": "Any positive increment of the trace distance between two reduced states is information backflow "
            "(Breuer-Laine-Piilo); a Markovian (CP-divisible) reduced dynamics would make D non-increasing.",
}
results["S2_nonmonotone_relaxation"] = S2

# ----------------------------------------------------------------------------- S3
def dyn_map(rho_bath, t, closed_obj=closed):
    def phi(opA):
        return ptrace(closed_obj.evolve(np.kron(opA, rho_bath), t), A, N)
    return map_matrix(phi)


def tcl_generator(rho_bath, t, h=1e-3):
    Mp = dyn_map(rho_bath, t + h)
    Mm = dyn_map(rho_bath, t - h)
    M = dyn_map(rho_bath, t)
    Md = (Mp - Mm) / (2 * h)
    return Md @ np.linalg.inv(M), M


S3 = {"sample_times": [0.5, 1.0, 2.0, 5.0, 10.0, 20.0], "cptp_check": {}, "tcl": {}, "cp_divisibility": {}}
min_choi = 1.0
max_tp_dev = 0.0
for t in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0]:
    M = dyn_map(rho_bath_th, t)
    C = choi_from_matrix(M)
    min_choi = min(min_choi, float(np.linalg.eigvalsh(herm(C)).min()))
    # trace preservation: first row of M must be (1,0,0,0)
    max_tp_dev = max(max_tp_dev, float(np.max(np.abs(M[0] - np.array([1, 0, 0, 0])))))
S3["cptp_check"] = {"min_choi_eigenvalue_over_samples": min_choi, "max_trace_preservation_deviation": max_tp_dev,
                    "note": "For a product preparation with a fixed bath state, rho_A(0) -> rho_A(t) is a linear CPTP map (standard)."}
gens = {}
for t in S3["sample_times"]:
    K_th, M_th = tcl_generator(rho_bath_th, t)
    K_gs, _ = tcl_generator(rho_bath_gs, t)
    K_inf, _ = tcl_generator(rho_bath_inf, t)
    gens[t] = (K_th, K_gs, K_inf, M_th)
K_th_ref = gens[0.5][0]
S3["tcl"] = {
    "generator_time_dependence_fro_norm_K(t)-K(0.5)_thermal": {str(t): float(np.linalg.norm(gens[t][0] - K_th_ref)) for t in S3["sample_times"]},
    "generator_bath_state_dependence_fro_norm_K_thermal-K_ground": {str(t): float(np.linalg.norm(gens[t][0] - gens[t][1])) for t in S3["sample_times"]},
    "generator_bath_state_dependence_fro_norm_K_thermal-K_infiniteT": {str(t): float(np.linalg.norm(gens[t][0] - gens[t][2])) for t in S3["sample_times"]},
    "map_condition_number_thermal": {str(t): float(np.linalg.cond(gens[t][3])) for t in S3["sample_times"]},
    "note": "K(t) = dPhi_t/dt Phi_t^{-1} is the exact time-convolutionless generator of the derived reduced dynamics; "
            "it is time dependent (no fixed Lindbladian) and depends on the bath preparation (kernel is downstream of the state).",
}
# CP-divisibility: intermediate map V(t+tau,t) = Phi_{t+tau} Phi_t^{-1}; Choi min eigenvalue
tau = 0.5
cpd = {}
worst = 1.0
for t in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0]:
    M1 = dyn_map(rho_bath_th, t)
    M2 = dyn_map(rho_bath_th, t + tau)
    V = M2 @ np.linalg.inv(M1)
    lam = float(np.linalg.eigvalsh(herm(choi_from_matrix(V))).min())
    cpd[str(t)] = lam
    worst = min(worst, lam)
S3["cp_divisibility"] = {"tau": tau, "min_choi_eigenvalue_of_intermediate_map": cpd, "most_negative": worst,
                         "note": "A negative Choi eigenvalue of an intermediate map means the derived reduced dynamics is not CP-divisible (non-Markovian)."}
results["S3_derived_reduced_dynamics"] = S3

# ----------------------------------------------------------------------------- S4
# Two global preparations with IDENTICAL marginals rho_A(0) = I/2 and rho_Abar(0), differing by a marginal-free
# correlation term chi (Tr_A chi = 0 = Tr_Abar chi). Variant (i): infinite-temperature bath, chi = Y0 X1 or Z0 Z1.
# Variant (ii): thermal bath at beta = 1, chi = Y0 (x) sqrt(rho_b)(X1 - <X1>)sqrt(rho_b), which lives in the range of rho_b.
rho_A_mm = np.eye(2) / 2


def eps_max(rho_base, chi):
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if np.linalg.eigvalsh(herm(rho_base + mid * chi)).min() >= -1e-13:
            lo = mid
        else:
            hi = mid
    return lo


def sqrtm_psd(rho):
    w, V = eigh(herm(rho))
    w = np.clip(w, 0, None)
    return (V * np.sqrt(w)) @ V.conj().T


Eb, Vb = eigh(herm(rho_bath_th))
sq_b = sqrtm_psd(rho_bath_th)
X1_b = site_op(X, 0, NB)  # X on bath qubit 1 (first bath qubit) in the 7-qubit bath space
c1 = float(np.real(np.trace(rho_bath_th @ X1_b)))
O_dressed = sq_b @ (X1_b - c1 * np.eye(2 ** NB)) @ sq_b
S4_cases = [
    ("inf_T_bath_quantum_correlation_Y0X1", np.kron(rho_A_mm, rho_bath_inf), two_site(Y, 0, X, 1, N), rho_bath_inf),
    ("inf_T_bath_classical_correlation_Z0Z1", np.kron(rho_A_mm, rho_bath_inf), two_site(Z, 0, Z, 1, N), rho_bath_inf),
    ("thermal_beta1_bath_dressed_quantum_correlation", np.kron(rho_A_mm, rho_bath_th), np.kron(Y, O_dressed), rho_bath_th),
]
S4 = {}
for name, rho1, chi, rb in S4_cases:
    em = eps_max(rho1, chi)
    eps = 0.5 * em
    rho2 = rho1 + eps * chi
    marg_dev = max(np.max(np.abs(ptrace(rho2, A, N) - rho_A_mm)), np.max(np.abs(ptrace(rho2, ABAR, N) - rb)))
    dists = []
    for t in t_grid:
        dists.append(trace_distance(ptrace(closed.evolve(rho1, t), A, N), ptrace(closed.evolve(rho2, t), A, N)))
    dists = np.array(dists)
    inhom0 = -1j * ptrace(H @ (eps * chi) - (eps * chi) @ H, A, N)
    S4[name] = {
        "eps_max_for_positivity": float(em), "eps_used": float(eps),
        "min_eigenvalue_rho2": float(np.linalg.eigvalsh(herm(rho2)).min()),
        "trace_norm_of_correlation_term": float(np.sum(np.abs(np.linalg.eigvalsh(herm(eps * chi))))),
        "max_marginal_deviation_between_rho1_and_rho2": float(marg_dev),
        "trace_distance_rhoA_t0": float(dists[0]),
        "max_trace_distance_rhoA_over_grid": float(dists.max()), "t_of_max": float(t_grid[int(np.argmax(dists))]),
        "trace_distance_rhoA_at_t1_t5_t20_t40": [float(dists[np.argmin(np.abs(t_grid - s))]) for s in (1, 5, 20, 40)],
        "norm_of_inhomogeneous_term_at_t0": float(np.linalg.norm(inhom0)),
    }
S4["note"] = ("Two global preparations with identical rho_A(0) AND identical rho_Abar(0) evolve to different rho_A(t). "
              "No map on rho_A(0) alone (nor on both marginals) reproduces both futures; the Nakajima-Zwanzig inhomogeneous term "
              "carries the difference. Standard initial-correlation result (Pechukas 1994; Stelmachovic-Buzek 2001), computed in the supplied model.")
results["S4_initial_correlations"] = S4

# ----------------------------------------------------------------------------- S5
# Same H_A, same rho_A(0) = |+>, same coupling form g X_0 X_1; different complements.
# (b) single bath spin, omega1 = 1.0 ; (c) clean XX chain (free-fermion, integrable), 7 spins, no fields.
H2 = omega0 / 2 * site_op(Z, 0, 2) + omega0 / 2 * site_op(Z, 1, 2) + g * two_site(X, 0, X, 1, 2)
closed2 = Closed(H2)
rho0_2 = np.kron(rho_A0, thermal(omega0 / 2 * Z, beta))
H_xx = H_A + chain_H(N, ABAR, J, 0.0, [0.0] * 7) + H_AB
closed_xx = Closed(H_xx)
rho0_xx = np.kron(rho_A0, thermal(chain_H(NB, list(range(NB)), J, 0.0, [0.0] * 7), beta))
S5 = {"models": {
    "a_nonintegrable_XXZ_disordered": DATUM["H_Abar"],
    "b_single_spin": "H_Abar = omega0/2 Z_1 (one bath qubit), thermal at beta=1",
    "c_clean_XX_chain": "J (XX+YY) chain on 7 qubits, Jz=0, no fields (free fermions), thermal at beta=1",
}, "sample_times": [1.0, 2.0, 5.0, 10.0, 20.0, 40.0], "trace_distance_rhoA": {}, "S_A": {}}
SA_a, SA_b, SA_c = [], [], []
for t in t_grid:
    SA_b.append(vn_entropy(ptrace(closed2.evolve(rho0_2, t), [0], 2)))
    SA_c.append(vn_entropy(ptrace(closed_xx.evolve(rho0_xx, t), A, N)))
SA_b, SA_c = np.array(SA_b), np.array(SA_c)
for t in S5["sample_times"]:
    i = int(np.argmin(np.abs(t_grid - t)))
    rA_a = rhoA_traj[i]
    rA_b = ptrace(closed2.evolve(rho0_2, t), [0], 2)
    rA_c = ptrace(closed_xx.evolve(rho0_xx, t), A, N)
    S5["trace_distance_rhoA"][str(t)] = {"a_vs_b": trace_distance(rA_a, rA_b), "a_vs_c": trace_distance(rA_a, rA_c), "b_vs_c": trace_distance(rA_b, rA_c)}
    S5["S_A"][str(t)] = {"a": float(SA[i]), "b": float(SA_b[i]), "c": float(SA_c[i])}
# recurrence in model b: S_A returns close to its initial value
S5["model_b_recurrence"] = {
    "S_A_t0": float(SA_b[0]), "S_A_max": float(SA_b.max()),
    "min_S_A_after_first_max": float(SA_b[int(np.argmax(SA_b)):].min()),
    "t_of_that_min": float(t_grid[int(np.argmax(SA_b)) + int(np.argmin(SA_b[int(np.argmax(SA_b)):]))]),
    "number_of_times_S_A_within_1e-6_of_S_A_t0_for_t>1": int(np.sum(np.abs(SA_b[t_grid > 1] - SA_b[0]) < 1e-6)),
}
S5["model_c_recurrence"] = {"S_A_max": float(SA_c.max()), "min_S_A_after_first_max": float(SA_c[int(np.argmax(SA_c)):].min())}
# (b') commensurate single-spin complement: g' = 0.75 makes all Bohr frequencies multiples of 0.5 -> exact period 4 pi
g2 = 0.75
H2c = omega0 / 2 * site_op(Z, 0, 2) + omega0 / 2 * site_op(Z, 1, 2) + g2 * two_site(X, 0, X, 1, 2)
closed2c = Closed(H2c)
rho0_2c = np.kron(rho_A0, thermal(omega0 / 2 * Z, beta))
Tper = 4 * np.pi
S5["model_b_prime_exact_recurrence"] = {
    "g_prime": g2, "period": float(Tper),
    "Bohr_frequencies": sorted(set(round(abs(a - b), 6) for a in closed2c.E for b in closed2c.E)),
    "trace_distance_rho_global_T_vs_0": trace_distance(closed2c.evolve(rho0_2c, Tper), rho0_2c),
    "S_A_at_half_period": vn_entropy(ptrace(closed2c.evolve(rho0_2c, Tper / 2), [0], 2)),
    "S_A_at_period": vn_entropy(ptrace(closed2c.evolve(rho0_2c, Tper), [0], 2)),
    "note": "A closed complement with commensurate spectrum returns A exactly to its prepared state: relaxation is not forced.",
}
S5["note"] = "Closure fixes nothing about the local law of A: the same subsystem datum with different (closed) complements has different reduced dynamics."
results["S5_nonunique_local_law"] = S5

# ----------------------------------------------------------------------------- S6
N9 = 9
R = 8
H9 = np.kron(H, I2)  # register idle: [H9, Z_R] = 0 -> record stable by construction of H, not by closure
rho0_9 = np.kron(rho0, proj(ZERO))
closed9 = Closed(H9)
t1, t2 = 3.0, 6.0
CNOT = kron_all([I2] * N9)  # build CNOT(A -> R): |0><0|_A x I + |1><1|_A x X_R
P0A = site_op(proj(ZERO), 0, N9)
P1A = site_op(proj(ONE), 0, N9)
CNOT = P0A + P1A @ site_op(X, R, N9)
rho_t1 = closed9.evolve(rho0_9, t1)
pz_plus_t1 = float(np.real(np.trace(rho_t1 @ P0A)))
rho_after_write = CNOT @ rho_t1 @ CNOT.conj().T
rR = ptrace(rho_after_write, [R], N9)
# record populations equal Born probabilities at t1
record_pop_0 = float(np.real(rR[0, 0]))
rho_t2_int = closed9.evolve(rho_after_write, t2 - t1)
rho_t2_noint = closed9.evolve(rho0_9, t2)
P0R = site_op(proj(ZERO), R, N9)
P1R = site_op(proj(ONE), R, N9)
joint = {}
for za, Pa in (("+1", P0A), ("-1", P1A)):
    for zr, Pr in (("+1", P0R), ("-1", P1R)):
        joint[f"zA={za},zR={zr}"] = float(np.real(np.trace(rho_t2_int @ Pa @ Pr)))
pA_int = float(np.real(np.trace(rho_t2_int @ P0A)))
pA_noint = float(np.real(np.trace(rho_t2_noint @ P0A)))
# decohered (classical-record) comparison: measure-and-forget = dephase A at t1 without a register
rho_t1_8 = closed.evolve(rho0, t1)
P0A8, P1A8 = site_op(proj(ZERO), 0, N), site_op(proj(ONE), 0, N)
rho_deph = P0A8 @ rho_t1_8 @ P0A8 + P1A8 @ rho_t1_8 @ P1A8
pA_deph = float(np.real(np.trace(closed.evolve(rho_deph, t2 - t1) @ P0A8)))
# record persistence under a Hamiltonian that acts on R
kappa = 0.3
H9k = H9 + kappa * site_op(X, R, N9)
closed9k = Closed(H9k)
agree = {}
for tt2 in (3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
    r_k = closed9k.evolve(rho_after_write, tt2 - t1)
    r_0 = closed9.evolve(rho_after_write, tt2 - t1)
    # record stability: probability that Z_R still equals the value written (conditional on the written value)
    # computed as Tr[rho * (P0A(t1-record) ... )]: compare register marginal to the one right after the write
    agree[str(tt2)] = {"kappa=0_trace_distance_register_vs_written": trace_distance(ptrace(r_0, [R], N9), rR),
                       "kappa=0.3_trace_distance_register_vs_written": trace_distance(ptrace(r_k, [R], N9), rR)}
S6 = {
    "construction": "Omega' = A u Abar u R (9 qubits); H' = H x I_R; record written at t1=3 by the unitary CNOT(A->R); readout at t2=6.",
    "born_pZ_plus_at_t1": pz_plus_t1, "register_population_0_after_write": record_pop_0,
    "abs_difference": abs(pz_plus_t1 - record_pop_0),
    "joint_distribution_zA_t2_zR": joint,
    "pZ_plus_A_at_t2_with_record_written": pA_int,
    "pZ_plus_A_at_t2_without_intervention": pA_noint,
    "pZ_plus_A_at_t2_with_measure_and_forget_dephasing": pA_deph,
    "note_two_time": "The one-time map Phi_t does not give the statistics after an intervention at t1: the process after the record is "
                     "written differs from the uninterrupted process (back-action by entanglement with the included register). "
                     "A measure-and-forget dephasing at t1 gives the same later statistics as the unitary record (register traced out).",
    "record_persistence": agree,
    "note_persistence": "With [H', Z_R] = 0 the record is exactly stable; adding kappa X_R to the supplied Hamiltonian makes it drift. "
                        "Record persistence is a property of the supplied dynamics, not of closure.",
}
results["S6_records_and_interventions"] = S6

# ----------------------------------------------------------------------------- S7
T = 40.0
rho_T = closed.evolve(rho0, T)
rho_back = closed.evolve(rho_T, -T)
# H is real in the computational basis -> complex conjugation is a time reversal; the conjugated final state is a legal preparation
is_real = float(np.max(np.abs(np.imag(H))))
rho_rev0 = rho_T.conj()
S_rev = []
for t in (0.0, 5.0, 10.0, 20.0, 30.0, 40.0):
    S_rev.append((t, vn_entropy(ptrace(closed.evolve(rho_rev0, t), A, N))))
S7 = {
    "trace_distance_after_forward_then_backward_evolution": trace_distance(rho_back, rho0),
    "H_max_imag_part": is_real,
    "entropy_of_A_along_the_conjugated_reversed_preparation": {str(t): s for t, s in S_rev},
    "S_A_original_t0_and_t40": [float(SA[0]), float(SA[-1])],
    "note": "The time-reversed final state is an admissible preparation of the same closed whole and its forward evolution lowers S_A "
            "back to the original initial value. Closure supplies no arrow; a low-entropy preparation does.",
}
results["S7_time_reversal"] = S7

# ----------------------------------------------------------------------------- S8
# what closure does force (standard theorems) -- checked in the pure-global-state version of S1
psi_bath = np.linalg.eigh(H_bath_only)[1][:, 3]  # an excited pure bath eigenstate (pure global state)
rho0_pure = np.kron(rho_A0, proj(psi_bath))
checks = {"S_A_equals_S_Abar_max_dev": 0.0, "araki_lieb_violations": 0, "subadditivity_violations": 0, "purity_const_dev": 0.0}
for t in (0.0, 1.0, 5.0, 10.0, 20.0, 40.0):
    rho = closed.evolve(rho0_pure, t)
    rA, rB = ptrace(rho, A, N), ptrace(rho, ABAR, N)
    sA, sB, sT = vn_entropy(rA), vn_entropy(rB), vn_entropy(rho)
    checks["S_A_equals_S_Abar_max_dev"] = max(checks["S_A_equals_S_Abar_max_dev"], abs(sA - sB))
    if sT > sA + sB + 1e-9:
        checks["subadditivity_violations"] += 1
    if sT < abs(sA - sB) - 1e-9:
        checks["araki_lieb_violations"] += 1
    checks["purity_const_dev"] = max(checks["purity_const_dev"], abs(purity(rho) - 1.0))
# Poincare recurrence is a theorem for finite dimension but the recurrence time is not short: report the closest return within the grid
dmin, tmin = 1.0, 0.0
for t in t_grid[t_grid >= 1.0]:
    d = trace_distance(closed.evolve(rho0_pure, t), rho0_pure)
    if d < dmin:
        dmin, tmin = d, float(t)
checks["closest_return_to_initial_pure_state_within_grid_t_ge_1"] = {"trace_distance": float(dmin), "t": tmin}
checks["note"] = ("Closure forces: unitarity (global purity and entropy constant), conservation of the fixed global H, "
                  "S_A = S_Abar for a pure global state, subadditivity and Araki-Lieb, and eventual quasi-recurrence in finite dimension. "
                  "All are theorems of standard quantum mechanics; none is a GRUT contribution.")
results["S8_what_closure_forces"] = checks

# ----------------------------------------------------------------------------- provenance
src = (HERE / "closed_whole_demo.py").read_bytes()
results["provenance"] = {
    "script_sha256": hashlib.sha256(src).hexdigest(),
    "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
    "platform": platform.platform(), "wall_seconds": round(time.time() - t0, 2),
}


def rnd(o):
    if isinstance(o, float):
        return float(f"{o:.12g}")
    if isinstance(o, dict):
        return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list):
        return [rnd(v) for v in o]
    return o


(OUT / "closed_whole_demo.json").write_text(json.dumps(rnd(results), indent=1))

# ----------------------------------------------------------------------------- summary
inv = S1["invariants"]
print("S1 invariants: max|dE_tot| = %.2e, max|dS_tot| = %.2e, max|d purity| = %.2e; Born local-vs-global max dev = %.1e" %
      (inv["max_abs_dE_tot"], inv["max_abs_dS_tot"], inv["max_abs_dpurity_tot"], S1["born_global_minus_local_max"]))
print("S1 ranges: E_A in [%.4f, %.4f], E_int in [%.4f, %.4f], S_A in [%.4f, %.4f] (ln2 = %.4f)" %
      (*inv["range_E_A"], *inv["range_E_int"], *inv["range_S_A"], np.log(2)))
rc = S1["reduced_equation_check_t3"]
print("S1 reduced-equation check at t=3: |lhs-rhs| = %.1e (full), %.1e (without H_Abar term); |Tr_Abar[H_Abar,rho]| = %.1e" %
      (rc["max_abs_lhs_minus_rhs_full"], rc["max_abs_lhs_minus_rhs_without_H_Abar"], rc["max_abs_Tr_Abar_commutator_H_Abar"]))
print("S2 S_A: t0 %.4f, max %.4f at t=%.1f, min after max %.4f at t=%.1f, sign changes of dS/dt: %d; BLP backflow sum %.4f (D: %.3f -> min %.3f)" %
      (S2["S_A_t0"], S2["S_A_max"], S2["t_of_S_A_max"], S2["S_A_min_after_max"], S2["t_of_min_after_max"],
       S2["number_of_sign_changes_of_dS_A_dt_on_grid"], S2["BLP"]["backflow_sum_of_positive_increments"], S2["BLP"]["D_t0"], S2["BLP"]["D_min"]))
print("S3 CPTP: min Choi eig %.2e, TP dev %.1e; CP-divisibility most negative intermediate Choi eig %.4f" %
      (S3["cptp_check"]["min_choi_eigenvalue_over_samples"], S3["cptp_check"]["max_trace_preservation_deviation"], S3["cp_divisibility"]["most_negative"]))
print("S3 TCL generator |K(t)-K(0.5)| thermal:", {k: round(v, 4) for k, v in S3["tcl"]["generator_time_dependence_fro_norm_K(t)-K(0.5)_thermal"].items()})
print("S3 TCL generator |K_thermal-K_ground|:", {k: round(v, 4) for k, v in S3["tcl"]["generator_bath_state_dependence_fro_norm_K_thermal-K_ground"].items()})
for k in [kk for kk in S4 if kk != "note"]:
    s = S4[k]
    print("S4 %s: eps=%.4g (max %.4g), marginal dev %.1e, D(rhoA) t0 %.1e, max %.4f at t=%.1f, |inhom(0)|=%.4f" %
          (k, s["eps_used"], s["eps_max_for_positivity"], s["max_marginal_deviation_between_rho1_and_rho2"], s["trace_distance_rhoA_t0"],
           s["max_trace_distance_rhoA_over_grid"], s["t_of_max"], s["norm_of_inhomogeneous_term_at_t0"]))
print("S5 b' exact recurrence: D(T)=%.1e, S_A(T/2)=%.4f, S_A(T)=%.1e" % (S5["model_b_prime_exact_recurrence"]["trace_distance_rho_global_T_vs_0"], S5["model_b_prime_exact_recurrence"]["S_A_at_half_period"], S5["model_b_prime_exact_recurrence"]["S_A_at_period"]))
print("S5 trace distances rhoA at t=5:", {k: round(v, 4) for k, v in S5["trace_distance_rhoA"]["5.0"].items()},
      "| model b recurrence min S_A after max = %.2e at t=%.1f" % (S5["model_b_recurrence"]["min_S_A_after_first_max"], S5["model_b_recurrence"]["t_of_that_min"]))
print("S6 Born p(Z+,t1)=%.6f register pop=%.6f; p(Z+,t2) with record %.6f, without %.6f, dephased %.6f" %
      (S6["born_pZ_plus_at_t1"], S6["register_population_0_after_write"], S6["pZ_plus_A_at_t2_with_record_written"],
       S6["pZ_plus_A_at_t2_without_intervention"], S6["pZ_plus_A_at_t2_with_measure_and_forget_dephasing"]))
print("S6 record drift at t=10: kappa=0 %.1e, kappa=0.3 %.4f" % (agree["10.0"]["kappa=0_trace_distance_register_vs_written"], agree["10.0"]["kappa=0.3_trace_distance_register_vs_written"]))
print("S7 forward-backward trace distance %.1e; reversed preparation S_A: %s" % (S7["trace_distance_after_forward_then_backward_evolution"], {k: round(v, 4) for k, v in S7["entropy_of_A_along_the_conjugated_reversed_preparation"].items()}))
print("S8 pure-global checks: |S_A-S_Abar| max %.1e, subadditivity violations %d, Araki-Lieb violations %d, closest return D=%.4f at t=%.1f" %
      (checks["S_A_equals_S_Abar_max_dev"], checks["subadditivity_violations"], checks["araki_lieb_violations"],
       checks["closest_return_to_initial_pure_state_within_grid_t_ge_1"]["trace_distance"], checks["closest_return_to_initial_pure_state_within_grid_t_ge_1"]["t"]))
print("wall %.1fs; sha256 %s" % (results["provenance"]["wall_seconds"], results["provenance"]["script_sha256"][:16]))
