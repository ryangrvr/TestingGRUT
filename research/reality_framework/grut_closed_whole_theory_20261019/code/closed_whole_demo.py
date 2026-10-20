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
  S9  unitarily equivariant state-selection rules F(H) are stationary (finite-class obstruction)
  S10 conditional equilibration on average (Linden-Popescu-Short-Winter 2009) for pure global states
  S11 pointwise local inequalities that closure does force: athermality budget (closure alone) and the
      second-law identity Sigma = I + D for a complement initially in a Gibbs state uncorrelated with A
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
              "longitudinal fields h = [0.3,-0.7,0.5,-0.2,0.8,-0.4,0.6] (fixed disorder; no known integrability; "
              "H_Abar conserves total Z of the bath, and the full H conserves the Z-parity of all 8 qubits; nothing below needs non-integrability)",
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
loc_max = [i for i in range(1, len(SA) - 1) if SA[i] > SA[i - 1] and SA[i] >= SA[i + 1]]
i1 = loc_max[0]
j1 = next(j for j in range(i1 + 1, len(SA) - 1) if SA[j] < SA[j - 1] and SA[j] <= SA[j + 1])
S2 = {
    "first_local_max": {"t": float(t_grid[i1]), "S_A": float(SA[i1])},
    "first_local_min_after_it": {"t": float(t_grid[j1]), "S_A": float(SA[j1])},
    "local_maxima_times": [float(t_grid[i]) for i in loc_max],
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
    "grid_note": "sum over the 0.1 grid: a lower bound on the continuous Breuer-Laine-Piilo measure",
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
ZPAR = kron_all([Z] * N)  # global Z-parity: [H, ZPAR] = 0, Tr_A ZPAR = 0 = Tr_Abar ZPAR
S4_cases = [
    ("inf_T_bath_quantum_correlation_Y0X1", np.kron(rho_A_mm, rho_bath_inf), two_site(Y, 0, X, 1, N), rho_bath_inf),
    ("inf_T_bath_conserved_correlation_Zparity", np.kron(rho_A_mm, rho_bath_inf), ZPAR, rho_bath_inf),
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
        "norm_of_inhomogeneous_term_at_t0_A_factor": float(np.linalg.norm(inhom0)),
        "frobenius_norm_of_commutator_H_chi": float(np.linalg.norm(H @ chi - chi @ H)),
    }
S4["note_converse"] = ("The Zparity case has Q rho(0) = eps Z^(x8) != 0 but [H, Z^(x8)] = 0, so the inhomogeneous term vanishes identically "
                       "and the futures coincide: the inhomogeneous term vanishes IF the preparation is rho_A (x) rho_ref, not only if.")
S4["note"] = ("Two global preparations with identical rho_A(0) AND identical rho_Abar(0) evolve to different rho_A(t). "
              "No map on rho_A(0) alone (nor on both marginals) reproduces both futures; the Nakajima-Zwanzig inhomogeneous term "
              "carries the difference (its t = 0 value is recorded; for the classical Z0Z1 correlation the first two orders vanish at t = 0 and the "
              "difference appears at third order in t). Standard initial-correlation result (Pechukas 1994; Stelmachovic-Buzek 2001).")
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
P0A = site_op(proj(ZERO), 0, N9)
P1A = site_op(proj(ONE), 0, N9)
CNOT = P0A + P1A @ site_op(X, R, N9)
rho_t1 = closed9.evolve(rho0_9, t1)
pz_plus_t1 = float(np.real(np.trace(rho_t1 @ P0A)))
rho_after_write = CNOT @ rho_t1 @ CNOT.conj().T
energy_injected_by_write = float(np.real(np.trace(H9 @ rho_after_write) - np.trace(H9 @ rho_t1)))
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
# composed one-time-map prediction: Phi_{t2-t1} o (Z-dephasing) o Phi_{t1}, bath re-prepared at t1 (what "one-time maps + instrument" would give)
def apply_map(M, rho):
    coeffs = np.array([np.trace(P.conj().T @ rho) / 2 for P in PAULI])
    outc = M @ coeffs
    return sum(c * P for c, P in zip(outc, PAULI))


rA_t1_map = apply_map(dyn_map(rho_bath_th, t1), rho_A0)
rA_t1_deph = PZp @ rA_t1_map @ PZp + (I2 - PZp) @ rA_t1_map @ (I2 - PZp)
pA_composed = float(np.real(apply_map(dyn_map(rho_bath_th, t2 - t1), rA_t1_deph)[0, 0]))
# record persistence: (i) kappa X_R acting on R alone (a basis rotation, no information loss); (ii) kappa X_R X_7 coupling R to the bath
kappa = 0.3
H9k = H9 + kappa * site_op(X, R, N9)
H9c = H9 + kappa * two_site(X, R, X, 7, N9)
closed9k, closed9c = Closed(H9k), Closed(H9c)
# branch states by recorded value: sigma_z = P_z rho(t1) P_z (x) |z><z|_R ; Helstrom success of reading the recorded bit from R alone
P0A8, P1A8 = site_op(proj(ZERO), 0, N), site_op(proj(ONE), 0, N)
rho_t1_8 = closed.evolve(rho0, t1)
sig_plus = np.kron(P0A8 @ rho_t1_8 @ P0A8, proj(ZERO))
sig_minus = np.kron(P1A8 @ rho_t1_8 @ P1A8, proj(ONE))


def helstrom_from_register(cl, tt):
    a = ptrace(cl.evolve(sig_plus, tt - t1), [R], N9)
    b = ptrace(cl.evolve(sig_minus, tt - t1), [R], N9)
    return float(0.5 * (1 + np.sum(np.abs(np.linalg.eigvalsh(herm(a - b))))))


agree = {}
for tt2 in (3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
    r_k = closed9k.evolve(rho_after_write, tt2 - t1)
    r_0 = closed9.evolve(rho_after_write, tt2 - t1)
    agree[str(tt2)] = {
        "idle_register_trace_distance_vs_written": trace_distance(ptrace(r_0, [R], N9), rR),
        "kappa_X_R_fixed_Zbasis_trace_distance_vs_written": trace_distance(ptrace(r_k, [R], N9), rR),
        "kappa_X_R_optimal_readout_success": helstrom_from_register(closed9k, tt2),
        "kappa_X_R_X_7_optimal_readout_success": helstrom_from_register(closed9c, tt2),
    }
S6 = {
    "construction": "Omega' = A u Abar u R (9 qubits); H' = H x I_R; record written at t1=3 by the unitary CNOT(A->R), applied as an "
                    "externally timed intervention (a supplied control and clock, not generated by H'); readout at t2=6.",
    "energy_injected_by_externally_timed_write": energy_injected_by_write,
    "born_pZ_plus_at_t1": pz_plus_t1, "register_population_0_after_write": record_pop_0,
    "abs_difference": abs(pz_plus_t1 - record_pop_0),
    "joint_distribution_zA_t2_zR": joint,
    "pZ_plus_A_at_t2_with_record_written": pA_int,
    "pZ_plus_A_at_t2_without_intervention": pA_noint,
    "pZ_plus_A_at_t2_with_measure_and_forget_dephasing": pA_deph,
    "pZ_plus_A_at_t2_composed_one_time_maps_bath_reprepared": pA_composed,
    "note_two_time": "The uninterrupted one-time map gives 0.449; composing one-time maps with the instrument (bath re-prepared at t1) gives "
                     "the composed value; the exact closed-whole result with the record written differs from both. Multi-time statistics are "
                     "properties of the process (process tensor), not of Phi_t. Two recorded equalities are identities, not evidence: "
                     "register population = Born probability (CNOT with a fresh |0> is a Stinespring dilation of Z-dephasing on A), and "
                     "unitary record = measure-and-forget (the idle register commutes with the partial trace).",
    "record_persistence": agree,
    "note_persistence": "With [H', Z_R] = 0 the register is unchanged. kappa X_R alone rotates the record out of the fixed Z_R readout basis but "
                        "loses no information (optimal readout success stays 1). kappa X_R X_7 couples the register to the bath and the recorded "
                        "bit becomes genuinely less retrievable from R. Record persistence is a property of the supplied dynamics, not of closure.",
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
rev_identity_dev = 0.0
for t in (0.0, 5.0, 10.0, 20.0, 30.0, 40.0):
    r_rev = closed.evolve(rho_rev0, t)
    S_rev.append((t, vn_entropy(ptrace(r_rev, A, N))))
    rev_identity_dev = max(rev_identity_dev, float(np.max(np.abs(r_rev - closed.evolve(rho0, T - t).conj()))))
SA_rev_grid = np.array([vn_entropy(ptrace(closed.evolve(rho_rev0, t), A, N)) for t in t_grid])
dSr = np.diff(SA_rev_grid)
S7 = {
    "trace_distance_after_forward_then_backward_evolution": trace_distance(rho_back, rho0),
    "H_max_imag_part": is_real,
    "entropy_of_A_along_the_conjugated_reversed_preparation": {str(t): s for t, s in S_rev},
    "S_A_original_t0_and_t40": [float(SA[0]), float(SA[-1])],
    "max_abs_rho_rev(t)_minus_rho(T-t)_conj": rev_identity_dev,
    "reversed_run_S_A_max_and_t": [float(SA_rev_grid.max()), float(t_grid[int(np.argmax(SA_rev_grid))])],
    "reversed_run_number_of_increasing_steps_of_400": int(np.sum(dSr > 0)),
    "note": "The conjugated final state is an admissible preparation of the same closed whole; rho_rev(t) = rho(T-t)^* exactly, so its forward "
            "evolution is the original run read backwards: S_A returns non-monotonically to 0 at t = 40. Closure supplies no arrow; the "
            "low-entropy preparation does.",
}
results["S7_time_reversal"] = S7

# ----------------------------------------------------------------------------- S8
# what closure does force (standard theorems) -- checked in the pure-global-state version of S1
psi_bath = np.linalg.eigh(H_bath_only)[1][:, 3]  # an excited pure bath eigenstate (pure global state)
rho0_pure = np.kron(rho_A0, proj(psi_bath))
checks = {"pure_global_state_S_A_equals_S_Abar_max_dev": 0.0, "pure_global_state_purity_const_dev": 0.0}
for t in (0.0, 1.0, 5.0, 10.0, 20.0, 40.0):
    rho = closed.evolve(rho0_pure, t)
    rA, rB = ptrace(rho, A, N), ptrace(rho, ABAR, N)
    checks["pure_global_state_S_A_equals_S_Abar_max_dev"] = max(checks["pure_global_state_S_A_equals_S_Abar_max_dev"], abs(vn_entropy(rA) - vn_entropy(rB)))
    checks["pure_global_state_purity_const_dev"] = max(checks["pure_global_state_purity_const_dev"], abs(purity(rho) - 1.0))
# subadditivity and Araki-Lieb are only non-trivial for a mixed global state: evaluate them on the S1 trajectory (S_tot = 1.26)
_SA, _SB, _ST = np.array(S1["S_A"]), np.array(S1["S_Abar"]), np.array(S1["S_tot"])
sub_slack = _SA + _SB - _ST          # = mutual information >= 0
al_slack = _ST - np.abs(_SA - _SB)   # Araki-Lieb >= 0
checks["mixed_S1_trajectory"] = {
    "S_tot": float(_ST[0]),
    "subadditivity_min_slack_MI": float(sub_slack.min()), "subadditivity_violations": int(np.sum(sub_slack < -1e-9)),
    "araki_lieb_min_slack": float(al_slack.min()), "araki_lieb_slack_at_t3_t20": [float(al_slack[np.argmin(np.abs(t_grid - 3))]), float(al_slack[np.argmin(np.abs(t_grid - 20))])],
    "araki_lieb_violations": int(np.sum(al_slack < -1e-9)),
    "note": "Both inequalities are saturated at t = 0 (product state with pure A) and strict afterwards; on a pure global state they reduce to S_A = S_Abar and test nothing more.",
}
# Poincare recurrence is a theorem for finite dimension but the recurrence time is not short: report the closest return within the grid
dmin, tmin = 1.0, 0.0
for t in t_grid[t_grid >= 1.0]:
    d = trace_distance(closed.evolve(rho0_pure, t), rho0_pure)
    if d < dmin:
        dmin, tmin = d, float(t)
checks["closest_return_to_initial_pure_state_within_grid_t_ge_1"] = {"trace_distance": float(dmin), "t": tmin}
checks["note"] = ("Closure forces: unitarity (global purity and entropy constant), conservation of the fixed global H, "
                  "S_A = S_Abar for a pure global state, subadditivity and Araki-Lieb (checked on the mixed trajectory), and eventual "
                  "quasi-recurrence in finite dimension. All are theorems of standard quantum mechanics; none is a GRUT contribution.")
results["S8_what_closure_forces"] = checks

# ----------------------------------------------------------------------------- S9
# Any rule F that picks a global state from H alone and is unitarily equivariant, F(U H U^+) = U F(H) U^+, satisfies
# [F(H), H] = 0 (take U = exp(iHs)); the selected state is stationary, so no subsystem has any dynamics.
S9 = {}
for name, F in (("gibbs_beta_1_of_global_H", thermal(H, 1.0)), ("ground_projector_of_global_H", ground(H))):
    devs = [trace_distance(ptrace(closed.evolve(F, t), A, N), ptrace(F, A, N)) for t in (1.0, 5.0, 10.0, 40.0)]
    S9[name] = {"frobenius_norm_commutator_F_H": float(np.linalg.norm(F @ H - H @ F)),
                "max_trace_distance_rhoA_t_vs_0": float(max(devs))}
S9["note"] = ("A unitarily equivariant selection rule F(H) yields a stationary state: every restriction is frozen. A selection rule "
              "with local dynamical content must depend on more than H (the access structure, a clock, or a non-equivariant datum).")
results["S9_equivariant_state_selection_is_stationary"] = S9

# ----------------------------------------------------------------------------- S10
# Linden-Popescu-Short-Winter, PRE 79, 061103 (2009), Theorem 1: pure global state, H with non-degenerate energy gaps =>
# <D(rho_A(t), omega_A)>_t <= (1/2) sqrt(d_A * <||rho_A(t) - omega_A||_2^2>_t) <= (1/2) sqrt(d_A^2 / d_eff),
# d_eff = 1 / sum_k p_k^2 with p_k = |<E_k|psi>|^2, omega = time-averaged (dephased) state.
E_all, V_all = closed.E, closed.V
spac = np.diff(np.sort(E_all))
gaps = np.array([E_all[k] - E_all[l] for k in range(len(E_all)) for l in range(k)])
gaps_sorted = np.sort(gaps)
min_gap_coincidence = float(np.min(np.diff(gaps_sorted)))
Vr = V_all.reshape(2, 2 ** NB, -1)                        # (a, b, k)
Mkl = np.einsum("abk,cbl->klac", Vr, Vr.conj())           # Tr_Abar |k><l|, shape (256,256,2,2)
Mkl_norm2 = np.einsum("klac,klac->kl", Mkl, Mkl.conj()).real


def lpsw(psi, label):
    c = V_all.conj().T @ psi
    pk = np.abs(c) ** 2
    d_eff = float(1.0 / np.sum(pk ** 2))
    omega = (V_all * pk) @ V_all.conj().T
    omega_A = ptrace(omega, A, N)
    hs2_avg = float(np.sum(np.outer(pk, pk) * Mkl_norm2) - np.sum(pk ** 2 * np.diag(Mkl_norm2)))  # exact infinite-time average, non-degenerate gaps
    rng = np.random.default_rng(11)
    ts = np.sort(rng.uniform(0, 1.0e4, 3000))
    Ds = []
    for t in ts:
        ph = c * np.exp(-1j * E_all * t)
        v = V_all @ ph
        rA = ptrace(np.outer(v, v.conj()), A, N)
        Ds.append(trace_distance(rA, omega_A))
    return {"preparation": label, "d_eff": d_eff,
            "lpsw_bound_half_sqrt_dA2_over_deff": float(0.5 * np.sqrt(4.0 / d_eff)),
            "intermediate_bound_half_sqrt_dA_times_exact_avg_HS2": float(0.5 * np.sqrt(2.0 * hs2_avg)),
            "sampled_time_average_D_3000_random_t_in_0_1e4": float(np.mean(Ds)),
            "sampled_max_D": float(np.max(Ds))}


psi_S8 = np.kron(PLUS, psi_bath)
neel = kron_all([np.array([[1.0], [0.0]]) if k % 2 == 0 else np.array([[0.0], [1.0]]) for k in range(NB)]).ravel().astype(complex)
psi_neel = np.kron(PLUS, neel)
S10 = {
    "spectrum_min_level_spacing": float(spac.min()),
    "min_separation_between_distinct_Bohr_gaps": min_gap_coincidence,
    "non_degenerate_gaps_numerically": bool(spac.min() > 1e-8 and min_gap_coincidence > 1e-8),
    "cases": [lpsw(psi_S8, "|+> (x) excited bath eigenstate (the S8 pure preparation)"),
              lpsw(psi_neel, "|+> (x) bath Neel product state |0101010>")],
    "note": ("Closure plus the generic hypothesis of non-degenerate energy gaps forces equilibration ON AVERAGE: the infinite-time average "
             "of D(rho_A(t), omega_A) is bounded by (1/2) sqrt(d_A^2/d_eff). The size of the bound is set by the supplied preparation "
             "through d_eff; it does not force pointwise or monotone relaxation, a thermal omega_A, or a time scale."),
}
results["S10_equilibration_on_average"] = S10

# ----------------------------------------------------------------------------- S11
def logm_psd(r, floor=1e-300):
    w, V = eigh(herm(r))
    return (V * np.log(np.clip(w, floor, None))) @ V.conj().T


def rel_entropy(r, sig):
    w, V = eigh(herm(r))
    w = np.clip(w, 0, None)
    lr = (V * np.log(np.where(w > 1e-15, w, 1.0))) @ V.conj().T
    return float(np.real(np.trace(r @ lr) - np.trace(r @ logm_psd(sig))))


t_coarse = np.arange(0.0, 40.0 + 1e-9, 1.0)
sigma_glob = thermal(H, 1.0)               # stationary: [sigma, H] = 0
sigma_A = ptrace(sigma_glob, A, N)
D_glob, D_loc = [], []
Sig, MIc, DB = [], [], []
SA0 = vn_entropy(ptrace(rho0, A, N))
EB0 = float(np.real(np.trace(H_B @ rho0)))
for t in t_coarse:
    rho = closed.evolve(rho0, t)
    rA, rB = ptrace(rho, A, N), ptrace(rho, ABAR, N)
    D_glob.append(rel_entropy(rho, sigma_glob))
    D_loc.append(rel_entropy(rA, sigma_A))
    sA, sB, sT = vn_entropy(rA), vn_entropy(rB), vn_entropy(rho)
    EB = float(np.real(np.trace(H_B @ rho)))
    Sig.append((sA - SA0) + 1.0 * (EB - EB0))
    MIc.append(sA + sB - sT)
    DB.append(rel_entropy(rB, rho_bath_th))
D_glob, D_loc, Sig, MIc, DB = map(np.array, (D_glob, D_loc, Sig, MIc, DB))
S11 = {
    "athermality_budget": {
        "sigma": "Gibbs state of the global H at beta = 1 (any stationary sigma works)",
        "D_global_t0": float(D_glob[0]), "max_abs_drift_D_global": float(np.max(np.abs(D_glob - D_glob[0]))),
        "max_D_local_over_t_in_0_40_step_1": float(D_loc.max()), "t_of_max": float(t_coarse[int(np.argmax(D_loc))]),
        "violations_of_D_local_le_D_global": int(np.sum(D_loc > D_glob[0] + 1e-9)),
        "note": "Closure alone: D(rho_Omega(t)||sigma) is conserved for stationary sigma, and monotonicity of relative entropy under "
                "partial trace gives D(rho_A(t)||Tr_Abar sigma) <= D(rho_Omega(0)||sigma) at every t, for every subsystem.",
    },
    "second_law_identity": {
        "hypothesis": "rho_Omega(0) = rho_A(0) (x) gamma_Abar with gamma_Abar the Gibbs state of H_Abar at beta = 1 (the S1 preparation)",
        "definition": "Sigma(t) = [S_A(t) - S_A(0)] + beta [E_Abar(t) - E_Abar(0)],  E_Abar = Tr(H_Abar rho)",
        "max_abs_Sigma_minus_I_plus_D": float(np.max(np.abs(Sig - (MIc + DB)))),
        "min_Sigma": float(Sig.min()),
        "Sigma_at_t1_t10_t20_t40": [float(Sig[int(k)]) for k in (1, 10, 20, 40)],
        "number_of_decreases_of_Sigma_on_unit_grid": int(np.sum(np.diff(Sig) < -1e-12)),
        "note": "Esposito-Lindenberg-Van den Broeck identity Sigma = I(A:Abar) + D(rho_Abar||gamma_Abar) >= 0 at every t (exact, from global "
                "unitarity). Pointwise non-negative but not monotone; Landauer's erasure bound follows from it.",
    },
}
results["S11_pointwise_inequalities_forced"] = S11

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
for k in [kk for kk in S4 if not kk.startswith("note")]:
    s = S4[k]
    print("S4 %s: eps=%.4g (max %.4g), marginal dev %.1e, D(rhoA) t0 %.1e, max %.4g at t=%.1f, |inhom(0) A-factor|=%.4f" %
          (k, s["eps_used"], s["eps_max_for_positivity"], s["max_marginal_deviation_between_rho1_and_rho2"], s["trace_distance_rhoA_t0"],
           s["max_trace_distance_rhoA_over_grid"], s["t_of_max"], s["norm_of_inhomogeneous_term_at_t0_A_factor"]))
print("S5 b' exact recurrence: D(T)=%.1e, S_A(T/2)=%.4f, S_A(T)=%.1e" % (S5["model_b_prime_exact_recurrence"]["trace_distance_rho_global_T_vs_0"], S5["model_b_prime_exact_recurrence"]["S_A_at_half_period"], S5["model_b_prime_exact_recurrence"]["S_A_at_period"]))
print("S5 trace distances rhoA at t=5:", {k: round(v, 4) for k, v in S5["trace_distance_rhoA"]["5.0"].items()},
      "| model b recurrence min S_A after max = %.2e at t=%.1f" % (S5["model_b_recurrence"]["min_S_A_after_first_max"], S5["model_b_recurrence"]["t_of_that_min"]))
print("S6 Born p(Z+,t1)=%.6f register pop=%.6f; p(Z+,t2) with record %.6f, without %.6f, dephased %.6f" %
      (S6["born_pZ_plus_at_t1"], S6["register_population_0_after_write"], S6["pZ_plus_A_at_t2_with_record_written"],
       S6["pZ_plus_A_at_t2_without_intervention"], S6["pZ_plus_A_at_t2_with_measure_and_forget_dephasing"]))
print("S6 composed one-time-map prediction p(Z+,t2) = %.6f (exact with record %.6f)" % (pA_composed, pA_int))
print("S6 record at t=10: idle %.1e; kappa X_R fixed-basis TD %.4f, optimal readout %.6f; kappa X_R X_7 optimal readout %.6f" % (
    agree["10.0"]["idle_register_trace_distance_vs_written"], agree["10.0"]["kappa_X_R_fixed_Zbasis_trace_distance_vs_written"],
    agree["10.0"]["kappa_X_R_optimal_readout_success"], agree["10.0"]["kappa_X_R_X_7_optimal_readout_success"]))
print("S7 forward-backward trace distance %.1e; reversed preparation S_A: %s" % (S7["trace_distance_after_forward_then_backward_evolution"], {k: round(v, 4) for k, v in S7["entropy_of_A_along_the_conjugated_reversed_preparation"].items()}))
print("S7 identity rho_rev(t)=rho(T-t)^*: %.1e; reversed run S_A max %.4f at t=%.1f, increasing steps %d/400" % (rev_identity_dev, *S7["reversed_run_S_A_max_and_t"], S7["reversed_run_number_of_increasing_steps_of_400"]))
print("S2 first local max S_A %.4f at t=%.1f -> min %.4f at t=%.1f; %d local maxima" % (S2["first_local_max"]["S_A"], S2["first_local_max"]["t"], S2["first_local_min_after_it"]["S_A"], S2["first_local_min_after_it"]["t"], len(S2["local_maxima_times"])))
print("S8 pure-global |S_A-S_Abar| max %.1e; mixed S1: min MI %.1e, min Araki-Lieb slack %.4f (t3 %.4f, t20 %.4f), violations %d/%d; closest return D=%.4f at t=%.1f" %
      (checks["pure_global_state_S_A_equals_S_Abar_max_dev"], checks["mixed_S1_trajectory"]["subadditivity_min_slack_MI"], checks["mixed_S1_trajectory"]["araki_lieb_min_slack"],
       *checks["mixed_S1_trajectory"]["araki_lieb_slack_at_t3_t20"], checks["mixed_S1_trajectory"]["subadditivity_violations"], checks["mixed_S1_trajectory"]["araki_lieb_violations"],
       checks["closest_return_to_initial_pure_state_within_grid_t_ge_1"]["trace_distance"], checks["closest_return_to_initial_pure_state_within_grid_t_ge_1"]["t"]))
print("S4 converse: Zparity case max D = %.1e with |[H,chi]| = %.1e" % (S4["inf_T_bath_conserved_correlation_Zparity"]["max_trace_distance_rhoA_over_grid"], S4["inf_T_bath_conserved_correlation_Zparity"]["frobenius_norm_of_commutator_H_chi"]))
print("S9 equivariant rules: %s" % {k: (round(v["frobenius_norm_commutator_F_H"], 16), round(v["max_trace_distance_rhoA_t_vs_0"], 16)) for k, v in S9.items() if k != "note"})
print("S10 gaps: min spacing %.3e, min gap coincidence %.3e, non-degenerate %s" % (S10["spectrum_min_level_spacing"], S10["min_separation_between_distinct_Bohr_gaps"], S10["non_degenerate_gaps_numerically"]))
for cse in S10["cases"]:
    print("S10 %s: d_eff %.2f, bound %.4f, intermediate %.4f, sampled <D> %.4f, max D %.4f" % (cse["preparation"][:30], cse["d_eff"], cse["lpsw_bound_half_sqrt_dA2_over_deff"], cse["intermediate_bound_half_sqrt_dA_times_exact_avg_HS2"], cse["sampled_time_average_D_3000_random_t_in_0_1e4"], cse["sampled_max_D"]))
a11, b11 = S11["athermality_budget"], S11["second_law_identity"]
print("S6 energy injected by the externally timed write: %.6f" % energy_injected_by_write)
print("S11 athermality: D_glob %.6f (drift %.1e), max D_loc %.6f at t=%.0f, violations %d" % (a11["D_global_t0"], a11["max_abs_drift_D_global"], a11["max_D_local_over_t_in_0_40_step_1"], a11["t_of_max"], a11["violations_of_D_local_le_D_global"]))
print("S11 second law: |Sigma-(I+D)| %.1e, min Sigma %.2e, Sigma(1,10,20,40) %s, decreases %d" % (b11["max_abs_Sigma_minus_I_plus_D"], b11["min_Sigma"], [round(x, 4) for x in b11["Sigma_at_t1_t10_t20_t40"]], b11["number_of_decreases_of_Sigma_on_unit_grid"]))
print("wall %.1fs; sha256 %s" % (results["provenance"]["wall_seconds"], results["provenance"]["script_sha256"][:16]))
