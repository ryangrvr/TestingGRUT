"""
Independent re-implementation of the closed-whole demo, written from the
specification only (hbar = 1).

Convention: qubit 0 is the leftmost (most significant) tensor factor.
Operator on qubit k = I x ... x P (position k) x ... x I.
"""
import json
import time
import numpy as np
import scipy.linalg as sla

np.set_printoptions(precision=12)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_all(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def op_on(n, pauli_dict):
    """pauli_dict: {qubit_index: 2x2 matrix}; identity elsewhere; n qubits total."""
    ops = [pauli_dict.get(k, I2) for k in range(n)]
    return kron_all(ops)


# ---------------------------------------------------------------- Hamiltonian
N = 8
h_field = [0.3, -0.7, 0.5, -0.2, 0.8, -0.4, 0.6]  # on qubits 1..7


def build_H_Abar_as_7qubit():
    """H_Abar on the 7-qubit chain (qubits 1..7 relabelled 0..6)."""
    n = 7
    H = np.zeros((2**n, 2**n), dtype=complex)
    for k in range(n - 1):  # k=0..5  <-> spec k=1..6
        H += op_on(n, {k: X, k + 1: X})
        H += op_on(n, {k: Y, k + 1: Y})
        H += 0.5 * op_on(n, {k: Z, k + 1: Z})
    for k in range(n):
        H += h_field[k] * op_on(n, {k: Z})
    return H


def build_H_global():
    H = np.zeros((2**N, 2**N), dtype=complex)
    # H_A
    H += 0.5 * op_on(N, {0: Z})
    # H_Abar on qubits 1..7
    for k in range(1, 7):  # k = 1..6, pairs (k, k+1)
        H += op_on(N, {k: X, k + 1: X})
        H += op_on(N, {k: Y, k + 1: Y})
        H += 0.5 * op_on(N, {k: Z, k + 1: Z})
    for k in range(1, 8):
        H += h_field[k - 1] * op_on(N, {k: Z})
    # H_AAbar
    H += 0.5 * op_on(N, {0: X, 1: X})
    return H


H_Abar7 = build_H_Abar_as_7qubit()
H = build_H_global()
assert np.allclose(H, H.conj().T)
# consistency: H == 0.5 Z_0 + I x H_Abar7 + 0.5 X0X1
assert np.allclose(H, 0.5 * op_on(N, {0: Z}) + np.kron(I2, H_Abar7) + 0.5 * op_on(N, {0: X, 1: X}))

# ---------------------------------------------------------------- initial states
plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
rho_plus = np.outer(plus, plus.conj())
rho_minus = np.outer(minus, minus.conj())

beta = 1.0
wB, VB = np.linalg.eigh(H_Abar7)
boltz = np.exp(-beta * (wB - wB.min()))
rho_bath = (VB * boltz) @ VB.conj().T
rho_bath /= np.trace(rho_bath).real
assert np.allclose(rho_bath, rho_bath.conj().T)
assert abs(np.trace(rho_bath) - 1) < 1e-12
# cross-check against expm
rb2 = sla.expm(-beta * H_Abar7)
rb2 /= np.trace(rb2)
assert np.allclose(rho_bath, rb2, atol=1e-12)

rho0_plus = np.kron(rho_plus, rho_bath)
rho0_minus = np.kron(rho_minus, rho_bath)

# ---------------------------------------------------------------- evolution
w, V = np.linalg.eigh(H)
Vd = V.conj().T


def evolve(rho0, t):
    """Exact unitary evolution via eigendecomposition."""
    phase = np.exp(-1j * w * t)
    U = (V * phase) @ Vd
    return U @ rho0 @ U.conj().T


def ptrace_to_A(rho, nA=1, n=N):
    dA = 2**nA
    dB = 2**(n - nA)
    r = rho.reshape(dA, dB, dA, dB)
    return np.einsum('ajbj->ab', r)


def vn_entropy(rho):
    lam = np.linalg.eigvalsh(rho)
    lam = lam[lam > 1e-15]
    return float(-np.sum(lam * np.log(lam)))


def trace_norm(M):
    return float(np.sum(np.linalg.svd(M, compute_uv=False)))


grid = np.round(np.arange(0, 401) * 0.1, 10)  # t in {0, 0.1, ..., 40}
assert len(grid) == 401 and abs(grid[-1] - 40.0) < 1e-12

results = {}

# ---------------------------------------------------------------- 1. energy conservation
E0 = float(np.trace(H @ rho0_plus).real)
E_t = []
S_A = []
D_t = []
t0 = time.time()
for t in grid:
    rp = evolve(rho0_plus, t)
    rm = evolve(rho0_minus, t)
    E_t.append(float(np.trace(H @ rp).real))
    rA_p = ptrace_to_A(rp)
    rA_m = ptrace_to_A(rm)
    S_A.append(vn_entropy(rA_p))
    D_t.append(0.5 * trace_norm(rA_p - rA_m))
E_t = np.array(E_t)
S_A = np.array(S_A)
D_t = np.array(D_t)
print(f"[evolution loop] {time.time()-t0:.2f} s")

results['1_energy'] = {
    'E_tot_t0': E0,
    'max_abs_drift': float(np.max(np.abs(E_t - E0))),
}

# ---------------------------------------------------------------- 2. S_A
i10 = int(np.argmin(np.abs(grid - 10.0)))
i40 = int(np.argmin(np.abs(grid - 40.0)))
imax = int(np.argmax(S_A))
results['2_entropy'] = {
    'S_A_t0': float(S_A[0]),
    'S_A_t10': float(S_A[i10]),
    'S_A_t40': float(S_A[i40]),
    'S_A_max': float(S_A[imax]),
    't_of_max': float(grid[imax]),
    'ln2': float(np.log(2)),
}

# ---------------------------------------------------------------- 3. BLP
dD = np.diff(D_t)
results['3_BLP'] = {
    'D_t0': float(D_t[0]),
    'sum_positive_increments': float(np.sum(dD[dD > 0])),
    'D_min': float(np.min(D_t)),
    't_of_min': float(grid[int(np.argmin(D_t))]),
    'n_positive_increments': int(np.sum(dD > 0)),
}

# ---------------------------------------------------------------- 4. record protocol (9 qubits)
N9 = 9
H9 = np.kron(H, I2)  # H' = H x I_R, R is qubit 8 (last factor)
w9, V9 = np.linalg.eigh(H9)
V9d = V9.conj().T


def evolve9(rho0, t):
    U = (V9 * np.exp(-1j * w9 * t)) @ V9d
    return U @ rho0 @ U.conj().T


ket0 = np.array([1, 0], dtype=complex)
rhoR0 = np.outer(ket0, ket0.conj())
rho9_0 = np.kron(rho0_plus, rhoR0)

# CNOT control qubit 0, target qubit 8
P0 = np.outer(ket0, ket0.conj())
P1 = I2 - P0
CNOT = op_on(N9, {0: P0}) + op_on(N9, {0: P1, 8: X})
assert np.allclose(CNOT @ CNOT, np.eye(2**N9))

t1, t2 = 3.0, 6.0
rho9_t1 = evolve9(rho9_0, t1)
rhoA_t1 = ptrace_to_A(rho9_t1, nA=1, n=N9)
pZ_t1 = float(rhoA_t1[0, 0].real)

rho9_after = CNOT @ rho9_t1 @ CNOT.conj().T
rho9_t2_rec = evolve9(rho9_after, t2 - t1)
rhoA_t2_rec = ptrace_to_A(rho9_t2_rec, nA=1, n=N9)
pZ_t2_rec = float(rhoA_t2_rec[0, 0].real)

rho9_t2_free = evolve9(rho9_0, t2)
rhoA_t2_free = ptrace_to_A(rho9_t2_free, nA=1, n=N9)
pZ_t2_free = float(rhoA_t2_free[0, 0].real)

# cross-check: no-intervention 9-qubit result == 8-qubit result at t=6
rhoA_t2_free8 = ptrace_to_A(evolve(rho0_plus, t2))
assert abs(rhoA_t2_free8[0, 0].real - pZ_t2_free) < 1e-10

# record readout: p(Z_R = +1) after CNOT should equal pZ_t1 and persist
rhoR_after = np.einsum('iaib->ab', rho9_t2_rec.reshape(256, 2, 256, 2))
pZR_t2 = float(rhoR_after[0, 0].real)

results['4_record'] = {
    'p_Z0_plus_t1_before_CNOT': pZ_t1,
    'p_Z0_plus_t2_with_record': pZ_t2_rec,
    'p_Z0_plus_t2_no_intervention': pZ_t2_free,
    'p_ZR_plus_t2_with_record (record persists)': pZR_t2,
    'Z0_t1_expectation': 2 * pZ_t1 - 1,
}

# ---------------------------------------------------------------- 5. initial correlations
eps = 1.0 / 512
rho1 = np.eye(256, dtype=complex) / 256
rho2 = rho1 + eps * op_on(N, {0: Y, 1: X})
assert np.allclose(rho2, rho2.conj().T)
assert np.linalg.eigvalsh(rho2).min() > 0
assert abs(np.trace(rho2) - 1) < 1e-12


def ptrace_to_Abar(rho):
    return np.einsum('iaib->ab', rho.reshape(2, 128, 2, 128))


marg_A_diff = trace_norm(ptrace_to_A(rho1) - ptrace_to_A(rho2))
marg_Abar_diff = trace_norm(ptrace_to_Abar(rho1) - ptrace_to_Abar(rho2))
TD = []
for t in grid:
    rA1 = ptrace_to_A(evolve(rho1, t))
    rA2 = ptrace_to_A(evolve(rho2, t))
    TD.append(0.5 * trace_norm(rA1 - rA2))
TD = np.array(TD)
iTD = int(np.argmax(TD))
results['5_initial_correlations'] = {
    'marginal_A_trace_norm_diff': marg_A_diff,
    'marginal_Abar_trace_norm_diff': marg_Abar_diff,
    'max_trace_distance_rhoA': float(TD[iTD]),
    't_of_max': float(grid[iTD]),
    'trace_distance_t0': float(TD[0]),
}

# ---------------------------------------------------------------- 6. Werner compatibility (SDP)
import cvxpy as cp

phi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
Phi = np.outer(phi, phi.conj()).real  # real symmetric


def werner(p):
    return p * Phi + (1 - p) * np.eye(4) / 4


# 3-qubit state, ordering A(0) B(1) C(2); real symmetric suffices (W real) but
# use Hermitian for generality.
rhoABC = cp.Variable((8, 8), hermitian=True)
p_var = cp.Variable()


def ptrace_C(R):
    # R (8x8) -> (4x4): index (a,b,c)
    blocks = []
    for i in range(4):
        row = []
        for j in range(4):
            row.append(cp.trace(R[2 * i:2 * i + 2, 2 * j:2 * j + 2]))
        blocks.append(row)
    return cp.bmat([[blocks[i][j] for j in range(4)] for i in range(4)]) if False else cp.vstack([cp.hstack([cp.reshape(blocks[i][j], (1, 1), order='C') for j in range(4)]) for i in range(4)])


def ptrace_A(R):
    # trace over first qubit: sum of the two diagonal 4x4 blocks
    return R[0:4, 0:4] + R[4:8, 4:8]


W_expr = p_var * Phi + (1 - p_var) * np.eye(4) / 4
cons = [rhoABC >> 0, cp.real(cp.trace(rhoABC)) == 1,
        ptrace_C(rhoABC) == W_expr,
        ptrace_A(rhoABC) == W_expr,
        p_var >= 0, p_var <= 1]
prob = cp.Problem(cp.Maximize(p_var), cons)
prob.solve(solver=cp.CLARABEL)
p_sdp = float(p_var.value)
rho_sol = rhoABC.value
# verify solution
rho_sol = (rho_sol + rho_sol.conj().T) / 2
min_eig = float(np.linalg.eigvalsh(rho_sol).min())
# numeric partial traces
r8 = rho_sol.reshape(2, 2, 2, 2, 2, 2)
rAB = np.einsum('abcdec->abde', r8).reshape(4, 4)
rBC = np.einsum('abcade->bcde', r8).reshape(4, 4)
res_AB = float(np.abs(rAB - werner(p_sdp)).max())
res_BC = float(np.abs(rBC - werner(p_sdp)).max())

# analytic cross-check: Werner state 2-extendible iff p <= 2/3
# (two-qubit symmetric extension criterion Tr rho_B^2 >= Tr rho_AB^2 - 4 sqrt(det rho_AB))
def sym_ext_margin(p):
    lam = np.array([(1 + 3 * p) / 4] + [(1 - p) / 4] * 3)
    return 0.5 - (np.sum(lam**2) - 4 * np.sqrt(np.prod(lam)))

from scipy.optimize import brentq
p_analytic = brentq(sym_ext_margin, 0.5, 0.99)

results['6_werner'] = {
    'p_max_sdp': p_sdp,
    'sdp_status': prob.status,
    'min_eig_solution': min_eig,
    'residual_AB_marginal': res_AB,
    'residual_BC_marginal': res_BC,
    'p_max_analytic_2_over_3': float(p_analytic),
}

out = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/reimpl_results.json'
with open(out, 'w') as f:
    json.dump(results, f, indent=2)
print(json.dumps(results, indent=2))

# also dump the time series for later comparison
np.savez('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/series.npz',
         grid=grid, E_t=E_t, S_A=S_A, D_t=D_t, TD=TD)
