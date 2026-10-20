"""
Independent re-implementation of the closed-whole finite-model checks,
written from the specification only (hbar = 1).

8 qubits: A = qubit 0, Abar = qubits 1..7 (qubit 0 is the leftmost kron factor).
H_A      = 0.5 Z_0
H_Abar   = sum_{k=1..6} [X_k X_{k+1} + Y_k Y_{k+1} + 0.5 Z_k Z_{k+1}] + sum_{k=1..7} h_k Z_k
H_AAbar  = 0.5 X_0 X_1
H        = H_A + H_Abar + H_AAbar
rho(0)   = |+><+|_0 (x) exp(-H_Abar)/Z     (beta = 1, 7-qubit chain operator)
"""
import json, sys, time
import numpy as np
from scipy.linalg import eigh, expm, logm

np.set_printoptions(precision=12, linewidth=140)

# ---------- Pauli algebra ----------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_all(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def op_on(n, placements):
    """n-qubit operator with single-qubit ops placed at given sites (dict site->op)."""
    return kron_all([placements.get(k, I2) for k in range(n)])


def ptrace_keep_first(rho, dA, dB):
    """Partial trace over the second factor of a dA*dB density matrix."""
    r = rho.reshape(dA, dB, dA, dB)
    return np.einsum('ijkj->ik', r)


def herm(M):
    return 0.5 * (M + M.conj().T)


def trace_norm(M):
    w = np.linalg.eigvalsh(herm(M))
    return float(np.sum(np.abs(w)))


def vn_entropy(rho):
    w = np.linalg.eigvalsh(herm(rho))
    w = w[w > 1e-15]
    return float(-np.sum(w * np.log(w)))


# ---------- Hamiltonian ----------
h = [0.3, -0.7, 0.5, -0.2, 0.8, -0.4, 0.6]  # fields on qubits 1..7


def chain_hamiltonian(n, offset):
    """XXZ-type chain with fields, on sites offset..offset+6 of an n-qubit register."""
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for j in range(6):  # bonds (k, k+1), k = 1..6  -> sites offset+j, offset+j+1
        a, b = offset + j, offset + j + 1
        H += op_on(n, {a: X, b: X}) + op_on(n, {a: Y, b: Y}) + 0.5 * op_on(n, {a: Z, b: Z})
    for j in range(7):
        H += h[j] * op_on(n, {offset + j: Z})
    return H


N = 8
H_A = 0.5 * op_on(N, {0: Z})
H_Abar_8 = chain_hamiltonian(N, 1)       # embedded in 8 qubits
H_Abar_7 = chain_hamiltonian(7, 0)       # 7-qubit chain operator (for the thermal state)
H_int = 0.5 * op_on(N, {0: X, 1: X})
H = H_A + H_Abar_8 + H_int
assert np.allclose(H, H.conj().T)
assert np.allclose(np.kron(I2, H_Abar_7), H_Abar_8)

# ---------- initial states ----------
plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
rho_plus = np.outer(plus, plus.conj())
rho_minus = np.outer(minus, minus.conj())

Eb, Vb = eigh(H_Abar_7)
w = np.exp(-Eb)
Zpart = w.sum()
rho_bath = (Vb * (w / Zpart)) @ Vb.conj().T
rho_bath = herm(rho_bath)
assert abs(np.trace(rho_bath).real - 1) < 1e-12
# cross-check thermal state against expm
assert np.allclose(rho_bath, expm(-H_Abar_7) / np.trace(expm(-H_Abar_7)), atol=1e-12)

rho0_plus = np.kron(rho_plus, rho_bath)
rho0_minus = np.kron(rho_minus, rho_bath)

# ---------- exact propagator via eigendecomposition ----------
E, V = eigh(H)
Vd = V.conj().T


def U_of(t):
    return (V * np.exp(-1j * E * t)) @ Vd


def evolve(rho0, t):
    U = U_of(t)
    return U @ rho0 @ U.conj().T


# spot-check propagator vs expm
for tt in (0.7, 3.0):
    assert np.allclose(U_of(tt), expm(-1j * H * tt), atol=1e-9), tt

ts = np.round(np.arange(0, 40 + 1e-9, 0.1), 10)
assert len(ts) == 401 and abs(ts[-1] - 40) < 1e-12

results = {}

# ========== 1. energy conservation ==========
E_tot0 = float(np.trace(H @ rho0_plus).real)
E_drift = []
S_A = []
D = []
rhoA_plus_list = []
t0 = time.time()
for t in ts:
    rt_p = evolve(rho0_plus, t)
    rt_m = evolve(rho0_minus, t)
    E_drift.append(abs(float(np.trace(H @ rt_p).real) - E_tot0))
    rA_p = ptrace_keep_first(rt_p, 2, 128)
    rA_m = ptrace_keep_first(rt_m, 2, 128)
    rhoA_plus_list.append(rA_p)
    S_A.append(vn_entropy(rA_p))
    D.append(0.5 * trace_norm(rA_p - rA_m))
E_drift = np.array(E_drift)
S_A = np.array(S_A)
D = np.array(D)
print(f"[timing] main grid loop: {time.time()-t0:.1f}s", file=sys.stderr)

results['E_tot'] = E_tot0
results['E_max_drift'] = float(E_drift.max())

# ========== 2. entropy of A ==========
i10 = int(np.argmin(np.abs(ts - 10.0)))
i40 = int(np.argmin(np.abs(ts - 40.0)))
imax = int(np.argmax(S_A))
results['S_A_t10'] = float(S_A[i10])
results['S_A_t40'] = float(S_A[i40])
results['S_A_max'] = float(S_A[imax])
results['S_A_argmax_t'] = float(ts[imax])
results['S_A_t0'] = float(S_A[0])
results['S_A_ln2'] = float(np.log(2))

# ========== 3. BLP backflow ==========
dD = np.diff(D)
results['BLP_sum_positive_increments'] = float(np.sum(dD[dD > 0]))
results['D_min'] = float(D.min())
results['D_argmin_t'] = float(ts[int(np.argmin(D))])
results['D_t0'] = float(D[0])
results['D_t40'] = float(D[-1])
results['BLP_n_positive_increments'] = int(np.sum(dD > 0))

# ========== 4. record protocol (9th qubit R) ==========
t1, t2 = 3.0, 6.0
P0_plus_8 = op_on(N, {0: (I2 + Z) / 2})


def p_up(rho8):
    return float(np.trace(P0_plus_8 @ rho8).real)


# 9-qubit implementation, R = qubit 8 (rightmost factor)
ket0 = np.array([[1, 0], [0, 0]], dtype=complex)
rho9_0 = np.kron(rho0_plus, ket0)
U8_3 = U_of(t1)
U9_3 = np.kron(U8_3, I2)             # H' = H (x) I_R  =>  U' = U (x) I
rho9_t1 = U9_3 @ rho9_0 @ U9_3.conj().T
rho8_t1 = ptrace_keep_first(rho9_t1, 256, 2)
p_t1 = p_up(rho8_t1)

# CNOT: control qubit 0, target qubit 8
P0 = (I2 + Z) / 2
P1 = (I2 - Z) / 2
CNOT = np.kron(np.kron(P0, np.eye(128)), I2) + np.kron(np.kron(P1, np.eye(128)), X)
assert np.allclose(CNOT @ CNOT, np.eye(512))
rho9_after = CNOT @ rho9_t1 @ CNOT.conj().T
rho9_t2 = U9_3 @ rho9_after @ U9_3.conj().T   # evolve another dt = 3
rho8_t2_rec = ptrace_keep_first(rho9_t2, 256, 2)
p_t2_record = p_up(rho8_t2_rec)

# no intervention
rho8_t2_free = evolve(rho0_plus, t2)
p_t2_free = p_up(rho8_t2_free)

# cross-check: record == Z_0 dephasing at t1 in the 8-qubit space
rho8_deph = P0_plus_8 @ rho8_t1 @ P0_plus_8 + (np.eye(256) - P0_plus_8) @ rho8_t1 @ (np.eye(256) - P0_plus_8)
rho8_t2_deph = U8_3 @ rho8_deph @ U8_3.conj().T
p_t2_deph = p_up(rho8_t2_deph)
# the record itself: p(R=0) should equal p(Z_0=+1) at t1, and R-A correlation
rho_R_t2 = np.einsum('ijij->j', rho9_t2.reshape(256, 2, 256, 2))  # wrong axes guard below
rho_R = np.trace(rho9_t2.reshape(256, 2, 256, 2), axis1=0, axis2=2)
p_R0 = float(rho_R[0, 0].real)

results['record_p_up_t1'] = p_t1
results['record_p_up_t2_with_record'] = p_t2_record
results['record_p_up_t2_no_intervention'] = p_t2_free
results['record_p_up_t2_dephasing_crosscheck'] = p_t2_deph
results['record_p_R0_at_t2'] = p_R0
results['record_p_up_t1_from_8q_grid'] = float(rhoA_plus_list[int(np.argmin(np.abs(ts - 3.0)))][0, 0].real)

# ========== 5. initial correlations ==========
eps = 1.0 / 512
rho1 = np.eye(256, dtype=complex) / 256
rho2 = rho1 + eps * op_on(N, {0: Y, 1: X})
assert np.allclose(rho2, rho2.conj().T)
w2 = np.linalg.eigvalsh(rho2)
assert w2.min() > 0, w2.min()
results['rho2_min_eig'] = float(w2.min())
results['rho2_max_eig'] = float(w2.max())
rA1, rA2 = ptrace_keep_first(rho1, 2, 128), ptrace_keep_first(rho2, 2, 128)
# marginal over Abar: trace out qubit 0 (the first factor)
rB1 = np.trace(rho1.reshape(2, 128, 2, 128), axis1=0, axis2=2)
rB2 = np.trace(rho2.reshape(2, 128, 2, 128), axis1=0, axis2=2)
results['corr_marginal_A_diff'] = float(np.abs(rA1 - rA2).max())
results['corr_marginal_Abar_diff'] = float(np.abs(rB1 - rB2).max())
TD = []
for t in ts:
    r1t = evolve(rho1, t) if t == 0 else rho1  # I/256 is stationary; still check once below
    r2t = evolve(rho2, t)
    TD.append(0.5 * trace_norm(ptrace_keep_first(r1t, 2, 128) - ptrace_keep_first(r2t, 2, 128)))
TD = np.array(TD)
# confirm stationarity of rho1 at a random time
assert np.allclose(evolve(rho1, 7.3), rho1)
results['corr_TD_max'] = float(TD.max())
results['corr_TD_argmax_t'] = float(ts[int(np.argmax(TD))])
results['corr_TD_t0'] = float(TD[0])

# ========== 6. Werner marginal compatibility (SDP) ==========
import cvxpy as cp

phi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
Phi = np.outer(phi, phi.conj()).real
I4 = np.eye(4)


def W(p):
    return p * Phi + (1 - p) * I4 / 4


def ptrace_np(rho, dims, keep):
    """Partial trace for a numpy array, dims list, keep = list of subsystem indices."""
    n = len(dims)
    r = rho.reshape(dims + dims)
    trace_axes = [k for k in range(n) if k not in keep]
    for k in sorted(trace_axes, reverse=True):
        r = np.trace(r, axis1=k, axis2=k + r.ndim // 2)
    d = int(np.prod([dims[k] for k in keep]))
    return r.reshape(d, d)


# linear partial-trace maps as matrices acting on vec(rho), to keep cvxpy expressions simple
def ptrace_matrix(dims, keep):
    D = int(np.prod(dims))
    d = int(np.prod([dims[k] for k in keep]))
    M = np.zeros((d * d, D * D))
    for idx in range(D * D):
        e = np.zeros(D * D); e[idx] = 1
        M[:, idx] = ptrace_np(e.reshape(D, D), dims, keep).reshape(-1)
    return M


M_AB = ptrace_matrix([2, 2, 2], [0, 1])   # trace out C
M_BC = ptrace_matrix([2, 2, 2], [1, 2])   # trace out A

# Werner states are real; a real symmetric feasible rho exists iff a Hermitian one does
# (take the real part), so a real SDP is exact. Also solve a complex (Hermitian) version as check.
rho_var = cp.Variable((8, 8), symmetric=True)
p_var = cp.Variable()
cons = [rho_var >> 0, cp.trace(rho_var) == 1,
        cp.reshape(M_AB @ cp.vec(rho_var, order='C'), (4, 4), order='C') == p_var * Phi + (1 - p_var) * I4 / 4,
        cp.reshape(M_BC @ cp.vec(rho_var, order='C'), (4, 4), order='C') == p_var * Phi + (1 - p_var) * I4 / 4]
prob = cp.Problem(cp.Maximize(p_var), cons)
prob.solve(solver=cp.CLARABEL)
p_star = float(p_var.value)
rho_star = np.array(rho_var.value)
# verify the certificate directly
chk_AB = np.abs(ptrace_np(rho_star, [2, 2, 2], [0, 1]) - W(p_star)).max()
chk_BC = np.abs(ptrace_np(rho_star, [2, 2, 2], [1, 2]) - W(p_star)).max()
chk_min_eig = float(np.linalg.eigvalsh(rho_star).min())

# Hermitian (complex) version
rho_c = cp.Variable((8, 8), hermitian=True)
p_c = cp.Variable()
cons_c = [rho_c >> 0, cp.real(cp.trace(rho_c)) == 1,
          cp.reshape(M_AB @ cp.vec(rho_c, order='C'), (4, 4), order='C') == p_c * Phi + (1 - p_c) * I4 / 4,
          cp.reshape(M_BC @ cp.vec(rho_c, order='C'), (4, 4), order='C') == p_c * Phi + (1 - p_c) * I4 / 4]
prob_c = cp.Problem(cp.Maximize(p_c), cons_c)
prob_c.solve(solver=cp.CLARABEL)
p_star_c = float(p_c.value)

# analytic reference: 2-extendibility of the 2-qubit isotropic state, F <= 3/4 <=> p <= 2/3
results['werner_p_max_sdp_real'] = p_star
results['werner_p_max_sdp_hermitian'] = p_star_c
results['werner_sdp_status'] = prob.status
results['werner_certificate_marginal_err_AB'] = float(chk_AB)
results['werner_certificate_marginal_err_BC'] = float(chk_BC)
results['werner_certificate_min_eig'] = chk_min_eig
results['werner_analytic_2_over_3'] = 2.0 / 3.0
# infeasibility check just above threshold: feasibility SDP at p = 2/3 + 0.01
rho_f = cp.Variable((8, 8), symmetric=True)
pf = 2.0 / 3.0 + 0.01
cons_f = [rho_f >> 0, cp.trace(rho_f) == 1,
          cp.reshape(M_AB @ cp.vec(rho_f, order='C'), (4, 4), order='C') == W(pf),
          cp.reshape(M_BC @ cp.vec(rho_f, order='C'), (4, 4), order='C') == W(pf)]
prob_f = cp.Problem(cp.Minimize(0), cons_f)
prob_f.solve(solver=cp.CLARABEL)
results['werner_feasibility_status_at_p_0.6767'] = prob_f.status

# ---------- save ----------
out = {
    'spec': 'closed-whole finite model, independent re-implementation (spec only)',
    'grid': {'t_min': 0.0, 't_max': 40.0, 'dt': 0.1, 'n': int(len(ts))},
    'results': results,
    'series': {
        't': ts.tolist(), 'S_A': S_A.tolist(), 'D': D.tolist(), 'E_drift': E_drift.tolist(), 'TD_corr': TD.tolist(),
    },
}
with open(sys.argv[1], 'w') as f:
    json.dump(out, f, indent=1)

for k, v in results.items():
    print(f"{k:45s} {v!r}")
