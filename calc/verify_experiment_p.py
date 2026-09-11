#!/usr/bin/env python3
"""REALITY_CHECK_04 — independent verification of Experiment P controls.

Re-computes Controls A, D, H from first principles. Does NOT modify any
Experiment P artifact. Read-only with respect to Experiment P files.

Run:  python3.12 calc/verify_experiment_p.py
"""
import json, numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
S2, A2 = 2, 2

def I(d): return np.eye(d)
def op(kron_list):
    m = np.array([[1.0+0j]])
    for x in kron_list: m = np.kron(m, x)
    return m

s0, s1 = np.array([1,0],complex), np.array([0,1],complex)
a_ready = s0

def cyclic_U1(dE):
    """EXPLICITLY UNITARY cyclic permutation on the FULL record space:
    U1|e_j> = |e_{(j+1) mod dE}>, i.e. U1[(j+1) mod dE, j] = 1 for all j.
    A cyclic permutation is a bijection of an orthonormal basis, hence unitary
    by construction. (Repairs the defective partial shift: the original U1
    left e_8..e_14 without images and e_0,e_9..e_14 unoccupied.)"""
    U1 = np.zeros((dE, dE), complex)
    for j in range(dE):
        U1[(j + 1) % dE, j] = 1
    return U1

def run_case(N, env0, alpha, beta):
    """Record-type measurement, H = H_sys(2) ⊗ H_env(2N).

    The apparatus stays in |A_ready> under U (factors out of every term),
    so it is dropped identically. Unitary: |0,e_m> -> |0,e_m>;
    |1,e_m> -> |1,e_{(m+1) mod 2N}> (cyclic shift, REPAIRED)."""
    dE = 2*N
    dim = 2*dE
    P0 = np.kron(np.outer(s0, s0), I(dE))
    P1 = np.kron(np.outer(s1, s1), I(dE))
    U1 = cyclic_U1(dE)
    U = P0 + P1 @ np.kron(I(2), U1)
    psi = alpha*np.kron(s0, env0) + beta*np.kron(s1, env0)
    out = U @ psi
    r = out.reshape(2, dE)
    rho_total = np.outer(r, r.conj())
    rho_S = np.einsum('ia,ka->ik', r, r.conj())
    c = rho_S[0,1]
    return dict(cross_term=complex(c),
                sa_purity=float(np.real(np.trace(rho_S @ rho_S))),
                branch_pop=[float(rho_S[0,0].real), float(rho_S[1,1].real)],
                global_purity=float(np.real(np.trace(rho_total @ rho_total))),
                unitary=float(np.max(np.abs(U.conj().T @ U - I(dim)))),
                unitary_T=float(np.max(np.abs(U @ U.conj().T - I(dim)))),
                E0=cyclic_U1(dE) @ np.eye(dE)[:, 0],
                E1=cyclic_U1(dE) @ env0)

def record_state(N):
    e = np.zeros(2*N, complex); e[0] = 1; return e
def generic_state(N, rng):
    v = rng.normal(size=2*N) + 1j*rng.normal(size=2*N)
    return v/np.linalg.norm(v)

def main():
    rng = np.random.default_rng(20260910)
    res = {}

    # ---- STEP 1/2: unitary proof of repaired U1 and full U -----------------
    dE = 16; dim = 32
    U1 = cyclic_U1(dE)
    e0 = I(dE)[:, 0]
    e1 = I(dE)[:, 1]
    res["U1_check"] = {
        "U1daggerU1_minus_I": float(np.max(np.abs(U1.conj().T @ U1 - I(dE)))),
        "U1U1dagger_minus_I": float(np.max(np.abs(U1 @ U1.conj().T - I(dE)))),
        "rank": int(np.linalg.matrix_rank(U1)),
        "all_column_norms_1": bool(np.allclose(
            np.linalg.norm(U1, axis=0), 1.0)),
        "all_row_norms_1": bool(np.allclose(
            np.linalg.norm(U1, axis=1), 1.0)),
        "permutation": [int(np.argmax(np.abs(U1[:, j]))) for j in range(dE)],
        "permutation_is_bijection": (
            sorted(int(np.argmax(np.abs(U1[:, j]))) for j in range(dE))
            == list(range(dE))),
        "analytic_argument": (
            "U1 maps the orthonormal basis bijectively: e_j -> e_{(j+1) mod dE}. "
            "A bijection of an orthonormal basis is unitary by construction: "
            "<U1 e_a | U1 e_b> = delta_ab, hence U1†U1 = I on the full space."),
        "e0_image": "e_1",
        "U1_e0_dot_U1_e0": complex(U1 @ e0 @ (U1 @ e0).conj()),
    }
    # full U on the full space
    P0f = np.kron(np.outer(s0, s0), I(dE))
    P1f = np.kron(np.outer(s1, s1), I(dE))
    Uf = P0f + P1f @ np.kron(I(2), U1)
    res["U_check"] = {
        "dimension": dim,
        "rank": int(np.linalg.matrix_rank(Uf)),
        "max_abs_UdaggerU_minus_I": float(np.max(np.abs(Uf.conj().T @ Uf - I(dim)))),
        "max_abs_UUdagger_minus_I": float(np.max(np.abs(Uf @ Uf.conj().T - I(dim)))),
    }

    # ---- STEP 3: record orthogonality (direct, from repaired operator) -----
    E0 = U1 @ e0
    E1 = U1 @ e1
    res["record_states"] = {
        "E0": "U1|e0> = |e1>",
        "E1": "U1|e1> = |e2>",
        "overlap": complex(E0 @ E1.conj()),
        "norm_E0": float(np.linalg.norm(E0)),
        "norm_E1": float(np.linalg.norm(E1)),
        "orthogonal": bool(abs(E0 @ E1.conj()) < 1e-12),
    }

    # ---- Control A recheck -------------------------------------------------
    res["A"] = run_case(8, record_state(8), 1/np.sqrt(2), 1/np.sqrt(2))
    res["A"]["analytic_note"] = ("REPAIRED cyclic U1: |1,e0> -> |1,e1>, "
                                 "<e0|e1>=0 -> decoherence exact; unitary on FULL space")

    # ---- Control D recheck (rerun under repaired model; original permanently
    #      INVALID_AS_EXECUTED — quarantined, not evidence) -------------------
    d = run_case(8, generic_state(8, rng), 1/np.sqrt(2), 1/np.sqrt(2))
    res["D_recheck"] = {
        "original_status": "INVALID_AS_EXECUTED (quarantined, preserved as history)",
        "second_state_cross_term": float(abs(d["cross_term"])),
        "identical_reduced_decoherence": abs(d["cross_term"]) < 1e-12,
        "conclusion": ("Under the repaired cyclic U1 a generic environment state still "
                       "does NOT decohere under U1 (generic nonzero diagonal support "
                       "survives a pure shift). The pair is still not an "
                       "identical-reduced-state pair. CONTROL D REMAINS "
                       "INVALID_AS_EXECUTED — supplies no evidence for "
                       "purificational non-identifiability.")
    }

    # ---- Control H recheck: same decoherence class, different weights -------
    h = []
    for a2 in (0.2, 0.5, 0.8):
        r = run_case(8, record_state(8), a2, np.sqrt(1-a2**2))
        h.append({"alpha": a2, "cross": float(abs(r["cross_term"])), "p0": r["branch_pop"][0],
                  "purity": r["sa_purity"]})
    res["H_recheck"] = h
    res["H_scoping"] = {
        "purity_matched": False,
        "note": ("REPAIRED RERUN. Only the pointer-basis decoherence class (cross=0) "
                 "is matched across amplitudes. p0=|alpha|^2 is itself a reduced-state "
                 "population element: the FULL reduced state does contain the weights "
                 "as INHERITED initial-state information. Permitted claim: decoherence "
                 "observables alone do not determine branch weights. Prohibited claim: "
                 "the full reduced state cannot determine the branch weights.")
    }
    # ---- STEP 8: weight-origin trace ---------------------------------------
    res["weight_origin"] = {
        "p0_p1_source": "initial amplitudes alpha,beta (inherited, not dynamically generated)",
        "classification": "WEIGHTS = STATE-DEPENDENT / INHERITED — NOT DYNAMICALLY DERIVED",
        "scope": ("within the repaired tested model class; NOT a universal proof "
                  "against the Born rule")
    }

    out = HERE.parent / "program" / "REALITY_CHECK_04_VERIFICATION.json"
    out.write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=2, default=str))

if __name__ == "__main__":
    main()
