"""Known controlled-unitary / NULL-PiT comparison; no candidate K evaluation."""
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
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1., -1.]).astype(complex)
P0 = np.diag([1., 0.]).astype(complex)
P1 = I - P0
E0 = np.diag([1., 0., 0., 0.]).astype(complex)
H0 = np.kron(Z, np.eye(4)) / 2  # hbar=omega=T=1
SWAP = np.eye(4, dtype=complex)[[0, 2, 1, 3]]
HAD = np.array([[1., 1.], [1., -1.]]) / math.sqrt(2)
PHASE = np.diag([1., 1j])
CNOT = np.eye(4, dtype=complex)[[0, 1, 3, 2]]
BELL_ENCODING = CNOT @ np.kron(HAD, I) @ np.kron(PHASE, I)


def reduce_to(x, keep, dims):
    """Partial trace using tensor axes; works also on non-Hermitian units."""
    out = np.asarray(x).reshape(tuple(dims) * 2)
    current = len(dims)
    for index in reversed(range(len(dims))):
        if index not in keep:
            out = np.trace(out, axis1=index, axis2=index + current)
            current -= 1
    size = math.prod(dims[i] for i in keep)
    return out.reshape(size, size)


def distinguishability(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum() / 2)


def interval_decision(v, d):
    """Exact supplied-box decision; no physical box certification claimed."""
    if (not all(isinstance(z, (tuple, list)) and len(z) == 2
                and all(isinstance(t, F) for t in z) for z in [v, d])
            or any(z[0] < 0 or z[0] > z[1] or z[1] > 1 for z in [v, d])):
        return {"status": "INVALID_INPUT"}
    lower = v[0]**2 + d[0]**2
    upper = v[1]**2 + d[1]**2
    return {"status": "EXCLUDES_STATED_QUANTUM_BASELINE" if lower > 1
            else "WITHIN_OR_UNRESOLVED", "sum_of_squares_lower": str(lower),
            "sum_of_squares_upper": str(upper)}


def run():
    maximum = {key: 0. for key in [
        "unitarity", "Hamiltonian_exponential", "nearest_neighbour_synthesis",
        "bare_energy_commutator", "matrix_unit_channel", "conditional_D",
        "pure_global_D", "Helstrom_guessing", "Bell_gate_synthesis",
        "encoded_channel_invariance", "encoded_single_fragment_D"]}
    minimum_choi = 0.
    channel_cells = 0
    records = []
    fractions = [F(1), F(4, 5), F(3, 5), F(0)]
    fixtures = [(f"grid_{i}_{j}", a, b) for i, a in enumerate(fractions)
                for j, b in enumerate(fractions)]
    fixtures += [("fixed_V_no_F_record", F(1), F(3, 5)),
                 ("fixed_V_partial_F_record", F(4, 5), F(3, 4)),
                 ("fixed_V_complete_F_record", F(3, 5), F(1))]

    def error(name, value):
        maximum[name] = max(maximum[name], float(value))

    for name, c_f, c_r in fixtures:
        theta_f, theta_r = math.acos(float(c_f)), math.acos(float(c_r))
        r_f, r_r = expm(-1j * theta_f * Y), expm(-1j * theta_r * Y)
        u1 = np.kron(r_f, r_r)
        u = np.kron(P0, np.eye(4)) + np.kron(P1, u1)
        h = np.kron(P1, theta_f*np.kron(Y, I) + theta_r*np.kron(I, Y))
        error("Hamiltonian_exponential", np.max(np.abs(expm(-1j*h)-u)))
        error("unitarity", np.max(np.abs(u.conj().T@u-np.eye(8))))
        error("bare_energy_commutator", np.max(np.abs(h@H0-H0@h)))
        # C_SF then SWAP_FR, C_SF(theta_R), SWAP_FR; four nearest-neighbour gates.
        def controlled_on_f(rotation):
            return np.kron(np.kron(P0, I) + np.kron(P1, rotation), I)
        swap_fr = np.kron(I, SWAP)
        chain = swap_fr @ controlled_on_f(r_r) @ swap_fr @ controlled_on_f(r_f)
        error("nearest_neighbour_synthesis", np.max(np.abs(chain-u)))
        q, d2 = c_f*c_r, 1-c_f**2
        slack = c_f**2 * (1-c_r**2)
        assert q*q + d2 + slack == 1  # exact rational algebra, including endpoints

        blocks = {}
        for i in range(2):
            for j in range(2):
                unit = np.zeros((2, 2), dtype=complex)
                unit[i, j] = 1
                total = u @ np.kron(unit, E0) @ u.conj().T
                reduced = reduce_to(total, [0], [2, 2, 2])
                target = unit * (1 if i == j else float(q))
                error("matrix_unit_channel", np.max(np.abs(reduced-target)))
                blocks[i, j] = reduced
                channel_cells += 1
        choi = np.block([[blocks[0, 0], blocks[0, 1]],
                         [blocks[1, 0], blocks[1, 1]]]) / 2
        minimum_choi = min(minimum_choi, float(np.linalg.eigvalsh(choi).min()))
        assert abs(np.trace(choi)-1) < 1e-13
        conditional = [E0, u1@E0@u1.conj().T]
        fs = [reduce_to(x, [0], [2, 2]) for x in conditional]
        d_f = distinguishability(*fs)
        d_fr = distinguishability(*conditional)
        error("conditional_D", abs(d_f-math.sqrt(float(d2))))
        error("pure_global_D", abs(d_fr-math.sqrt(float(1-q*q))))
        eigenvalues, eigenvectors = np.linalg.eigh(fs[0]-fs[1])
        positive = eigenvectors[:, eigenvalues > 1e-13]
        effect = positive@positive.conj().T
        guessing = (np.trace(effect@fs[0]) + np.trace((I-effect)@fs[1])).real/2
        error("Helstrom_guessing", abs(guessing-(1+d_f)/2))
        records.append({"name": name, "cos_theta_F": str(c_f),
                        "cos_theta_R": str(c_r), "V_exact": str(q),
                        "D_F_squared_exact": str(d2), "slack_exact": str(slack),
                        "D_F_matrix": d_f, "D_FR_matrix": d_fr,
                        "Helstrom_guessing_matrix": float(guessing),
                        "normalized_Choi_min_eigenvalue": float(np.linalg.eigvalsh(choi).min())})

    # Hostile control: a single fixed environment-only encoding, no changed F readout.
    phi_plus = np.array([1, 0, 0, 1]) / math.sqrt(2)
    phi_minus = np.array([1, 0, 0, -1]) / math.sqrt(2)
    psi_plus = np.array([0, 1, 1, 0]) / math.sqrt(2)
    psi_minus = np.array([0, 1, -1, 0]) / math.sqrt(2)
    specified_w = np.column_stack([phi_plus, psi_plus, 1j*phi_minus, 1j*psi_minus])
    error("Bell_gate_synthesis", np.max(np.abs(BELL_ENCODING-specified_w)))
    error("unitarity", np.max(np.abs(BELL_ENCODING.conj().T@BELL_ENCODING-np.eye(4))))
    encoded = []
    for q in [F(0), F(3, 5), F(1)]:
        u1 = np.kron(expm(-1j*math.acos(float(q))*Y), I)
        u = np.kron(P0, np.eye(4)) + np.kron(P1, u1)
        transformed = np.kron(I, BELL_ENCODING)@u
        branch = [BELL_ENCODING@E0@BELL_ENCODING.conj().T,
                  BELL_ENCODING@u1@E0@u1.conj().T@BELL_ENCODING.conj().T]
        d_f = distinguishability(*(reduce_to(x, [0], [2, 2]) for x in branch))
        d_r = distinguishability(*(reduce_to(x, [1], [2, 2]) for x in branch))
        d_fr = distinguishability(*branch)
        error("encoded_single_fragment_D", max(d_f, d_r))
        error("pure_global_D", abs(d_fr-math.sqrt(float(1-q*q))))
        error("bare_energy_commutator", np.max(np.abs(transformed@H0-H0@transformed)))
        for i in range(2):
            for j in range(2):
                unit = np.zeros((2, 2), dtype=complex)
                unit[i, j] = 1
                initial = np.kron(unit, E0)
                before = reduce_to(u@initial@u.conj().T, [0], [2, 2, 2])
                after = reduce_to(transformed@initial@transformed.conj().T, [0], [2, 2, 2])
                error("encoded_channel_invariance", np.max(np.abs(after-before)))
                channel_cells += 1
        encoded.append({"V_exact": str(q), "D_F": d_f, "D_R": d_r,
                        "D_FR": d_fr, "D_FR_squared_exact": str(1-q*q)})

    decisions = [
        interval_decision((F(79, 100), F(81, 100)), (F(69, 100), F(71, 100))),
        interval_decision((F(59, 100), F(61, 100)), (F(79, 100), F(81, 100))),
        interval_decision((F(59, 100), F(61, 100)), (F(0), F(1, 100))),
        interval_decision((F(-1, 100), F(1, 100)), (F(0), F(1)))
    ]
    assert [x["status"] for x in decisions] == ["EXCLUDES_STATED_QUANTUM_BASELINE",
        "WITHIN_OR_UNRESOLVED", "WITHIN_OR_UNRESOLVED", "INVALID_INPUT"]
    assert all(x < 1e-12 for x in maximum.values()) and minimum_choi > -1e-12
    ledger = hashlib.sha256((ROOT.parent/'card01'/'ATTEMPT_LEDGER.json').read_bytes()).hexdigest()
    assert ledger == 'ef4fc46ec5ae752ec3defd8aba9c8d6521dd23298b61ca966008e64ea064c5fb'
    result = {
        "status": "CONVENTIONAL_JOINT_IMAGE_AUDIT_PASS_NO_GRUT_SELECTOR",
        "scope": "Known controlled-unitary / NULL-PiT baseline; no generative candidate",
        "new_GRUT_candidate_formulations": 0, "new_GRUT_candidate_evaluations": 0,
        "frozen_slots_consumed_total": 1, "remaining_slots": 2,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "counts": {"product_fixtures": len(fixtures), "encoded_fixtures": len(encoded),
                   "matrix_unit_channel_cells": channel_cells, "interval_fixtures": len(decisions)},
        "maximum_errors": maximum, "minimum_normalized_Choi_eigenvalue": minimum_choi,
        "product_fixtures": records, "encoded_fixtures": encoded, "interval_decisions": decisions,
        "Card1_ledger_sha256": ledger,
        "limitations": "Synthetic finite controls; floating tolerances are not rigorous enclosures; no experimental lock, new GRUT restriction, external review or QFT/KMS completion."
    }
    (ROOT/'results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    names = ['audit.py', 'REPORT.md', 'SOURCES.json', 'SCOPE.md', 'reproduce.sh', 'requirements.txt', 'results.json']
    hashes = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    (ROOT/'MANIFEST.json').write_text(json.dumps({"sha256": hashes}, indent=2)+'\n')
    print(result['status'])
    print(json.dumps({"counts": result['counts'], "maximum_errors": maximum,
                      "minimum_normalized_Choi_eigenvalue": minimum_choi}, indent=2))


if __name__ == '__main__':
    run()
