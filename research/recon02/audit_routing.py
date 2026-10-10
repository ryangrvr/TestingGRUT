"""Conventional control-closure audit; no candidate law or selector tested."""
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

from audit import E0, H0, I, P0, P1, Y, Z, distinguishability, reduce_to

ROOT = Path(__file__).resolve().parent
G = np.zeros((4, 4), dtype=complex)
G[1, 2], G[2, 1] = 1j, -1j
PERTURBATION = np.kron(Z, Y)


def router(phi):
    c, s = math.cos(phi), math.sin(phi)
    return np.array([[1, 0, 0, 0], [0, c, s, 0],
                     [0, -s, c, 0], [0, 0, 0, 1]], dtype=complex)


def canonicalizer(e0, e1, rank_one):
    """Gram-preserving span isometry extended to one 4x4 unitary."""
    q = np.vdot(e1, e0)
    vectors = [e0/np.linalg.norm(e0)]
    if not rank_one:
        residual = e1 - q.conjugate()*e0
        vectors.append(residual/np.linalg.norm(residual))
    for e in np.eye(4, dtype=complex):
        residual = e.copy()
        for f in vectors:
            residual -= f*np.vdot(f, residual)
        norm = np.linalg.norm(residual)
        if norm > 1e-10:
            vectors.append(residual/norm)
        if len(vectors) == 4:
            break
    source = np.column_stack(vectors)
    target = np.eye(4, dtype=complex)[:, [0, 2, 1, 3]]
    return target@source.conj().T


def channel(unitary):
    out = []
    for i in range(2):
        for j in range(2):
            element = np.zeros((2, 2), dtype=complex)
            element[i, j] = 1
            total = unitary@np.kron(element, E0)@unitary.conj().T
            out.append(reduce_to(total, [0], [2, 2, 2]))
    return out


def local_d(e0, e1):
    states = [np.outer(e, e.conj()) for e in [e0, e1]]
    return distinguishability(*(reduce_to(r, [0], [2, 2]) for r in states))


def run():
    maximum = {key: 0. for key in ["router_exponential", "unitarity",
        "bare_energy_commutator", "excitation_commutator", "canonical_pair",
        "channel_invariance", "D_formula", "inversion_attainment",
        "perturbed_channel_invariance", "error_bound_violation"]}

    def error(key, value):
        maximum[key] = max(maximum[key], float(value))

    rational = []
    attainments = []
    for v in [F(0), F(3, 5), F(4, 5), F(1)]:
        s2 = 1-v*v
        for x in [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]:
            d2 = s2*x*(v*v+s2*x)
            assert 0 <= d2 <= s2
            assert (x != 0 or d2 == 0) and (x != 1 or d2 == s2)
            if s2:
                # Exact identity for the positive root used in the proof.
                assert (v*v+2*s2*x)**2 == v**4+4*d2
            rational.append({"v": str(v), "x": str(x), "D_squared": str(d2)})
        for portion in [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]:
            desired = float(portion)*math.sqrt(float(s2))
            if s2:
                x = (math.sqrt(float(v**4)+4*desired**2)-float(v*v))/(2*float(s2))
                got = math.sqrt(float(s2)*x*(float(v*v)+float(s2)*x))
            else:
                x, got = 0., 0.
            error("inversion_attainment", abs(got-desired))
            attainments.append({"v": str(v), "fraction_of_max_D": str(portion),
                                "desired_D": desired, "x": x, "attained_D": got})

    rng = np.random.default_rng(20261010)
    routes, calibrations = [], []
    excitation = np.kron(np.diag([0., 1.]), I)+np.kron(I, np.diag([0., 1.]))
    for v in [F(0), F(3, 5), F(1)]:
        for phase in [0., math.pi/3]:
            raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
            mixing, _ = np.linalg.qr(raw)
            rotation = np.kron(expm(-1j*math.acos(float(v))*Y), I)
            u0, u1 = mixing, np.exp(1j*phase)*mixing@rotation
            source_u = np.kron(P0, u0)+np.kron(P1, u1)
            e0, e1 = u0[:, 0], u1[:, 0]
            q = np.vdot(e1, e0)
            w = canonicalizer(e0, e1, rank_one=(v == 1))
            target0 = np.array([1, 0, 0, 0], dtype=complex)
            target1 = np.array([q.conjugate(), 0, math.sqrt(float(1-v*v)), 0])
            error("canonical_pair", max(np.max(np.abs(w@e0-target0)),
                                          np.max(np.abs(w@e1-target1))))
            original_channel = channel(source_u)
            for phi in [0., math.pi/6, math.pi/4, math.pi/3, math.pi/2]:
                b = router(phi)
                error("router_exponential", np.max(np.abs(b-expm(-1j*phi*G))))
                error("excitation_commutator", np.max(np.abs(b@excitation-excitation@b)))
                combined = b@w
                total_u = np.kron(I, combined)@source_u
                error("unitarity", np.max(np.abs(total_u.conj().T@total_u-np.eye(8))))
                error("bare_energy_commutator", np.max(np.abs(total_u@H0-H0@total_u)))
                reduced = channel(total_u)
                error("channel_invariance", max(np.max(np.abs(a-z))
                                                for a, z in zip(reduced, original_channel)))
                d = local_d(combined@e0, combined@e1)
                s2, x = float(1-v*v), math.cos(phi)**2
                expected = math.sqrt(s2*x*(float(v*v)+s2*x))
                error("D_formula", abs(d-expected))
                routes.append({"v": str(v), "seed_relative_phase": phase, "phi": phi,
                               "q_real": float(q.real), "q_imag": float(q.imag),
                               "D_F_matrix": d, "D_F_formula": expected,
                               "two_qubit_pulse_slots": 4,
                               "additional_router_norm_hbar_over_T": 16*phi})
                for epsilon in [F(1, 50), F(1, 10)]:
                    perturbed = expm(-1j*float(epsilon)*PERTURBATION)@combined
                    eta = float(np.linalg.norm(perturbed-combined, ord=2))
                    d_tilde = local_d(perturbed@e0, perturbed@e1)
                    delta = abs(d_tilde-d)
                    error("error_bound_violation", max(0., delta-2*eta))
                    perturbed_channel = channel(np.kron(I, perturbed)@source_u)
                    error("perturbed_channel_invariance", max(np.max(np.abs(a-z))
                        for a, z in zip(perturbed_channel, original_channel)))
                    calibrations.append({"v": str(v), "seed_relative_phase": phase,
                        "phi": phi, "epsilon": str(epsilon), "operator_error_eta": eta,
                        "D_F_change": delta, "certified_algebraic_bound_2_eta": 2*eta})

    assert all(e < 1e-12 for e in maximum.values())
    assert np.max(np.abs(G-G.conj().T)) == 0
    assert abs(np.linalg.norm(G, ord=2)-1) < 1e-13
    ledger = hashlib.sha256((ROOT.parent/'card01'/'ATTEMPT_LEDGER.json').read_bytes()).hexdigest()
    assert ledger == 'ef4fc46ec5ae752ec3defd8aba9c8d6521dd23298b61ca966008e64ea064c5fb'
    result = {"status": "CONVENTIONAL_ROUTING_CLOSURE_CONTROLS_PASS_NO_GRUT_LAW",
        "scope": "Conditional control-closure theorem; no universal no-go or generating K",
        "new_candidate_formulations": 0, "new_candidate_evaluations": 0,
        "frozen_slots_consumed_total": 1, "remaining_slots": 2,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "counts": {"exact_rational_controls": len(rational), "inversion_controls": len(attainments),
                   "routing_fixtures": len(routes), "routing_channel_matrix_units": 4*len(routes),
                   "perturbation_fixtures": len(calibrations),
                   "perturbation_channel_matrix_units": 4*len(calibrations)},
        "maximum_errors": maximum, "rational_controls": rational, "inversion_controls": attainments,
        "routing_fixtures": routes, "unitary_error_fixtures": calibrations,
        "Card1_ledger_sha256": ledger,
        "limits": "Synthetic finite controls. Closure, spare resources and arbitrary neighbouring two-qubit pulse capability are explicit hypotheses; no real calibration, certified numerical enclosures, external review or GRUT selector."}
    (ROOT/'routing_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    names = ['OPERATION_CLOSURE.md', 'audit_routing.py', 'reproduce_routing.sh',
             'ROUTING_SOURCES.json', 'routing_results.json', 'audit.py', 'requirements.txt']
    (ROOT/'ROUTING_MANIFEST.json').write_text(json.dumps({"sha256": {
        p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names}}, indent=2)+'\n')
    print(result['status'])
    print(json.dumps({"counts": result['counts'], "maximum_errors": maximum}, indent=2))


if __name__ == '__main__':
    run()
