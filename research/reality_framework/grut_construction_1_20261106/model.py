"""Frozen microscopic model, state preparation, and projective readouts."""

import hashlib
import json
from pathlib import Path

import numpy as np


DIRECTORY = Path(__file__).resolve().parent
CARD_PATH = DIRECTORY / "parameters.json"
CARD = json.loads(CARD_PATH.read_text())
CARD_SHA = hashlib.sha256(CARD_PATH.read_bytes()).hexdigest()
SPEC_SHA = hashlib.sha256((DIRECTORY / "SPECIFICATION.md").read_bytes()).hexdigest()
REGISTERS = ("S", "M1", "M2", "A", "D1", "D2")
DIM = 2 ** len(REGISTERS)
LINKS = CARD["clock_links"]
OMEGA = CARD["omega"]
TOL = CARD["absolute_tolerance"]
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
ONE = np.eye(DIM, dtype=complex)
BITS = np.array(list(np.ndindex((2,) * len(REGISTERS))))
S_BIT, M1_BIT, M2_BIT, A_BIT, D1_BIT, D2_BIT = BITS.T
PARITY_BIT = M1_BIT ^ M2_BIT
READOUT_BITS = (S_BIT, PARITY_BIT, A_BIT)


def tensor(factors):
    product = np.ones((1, 1), dtype=complex)
    for factor in factors:
        product = np.kron(product, factor)
    return product


def work_operator(**factors):
    return tensor(factors.get(name, I) for name in REGISTERS)


P = work_operator(M1=Z, M2=Z)
TAU_X = work_operator(M2=X)
TAU_Y = work_operator(M1=Z, M2=Y)
GAUGE_X = work_operator(M1=X, M2=X)
GAUGE_Y = work_operator(M1=Y, M2=X)
GAUGE_Z = work_operator(M1=Z)


def controlled(observable, unitary):
    """Apply unitary in the negative eigenspace of the control."""
    return ((ONE + observable) + (ONE - observable) @ unitary) / 2


def gates(*, no_writer=False, fixed_action=False):
    lam = np.pi * CARD["writer_hold_angle_over_pi"]
    alpha = np.pi * CARD["disturbance_angle_over_pi"]
    beta = np.pi * CARD["parity_phase_angle_over_pi"]
    pulse = np.cos(alpha / 2) * ONE - 1j * np.sin(alpha / 2) * TAU_X
    return (
        ONE if no_writer else controlled(work_operator(S=Z), TAU_X),
        np.diag(np.exp(-0.5j * lam * np.diag(GAUGE_Z))),
        controlled(work_operator(D1=Z), pulse),
        np.diag(np.exp(-0.5j * beta * np.diag(P))),
        controlled(work_operator(D2=Z), pulse),
        ONE if fixed_action else controlled(P, work_operator(A=X)),
        ONE,
    )


def hamiltonian(links):
    """Assemble the actual full Hamiltonian from the off-diagonal blocks."""
    h = np.zeros(((LINKS + 1) * DIM,) * 2, dtype=complex)
    for k, unitary in enumerate(links):
        left = slice(k * DIM, (k + 1) * DIM)
        right = slice((k + 1) * DIM, (k + 2) * DIM)
        g = OMEGA * np.sqrt((k + 1) * (LINKS - k)) / 2
        h[right, left] = g * unitary
        h[left, right] = g * unitary.conj().T
    return h


def program_mask(program):
    return (D1_BIT == program[0]) & (D2_BIT == program[1])


def preparation(program, *, scrambled=False):
    # A factor matrix B represents rho = B B^dagger, here with one column.
    basis = np.eye(2)
    signal = np.array([1, 1]) / np.sqrt(2)
    blank = np.array([1, 0, 0, 1]) / np.sqrt(2)
    work = tensor(x[:, None] for x in (signal, blank, basis[0], *basis[list(program)]))
    if scrambled:
        writer = gates()[0]
        work = writer.conj().T @ ((TAU_X + P) / np.sqrt(2)) @ writer @ work
    return work


def mixed_preparation(program):
    rng = np.random.default_rng(CARD["mixed_reference_seed"])
    shape = (DIM // 4, CARD["mixed_reference_rank"])
    matrix = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    work = np.zeros((DIM, shape[1]), dtype=complex)
    work[program_mask(program)] = matrix / np.linalg.norm(matrix)
    return work


def at_tick(tick, work):
    state = np.zeros((LINKS + 1, DIM, work.shape[1]), dtype=complex)
    state[tick] = work
    return state.reshape(-1, work.shape[1])


def clock_distribution(state):
    blocks = state.reshape(LINKS + 1, DIM, -1)
    return np.sum(abs(blocks) ** 2, axis=(1, 2))


def conditioned_work(state, tick):
    block = state.reshape(LINKS + 1, DIM, -1)[tick]
    return block / np.linalg.norm(block)


def unobserved_read(work, readout_bit):
    """Lüders measurement in the declared basis, with outcome ignored."""
    return np.concatenate([work * (readout_bit == b)[:, None] for b in (0, 1)], axis=1)


def probability(work):
    result = np.zeros((2, 2, 2))
    np.add.at(result, READOUT_BITS, np.sum(abs(work) ** 2, axis=1))
    return result


def expectation(work, operator):
    return float(np.vdot(work, operator @ work).real)


def marginal(work, registers):
    kept = [REGISTERS.index(name) for name in registers]
    axes = kept + [i for i in range(len(REGISTERS) + 1) if i not in kept]
    matrix = work.reshape((2,) * len(REGISTERS) + (-1,)).transpose(axes)
    matrix = matrix.reshape(2 ** len(kept), -1)
    return matrix @ matrix.conj().T
