"""Exact parity-block solver and population-only comparator for the readout."""

import numpy as np

import model as m


def parity_blocks(work, program):
    """Trace gauge and work factor columns in the (g, parity=g XOR M2) basis."""
    state = work.reshape((2,) * len(m.REGISTERS) + (-1,))
    weights = np.zeros((2, 2))
    bloch = np.zeros((2, 2, 3))
    for s, a in np.ndindex(2, 2):
        amplitudes = np.array([
            np.concatenate([state[s, g, g ^ p, a, *program] for g in (0, 1)])
            for p in (0, 1)
        ])
        block = amplitudes @ amplitudes.conj().T
        weights[s, a] = np.trace(block).real
        bloch[s, a] = (2 * block[0, 1].real, -2 * block[0, 1].imag,
                       (block[0, 0] - block[1, 1]).real)
    return weights, bloch


def joint(weights, bloch):
    p_plus = (weights + bloch[..., 2]) / 2
    p_minus = (weights - bloch[..., 2]) / 2
    return np.stack((p_plus, p_minus), axis=1)


def act(probability, fixed_action):
    result = probability.copy()
    if not fixed_action:
        result[:, 1, :] = probability[:, 1, ::-1]
    return result


def rotate(vectors, first, second, angle):
    cosine, sine = np.cos(angle), np.sin(angle)
    a, b = vectors[..., first].copy(), vectors[..., second].copy()
    vectors[..., first] = cosine * a - sine * b
    vectors[..., second] = sine * a + cosine * b


def solve(work, program, *, no_writer=False, fixed_action=False):
    weights, bloch = parity_blocks(work, program)
    alpha = np.pi * m.CARD["disturbance_angle_over_pi"]
    beta = np.pi * m.CARD["parity_phase_angle_over_pi"]
    ticks = [joint(weights, bloch)]
    if not no_writer:
        bloch[1, :, 1:] *= -1
    ticks.extend((joint(weights, bloch), joint(weights, bloch)))
    rotate(bloch, 1, 2, alpha * program[0])
    ticks.append(joint(weights, bloch))
    rotate(bloch, 0, 1, beta)
    ticks.append(joint(weights, bloch))
    rotate(bloch, 1, 2, alpha * program[1])
    ticks.append(joint(weights, bloch))
    post_action = act(ticks[-1], fixed_action)
    ticks.extend((post_action, post_action.copy()))
    return np.array(ticks)


def populations(work, program):
    """Drop parity coherences and apply each enabled pulse as a stochastic flip."""
    initial = m.probability(work)
    after_write = initial.copy()
    after_write[1] = initial[1, ::-1]
    flip = np.sin(np.pi * m.CARD["disturbance_angle_over_pi"] / 2) ** 2
    def pulse(probability, enabled):
        q = flip * enabled
        return (1 - q) * probability + q * probability[:, ::-1]

    first = pulse(after_write, program[0])
    second = pulse(first, program[1])
    post_action = act(second, False)
    return np.array((initial, after_write, after_write, first, first, second,
                     post_action, post_action))
