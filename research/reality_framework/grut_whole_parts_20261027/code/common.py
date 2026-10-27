"""Small standalone helpers. No imports of legacy GRUT numerical engines."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def save(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def kron(*operators):
    out = np.array([1.0 + 0j])
    for operator in operators:
        out = np.kron(out, operator)
    return out


def density(psi):
    return np.outer(psi, psi.conj())


def reduced(rho, dims, keep):
    """Partial trace by explicitly pairing indices in the discarded factors."""
    keep = tuple(keep)
    traced = tuple(i for i in range(len(dims)) if i not in keep)
    perm = keep + traced + tuple(i + len(dims) for i in keep + traced)
    n_keep = int(np.prod([dims[i] for i in keep]))
    n_trace = int(np.prod([dims[i] for i in traced]))
    tensor = rho.reshape(tuple(dims) * 2).transpose(perm)
    return np.einsum('aibi->ab', tensor.reshape(n_keep, n_trace, n_keep, n_trace))


def entropy(rho):
    eigs = np.linalg.eigvalsh(rho)
    eigs = eigs[eigs > 1e-14]
    return float(-np.sum(eigs * np.log2(eigs)))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def expectation(rho, operator):
    return float(np.trace(rho @ operator).real)


def error(value):
    return float(np.max(np.abs(value)))
