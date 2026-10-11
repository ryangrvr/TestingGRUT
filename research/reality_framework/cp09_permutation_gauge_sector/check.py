"""Exact permutation-invariant CP-09 Hilbert space and relational response.

Keeps frozen CP-09 unchanged. All labelled N=2..5 causal orders retained
before projection. Independent full-H assembler imported from the already
published five-event geometry control, NOT from the frozen CP-09 author
builder. Requires NumPy and SciPy.

Run from any directory: python check.py
"""
from pathlib import Path
import itertools
import json
import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh, expm_multiply

here = Path(__file__).resolve()
parent = here.parents[1]
source = parent / "cp09_five_event_geometry" / "check.py"
src = source.read_text()
ns = {"__file__": str(source)}
# Load only the independently implemented definitions, without executing
# its full test suite and without rewriting its results.json.
exec(compile(src.split("\nout=[]", 1)[0], str(source), "exec"), ns)
build = ns["build"]


def invariant_orbits(n, states, pairs):
    """Orbit-isometry T from unlabelled geometry+matter orbits into full space."""
    m = 1 << n
    lookup = {v: i for i, v in enumerate(states)}
    pairidx = {v: k for k, v in enumerate(pairs)}
    maps = []
    for a in range(n-1):
        p = list(range(n))
        p[a], p[a+1] = p[a+1], p[a]
        geom = np.empty(len(states), dtype=np.int32)
        for i, mask in enumerate(states):
            new = 0
            for k, (u, v) in enumerate(pairs):
                if mask >> k & 1:
                    new |= 1 << pairidx[p[u], p[v]]
            geom[i] = lookup[new]
        matter = np.empty(m, dtype=np.int32)
        for bits in range(m):
            new = 0
            for i in range(n):
                new |= ((bits >> (n-1-i)) & 1) << (n-1-p[i])
            matter[bits] = new
        maps.append((geom, matter))
    seen = np.zeros(len(states)*m, dtype=bool)
    rows, cols, vals, orbits = [], [], [], []
    for first in range(len(seen)):
        if seen[first]:
            continue
        orbit = [first]
        seen[first] = True
        for old in orbit:
            poset, bits = divmod(old, m)
            for geom, matter in maps:
                dest = int(geom[poset])*m+int(matter[bits])
                if not seen[dest]:
                    seen[dest] = True
                    orbit.append(dest)
        orbit.sort()
        k = len(orbits)
        orbits.append(orbit)
        for item in orbit:
            rows.append(item)
            cols.append(k)
            vals.append(1/math.sqrt(len(orbit)))
    T = csr_matrix((vals, (rows, cols)), shape=(len(seen), len(orbits)))
    assert np.linalg.norm((T.conj().T @ T - csr_matrix(np.eye(len(orbits)))).data) < 1e-11
    return T, orbits, maps


def causal_alignment(n, states, pairs):
    """Predeclared permutation-invariant joint geometry-matter operator.

    G(C,b) = sum_i (number of successors - number of predecessors)_i Z_i
             / [N(N-1)].
    This is scalar under simultaneous permutations, NOT a local detector
    and NOT a spacetime metric.
    """
    m = 1 << n
    bits = np.arange(m)
    z = np.array([1-2*((bits >> (n-1-i)) & 1) for i in range(n)])
    raw = np.zeros(len(states)*m)
    for a, mask in enumerate(states):
        net = np.zeros(n)
        for k, (i, j) in enumerate(pairs):
            if mask >> k & 1:
                net[i] += 1
                net[j] -= 1
        raw[a*m:(a+1)*m] = (net[:, None]*z).sum(axis=0)/(n*(n-1))
    return raw


def analyze(n):
    H, states, pairs = build(n)
    expected_counts = {2:3, 3:19, 4:219, 5:4231}
    assert len(states) == expected_counts[n]
    T, orbits, maps = invariant_orbits(n, states, pairs)
    R = (T.conj().T @ H @ T).tocsr()
    err = float(np.linalg.norm((H@T-T@R).data))
    assert err < 1e-9
    assert np.linalg.norm((R-R.conj().T).data) < 1e-9
    k = min(8, R.shape[0]-2)
    e, v = eigsh(R, k=k, which="SA", tol=1e-11,
                 v0=np.ones(R.shape[0]))
    ix = np.argsort(e)
    e, v = e[ix], v[:, ix]
    residual = max(np.linalg.norm(R@v[:,j]-e[j]*v[:,j])
                   for j in range(k))
    assert residual < 1e-8
    ground, phys_gap = v[:,0], float(e[1]-e[0])
    expected_phys = {
        2: .38819162539638,
        3: .28878172572616,
        4: .20924487505959,
        5: .14941863234445,
    }
    assert abs(phys_gap-expected_phys[n]) < 2e-8
    raw = causal_alignment(n, states, pairs)
    G = np.array([raw[orbit[0]] for orbit in orbits])
    assert all(np.max(np.abs(raw[orbit]-G[j])) < 1e-12
               for j, orbit in enumerate(orbits))
    # Nonzero gauge-invariant spectral coupling and retarded susceptibility.
    mode_weight = float(abs(np.vdot(v[:,1], G*ground))**2)
    t, eps = .4, 1e-5
    evolution = expm_multiply(-1j*R*t, ground)
    response = float(-2*np.vdot(evolution,
        G*expm_multiply(-1j*R*t, G*ground)).imag)
    plus = expm_multiply(-1j*R*t, np.exp(1j*eps*G)*ground)
    minus = expm_multiply(-1j*R*t, np.exp(-1j*eps*G)*ground)
    derivative = float((np.vdot(plus, G*plus).real-
                         np.vdot(minus, G*minus).real)/(2*eps))
    assert abs(derivative-response) < 1e-8
    assert mode_weight > .01 and response > 0
    expected_response = {2:.377611371755, 3:.133426761445,
                         4:.070987431709, 5:.043643906634}
    assert abs(response-expected_response[n]) < 2e-7

    # Full labelled spectrum and orthonormalized FIRST shell, for comparison.
    kc = 5 if n==5 else min(16, H.shape[0]-2)
    ef, vf = eigsh(H, k=kc, which="SA", tol=1e-9,
                   v0=np.ones(H.shape[0]))
    ix = np.argsort(ef)
    ef, vf = ef[ix], vf[:, ix]
    full_gap = float(ef[1]-ef[0])
    expected_full = {2:.318298541788, 3:.185712205054,
                     4:.096503863971, 5:.041881529674}
    assert abs(full_gap-expected_full[n]) < 2e-7
    shell_indices = np.flatnonzero(abs(ef-ef[1]) < 1e-7)
    Q, _ = np.linalg.qr(vf[:,shell_indices])
    # Check shell orthogonality *after* QR (important in degenerate ARPACK).
    assert np.linalg.norm(Q.conj().T@Q - np.eye(len(shell_indices))) < 1e-9
    invariant_overlap = float(np.linalg.norm(T.conj().T@Q))
    assert invariant_overlap < 1e-6
    character = None
    if n in (3,4):
        gmap, bmap = maps[0]
        mm = 1 << n
        destination = (gmap[:,None]*mm+bmap[None,:]).reshape(-1)
        image = np.zeros_like(Q)
        image[destination,:] = Q
        rep = Q.conj().T@image
        character = float(np.trace(rep).real)
        assert abs(character-(n-3)) < 1e-7
        assert np.linalg.norm(image-Q@rep) < 1e-7
    # n=5 sparse k=5 can MISS full degeneracy; invariant exclusion is
    # separately exact from the projected spectrum and T^dagger Q.
    row = dict(N=n, posets=len(states), full_dim=H.shape[0],
        invariant_dim=R.shape[0],
        full_gap=full_gap, invariant_gap=phys_gap,
        initial_shell_modes_computed=len(shell_indices),
        first_shell_invariant_projection_norm=invariant_overlap,
        first_shell_adjacent_swap_character=character,
        scalar_first_invariant_mode_weight=mode_weight,
        scalar_retarded_response_t04=response,
        scalar_independent_impulse_derivative_t04=derivative,
        subspace_invariance_residual=err,
        solver_residual=float(residual))
    print(json.dumps(row, indent=2), flush=True)
    return row


if __name__ == "__main__":
    results = [analyze(n) for n in (2,3,4,5)]
    here.with_name("results.json").write_text(json.dumps(results, indent=2)+"\n")
    print("CP09_PERMUTATION_GAUGE_SECTOR_PASS")
