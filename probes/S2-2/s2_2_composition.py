"""SCOUT-2 S2-2: which composition rule is forced? (Local tomography NOT used as an input.)

Criteria tested:
  (i)  independent preparability (product states exist)       -> holds in every composite: selects nothing
  (ii) existence of a CONTINUOUS reversible interaction        -> computed below
  (iii) copying / Cartesian structure                           -> classical simplex only

A. Two gbits (square state spaces):
   max tensor = no-signalling polytope (16 local deterministic + 8 PR vertices, 16 positivity facets);
   min tensor = local polytope (16 vertices; 16 positivity + 8 CHSH facets).
   Combinatorial automorphism group of each polytope (vertex-facet incidence graph). Compare with
   local relabellings x swap (order 8*8*2 = 128). Does any automorphism map a product vertex to a non-product one?
B. Classical bits: composite = simplex on 4 points; automorphisms = S_4 (finite): interactions exist (e.g. CNOT)
   but none is continuous (Lie algebra dimension 0).
C. Qubits (complex) and rebits (real): continuous entangling one-parameter groups exist:
   complex exp(-i t Z x Z); real exp(t G), G = i Y x X (real antisymmetric). Entanglement generated from a product state.
"""
import itertools
import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism

np.set_printoptions(precision=4, suppress=True)
IDX = {(a, b, x, y): 8 * x + 4 * y + 2 * a + b for a in (0, 1) for b in (0, 1) for x in (0, 1) for y in (0, 1)}


def box(fn):
    v = np.zeros(16)
    for (a, b, x, y), i in IDX.items():
        v[i] = fn(a, b, x, y)
    return v


def local_deterministic():
    out = []
    for fa in itertools.product((0, 1), repeat=2):        # a = fa[x]
        for fb in itertools.product((0, 1), repeat=2):    # b = fb[y]
            out.append(box(lambda a, b, x, y: float(a == fa[x] and b == fb[y])))
    return out


def pr_boxes():
    out = []
    for al, be, ga in itertools.product((0, 1), repeat=3):
        out.append(box(lambda a, b, x, y: 0.5 * float((a ^ b) == ((x & y) ^ (al & x) ^ (be & y) ^ ga))))
    return out


def chsh_facets():
    fs = []
    for al, be, ga in itertools.product((0, 1), repeat=3):
        # S = sum_{xy} s_xy E_xy <= 2,  E_xy = sum_ab (-1)^{a^b} P(ab|xy),  s_xy = (-1)^{xy ^ al x ^ be y ^ ga}
        w = box(lambda a, b, x, y: (-1) ** ((x & y) ^ (al & x) ^ (be & y) ^ ga) * (-1) ** (a ^ b))
        fs.append(w)
    return fs


def incidence_graph(vertices, facets):
    """facets: list of (vector w, bound c) meaning w.v <= c tight when w.v == c"""
    G = nx.Graph()
    for i, v in enumerate(vertices):
        G.add_node(("v", i), kind="v")
    for j, (w, c) in enumerate(facets):
        G.add_node(("f", j), kind="f")
        for i, v in enumerate(vertices):
            if abs(w @ v - c) < 1e-9:
                G.add_edge(("v", i), ("f", j))
    return G


def automorphisms(G):
    gm = isomorphism.GraphMatcher(G, G, node_match=lambda a, b: a["kind"] == b["kind"])
    return list(gm.isomorphisms_iter())


if __name__ == "__main__":
    L = local_deterministic(); PR = pr_boxes()
    pos = [(-np.eye(16)[k], 0.0) for k in range(16)]                 # -P_k <= 0
    print("=== A. two gbits ===")
    # facet-incidence counts are invariants of every (combinatorial or linear) polytope automorphism
    V_ns = L + PR
    G_ns = incidence_graph(V_ns, pos)
    deg_L = sorted({G_ns.degree(("v", i)) for i in range(16)})
    deg_PR = sorted({G_ns.degree(("v", i)) for i in range(16, 24)})
    print(f"  max tensor (NS polytope): {len(V_ns)} vertices, {len(pos)} facets")
    print(f"     facets through a product vertex: {deg_L}; through a PR vertex: {deg_PR}")
    print("     -> no automorphism maps a product vertex to a PR vertex (incidence count is invariant): no entangling map")
    fac_loc = pos + [(w, 2.0) for w in chsh_facets()]
    G_loc = incidence_graph(L, fac_loc)
    print(f"  min tensor (local polytope): {len(L)} vertices (all product), {len(fac_loc)} facets;"
          f" facets per vertex: {sorted({G_loc.degree(('v', i)) for i in range(16)})}")
    print("     every vertex is a product state, so every automorphism permutes product states")
    print("  any polytope has a FINITE automorphism group (a subgroup of the vertex permutations):")
    print("  no continuous reversible dynamics, hence no continuous interaction, in either gbit composite")
    print("\n=== B. classical bits ===")
    print("  composite = simplex on 4 pure states; reversible maps = S_4 (24 permutations, includes CNOT);")
    print("  Lie algebra of the automorphism group = 0 -> no CONTINUOUS reversible interaction")

    print("\n=== C. continuous entangling one-parameter groups ===")
    X = np.array([[0, 1], [1, 0]]); Z = np.diag([1., -1]); Y = np.array([[0, -1j], [1j, 0]])
    plus = np.array([1, 1]) / np.sqrt(2)
    psi0 = np.kron(plus, plus)
    from scipy.linalg import expm
    zero2 = np.array([1., 0, 0, 0])
    # NOTE: |++> is a bad probe for G = iY x X (it maps to a product state); a first run used it and reported 0.
    for name, Gen, real, psi_init in (("complex qubits, exp(-i t ZZ), from |++>", -1j * np.kron(Z, Z), False, psi0),
                                      ("real rebits,   exp(t iYX),  from |00> ", np.real(1j * np.kron(Y, X)), True, zero2)):
        row = []
        for t in (0.0, 0.2, 0.4, np.pi / 4):
            U = expm(t * Gen)
            if real:
                assert np.allclose(U.imag if np.iscomplexobj(U) else 0, 0) and np.allclose(U @ U.T, np.eye(4))
            psi = U @ psi_init
            s = np.linalg.svd(psi.reshape(2, 2), compute_uv=False)
            row.append(2 * s[0] * s[1])          # concurrence of a pure two-qubit state = 2 s0 s1
        print(f"  {name}: entanglement (2 s0 s1) at t = 0, 0.2, 0.4, pi/4: {np.round(row, 4)}")
    print("  Lie algebras: classical 0; gbit (square) 0; rebit (disc) so(2) dim 1; qubit (ball) so(3) dim 3;")
    print("  composites: so(4) (dim 6) > local so(2)+so(2) (dim 2); su(4) (dim 15) > local su(2)+su(2) (dim 6)")
