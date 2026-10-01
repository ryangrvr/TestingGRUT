"""SCOUT-1 W1-G: do Gleason (d>=3), Busch (POVMs) or composition + no-signalling force h(p) = p?

Test family: outcome weight w = h(p), p = |<e|psi>|^2, with h(1-p) = 1 - h(p) (so d = 2 projective
normalization holds automatically):   h_eps(p) = p + eps * sin(2 pi p) / (2 pi).
 1. d = 2 projective: every h_eps is a valid frame function (normalized on every orthonormal basis).
 2. d = 3 projective (Gleason): normalization fails for eps != 0 on random bases; the functional equation
    h(x) + h(y) + h(1-x-y) = 1 forces h linear.
 3. d = 2 POVMs (Busch): trine/tetrahedral POVMs break normalization for eps != 0.
 4. composition + no-signalling: Bell pair, joint weights ~ h(p_ab) normalized; Bob's marginal depends on
    Alice's basis choice iff h nonlinear.
 5. premise audit: noncontextuality = the frame function depends only on the projector.
"""
import numpy as np

rng = np.random.default_rng(11)


def h(p, eps):
    return p + eps * np.sin(2 * np.pi * p) / (2 * np.pi)


def haar(d):
    z = (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def rand_state(d):
    v = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    return v / np.linalg.norm(v)


eps_list = (0.0, 0.05, 0.2, 0.5)
print("h_eps monotone on [0,1] for |eps| <= 1:", all(np.all(np.diff(h(np.linspace(0, 1, 1001), e)) >= 0) for e in eps_list))

print("\n=== 1. d = 2 projective: sum_i h(p_i) over random orthonormal bases ===")
for eps in eps_list:
    dev = max(abs(sum(h(abs(np.vdot(U[:, i], psi)) ** 2, eps) for i in range(2)) - 1)
              for U, psi in ((haar(2), rand_state(2)) for _ in range(2000)))
    print(f"  eps={eps:4.2f}: max |sum - 1| = {dev:.2e}   -> valid non-Born frame function" if eps else
          f"  eps={eps:4.2f}: max |sum - 1| = {dev:.2e}   (Born)")

print("\n=== 2. d = 3 projective (Gleason regime) ===")
for eps in eps_list:
    devs = [abs(sum(h(abs(np.vdot(U[:, i], psi)) ** 2, eps) for i in range(3)) - 1)
            for U, psi in ((haar(3), rand_state(3)) for _ in range(5000))]
    print(f"  eps={eps:4.2f}: max |sum - 1| = {max(devs):.3e}, median {np.median(devs):.3e}")
print("  functional equation: h(x)+h(y)+h(1-x-y)=1 on the simplex, with (1,0,0): h(1)+2h(0)=1, d=2 embedding h(0)=0")
print("  => h(x+y) = h(x)+h(y) on [0,1] => (monotone) h(p) = p. Check residual of the equation for h_eps at random (x,y):")
xs = rng.dirichlet([1, 1, 1], 5)
for eps in eps_list[1:]:
    print(f"   eps={eps}: residuals {np.round([h(a, eps) + h(b, eps) + h(c, eps) - 1 for a, b, c in xs], 4)}")

print("\n=== 3. d = 2 POVMs (Busch regime): effects E_k = (2/n)|phi_k><phi_k|, value v(E) = h(Tr rho E) ===")
def trine(n):
    out = []
    for k in range(n):
        th = 2 * np.pi * k / n
        phi = np.array([np.cos(th / 2), np.sin(th / 2)])
        out.append((2 / n) * np.outer(phi, phi))
    return out
tet_dirs = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
tet = [0.25 * (np.eye(2) + n[0] * sx + n[1] * sy + n[2] * sz) for n in tet_dirs]
for name, povm in (("trine", trine(3)), ("5-star", trine(5)), ("tetrahedron", tet)):
    assert np.allclose(sum(povm), np.eye(2))
    for eps in eps_list:
        devs = []
        for _ in range(2000):
            psi = rand_state(2); rho = np.outer(psi, psi.conj())
            devs.append(abs(sum(h(np.real(np.trace(rho @ E)), eps) for E in povm) - 1))
        if eps in (0.0, 0.2):
            print(f"  {name:11s} eps={eps:4.2f}: max |sum - 1| = {max(devs):.3e}")
print("  Busch: additivity over all effect decompositions => v linear in E => v(E) = Tr(rho E), already for d = 2")

print("\n=== 4. composition + no-signalling: Bell pair, joint weights h(p_ab), normalized ===")
bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
non_max = np.array([np.cos(0.4), 0, 0, np.sin(0.4)])
def basis(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
for name, psi in (("Bell", bell), ("cos/sin 0.4", non_max)):
    for eps in eps_list:
        marg = []
        for thA in (0.0, 0.3, 0.7, np.pi / 4):
            A = basis(thA); B = basis(0.2)
            w = np.zeros((2, 2))
            for a in range(2):
                for b in range(2):
                    v = np.kron(A[:, a], B[:, b])
                    w[a, b] = h(abs(np.vdot(v, psi)) ** 2, eps)
            w /= w.sum()
            marg.append(w[:, 0].sum())
        spread = max(marg) - min(marg)
        print(f"  {name:11s} eps={eps:4.2f}: Bob P(b=0) across Alice's 4 basis choices = {np.round(marg, 5)}  spread {spread:.2e}"
              + ("  SIGNALLING" if spread > 1e-12 else ""))

print("\n=== 5. premise audit ===")
print("  Gleason: frame function on projectors of a SUPPLIED Hilbert space, dim >= 3, noncontextual (depends on P only)")
print("  Busch: additive probability on the SUPPLIED effect algebra (all POVMs admitted as measurements)")
print("  no-signalling: SUPPLIED tensor-product composition + the joint-weight rule applied to the composite")
print("  contextual escape: Kochen-Specker forbids noncontextual value assignments, not contextual ones;")
print("  Bohmian mechanics gets h(p)=p only via the quantum-equilibrium initial distribution (supplied preparation)")
