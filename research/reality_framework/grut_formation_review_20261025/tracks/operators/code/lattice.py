"""Independent (Claude, operators track) assembly of the staggered fourth-order
radial operators and the causal CFS linear system.

Written from the equations, not from Codex's code. Conventions:
  nodes j = 0..n with U_0 = U_n = 0 (r = j h, wall at r = n h);
  interior unknowns j = 1..n-1 (index j-1);
  faces f = 0..n-1 sit at r = (f + 1/2) h;
  node ghosts are odd about j=0 and j=n, face ghosts even about the same points.
"""
import numpy as np

C1, C3 = 9.0 / 8.0, -1.0 / 24.0          # staggered fourth-order weights
A0, A1, A2, A3 = -365 / 144, 87 / 64, -3 / 32, 1 / 576


def node_fold(j, n):
    """Return (index, sign) of node j after odd continuation; index None for a zero node."""
    s = 1.0
    if j < 0:
        j, s = -j, -1.0
    if j > n:
        j, s = 2 * n - j, -s
    if j == 0 or j == n:
        return None, 0.0
    return j, s


def face_fold(f, n):
    if f < 0:
        f = -1 - f
    if f >= n:
        f = 2 * n - 1 - f
    return f


def G_matrix(n, h):
    """Face gradient: (GU)_{f+1/2} = [C1 (U_{f+1}-U_f) + C3 (U_{f+2}-U_{f-1})]/h."""
    G = np.zeros((n, n - 1))
    for f in range(n):
        for j, c in ((f + 1, C1), (f, -C1), (f + 2, C3), (f - 1, -C3)):
            idx, s = node_fold(j, n)
            if idx is not None:
                G[f, idx - 1] += s * c / h
    return G


def D_matrix(n, h):
    """Node divergence with even face ghosts."""
    D = np.zeros((n - 1, n))
    for j in range(1, n):
        for f, c in ((j, C1), (j - 1, -C1), (j + 1, C3), (j - 2, -C3)):
            D[j - 1, face_fold(f, n)] += c / h
    return D


def L_matrix(n, h):
    """Seven-point stencil quoted by Codex with odd node ghosts (built independently)."""
    L = np.zeros((n - 1, n - 1))
    for j in range(1, n):
        for d, c in ((0, A0), (1, A1), (-1, A1), (2, A2), (-2, A2), (3, A3), (-3, A3)):
            idx, s = node_fold(j + d, n)
            if idx is not None:
                L[j - 1, idx - 1] += s * c / h ** 2
    return L


def profiles(n, h, R, sigma_max, power=3):
    """sigma at nodes and faces; zero for r <= R, cubic ramp, clipped at the wall."""
    W = n * h - R
    rn = np.arange(1, n) * h
    rf = (np.arange(n) + 0.5) * h
    sn = sigma_max * np.clip((rn - R) / W, 0, 1) ** power
    sf = sigma_max * np.clip((rf - R) / W, 0, 1) ** power
    return sn, sf


def cfs_system(n, h, m, R, sigma_max=20.0, alpha=0.1, sn=None, sf=None):
    """Linear first-order system Y' = A Y, Y = (U, V, psi, q) with my own derivation:
        U' = V
        V' = L U - m^2 U + D psi + q
        psi' = -(alpha+sf) psi - sf G U
        q'   = -(alpha+sn) q - sn (L U + D psi)
    """
    G, D = G_matrix(n, h), D_matrix(n, h)
    L = D @ G
    if sn is None:
        sn, sf = profiles(n, h, R, sigma_max)
    N, F = n - 1, n
    I = np.eye(N)
    A = np.zeros((3 * N + F, 3 * N + F))
    u, v, p, q = slice(0, N), slice(N, 2 * N), slice(2 * N, 2 * N + F), slice(2 * N + F, 3 * N + F)
    A[u, v] = I
    A[v, u] = L - m * m * I
    A[v, p] = D
    A[v, q] = I
    A[p, p] = -np.diag(alpha + sf)
    A[p, u] = -np.diag(sf) @ G
    A[q, q] = -np.diag(alpha + sn)
    A[q, u] = -np.diag(sn) @ L
    A[q, p] = -np.diag(sn) @ D
    return A, dict(G=G, D=D, L=L, sn=sn, sf=sf, slices=(u, v, p, q))


def stretched_operator(s, G, D, sn, sf, alpha):
    """Frequency-domain C_n D C_f G for Laplace variable s."""
    Cn = (s + alpha) / (s + alpha + sn)
    Cf = (s + alpha) / (s + alpha + sf)
    return (Cn[:, None] * D) @ (Cf[:, None] * G)


def lattice_momentum(k, h):
    t = np.sin(k * h / 2)
    return 2 * t / h * (1 + t * t / 6)
