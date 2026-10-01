"""SCOUT-1 W2-QE: de Broglie-Bohm relaxation in a 2D box (hbar = m = 1, box [0, pi]^2).
psi(x,y,t) = sum_{m,n<=M} c_mn (2/pi) sin(m x) sin(n y) exp(-i E_mn t),  E_mn = (m^2 + n^2)/2.
Particles move with v = Im(grad psi / psi). Coarse-grained H-function on a C x C grid:
  Hbar(t) = sum_cells p_c ln(p_c / q_c),  p_c = particle fraction, q_c = int_cell |psi|^2  (KL >= 0)
Runs: many-mode (16 modes, random phases), two-mode, single-mode, equilibrium start; dt convergence.
"""
import sys
import numpy as np

rng = np.random.default_rng(2024)


class Box:
    def __init__(self, M, amps, phases):
        self.M = M
        c = amps * np.exp(1j * phases)
        self.c = c / np.sqrt(np.sum(np.abs(c) ** 2))
        m = np.arange(1, M + 1)
        self.E = 0.5 * (m[:, None] ** 2 + m[None, :] ** 2)
        self.m = m

    def A(self, t):
        return (2 / np.pi) * self.c * np.exp(-1j * self.E * t)

    def fields(self, x, y, t):
        m = self.m
        Sx, Sy = np.sin(np.outer(x, m)), np.sin(np.outer(y, m))
        Cx, Cy = np.cos(np.outer(x, m)) * m, np.cos(np.outer(y, m)) * m
        A = self.A(t)
        SA = Sx @ A
        psi = np.sum(SA * Sy, axis=1)
        px = np.sum((Cx @ A) * Sy, axis=1)
        py = np.sum(SA * Cy, axis=1)
        return psi, px, py

    def vel(self, x, y, t, vmax=200.0):
        psi, px, py = self.fields(x, y, t)
        d = np.abs(psi) ** 2 + 1e-14
        vx = np.imag(np.conj(psi) * px) / d
        vy = np.imag(np.conj(psi) * py) / d
        return np.clip(vx, -vmax, vmax), np.clip(vy, -vmax, vmax)

    def density_cells(self, t, C, sub=8):
        g = (np.arange(C * sub) + 0.5) * np.pi / (C * sub)
        X, Y = np.meshgrid(g, g, indexing="ij")
        psi, _, _ = self.fields(X.ravel(), Y.ravel(), t)
        rho = (np.abs(psi) ** 2).reshape(C * sub, C * sub)
        q = rho.reshape(C, sub, C, sub).sum(axis=(1, 3))
        return q / q.sum()


def sample(f, N):
    """rejection sampling of a density f(x, y) on the box"""
    out = np.empty((0, 2))
    fmax = None
    while len(out) < N:
        P = rng.uniform(0, np.pi, (4 * N, 2))
        v = f(P[:, 0], P[:, 1])
        if fmax is None:
            fmax = 1.2 * v.max()
        keep = rng.uniform(0, fmax, len(v)) < v
        out = np.vstack([out, P[keep]])
    return out[:N]


def hbar(box, X, t, C):
    q = box.density_cells(t, C)
    H, _, _ = np.histogram2d(X[:, 0], X[:, 1], bins=C, range=[[0, np.pi], [0, np.pi]])
    p = H / H.sum()
    m = p > 0
    return np.sum(p[m] * np.log(p[m] / q[m]))


def evolve(box, X0, T, dt, C, nrec=9):
    X = X0.copy()
    rec_t = np.linspace(0, T, nrec)
    out = [(0.0, hbar(box, X, 0.0, C))]
    t = 0.0
    k = 1
    nsteps = int(round(T / dt))
    for s in range(nsteps):
        x, y = X[:, 0], X[:, 1]
        k1 = box.vel(x, y, t)
        k2 = box.vel(x + 0.5 * dt * k1[0], y + 0.5 * dt * k1[1], t + 0.5 * dt)
        k3 = box.vel(x + 0.5 * dt * k2[0], y + 0.5 * dt * k2[1], t + 0.5 * dt)
        k4 = box.vel(x + dt * k3[0], y + dt * k3[1], t + dt)
        X[:, 0] = x + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        X[:, 1] = y + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        X = np.clip(X, 1e-9, np.pi - 1e-9)
        t += dt
        if k < nrec and t >= rec_t[k] - 1e-9:
            out.append((t, hbar(box, X, t, C)))
            k += 1
    return out, X


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    C = 12
    T = 4 * np.pi
    ground = lambda x, y: np.sin(x) ** 2 * np.sin(y) ** 2
    print(f"N = {N} particles, coarse cells {C}x{C}, T = 4 pi (|psi|^2 recurs at 4 pi); sampling bias ~ {C*C/(2*N):.4f}")

    runs = {}
    M = 4
    phases = rng.uniform(0, 2 * np.pi, (M, M))
    many = Box(M, np.ones((M, M)), phases)
    X0 = sample(ground, N)
    runs["16 modes, ground-state start, dt=0.002"] = evolve(many, X0, T, 0.002, C)[0]
    runs["16 modes, ground-state start, dt=0.001 (convergence, to 2 pi)"] = evolve(many, X0, T / 2, 0.001, C, nrec=5)[0]

    eq_start = sample(lambda x, y: np.abs(many.fields(x, y, 0.0)[0]) ** 2, N)
    runs["16 modes, EQUILIBRIUM start (control)"] = evolve(many, eq_start, T, 0.002, C)[0]

    amps2 = np.zeros((M, M)); amps2[0, 1] = amps2[1, 0] = 1.0
    two = Box(M, amps2, phases)
    runs["2 modes (1,2)+(2,1), ground-state start"] = evolve(two, X0, T, 0.002, C)[0]

    amps1 = np.zeros((M, M)); amps1[1, 1] = 1.0
    one = Box(M, amps1, phases)
    runs["1 mode (2,2) stationary, ground-state start"] = evolve(one, X0, T, 0.002, C)[0]

    amps_r = rng.uniform(0.2, 1.0, (M, M))
    manyr = Box(M, amps_r, rng.uniform(0, 2 * np.pi, (M, M)))
    runs["16 modes, random amplitudes, ground-state start"] = evolve(manyr, X0, T, 0.002, C)[0]

    for name, out in runs.items():
        print(f"\n{name}")
        print("   t     : " + " ".join(f"{t:6.2f}" for t, _ in out))
        print("   Hbar  : " + " ".join(f"{h:6.4f}" for _, h in out))
