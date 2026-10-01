"""SCOUT-1 W2-QE firewall: exact (backward-trajectory) evaluation of the transported ratio f = rho/|psi|^2.

Since f is advected by the Bohm flow, rho(x,t) = |psi(x,t)|^2 * f0(X0(x,t)), where X0 is the backward trajectory
of the grid point x from t to 0. On a fine G x G grid this gives the fine-grained density exactly, up to
trajectory error. Reports:
  1. fine-grained H_fine(t) = int rho ln f (should stay = H_fine(0): conserved);
  2. coarse-grained Hbar(t) for cell counts C = 4, 8, 16, 32 (refinement dependence);
  3. mode-count dependence (1, 2, 4, 9, 16, 25 modes);
  4. initial microstructure (smooth start vs a start with fine-scale structure);
  5. dt convergence of one case.
usage: python3 w2qe_firewall.py <config>
"""
import sys
import numpy as np
from w2qe_bohm_relaxation import Box

G = 96
Cs = (4, 8, 16, 32)


def back_trace(box, t_end, dt):
    g = (np.arange(G) + 0.5) * np.pi / G
    X, Y = np.meshgrid(g, g, indexing="ij")
    x, y = X.ravel().copy(), Y.ravel().copy()
    n = int(round(t_end / dt))
    h = -t_end / n
    t = t_end
    for _ in range(n):
        k1 = box.vel(x, y, t)
        k2 = box.vel(x + 0.5 * h * k1[0], y + 0.5 * h * k1[1], t + 0.5 * h)
        k3 = box.vel(x + 0.5 * h * k2[0], y + 0.5 * h * k2[1], t + 0.5 * h)
        k4 = box.vel(x + h * k3[0], y + h * k3[1], t + h)
        x = np.clip(x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), 1e-9, np.pi - 1e-9)
        y = np.clip(y + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]), 1e-9, np.pi - 1e-9)
        t += h
    return X.ravel(), Y.ravel(), x, y


def measures(box, rho0_fn, t_end, dt):
    Xg, Yg, x0, y0 = back_trace(box, t_end, dt) if t_end > 0 else (None, None, None, None)
    g = (np.arange(G) + 0.5) * np.pi / G
    X, Y = np.meshgrid(g, g, indexing="ij")
    psi_t = box.fields(X.ravel(), Y.ravel(), t_end)[0]
    q = np.abs(psi_t) ** 2
    if t_end > 0:
        psi0 = box.fields(x0, y0, 0.0)[0]
        f = rho0_fn(x0, y0) / (np.abs(psi0) ** 2 + 1e-300)
    else:
        f = rho0_fn(X.ravel(), Y.ravel()) / (q + 1e-300)
    rho = q * f
    dA = (np.pi / G) ** 2
    Z = rho.sum() * dA
    rho_n = rho / Z
    q_n = q / (q.sum() * dA)
    m = rho_n > 0
    H_fine = np.sum(rho_n[m] * np.log(rho_n[m] / q_n[m])) * dA
    Hbar = {}
    for C in Cs:
        s = G // C
        P = rho_n.reshape(C, s, C, s).sum(axis=(1, 3)) * dA
        Qc = q_n.reshape(C, s, C, s).sum(axis=(1, 3)) * dA
        mm = P > 0
        Hbar[C] = np.sum(P[mm] * np.log(P[mm] / Qc[mm]))
    return H_fine, Hbar, Z


def run(label, box, rho0_fn, times, dt=0.002):
    print(f"\n{label}")
    print("   t      norm(rho)   H_fine    " + "  ".join(f"Hbar(C={C:2d})" for C in Cs))
    for t in times:
        Hf, Hb, Z = measures(box, rho0_fn, t, dt)
        print(f"  {t:5.2f}   {Z:.4f}     {Hf:.4f}    " + "  ".join(f"{Hb[C]:11.4f}" for C in Cs), flush=True)


if __name__ == "__main__":
    cfg = sys.argv[1]
    rng = np.random.default_rng(2024)
    M = 5
    phases = rng.uniform(0, 2 * np.pi, (M, M))
    smooth = lambda x, y: (4 / np.pi ** 2) * np.sin(x) ** 2 * np.sin(y) ** 2   # ground-state density
    times = (0.0, np.pi, 4 * np.pi)

    def box_with(nm):
        amps = np.zeros((M, M))
        for k in range(nm):
            amps[k // M, k % M] = 1.0 if nm != 2 else 0.0
        if nm == 2:
            amps[0, 1] = amps[1, 0] = 1.0
        if nm == 4:
            amps[:] = 0; amps[:2, :2] = 1.0
        if nm == 9:
            amps[:] = 0; amps[:3, :3] = 1.0
        if nm == 16:
            amps[:] = 0; amps[:4, :4] = 1.0
        if nm == 25:
            amps[:] = 1.0
        if nm == 1:
            amps[:] = 0; amps[1, 1] = 1.0
        return Box(M, amps, phases)

    if cfg.startswith("modes"):
        for nm in [int(v) for v in cfg.split("_")[1:]]:
            run(f"{nm} mode(s), smooth ground-state start", box_with(nm), smooth, times)
    elif cfg == "micro":
        b = box_with(16)
        micro = lambda x, y: smooth(x, y) * (1 + 0.9 * np.sin(24 * x) * np.sin(24 * y))
        run("16 modes, start WITH fine microstructure (k = 24)", b, micro, times)
        eq_micro = lambda x, y: np.abs(b.fields(x, y, 0.0)[0]) ** 2 * (1 + 0.9 * np.sin(24 * x) * np.sin(24 * y))
        run("16 modes, |psi0|^2 start with fine microstructure only (f0 = 1 + 0.9 sin sin)", b, eq_micro, times)
    elif cfg == "dt":
        run("16 modes, smooth start, dt = 0.001 (convergence)", box_with(16), smooth, (np.pi, 4 * np.pi), dt=0.001)
