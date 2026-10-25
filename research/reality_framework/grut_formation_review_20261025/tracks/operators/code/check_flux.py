#!/usr/bin/env python3
"""Item 3: are the extraction currents exact discrete conservation laws? Own code.

Own vectorised implementation of the radial three-field system (U = r x, r y, r(S-1))
with the quartic force, the staggered fourth-order Hamiltonian and the causal CFS memory.
(a) instantaneous (semi-discrete) identity, by an exact complex-step derivative of the
    polynomial E_inside / Q_inside along the RHS:  dE_in/dt + J_E = 0, dQ_in/dt + J_Q = 0;
    global ledgers dH/dt = sum V.B and dQ/dt = sum (y B_x - x B_y);
(b) fully discrete RK4 closure of E_in + int J_E (same RK stages) versus dt.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from lattice import G_matrix, D_matrix, profiles

OUT = Path(__file__).resolve().parents[1] / 'results'


class Radial:
    def __init__(self, n, h, R, smax=20., alpha=.1, cut=None):
        self.n, self.h, self.alpha = n, h, alpha
        self.G = G_matrix(n, h); self.D = D_matrix(n, h); self.L = self.D @ self.G
        self.sn, self.sf = profiles(n, h, R, smax)
        self.r = np.arange(1, n) * h
        self.w = 4 * np.pi * h
        self.cut = cut                      # last node index (1-based) inside the cut
        N = n - 1
        self.N = N

    # state: dict of arrays U(3,N) V(3,N) P(3,n) Qm(3,N)
    def force(self, U):
        r = self.r
        x, y, z = U
        rho = (x * x + y * y) / (r * r)
        S = 1 + z / r
        g = S * S + .1 * rho
        return np.array([g * x, g * y, r * (3 * (S * S - 1) + rho) * S])

    def potential(self, U):
        r = self.r
        x, y, z = U
        rho = (x * x + y * y) / (r * r)
        S = 1 + z / r
        return r * r * (.75 * (S * S - 1) ** 2 + .5 * S * S * rho + .025 * rho * rho)

    def rhs(self, Y):
        U, V, P, Qm = Y
        LU = U @ self.L.T
        DP = P @ self.D.T
        B = DP + Qm
        dU = V.copy()
        dV = LU - self.force(U) + B
        dP = -(self.alpha + self.sf) * P - self.sf * (U @ self.G.T)
        dQ = -(self.alpha + self.sn) * Qm - self.sn * (LU + DP)
        return [dU, dV, dP, dQ], B

    def node_energy(self, U, V):
        LU = U @ self.L.T
        return np.sum(.5 * V * V - .5 * U * LU, axis=0) + self.potential(U)

    def charge_density(self, U, V):
        return U[1] * V[0] - U[0] * V[1]

    def E_in(self, U, V):
        return self.w * np.sum(self.node_energy(U, V)[:self.cut])

    def Q_in(self, U, V):
        return self.w * np.sum(self.charge_density(U, V)[:self.cut])

    def H(self, U, V):
        return self.w * np.sum(self.node_energy(U, V))

    def Q(self, U, V):
        return self.w * np.sum(self.charge_density(U, V))

    def currents(self, U, V):
        """Outward signed currents across the cut between node c and c+1 (1-based)."""
        c, Lm = self.cut, self.L
        JE = JQ = 0.0
        for i in range(max(1, c - 2), c + 1):
            for j in range(c + 1, min(self.n - 1, i + 3) + 1):
                Lij = Lm[i - 1, j - 1]
                JE += .5 * self.w * Lij * np.sum(U[:, i - 1] * V[:, j - 1] - V[:, i - 1] * U[:, j - 1])
                JQ += self.w * Lij * (U[0, i - 1] * U[1, j - 1] - U[1, i - 1] * U[0, j - 1])
        return JE, JQ

    def ledger_rates(self, Y, B):
        U, V = Y[0], Y[1]
        return -self.w * np.sum(V * B), self.w * np.sum(U[0] * B[1] - U[1] * B[0])


def initial(sys_, rng=None, layer_memory=True):
    r = sys_.r
    env = np.exp(-(r - 8) ** 2 / 18)
    U = np.array([r * .9 * env * np.cos(.3 * r), r * .8 * env * np.sin(.3 * r), r * .05 * np.exp(-(r - 6) ** 2 / 8)])
    V = np.array([U[1] * 0.9, -U[0] * 0.9, r * .02 * np.exp(-(r - 10) ** 2 / 10)])
    P = np.zeros((3, sys_.n)); Qm = np.zeros((3, sys_.N))
    if layer_memory:   # arbitrary nonzero memory in the layer to exercise the ledger
        P[:, sys_.sf > 0] = .01 * np.sin(np.arange((sys_.sf > 0).sum()))[None, :]
        Qm[:, sys_.sn > 0] = .01 * np.cos(np.arange((sys_.sn > 0).sum()))[None, :]
    return [U, V, P, Qm]


def instantaneous(sys_, seed=3):
    # random smooth fields spanning the cut, plus nonzero memory only where sigma > 0
    rng = np.random.default_rng(seed)
    Y = initial(sys_)
    r = sys_.r
    for a in (0, 1):
        for f in range(3):
            c = rng.normal(size=12)
            Y[a][f] = Y[a][f] + r * sum(c[m] * np.sin((m + 1) * np.pi * r / r[-1] + c[-1]) for m in range(10)) * .05
    Y[0][2] *= .2
    d, B = sys_.rhs(Y)
    eps = 1e-30
    Uc = Y[0] + 1j * eps * d[0]; Vc = Y[1] + 1j * eps * d[1]
    dEin = np.imag(sys_.E_in(Uc, Vc)) / eps
    dQin = np.imag(sys_.Q_in(Uc, Vc)) / eps
    dH = np.imag(sys_.H(Uc, Vc)) / eps
    dQ = np.imag(sys_.Q(Uc, Vc)) / eps
    JE, JQ = sys_.currents(Y[0], Y[1])
    Er, Qr = sys_.ledger_rates(Y, B)
    scaleE = abs(sys_.w * np.sum(Y[1] * (Y[0] @ sys_.L.T)))
    scaleQ = abs(sys_.w * np.sum(Y[0][0] * (Y[0][1] @ sys_.L.T)))
    return dict(dEin_dt=float(dEin), J_E=float(JE), cut_energy_identity_residual_rel=float(abs(dEin + JE) / scaleE),
                dQin_dt=float(dQin), J_Q=float(JQ), cut_charge_identity_residual_rel=float(abs(dQin + JQ) / scaleQ),
                global_energy_ledger_residual_rel=float(abs(dH + Er) / scaleE),
                global_charge_ledger_residual_rel=float(abs(dQ + Qr) / scaleQ))


def rk4_run(sys_, dt, T):
    Y = initial(sys_, layer_memory=False)
    E0, Q0 = sys_.E_in(Y[0], Y[1]), sys_.Q_in(Y[0], Y[1])
    H0, QQ0 = sys_.H(Y[0], Y[1]), sys_.Q(Y[0], Y[1])
    acc = np.zeros(4)   # int J_E, int J_Q, E_rec, Q_rec
    steps = int(round(T / dt))
    worst = np.zeros(4)

    def aux(Y):
        d, B = sys_.rhs(Y)
        JE, JQ = sys_.currents(Y[0], Y[1])
        Er, Qr = sys_.ledger_rates(Y, B)
        return d, np.array([JE, JQ, Er, Qr])

    for _ in range(steps):
        k1, a1 = aux(Y)
        Y2 = [y + .5 * dt * k for y, k in zip(Y, k1)]; k2, a2 = aux(Y2)
        Y3 = [y + .5 * dt * k for y, k in zip(Y, k2)]; k3, a3 = aux(Y3)
        Y4 = [y + dt * k for y, k in zip(Y, k3)]; k4, a4 = aux(Y4)
        Y = [y + dt / 6 * (p + 2 * q + 2 * s + t) for y, p, q, s, t in zip(Y, k1, k2, k3, k4)]
        acc += dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        e = np.array([sys_.E_in(Y[0], Y[1]) + acc[0] - E0, sys_.Q_in(Y[0], Y[1]) + acc[1] - Q0,
                      sys_.H(Y[0], Y[1]) + acc[2] - H0, sys_.Q(Y[0], Y[1]) + acc[3] - QQ0])
        worst = np.maximum(worst, abs(e))
    return dict(dt=dt, T=T, max_cut_energy_closure_over_E0=float(worst[0] / abs(E0)),
                max_cut_charge_closure_over_Q0=float(worst[1] / abs(Q0)),
                max_global_energy_ledger_over_H0=float(worst[2] / abs(H0)),
                max_global_charge_ledger_over_Q0=float(worst[3] / abs(QQ0)),
                net_E_through_cut=float(acc[0]), net_Q_through_cut=float(acc[1]))


def main():
    out = {}
    s1 = Radial(n=240, h=.25, R=45., cut=120)       # cut at r = 30.125 (between nodes 120,121)
    out['instantaneous_semidiscrete'] = [instantaneous(s1, seed) for seed in (3, 4, 5)]
    print(out['instantaneous_semidiscrete'], flush=True)
    rows = []
    for dt in (.1, .05, .025, .0125):
        rows.append(rk4_run(s1, dt, 60.))
        print(rows[-1], flush=True)
    for a, b in zip(rows, rows[1:]):
        for key in ('max_cut_energy_closure_over_E0', 'max_cut_charge_closure_over_Q0', 'max_global_energy_ledger_over_H0'):
            b['order_' + key] = float(np.log(a[key] / b[key]) / np.log(2))
    out['rk4_closure_vs_dt'] = rows
    (OUT / 'flux_checks.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
