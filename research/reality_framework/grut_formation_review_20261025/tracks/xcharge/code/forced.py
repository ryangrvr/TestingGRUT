#!/usr/bin/env python3
"""Independent forced-radiation calculation for delta V = lam X^d + h.c. (track xcharge).

Linearize around the rotating ground state X0 = (v/sqrt2) F exp(-i Omega tau), phi0 = v S:
  dX   = v exp(-i Omega tau) [a exp(-i nu tau) + conj(b) exp(+i nu tau)],   nu = d Omega
  dphi = v [c exp(-i nu tau) + c.c.]
With U = rho (a, b, c) (s-wave) the sourced system is
  U'' + (K^2 - V) U = eps * rho F^(d-1) e_b,   K^2 = (W_i^2 - (1,1,6)),  W = ((d+1), (d-1), d) Omega
  V = [[S^2-1+2AF^2, AF^2, sqrt2 SF], [AF^2, S^2-1+2AF^2, sqrt2 SF], [sqrt2 SF, sqrt2 SF, B(3S^2-1)+F^2-6]]
  eps = d 2^((1-d)/2) lam v^(d-2) / m_X^2      (my normalization; Codex's eps is sqrt2 larger)
Physical number rate in channel i (derived in NOTES.md):  Ndot_i = (8 pi m_X / (g^2 hbar)) K_i |eps A_i|^2.

Two numerically distinct solutions of the same boundary problem:
  (A) particular + homogeneous shooting with DOP853 and outgoing matching at R,
  (B) Wronskian/Green formula A = W^{-1} int Phi^T S drho using only regular solutions.
"""
import numpy as np
from scipy.integrate import solve_ivp, simpson

from profile import G, VPHI, MX, A as AQ, B as BQ, ground_state

HBAR = 6.582119569e-25  # GeV s
MRED = 2.435e18          # GeV (Codex's value; 2.4353e18 is the CODATA-derived reduced Planck mass)
YEAR = 365.25 * 86400.0
T_REF = 13.8e9 * YEAR


class Bg:
    def __init__(self, sol, Om):
        self.sol, self.Om = sol, Om
        self.R = sol.x[-1]

    def FS(self, r):
        r = np.atleast_1d(r)
        y = self.sol.sol(np.minimum(r, self.R))
        F = np.where(r <= self.R, y[0], 0.0)
        S = np.where(r <= self.R, y[2], 1.0)
        return F, S

    def Vmat(self, r, scale=1.0):
        F, S = self.FS(r)
        F = F[0]; S = S[0]
        d = S * S - 1 + 2 * AQ * F * F
        o = np.sqrt(2) * S * F
        return scale * np.array([[d, AQ * F * F, o], [AQ * F * F, d, o], [o, o, BQ * (3 * S * S - 1) + F * F - 6]])


def channel_data(Om, d):
    W = np.array([(d + 1) * Om, (d - 1) * Om, d * Om])
    K2 = W**2 - np.array([1.0, 1.0, 6.0])
    return W, K2


def solve(bg, d, R=40.0, scale=1.0, rtol=1e-12, atol=1e-15):
    Om = bg.Om
    W, K2 = channel_data(Om, d)
    assert np.all(K2 > 0)
    K = np.sqrt(K2)

    def rhs(r, y, with_source):
        U = y[:3]; dU = y[3:]
        M = bg.Vmat(r, scale) - np.diag(K2)
        ddU = M @ U
        if with_source:
            F, _ = bg.FS(r)
            ddU = ddU.copy(); ddU[1] += r * F[0]**(d - 1)
        return np.concatenate([dU, ddU])

    grid = np.linspace(0, R, 8001)
    Phi = np.zeros((3, 3, grid.size)); dPhi = np.zeros((3, 3, grid.size))
    for j in range(3):
        y0 = np.zeros(6); y0[3 + j] = 1.0
        s = solve_ivp(rhs, (0, R), y0, args=(False,), method='DOP853', rtol=rtol, atol=atol, t_eval=grid)
        Phi[:, j] = s.y[:3]; dPhi[:, j] = s.y[3:]
    sp = solve_ivp(rhs, (0, R), np.zeros(6), args=(True,), method='DOP853', rtol=rtol, atol=atol, t_eval=grid)
    Up, dUp = sp.y[:3], sp.y[3:]

    E = np.exp(1j * K * R)
    # (A) shooting: U = Phi alpha + Up, impose U' = iK U at R
    M = dPhi[:, :, -1] - 1j * K[:, None] * Phi[:, :, -1]
    alpha = np.linalg.solve(M, -(dUp[:, -1] - 1j * K * Up[:, -1]))
    U = np.einsum('ijn,j->in', Phi.astype(complex), alpha) + Up
    A_shoot = U[:, -1] / E

    # (B) Green/Wronskian formula
    Psi = np.diag(E); dPsi = np.diag(1j * K * E)
    Wr = Phi[:, :, -1].T @ dPsi - dPhi[:, :, -1].T @ Psi
    F, _ = bg.FS(grid)
    src = grid * F**(d - 1)
    J = simpson(Phi[1, :, :] * src[None, :], x=grid, axis=1)
    A_green = np.linalg.solve(Wr, J)

    flux = K * np.abs(A_shoot)**2
    work = simpson(np.imag(np.conj(U[1]) * src), x=grid)
    return dict(d=d, R=R, scale=scale, Omega=Om, W=W, K=K, A_shoot=A_shoot, A_green=A_green,
                flux=flux, flux_green=K * np.abs(A_green)**2, source_work=work,
                rel_shoot_vs_green=float(np.max(np.abs(A_shoot - A_green)) / np.max(np.abs(A_green))),
                rel_work=float(abs(flux.sum() - work) / flux.sum()),
                max_field=float(np.max(np.abs(np.divide(U, grid[None, :], out=np.zeros_like(U), where=grid[None, :] > 0)))))


def eps_mine(d, c=1.0, M=MRED):
    return d * 2**((1 - d) / 2) * abs(c) * VPHI**(d - 2) / (MX**2 * M**(d - 4))


def rates(res, Q, c=1.0, M=MRED):
    d = res['d']; Om = res['Omega']
    e = eps_mine(d, c, M)
    Ndot = 8 * np.pi * MX / (G**2 * HBAR) * e**2 * res['flux']  # per second
    core = float(np.dot(res['W'] / Om, Ndot))
    total_violation = float(d * Ndot.sum())
    frac = core / Q
    return dict(eps=e, Ndot=Ndot.tolist(), core_loss_per_s=core, total_X_violation_per_s=total_violation,
                frac_rate_per_s=frac, frac_times_age=frac * T_REF,
                c_max_1pct=float(np.sqrt(0.01 / (frac * T_REF)) * abs(c)),
                M_unit_coeff_1pct=float(M * (frac * T_REF / 0.01)**(1 / (2 * (d - 4)))))


def first_order_excursion(bg, d, Q, c=1.0, M=MRED):
    """Peak-to-peak first-order core-charge excursion fraction: 4|lam| int |X0|^d d^3x / (omega Q)."""
    x = np.linspace(0, bg.R, 40001)
    F, _ = bg.FS(x)
    integral = 4 * np.pi * simpson(x * x * (VPHI * F / np.sqrt(2))**d, x=x) / MX**3  # GeV^(d-3)
    lam = abs(c) / M**(d - 4)
    return float(4 * lam * integral / (bg.Om * MX * Q))


if __name__ == '__main__':
    import json, sys
    Q = 200000.
    Om, sol = ground_state(Q)
    bg = Bg(sol, Om)
    out = dict(Q=Q, Omega=Om, rows=[])
    for d in [5, 6]:
        for R in [30., 40., 55.]:
            r = solve(bg, d, R=R)
            rr = rates(r, Q)
            row = dict(d=d, R=R, K=r['K'].tolist(), flux=r['flux'].tolist(), flux_sum=float(r['flux'].sum()),
                       source_work=r['source_work'], rel_work=r['rel_work'], rel_shoot_vs_green=r['rel_shoot_vs_green'],
                       A_re=r['A_shoot'].real.tolist(), A_im=r['A_shoot'].imag.tolist(), max_unit_field=r['max_field'], **rr)
            out['rows'].append(row)
            print(d, R, 'flux', r['flux'], 'frac', rr['frac_rate_per_s'], 'cmax', rr['c_max_1pct'], 'shoot/green', r['rel_shoot_vs_green'], 'work', r['rel_work'], flush=True)
        f = solve(bg, d, R=40., scale=0.0)
        rf = rates(f, Q)
        out['rows'].append(dict(d=d, R=40., free=True, flux=f['flux'].tolist(), frac_rate_per_s=rf['frac_rate_per_s'], A_free=f['A_shoot'].real.tolist()))
        print(d, 'free flux', f['flux'], 'frac', rf['frac_rate_per_s'])
        out['first_order_pp_excursion_d%d' % d] = first_order_excursion(bg, d, Q)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else '/dev/stdout', 'w'), indent=1, default=float)
