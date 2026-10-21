"""Verifier v_reach: Q-ball family by DIRECT CONSTRAINED ENERGY MINIMISATION at fixed charge (independent of the
reach track's solve_bvp continuation, and of the archived profile: own tanh initial guess).

Units (Codex/charge_phase): S = phi/v, F^2 = 2|X|^2/v^2, ell = 1/(sqrt(lam_phi) v), a = 1, b = 0.1.
At fixed Qbar the X rotation energy is Qbar^2/(2 I), I = int F^2 d^3x, Omega = Qbar/I, so
  E[F,S] = int [ S'^2/2 + F'^2/2 + (S^2-1)^2/4 + a S^2 F^2/2 + b F^4/4 ] d^3x + Qbar^2/(2 I).
Cell-centred radial grid r_j = (j+1/2) h, Dirichlet F = 0, S = 1 at r = (N+1/2) h. L-BFGS with analytic gradient,
two resolutions (h = 0.08, 0.04) and Richardson extrapolation (second-order scheme). Q_phys = 100 Qbar,
e = E/(m_X Q) = Ebar/(sqrt(a) Qbar).
"""
import json, sys, time
import numpy as np
from scipy.optimize import minimize

a, b = 1.0, 0.1


class Ball:
    def __init__(self, h, L):
        self.h = h
        self.N = int(round(L / h))
        j = np.arange(self.N)
        self.r = (j + 0.5) * h
        self.wv = 4 * np.pi * self.r**2 * h               # volume weights
        self.wl = 4 * np.pi * ((j + 1) * h) ** 2 / h       # link weights (gradient), link j -> j+1

    def energy(self, z, Qb):
        N, h = self.N, self.h
        F, S = z[:N], z[N:]
        Fe = np.append(F, 0.0); Se = np.append(S, 1.0)
        dF = np.diff(Fe); dS = np.diff(Se)
        Egrad = 0.5 * np.sum(self.wl * (dF**2 + dS**2))
        V = 0.25 * (S**2 - 1) ** 2 + 0.5 * a * S**2 * F**2 + 0.25 * b * F**4
        Epot = np.sum(self.wv * V)
        I = np.sum(self.wv * F**2)
        Ekin = Qb**2 / (2 * I)
        E = Egrad + Epot + Ekin
        # gradient
        gF = np.zeros(N); gS = np.zeros(N)
        tF = self.wl * dF; tS = self.wl * dS
        gF += -tF; gF[1:] += tF[:-1]
        gS += -tS; gS[1:] += tS[:-1]
        gF += self.wv * (a * S**2 * F + b * F**3) - Qb**2 / I**2 * self.wv * F
        gS += self.wv * (S * (S**2 - 1) + a * S * F**2)
        return E, np.concatenate([gF, gS]), dict(Egrad=Egrad, Epot=Epot, Ekin=Ekin, I=I)

    def solve(self, Qb, z0, gtol=1e-11):
        f = lambda z: self.energy(z, Qb)[:2]
        res = minimize(f, z0, jac=True, method='L-BFGS-B',
                       options=dict(maxiter=200000, maxcor=50, ftol=1e-16, gtol=gtol, maxfun=400000))
        E, g, parts = self.energy(res.x, Qb)
        return res, E, parts

    def guess(self, Qb, R=None, w=1.0):
        nst = b ** -0.25
        if R is None:
            R = (3 * Qb / (4 * np.pi * nst)) ** (1 / 3)
        prof = 0.5 * (1 - np.tanh((self.r - R) / w))
        F = 1.6 * prof
        S = 1 - 0.95 * prof
        return np.concatenate([F, S])

    def interp(self, other_r, z_other, N_other):
        F = np.interp(self.r, other_r, z_other[:N_other], right=0.0)
        S = np.interp(self.r, other_r, z_other[N_other:], right=1.0)
        return np.concatenate([F, S])

    def observables(self, z, Qb, E, parts):
        N = self.N
        F, S = z[:N], z[N:]
        Om = Qb / parts['I']
        # energy density on cells (gradient split to cells approximately) for R90
        Fe = np.append(F, 0.0); Se = np.append(S, 1.0)
        dF = np.diff(Fe) / self.h; dS = np.diff(Se) / self.h
        edens = 0.5 * (dF**2 + dS**2) + 0.25 * (S**2 - 1) ** 2 + 0.5 * a * S**2 * F**2 + 0.25 * b * F**4 + 0.5 * Om**2 * F**2
        cum = np.cumsum(self.wv * edens)
        r90 = float(np.interp(0.9 * cum[-1], cum, self.r + 0.5 * self.h))
        return dict(Qbar=Qb, Q_phys=100 * Qb, E_bar=E, e=E / (np.sqrt(a) * Qb), Omega=Om, R90E=r90, S0=float(S[0]),
                    F0=float(F[0]), Fedge=float(np.max(np.abs(F[-20:]))))


def run_family(Qlist, h, L_of_Q, prev=None):
    out = []
    for Qb in Qlist:
        L = L_of_Q(Qb)
        B = Ball(h, L)
        z0 = B.guess(Qb) if prev is None else B.interp(prev[0], prev[1], prev[2])
        t = time.time()
        res, E, parts = B.solve(Qb, z0)
        ob = B.observables(res.x, Qb, E, parts)
        ob.update(h=h, L=L, converged=bool(res.success), nit=int(res.nit), msg=str(res.message)[:60],
                  gnorm=float(np.max(np.abs(B.energy(res.x, Qb)[1]))), sec=round(time.time() - t, 1))
        print(f"h={h} Q={100*Qb:.4g} e={ob['e']:.7f} Om={ob['Omega']:.5f} R90={ob['R90E']:.4f} S0={ob['S0']:.3f} "
              f"F0={ob['F0']:.3f} edge={ob['Fedge']:.1e} it={res.nit} g={ob['gnorm']:.1e} {ob['sec']}s", flush=True)
        out.append(ob)
        prev = (B.r, res.x, B.N)
    return out


if __name__ == '__main__':
    outpath = sys.argv[1]
    L_of_Q = lambda Qb: max(30.0, 3.0 * (3 * Qb / (4 * np.pi * 1.778)) ** (1 / 3) + 18.0)
    # archived point first, then up and down from it (continuation in Q from a fresh tanh guess)
    Q_down = [399.44, 330.0, 280.0, 240.0, 210.0, 195.0, 185.0, 180.0, 177.0, 174.0, 170.0, 165.0, 160.0, 155.0, 150.0]
    Q_up = [399.44, 600.0, 1000.0, 1500.0, 2500.0, 4000.0, 6185.3]
    fam = {}
    for h in (0.08, 0.04):
        fam[h] = dict(down=run_family(Q_down, h, L_of_Q), up=run_family(Q_up, h, L_of_Q))
    json.dump({str(k): v for k, v in fam.items()}, open(outpath, 'w'), indent=1)
