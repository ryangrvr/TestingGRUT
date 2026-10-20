"""Independent fixed-charge Hessian check for the isolated Candidate-1 Q-ball (D=0, a=1, b=.1).
Own BVP solve from a generic tanh guess (does not read profile.dat), own FD Hessians.
F'' + 2F'/r = (S^2 + b F^2 - W^2) F ;  S'' + 2S'/r = (S^2 - 1 + F^2) S
Usage: python -I independent_spectral.py
"""
import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import eigh_tridiagonal, eigh
from scipy.interpolate import CubicSpline

b = 0.1


def solve(W, R=50.0, guess=None):
    r = np.linspace(1e-6, R, 2001)

    def f(r, y):
        F, Fp, S, Sp = y
        return np.vstack([Fp, (S**2 + b*F**2 - W**2)*F - 2*Fp/r, Sp, (S**2 - 1 + F**2)*S - 2*Sp/r])

    def bc(ya, yb):
        return np.array([ya[1], ya[3], yb[0], yb[2] - 1])
    if guess is None:
        F0 = 2.0/(1 + np.exp(r - 6.0)); S0 = 1 - 0.95/(1 + np.exp(r - 6.0))
        y0 = np.vstack([F0, np.gradient(F0, r), S0, np.gradient(S0, r)])
    else:
        y0 = guess.sol(r)
    s = solve_bvp(f, bc, r, y0, tol=1e-9, max_nodes=200000)
    return s


def charges(s, W, R=50.0):
    r = np.linspace(0, R, 200001); F, Fp, S, Sp = s.sol(np.maximum(r, 1e-6))
    I = 4*np.pi*np.trapezoid(r**2*F**2, r)
    E = 4*np.pi*np.trapezoid(r**2*(.5*W**2*F**2 + .5*Fp**2 + .5*Sp**2 + .25*(S**2-1)**2 + .5*S**2*F**2 + .25*b*F**4), r)
    return W*I, E, I


W = 0.8000000019763136
s = solve(W)
print('BVP status', s.status, s.message, 'nodes', s.x.size)
Q, E, I = charges(s, W)
print('Q=%.6f E=%.6f E/Q=%.7f F(0)=%.6f S(0)=%.6f' % (Q, E, E/Q, s.sol(1e-6)[0], s.sol(1e-6)[2]))
dW = 2e-4
sp = solve(W + dW, guess=s); sm = solve(W - dW, guess=s)
Qp = charges(sp, W + dW)[0]; Qm = charges(sm, W - dW)[0]
print('dQ/dOmega (central, dW=%g) = %.3f' % (dW, (Qp - Qm)/(2*dW)))


def spectra(R, dr, k=3):
    r = np.arange(1, int(round(R/dr)))*dr; n = r.size
    F, Fp, S, Sp = s.sol(r) if R <= 50 else ext(r)
    lap_d = 2/dr**2*np.ones(n); lap_o = -1/dr**2*np.ones(n-1)
    out = {}
    # phase block l=0
    dm = lap_d + S**2 + b*F**2 - W**2
    out['phase_l0'] = eigh_tridiagonal(dm, lap_o, select='i', select_range=(0, k-1))[0]

    def plus(l, fixedQ=False):
        c = l*(l+1)/r**2
        A = np.diag(lap_d + c + S**2 + 3*b*F**2 - W**2) + np.diag(lap_o, 1) + np.diag(lap_o, -1)
        B = np.diag(lap_d + c + 3*S**2 - 1 + F**2) + np.diag(lap_o, 1) + np.diag(lap_o, -1)
        C = np.diag(2*S*F)
        H = np.block([[A, C], [C, B]])
        if fixedQ:
            a = r*F; norm = np.sum(a*a)*dr
            H[:n, :n] += 4*W**2/norm*np.outer(a, a)*dr
        return eigh(H, subset_by_index=[0, k-1], eigvals_only=True)
    out['unconstrained_l0'] = plus(0); out['fixed_charge_l0'] = plus(0, True)
    out['translation_l1'] = plus(1); out['shape_l2'] = plus(2); out['l3'] = plus(3)
    return out


def ext(r):
    # extend beyond R=50 with vacuum tails (F~0,S~1): profile is ~1e-14 there
    y = np.zeros((4, r.size)); m = r <= 50
    y[:, m] = s.sol(r[m]); y[2, ~m] = 1.0
    return y


for R, dr in [(50, .1), (50, .05), (50, .025), (60, .05)]:
    o = spectra(R, dr)
    print('R=%g dr=%g ' % (R, dr) + ' '.join('%s=%s' % (k, np.array2string(v, precision=9)) for k, v in o.items()))
