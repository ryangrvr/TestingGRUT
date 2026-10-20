"""Q(Omega) along the nodeless branch and fixed-Q Hessian at the charges of the outgoing products (~260, ~519).
Usage: python -I products_branch.py"""
import numpy as np, warnings
from scipy.integrate import solve_bvp
from scipy.linalg import eigh, eigh_tridiagonal
warnings.filterwarnings('ignore')
b = 0.1; R = 60.0
def solve(W, guess=None):
    r = np.linspace(1e-6, R, 3001)
    f = lambda r, y: np.vstack([y[1], (y[2]**2 + b*y[0]**2 - W**2)*y[0] - 2*y[1]/r, y[3], (y[2]**2 - 1 + y[0]**2)*y[2] - 2*y[3]/r])
    bc = lambda ya, yb: np.array([ya[1], ya[3], yb[0], yb[2]-1])
    if guess is None:
        F0 = 2.0/(1+np.exp(r-6.0)); S0 = 1-0.95/(1+np.exp(r-6.0)); y0 = np.vstack([F0, np.gradient(F0, r), S0, np.gradient(S0, r)])
    else:
        y0 = guess.sol(r)
    return solve_bvp(f, bc, r, y0, tol=1e-7, max_nodes=300000)
def QE(s, W):
    r = np.linspace(1e-6, R, 120001); F, Fp, S, Sp = s.sol(r)
    I = 4*np.pi*np.trapezoid(r**2*F**2, r)
    E = 4*np.pi*np.trapezoid(r**2*(.5*W**2*F**2+.5*Fp**2+.5*Sp**2+.25*(S**2-1)**2+.5*S**2*F**2+.25*b*F**4), r)
    return W*I, E
def spec(s, W, dr=0.05, k=2):
    r = np.arange(1, int(round(R/dr)))*dr; n = r.size; F, Fp, S, Sp = s.sol(r)
    ld = 2/dr**2*np.ones(n); lo = -1/dr**2*np.ones(n-1); out = {}
    out['phase0'] = eigh_tridiagonal(ld+S**2+b*F**2-W**2, lo, select='i', select_range=(0, 0))[0][0]
    for l, fq in [(0, True), (1, False), (2, False)]:
        c = l*(l+1)/r**2
        A = np.diag(ld+c+S**2+3*b*F**2-W**2)+np.diag(lo, 1)+np.diag(lo, -1); B = np.diag(ld+c+3*S**2-1+F**2)+np.diag(lo, 1)+np.diag(lo, -1)
        H = np.block([[A, np.diag(2*S*F)], [np.diag(2*S*F), B]])
        if fq:
            a = r*F; H[:n, :n] += 4*W**2/(np.sum(a*a)*dr)*np.outer(a, a)*dr
        out[f'l{l}'+('fixedQ' if fq else '')] = eigh(H, subset_by_index=[0, k-1], eigvals_only=True)
    out['threshold 1-W^2'] = 1-W**2
    return out
Ws = [0.72, 0.74, 0.76, 0.77, 0.78, 0.80, 0.82, 0.83, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94]
sols = {}; prev = None
for W in sorted(Ws, key=lambda w: abs(w-0.8)):
    pass
# continuation outward from 0.8
order = [0.80, 0.78, 0.77, 0.76, 0.74, 0.72, 0.82, 0.83, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94]
for W in order:
    g = None
    if sols:
        g = sols[min(sols, key=lambda w: abs(w-W))]
    s = solve(W, g); sols[W] = s
rows = []
for W in sorted(sols):
    Q, E = QE(sols[W], W); rows.append((W, Q, E))
    print('W=%.3f status=%d Q=%10.3f E=%10.3f E/Q=%.5f F0=%.4f S0=%.4f' % (W, sols[W].status, Q, E, E/Q, sols[W].sol(1e-6)[0], sols[W].sol(1e-6)[2]))
rows = np.array(rows)
for target in [260.0, 519.0]:
    # interpolate W at target Q (Q monotone decreasing in W expected)
    Wt = np.interp(target, rows[::-1, 1], rows[::-1, 0])
    s = solve(Wt, sols[min(sols, key=lambda w: abs(w-Wt))])
    Q, E = QE(s, Wt); h = 1e-3
    sp = solve(Wt+h, s); sm = solve(Wt-h, s)
    dQ = (QE(sp, Wt+h)[0]-QE(sm, Wt-h)[0])/(2*h)
    o = spec(s, Wt)
    print('target Q=%.0f -> W=%.5f Q=%.3f E/Q=%.5f dQ/dW=%.1f ; ' % (target, Wt, Q, E/Q, dQ) + ' '.join('%s=%s' % (k, np.round(v, 6)) for k, v in o.items()))
