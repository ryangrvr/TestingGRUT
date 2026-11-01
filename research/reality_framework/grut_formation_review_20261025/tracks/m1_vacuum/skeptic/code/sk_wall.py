"""m1_vacuum SKEPTIC, part 5 [NEW]: planar interface (surface tension) between the zero-pressure bulk Q-matter phase
(S_b, F_b, Omega_b) and the vacuum (1, 0), and the thin-wall coefficient of E/Q - Omega_b = sigma (36 pi)^{1/3} n_b^{-2/3} Q^{-1/3}.
Independent code (Claude skeptic, 1 Nov 2026). Purpose: test the track's two-point fit E/Q ~ 0.77586 + 8.2 Q^{-1/3} + 100 Q^{-2/3}
against a direct surface tension, and give the asymptotic coefficient with its own error estimate.

1D static equations at Omega = Omega_b:  F'' = (S^2 + a F^2 - Omega_b^2) F,  S'' = U_phi'(S) + F^2 S,
with F -> F_b, S -> S_b at x -> -inf and F -> 0, S -> 1 at x -> +inf.  First integral: (F'^2 + S'^2)/2 = W(F, S)
= U - Omega_b^2 F^2/2 (which vanishes in both phases at the zero-pressure point), so sigma = int (F'^2 + S'^2) dx.
Second-order FD on a uniform grid, Dirichlet ends, damped Newton; Richardson over h.
"""
import json, math, sys
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import splu

a, b, c = 0.1, 3.0, 0.2
G2 = 0.01
Sb, Fb, Omb = 0.1511721078446304, 2.406476756777611, 0.7758647074595397
nb = Omb * Fb * Fb / G2


def Uphi(S):   return b / 4 * (S * S - 1)**2 - c * (S - 1)**3
def Uphi_p(S): return b * (S * S - 1) * S - 3 * c * (S - 1)**2
def Uphi_pp(S): return b * (3 * S * S - 1) - 6 * c * (S - 1)


def solve(Lx, h):
    N = int(round(2 * Lx / h)) - 1
    x = -Lx + (np.arange(N) + 1) * h
    lap = sps.diags([np.ones(N - 1), -2 * np.ones(N), np.ones(N - 1)], [-1, 0, 1]) / h**2
    bF = np.zeros(N); bF[0] = Fb / h**2            # left boundary F = Fb, right F = 0
    bS = np.zeros(N); bS[0] = Sb / h**2; bS[-1] = 1.0 / h**2
    F = Fb / (1 + np.exp(x / 1.0)); S = 1 - (1 - Sb) / (1 + np.exp(x / 1.0))
    for it in range(100):
        GF = (S * S + a * F * F - Omb * Omb) * F; GS = Uphi_p(S) + F * F * S
        R = np.r_[lap @ F + bF - GF, lap @ S + bS - GS]
        nr = np.max(np.abs(R))
        if nr < 1e-11:
            break
        J = sps.bmat([[lap - sps.diags(S * S + 3 * a * F * F - Omb * Omb), -sps.diags(2 * S * F)],
                      [-sps.diags(2 * F * S), lap - sps.diags(Uphi_pp(S) + F * F)]], format='csc')
        dx = splu(J).solve(-R)
        lam = 1.0
        while lam > 1e-4:
            Fn, Sn = F + lam * dx[:N], S + lam * dx[N:]
            Rn = np.r_[lap @ Fn + bF - (Sn * Sn + a * Fn * Fn - Omb * Omb) * Fn, lap @ Sn + bS - Uphi_p(Sn) - Fn * Fn * Sn]
            if np.linalg.norm(Rn) < np.linalg.norm(R):
                break
            lam *= 0.5
        F, S = Fn, Sn
    Fe = np.r_[Fb, F, 0.0]; Se = np.r_[Sb, S, 1.0]
    dF = (Fe[1:] - Fe[:-1]) / h; dS = (Se[1:] - Se[:-1]) / h          # face derivatives
    sigma_grad = np.sum(dF * dF + dS * dS) * h                        # int (F'^2 + S'^2) dx
    W = Uphi(S) + 0.5 * S * S * F * F + a / 4 * F**4 - 0.5 * Omb * Omb * F * F
    sigma_pot = 2 * np.sum(W) * h                                     # int 2 W dx (equal by the first integral)
    sigma_mixed = 0.5 * np.sum(dF * dF + dS * dS) * h + np.sum(W) * h
    # wall position (F = Fb/2) and width
    xm = float(np.interp(Fb / 2, F[::-1], x[::-1]))
    return dict(Lx=Lx, h=h, N=N, newton_res=float(nr), iters=it, sigma_grad=float(sigma_grad), sigma_pot=float(sigma_pot),
                sigma=float(sigma_mixed), first_integral_mismatch=float(sigma_grad - sigma_pot), x_half=xm,
                F_left=float(F[0]), S_left=float(S[0]), F_right=float(F[-1]), S_right=float(S[-1]))


out = {'runs': []}
for (Lx, h) in [(30.0, 0.02), (30.0, 0.01), (30.0, 0.005), (45.0, 0.01)]:
    r = solve(Lx, h); out['runs'].append(r)
    print('Lx=%g h=%.3f sigma=%.8f (grad %.8f pot %.8f) res=%.1e it=%d x_half=%.3f ends F %.1e,%.1e S %.6f,%.10f' % (
        Lx, h, r['sigma'], r['sigma_grad'], r['sigma_pot'], r['newton_res'], r['iters'], r['x_half'], r['F_left'] - Fb, r['F_right'], r['S_left'], r['S_right']), flush=True)
s1, s2, s3 = out['runs'][0]['sigma'], out['runs'][1]['sigma'], out['runs'][2]['sigma']
sig_dimless = (16 * (4 * s3 - s2) / 3 - (4 * s2 - s1) / 3) / 15
# UNITS: the dimensionless energy density is in g^2 v^4 = 100 m_X^4 and lengths in 1/m_X, so sigma is in g^2 v^4 / m_X^2 = 100 m_X^3
UNIT_SIGMA = (1 / 0.1)**4 * G2          # (v/m_X)^4 g^2 = 100
sig = sig_dimless * UNIT_SIGMA
coef = sig * (36 * math.pi)**(1 / 3) * nb**(-2 / 3)
out['sigma_dimensionless_g2v4_over_mX2'] = sig_dimless
out['sigma_mX3'] = sig
out['error_ratio_expect_4'] = (s2 - s1) / (s3 - s2)
out['n_bulk_mX3'] = nb
out['thin_wall_coefficient_of_Q_to_minus_third'] = coef
out['track_fit_coefficient'] = 8.2
out['track_surface_coefficient_sequence'] = dict(Q1e6=9.12, Q1e7=8.66, Q3e7=8.52)
# what the thin-wall law predicts at the track's charges, and the implied Q^{-2/3} term
pred = {}
for Q, EQ in ((2e5, 0.9373999905), (1e6, 0.8670265118), (1e7, 0.8160776689), (3e7, 0.8032854264)):
    pred[str(int(Q))] = dict(E_over_Q=EQ, thin_wall=Omb + coef * Q**(-1 / 3), residual=EQ - Omb - coef * Q**(-1 / 3),
                             residual_times_Q23=(EQ - Omb - coef * Q**(-1 / 3)) * Q**(2 / 3))
out['thin_wall_vs_track_EQ'] = pred
# old symmetric (eps = 0) model for comparison: bulk S = 0, F^2 = 2 sqrt(U_phi(0)/a), Omega^2 = 2 sqrt(a b/4)
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
print('sigma = %.6f m_X^3 (ratio %.2f); n_b = %.4f; coefficient of Q^{-1/3} = %.4f (track fit 8.2)' % (sig, out['error_ratio_expect_4'], nb, coef))
for k, v in pred.items():
    print('  Q=%s E/Q=%.7f thin-wall %.7f residual %.5f -> residual*Q^{2/3} = %.1f' % (k, v['E_over_Q'], v['thin_wall'], v['residual'], v['residual_times_Q23']))
