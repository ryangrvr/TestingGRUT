#!/usr/bin/env python3
"""Independent one-loop estimate of the homogeneous charge-energy floor m*^2.

Candidate 1 benchmark, units v=1:
  V0 = (lp/4)(phi^2-1)^2 + (c/2) phi^2 R^2 + (lX/4) R^4 + (1/2)(mD0^2 + d phi^2) D^2 + (lD/4) D^4
  X = (x+iy)/sqrt2, R^2 = x^2+y^2 = 2|X|^2, c=gX^2, d=hD^2.
MSbar, mu_R = mX = gX v = 0.1.

(1) Static Coleman-Weinberg (derivative-expansion / constant-field) estimate
    m*_eff^2 = min_{phi,R} 2[Veff(phi,R) - Veff(vac)]/R^2.
(2) Homogeneous rotating Q-matter at chemical potential mu (no derivative expansion needed
    for a homogeneous phase): one-loop grand potential with gyroscopic (x,y) mixing,
    m*_eff = mu at which the one-loop pressure vanishes.
Neither is a quantum energy inequality.
"""
import json, sys
import numpy as np
from scipy.optimize import minimize

gX, lp, lX = 0.1, 0.01, 0.001
c = gX**2
mD0sq, d, lD = gX**2, 0.01**2, 0.001
muR = gX           # MSbar scale = mX = mu0
K = 1.0/(64*np.pi**2)

def V0(phi, R, D=0.0):
    return lp/4*(phi**2-1)**2 + c/2*phi**2*R**2 + lX/4*R**4 + 0.5*(mD0sq + d*phi**2)*D**2 + lD/4*D**4

def hess(phi, R, D=0.0):
    """Hessian in (phi, x, y, D) at x=R, y=0."""
    H = np.zeros((4, 4))
    H[0, 0] = lp*(3*phi**2-1) + c*R**2 + d*D**2
    H[1, 1] = c*phi**2 + 3*lX*R**2
    H[2, 2] = c*phi**2 + lX*R**2
    H[3, 3] = mD0sq + d*phi**2 + 3*lD*D**2
    H[0, 1] = H[1, 0] = 2*c*phi*R
    H[0, 3] = H[3, 0] = 2*d*phi*D
    return H

def cw(h):
    h = np.asarray(h, dtype=float)
    out = 0.0
    for hi in h:
        if abs(hi) < 1e-300:
            continue
        out += K*hi**2*(np.log(abs(hi)/muR**2) - 1.5)   # real part; negative h flagged separately
    return out

def V1(phi, R, D=0.0):
    return cw(np.linalg.eigvalsh(hess(phi, R, D)))

def neg_eigs(phi, R):
    return np.linalg.eigvalsh(hess(phi, R))[0] < 0

res = {}
# --- vacuum
h_vac = np.linalg.eigvalsh(hess(1.0, 0.0))
V1vac_tree = V1(1.0, 0.0)
res['vacuum_hessian_eigs_over_v2'] = h_vac.tolist()
res['strict_one_loop_vacuum_energy_over_v4'] = V1vac_tree
res['quantum_REPORT_value_over_v4'] = -1.22652632365e-6
# minimized vacuum of V0+V1 along phi (X=D=0)
from scipy.optimize import minimize_scalar
mv = minimize_scalar(lambda p: V0(p, 0.0) + V1(p, 0.0), bracket=(0.99, 1.0, 1.01), tol=1e-14)
res['truncated_vacuum_phi_over_v'] = mv.x
res['truncated_vacuum_energy_over_v4'] = mv.fun

# --- tree floor
mstar2_tree = np.sqrt(lp*lX)
R2star = np.sqrt(lp/lX)    # F^2 = sqrt(a/b) with a=lp/gX^2, b=lX/gX^2 ; F=R/v
Rstar = np.sqrt(R2star)
res['tree_mstar2_over_v2'] = mstar2_tree
res['tree_minimizer'] = {'phi': 0.0, 'R_over_v': Rstar}
res['hessian_eigs_at_tree_minimizer'] = np.linalg.eigvalsh(hess(0.0, Rstar)).tolist()

# --- (1a) strict first order (envelope theorem): delta m*^2 = 2[V1(p*) - V1(vac)]/R*^2
dm2_strict = 2*(V1(0.0, Rstar) - V1vac_tree)/R2star
res['static_CW_strict_delta_mstar2_over_v2'] = dm2_strict
res['static_CW_strict_relative_delta_mstar2'] = dm2_strict/mstar2_tree
res['quantum_REPORT_strict_delta_over_v2'] = 1.17872091832e-7

# --- (1b) full numerical minimization of the truncated ratio with the truncated vacuum
def ratio(p, Vvac):
    phi, R = p
    return 2*(V0(phi, R) + V1(phi, R) - Vvac)/R**2
best = None
for start in [(0.0, Rstar), (0.05, 1.7), (0.1, 1.9), (0.0, 1.5), (0.0, 2.1)]:
    r = minimize(ratio, start, args=(mv.fun,), method='Nelder-Mead',
                 options={'xatol': 1e-12, 'fatol': 1e-18, 'maxiter': 20000})
    if best is None or r.fun < best.fun:
        best = r
res['static_CW_numerical_min_mstar2_over_v2'] = best.fun
res['static_CW_numerical_minimizer'] = {'phi': abs(best.x[0]), 'R_over_v': best.x[1]}
res['static_CW_numerical_relative_delta'] = best.fun/mstar2_tree - 1
res['negative_eig_at_numerical_minimizer'] = bool(neg_eigs(*best.x))

# --- negative mass-squared map along the tree ratio-minimizing path and on a grid
path = []
for phi in np.linspace(0, 0.999, 1000):
    R2 = np.sqrt(lp/lX)*(1-phi**2)  # optimal y at fixed s for s<1
    R = np.sqrt(R2)
    path.append((phi, R, np.linalg.eigvalsh(hess(phi, R))[0]))
negpath = [p for p in path if p[2] < 0]
res['tree_path_min_eig_over_v2'] = min(p[2] for p in path)
res['tree_path_negative_points'] = len(negpath)
res['tree_path_negative_phi_range'] = [negpath[0][0], negpath[-1][0]] if negpath else None
grid_phi = np.linspace(0, 1.2, 241); grid_R = np.linspace(0.0, 3.0, 301)
neg = [(p, R) for p in grid_phi for R in grid_R if neg_eigs(p, R)]
res['grid_negative_fraction'] = len(neg)/(len(grid_phi)*len(grid_R))
res['grid_negative_R_max'] = max(R for p, R in neg) if neg else None
res['grid_negative_phi_max'] = max(p for p, R in neg) if neg else None
# ratio value inside the negative region (to show it is far above the minimum)
if neg:
    res['min_tree_ratio_inside_negative_region_over_v2'] = min(2*V0(p, R)/R**2 for p, R in neg if R > 0)

# --- (2) rotating homogeneous Q-matter at chemical potential mu = m* (tree point phi=0, R*)
from scipy.integrate import quad
mu = np.sqrt(mstar2_tree)
hx = 3*lX*R2star          # Hess V0 eigenvalues at (phi=0,R*) in the lab basis
hy = lX*R2star
a = hx - mu**2            # rotating-frame curvature M = HessV0 - mu^2 P_xy
b = max(hy - mu**2, 0.0)   # = 0 (Goldstone); clip -4e-19 roundoff
def diff_stable(k, a, b, mu2):
    """omega_+ + omega_- - sqrt(k^2+a+mu2) - sqrt(k^2+b+mu2), cancellation-free.
    Uses (w+ + w-)^2 = A+B+4mu^2+2sqrt(AB), A=k^2+a, B=k^2+b."""
    A, B = k*k + a, k*k + b
    Ap, Bp = A + mu2, B + mu2
    sw = np.sqrt(A + B + 4*mu2 + 2*np.sqrt(A*B))
    T = np.sqrt(Ap) + np.sqrt(Bp)
    N = -(a-b)**2/2*(1/(np.sqrt(A)+np.sqrt(B))**2 + 1/(np.sqrt(Ap)+np.sqrt(Bp))**2)
    return 2*mu2*N/(np.sqrt(A*B)+np.sqrt(Ap*Bp))/(sw + T)
def diff_direct(k, a, b, mu2):
    A, B = k*k + a, k*k + b
    return np.sqrt(A + B + 4*mu2 + 2*np.sqrt(A*B)) - np.sqrt(k*k+a+mu2) - np.sqrt(k*k+b+mu2)
chk = [(kk, diff_stable(kk, a, b, mu**2), diff_direct(kk, a, b, mu**2)) for kk in (1e-6, 0.01, 0.05, 0.2, 1.0)]
uv = [kk**5*diff_stable(kk, a, b, mu**2) for kk in (1e1, 1e2, 1e3, 1e4)]
f = lambda k: k*k*diff_stable(k, a, b, mu**2)
pieces = [quad(f, 0, mu, epsabs=0, epsrel=1e-13, limit=400)[0],
          quad(f, mu, 1, epsabs=0, epsrel=1e-13, limit=400)[0],
          quad(f, 1, np.inf, epsabs=0, epsrel=1e-13, limit=400)[0]]
rem = sum(pieces)/(4*np.pi**2)
hphi = lp*(-1) + c*R2star
hD = mD0sq
Omega1 = cw([hphi, hx, hy, hD]) + rem
dOmega = Omega1 - V1vac_tree
dm2_rot = 2*dOmega/R2star     # delta mu^2 = 2 m* delta mu, delta mu = dOmega/(m* R*^2)
res['rotating_Qmatter'] = {
    'mu_over_v': mu, 'a_over_v2': a, 'b_over_v2': b,
    'stable_vs_direct_difference_check': chk,
    'UV_check_k^5*diff_at_k=1e1..1e4 (-> -mu^2(a-b)^2/16; integrand k^2*diff ~ k^-3)': uv,
    'UV_expected': -mu**2*(a-b)**2/16,
    'finite_remainder_over_v4': rem,
    'static_V1_at_point_over_v4': V1(0.0, Rstar),
    'Omega1_over_v4': Omega1,
    'delta_mstar2_over_v2': dm2_rot,
    'relative_delta_mstar2': dm2_rot/mstar2_tree,
    'relative_delta_mstar': dm2_rot/mstar2_tree/2,
}
# free-particle check of the method: a=b=m^2-mu^2 should give exactly zero remainder
res['method_check_free_complex_scalar_diff'] = [diff_stable(kk, 0.0075, 0.0075, 0.0025) for kk in (0.0, 0.1, 1.0)]

# --- margin comparison (from COMPARISON.json @2ee7f82: product limits at fixed mu)
stress_limit, primary_limit = 74.64290523646463, 63.93439835872732
res['required_fractional_reduction_of_mstar'] = {
    'stress_NGC5765b': 1 - stress_limit/100.0,
    'primary_all_rows(NGC5765b)': 1 - primary_limit/100.0,
    'primary_UGC3789': 1 - 79.05068738513258/100.0}
res['margin_over_static_loop_shift_in_mstar'] = (1 - stress_limit/100.0)/abs(dm2_strict/mstar2_tree/2)
res['margin_over_rotating_loop_shift_in_mstar'] = (1 - stress_limit/100.0)/abs(dm2_rot/mstar2_tree/2)
# generated running m_X^2 at mu/mu0 = 4 (quantum REPORT table), treated (incorrectly) as an uncancelled shift
res['running_mX2_at_4mu0_relative_shift_mstar2'] = -1.75743e-6/mstar2_tree
print(json.dumps(res, indent=1, default=float))
