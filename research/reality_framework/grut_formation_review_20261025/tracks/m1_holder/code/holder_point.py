#!/usr/bin/env python3
"""m1_holder track (Claude, own code) -- parts 1 and 2 of the review of Codex f8186fa (Option 1A curvature holder).

Independent of Codex's model.py: own potential, own starting guess, own Newton/BFGS, Hessian by exact symbolic
differentiation (sympy) of V(r_1..r_6) and by a fourth-order finite-difference stencil, own multi-start global
minimisation, own re-derivation of the fixed-R lower bound.  Codex's numbers are quoted [Q] only for comparison.
Units: fields and H_i = 1, R in H_i^2, V in H_i^4.  Pinned point [Q CANDIDATE_M1.json, d4747cb].
"""
import json, sys, time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import minimize, root

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
RES.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261101)

H_GEV = 134284273491.30815
MY_GEV, H_VEV_GEV, MU_GEV = 4028528204739.244, 5371370939652.326, 961402368973.181
KAP, LX, LY, LS, XIX, MPL = 0.0034586847211042523, 1e-3, 1.0, 1.0, -5.0 / 6.0, 2.435e18
M, HV, MU = MY_GEV / H_GEV, H_VEV_GEV / H_GEV, MU_GEV / H_GEV
ETA2 = (H_GEV / MPL) ** 2
LAMBDA6_GEV, XI6 = 1e15, 4e8          # the 94f660d sextic holder [Q]
CODEX = dict(xi=30.746728738599167, alpha=2.192144074179184, amps=[87.19556213287578, 5.75424796328559, 28.387181944361622],
             V0=-32480.585311406205, hess=[13.091408829453263, 14.589249978666821, 538.1466434275541, 1546.6597219731073, 3239.885251158698],
             bound=-285369.6842985777, B=0.5283522987365723, planck=1.0000000000480667, ddelta_dlnR=0.715945616,
             sextic=dict(amps=[86.657, 3.823, 28.338], delta=-2.09993, lightest=10.669))

# ----------------------------------------------------------------------------------------------------------------
# part 1: the operator
# ----------------------------------------------------------------------------------------------------------------
LOCK = (10j - MU) / 12.0                      # requirement mu_eff(R = 12) = 10 i  ->  mu + 12 lock = 10 i
XI_MUR = abs(LOCK) * M                        # V ⊃ (xi R / M_Y)(e^{i alpha} Sigma* Y^2 + h.c.) = 2 Re[(xi e^{i alpha}/M_Y) R Sigma* Y^2]
ALPHA = float(np.angle(LOCK))
assert abs(XI_MUR * np.exp(1j * ALPHA) / M - LOCK) < 1e-14


def u_of_R(R, lock=LOCK):
    return MU + lock * R


# ----------------------------------------------------------------------------------------------------------------
# the potential (own coding; six real canonical fields r = sqrt2 (Re X, Im X, Re Y, Im Y, Re S, Im S))
# ----------------------------------------------------------------------------------------------------------------
def z_of_r(r):
    return np.array([(r[0] + 1j * r[1]), (r[2] + 1j * r[3]), (r[4] + 1j * r[5])]) / np.sqrt(2)


def r_of_z(z):
    return np.sqrt(2) * np.array([z[0].real, z[0].imag, z[1].real, z[1].imag, z[2].real, z[2].imag])


def V_complex(z, R=12.0, lock=LOCK, c6=0.0):
    x, y, s = z
    u = u_of_R(R, lock)
    V = (LX * abs(x) ** 4 + LY * abs(y) ** 4 + M * M * abs(y) ** 2 + LS * (abs(s) ** 2 - HV * HV / 2) ** 2
         + 2 * np.real(KAP * np.conj(y) * x ** 3 + u * np.conj(s) * y * y) + XIX * R * abs(x) ** 2)
    if c6:                                                     # sextic holder of 94f660d: c6 (i S* X^6 + h.c.)
        V += 2 * np.real(1j * c6 * np.conj(s) * x ** 6)
    return float(V)


def V_real(r, **kw):
    return V_complex(z_of_r(r), **kw)


def grad_real(r, step=1e-4, **kw):
    """fourth-order central differences of V itself (no analytic gradient used)."""
    g = np.zeros(6)
    for i in range(6):
        e = np.zeros(6); e[i] = step
        g[i] = (-V_real(r + 2 * e, **kw) + 8 * V_real(r + e, **kw) - 8 * V_real(r - e, **kw) + V_real(r - 2 * e, **kw)) / (12 * step)
    return g


def hessian_fd(r, step=2e-3, **kw):
    Hm = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            ei = np.zeros(6); ei[i] = step; ej = np.zeros(6); ej[j] = step
            f = lambda a, b: V_real(r + a * ei + b * ej, **kw)
            if i == j:
                d = (-f(2, 0) + 16 * f(1, 0) - 30 * f(0, 0) + 16 * f(-1, 0) - f(-2, 0)) / (12 * step ** 2)
            else:
                d = (f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1)) / (4 * step ** 2)
                d2 = (f(2, 2) - f(2, -2) - f(-2, 2) + f(-2, -2)) / (16 * step ** 2)
                d = (4 * d - d2) / 3          # Richardson
            Hm[i, j] = Hm[j, i] = d
    return Hm


# symbolic Hessian (exact), independent of the A/B block formula used by Codex
_r = sp.symbols('r1:7', real=True)
_x = (_r[0] + sp.I * _r[1]) / sp.sqrt(2); _y = (_r[2] + sp.I * _r[3]) / sp.sqrt(2); _s = (_r[4] + sp.I * _r[5]) / sp.sqrt(2)
_R, _ur, _ui, _c6 = sp.symbols('R u_r u_i c6', real=True)
_u = _ur + sp.I * _ui
_V = (LX * (_x * sp.conjugate(_x)) ** 2 + LY * (_y * sp.conjugate(_y)) ** 2 + M ** 2 * _y * sp.conjugate(_y)
      + LS * (_s * sp.conjugate(_s) - sp.Rational(1, 2) * HV ** 2) ** 2
      + (KAP * sp.conjugate(_y) * _x ** 3 + _u * sp.conjugate(_s) * _y ** 2) + sp.conjugate(KAP * sp.conjugate(_y) * _x ** 3 + _u * sp.conjugate(_s) * _y ** 2)
      + XIX * _R * _x * sp.conjugate(_x)
      + sp.I * _c6 * sp.conjugate(_s) * _x ** 6 - sp.I * _c6 * _s * sp.conjugate(_x) ** 6)
_V = sp.expand(_V)
_hess = sp.hessian(_V, _r)
hess_sym = sp.lambdify((_r, _R, _ur, _ui, _c6), _hess, 'numpy')
grad_sym = sp.lambdify((_r, _R, _ur, _ui, _c6), [sp.diff(_V, v) for v in _r], 'numpy')


def hessian_exact(r, R=12.0, lock=LOCK, c6=0.0):
    u = u_of_R(R, lock)
    return np.array(hess_sym(r, R, u.real, u.imag, c6), float)


def gradient_exact(r, R=12.0, lock=LOCK, c6=0.0):
    u = u_of_R(R, lock)
    return np.array(grad_sym(r, R, u.real, u.imag, c6), float)


# ----------------------------------------------------------------------------------------------------------------
# part 2: stationary point, Hessian, bound
# ----------------------------------------------------------------------------------------------------------------
def stationary(R=12.0, lock=LOCK, c6=0.0, guess=None):
    """Own starting guess: x^2 = -xi R/(2 lx) (curvature mass balanced by the quartic), Y = -kap x^3/M^2 (forced),
    Sigma = -i h/sqrt2 rotated so that u Sigma* Y^2 is real negative.  Then BFGS on V, then Levenberg-Marquardt root
    on the exact gradient (the gauge zero mode makes the Jacobian singular; lm handles it)."""
    if guess is None:
        u = u_of_R(R, lock)
        x0 = np.sqrt(max(-XIX * R, 1.0) / (2 * LX))
        y0 = -KAP * x0 ** 3 / M ** 2
        s0 = -(u / abs(u)) * HV / np.sqrt(2)                   # makes u s* y^2 = -|u| s y^2 < 0
        guess = np.array([x0, y0, s0])
    r0 = r_of_z(np.asarray(guess, complex))
    kw = dict(R=R, lock=lock, c6=c6)
    m = minimize(lambda r: V_real(r, **kw), r0, jac=lambda r: gradient_exact(r, **kw), method='BFGS', options=dict(gtol=1e-9, maxiter=2000))
    sol = root(lambda r: gradient_exact(r, **kw), m.x, method='lm', options=dict(xtol=1e-15, ftol=1e-15, maxiter=20000))
    r = sol.x
    # gauge-fix: rotate so that X is real positive (charges 1,3,6)
    z = z_of_r(r); th = np.angle(z[0]); z = z * np.exp(-1j * th * np.array([1, 3, 6]))
    r = r_of_z(z)
    return z, r, float(np.linalg.norm(gradient_exact(r, **kw))), float(np.linalg.norm(grad_real(r, **kw)))


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def lower_bound_codex_chain(R, lock=LOCK, reserve=0.2):
    """Codex's inequality chain re-derived (NOTES.md section 2); returns the bound and the pieces."""
    A = reserve * LX
    c = (1 - reserve) * LX
    B = LY - 27 * KAP ** 4 / (16 * c ** 3)                    # min_x [c x^4 - 2|kap| y x^3] = -27 kap^4 y^4 /(16 c^3)
    u = abs(u_of_R(R, lock)); ax = max(-XIX * R, 0.0)
    bound = LS * HV ** 4 / 4 - ax ** 2 / (4 * A) - (LS * HV ** 2 + u * u / B) ** 2 / (4 * LS)
    return float(bound), float(A), float(B), float(u), float(ax)


def global_search(R=12.0, lock=LOCK, nstart=600, box=200.0):
    """multi-start L-BFGS-B over the six real fields; log-uniform radii up to `box`, random phases."""
    best = []
    t0 = time.time()
    kw = dict(R=R, lock=lock)
    for n in range(nstart):
        amps = np.exp(rng.uniform(np.log(0.3), np.log(box), 3)); ph = rng.uniform(-np.pi, np.pi, 3)
        z0 = amps * np.exp(1j * ph)
        r0 = r_of_z(z0)
        m = minimize(lambda r: V_real(r, **kw), r0, jac=lambda r: gradient_exact(r, **kw), method='L-BFGS-B',
                     options=dict(maxiter=4000, ftol=1e-15, gtol=1e-8))
        z = z_of_r(m.x); th = np.angle(z[0]); z = z * np.exp(-1j * th * np.array([1, 3, 6]))
        best.append((float(m.fun), [abs(v) for v in z], float(wrap(np.angle(z[2]) - 6 * np.angle(z[0]))), float(wrap(np.angle(z[1]) - 3 * np.angle(z[0])))))
    best.sort(key=lambda t: t[0])
    # distinct minima (by V to 1e-6 relative)
    distinct = []
    for b in best:
        if not any(abs(b[0] - d[0]) <= 1e-6 * max(1, abs(d[0])) for d in distinct):
            distinct.append(b)
    return best[0], distinct[:12], time.time() - t0


def main():
    out = {}
    # ---------------- part 1
    out['part1_operator'] = dict(
        mu_over_H=MU, M_Y_over_H=M, h_over_H=HV, kappa=KAP, eta2=ETA2,
        lock_over_H=[LOCK.real, LOCK.imag], xi_muR_own=XI_MUR, alpha_R_own=ALPHA,
        xi_codex=CODEX['xi'], alpha_codex=CODEX['alpha'],
        xi_rel_diff=XI_MUR / CODEX['xi'] - 1, alpha_diff=ALPHA - CODEX['alpha'],
        mu_eff_at_R12=[u_of_R(12).real, u_of_R(12).imag], mu_eff_at_R0=[u_of_R(0).real, u_of_R(0).imag],
        mu_eff_at_R3_matter=[u_of_R(3).real, u_of_R(3).imag],
        holder_coefficient_at_R12_over_H=abs(LOCK) * 12, ratio_holder_to_tree_mu_at_R12=abs(LOCK) * 12 / MU,
        xi_R_over_MY2_at_R12=XI_MUR * 12 / M ** 2, R_over_MY2=12 / M ** 2,
        dim7_estimate_same_xi_over_dim5=12 / M ** 2,
        loop_size_estimates=dict(
            one_loop_Y_loop_xi=LY * MU / (16 * np.pi ** 2 * M),       # c * lambda_Y mu/(16 pi^2 M_Y) with c = O(1)
            graviton_MY_over_MPl_sq=(MY_GEV / MPL) ** 2,
            required_over_one_loop=XI_MUR / (LY * MU / (16 * np.pi ** 2 * M)),
            tree_singlet_exchange='xi = g1 g2 M_Y/M_S for L ⊃ g1 M_S S R + g2 S Sigma* Y^2 (+h.c.)'),
        d_delta_eq_d_lnR_own=float(12 * np.imag(LOCK / u_of_R(12))), d_delta_eq_d_lnR_codex=CODEX['ddelta_dlnR'],
        delta_eq_R12=float(wrap(np.angle(u_of_R(12)) + np.pi)), delta_eq_R0=float(wrap(np.angle(u_of_R(0)) + np.pi)),
        delta_eq_R3=float(wrap(np.angle(u_of_R(3)) + np.pi)))

    # ---------------- part 2: stationary point
    z, r, gnorm, gnorm_fd = stationary()
    amps = np.abs(z); delta = float(wrap(np.angle(z[2]) - 6 * np.angle(z[0]))); dy = float(wrap(np.angle(z[1]) - 3 * np.angle(z[0])))
    V0 = V_real(r)
    Hs = hessian_exact(r); Hf = hessian_fd(r)
    ev, U = np.linalg.eigh(Hs); evf = np.linalg.eigvalsh(Hf)
    # gauge direction check: T r should be the null vector
    T = np.zeros((6, 6))
    for i, q in enumerate([1, 3, 6]):
        T[2 * i, 2 * i + 1] = -q; T[2 * i + 1, 2 * i] = q
    gvec = T @ r; gvec /= np.linalg.norm(gvec)
    gauge_residual = float(np.linalg.norm(Hs @ gvec))
    # name the modes
    names = ['ReX', 'ImX', 'ReY', 'ImY', 'ReS', 'ImS']
    modes = []
    for k in range(6):
        comp = U[:, k]; idx = np.argsort(-abs(comp))[:3]
        modes.append(dict(m2=float(ev[k]), dominant=[(names[i], round(float(comp[i]), 4)) for i in idx]))
    VR = XIX * amps[0] ** 2 + 2 * np.real(LOCK * np.conj(z[2]) * z[1] ** 2)
    planck = 1 - 2 * ETA2 * VR
    # local minimum? all physical eigenvalues positive
    phys = ev[np.abs(ev) > 1e-6]
    out['part2_stationary_point'] = dict(
        fields_over_H=[[v.real, v.imag] for v in z], amplitudes_own=amps.tolist(), amplitudes_codex=CODEX['amps'],
        amplitude_rel_diff=(amps / np.array(CODEX['amps']) - 1).tolist(),
        phase_sigma_minus_6x=delta, phase_y_minus_3x=dy, phase_minus_pi_over_2_residual=delta + np.pi / 2,
        gradient_norm_exact=gnorm, gradient_norm_fd=gnorm_fd, V_own=V0, V_codex=CODEX['V0'], V_rel_diff=V0 / CODEX['V0'] - 1,
        hessian_eigenvalues_sympy=ev.tolist(), hessian_eigenvalues_fd=evf.tolist(),
        physical_eigenvalues_own=phys.tolist(), physical_eigenvalues_codex=CODEX['hess'],
        physical_rel_diff=(phys / np.array(CODEX['hess']) - 1).tolist(), gauge_null_vector_residual=gauge_residual,
        modes=modes, all_physical_positive=bool(np.all(phys > 0)),
        dV_dR_own=float(VR), effective_Planck_ratio_own=float(planck), effective_Planck_ratio_codex=CODEX['planck'],
        X_over_MY=float(amps[0] / M), Y_over_MY=float(amps[1] / M), S_over_MY=float(amps[2] / M),
        kappa_X6_over_MY2_vs_lambdaX_X4=float(KAP ** 2 * amps[0] ** 6 / M ** 2 / (LX * amps[0] ** 4)))
    # hand check of the three amplitude equations at delta = -pi/2 (X real, Y real negative, Sigma = -i s)
    x, y, s = amps
    hand = dict(X_eq=2 * LX * x ** 3 - 3 * KAP * y * x ** 2 + XIX * 12 * x,
                Y_eq=-(M * M + 2 * LY * y * y) * y + KAP * x ** 3 + 2 * abs(u_of_R(12)) * s * y,
                S_eq=2 * LS * (s * s - HV * HV / 2) * s - abs(u_of_R(12)) * y * y)
    out['part2_stationary_point']['hand_equations_residual'] = hand
    np.save(RES / 'physical_basis.npy', np.c_[ev, U.T])         # rows: eigenvalue, eigenvector
    json.dump(dict(eigenvalues=ev.tolist(), eigenvectors_columns=U.tolist(), r_stationary=r.tolist(), z_stationary=[[v.real, v.imag] for v in z]),
              open(RES / 'stationary_point.json', 'w'), indent=1)

    # ---------------- part 2: bound
    bound, A, B, uabs, ax = lower_bound_codex_chain(12.0)
    # the pieces of the bound: check each inequality numerically at random points
    worst = 0.0
    for _ in range(20000):
        yv = rng.uniform(0, 60); xv = rng.uniform(0, 300)
        c = 0.8 * LX
        worst = min(worst, c * xv ** 4 - 2 * KAP * yv * xv ** 3 + 27 * KAP ** 4 * yv ** 4 / (16 * c ** 3))
    # direct numerical check of V >= bound on random points (large box)
    vmin_random = np.inf
    for _ in range(200000):
        amps_r = np.exp(rng.uniform(np.log(0.1), np.log(400), 3)); ph = rng.uniform(-np.pi, np.pi, 3)
        vmin_random = min(vmin_random, V_complex(amps_r * np.exp(1j * ph)))
    best, distinct, secs = global_search()
    out['part2_bound'] = dict(
        codex_chain_bound_own=bound, codex_bound=CODEX['bound'], A=A, B_own=B, B_codex=CODEX['B'], abs_u=uabs, a_x=ax,
        mixed_quartic_inequality_min_violation=worst, random_200k_min_V=float(vmin_random),
        multistart_best=dict(V=best[0], amplitudes=best[1], delta=best[2], delta_y=best[3]),
        multistart_distinct_minima=[dict(V=d[0], amplitudes=d[1], delta=d[2], delta_y=d[3]) for d in distinct],
        multistart_seconds=secs, bounds_at_other_R={str(R): lower_bound_codex_chain(R)[0] for R in (-1.0633, 0.0, 3.0, 6.0, 12.0)},
        comment='bound is valid but loose: the prepared point is the lowest minimum found by the multi-start search')

    # ---------------- 94f660d sextic holder cross-check [C]
    lam6 = LAMBDA6_GEV / H_GEV
    c6 = XI6 * 12.0 / lam6 ** 5
    z6, r6, g6, _ = stationary(lock=0j, c6=c6, guess=np.array([80.0, -3.4, -28.3 * 1j]))
    ev6 = np.linalg.eigvalsh(hessian_exact(r6, lock=0j, c6=c6))
    out['sextic_94f660d'] = dict(c6_over_H5=float(c6), amplitudes_own=np.abs(z6).tolist(), amplitudes_codex=CODEX['sextic']['amps'],
                                 delta_own=float(wrap(np.angle(z6[2]) - 6 * np.angle(z6[0]))), delta_codex=CODEX['sextic']['delta'],
                                 lightest_physical_own=float(np.sort(ev6)[1]), lightest_codex=CODEX['sextic']['lightest'],
                                 eigenvalues_own=ev6.tolist(), gradient_norm=g6,
                                 holder_term_at_point_over_H4=float(2 * np.real(1j * c6 * np.conj(z6[2]) * z6[0] ** 6)))
    json.dump(out, open(RES / 'holder_point.json', 'w'), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
