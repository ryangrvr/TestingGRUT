#!/usr/bin/env python3
"""hw_gates gate 4 (Claude, own code): one-loop running of the scalar couplings at the spot point and the
large-field boundedness bound on kappa.

Method (own implementation of the standard result): write V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d in the
8 real fields (X, Y, Sigma complex; phi, D real), then
  16 pi^2 d lambda_abcd / d ln mu = lambda_abef lambda_efcd + lambda_acef lambda_efbd + lambda_adef lambda_efbc
(pure scalar sector; the U(1) with e <= 1e-2 contributes O(e^2 lambda, e^4): quoted as a bound, not run). The beta
tensor is projected on the closed basis of gauge-invariant, Z2_phi-even, real-coupling quartics (16 operators).
Control: a lone complex field gives beta(lambda) = 20 lambda^2/(16 pi^2). Flow from mu0 = M_Y up and down.
Bound [D, ls_scaling]: along the aligned large-field direction (Y held by its quartic, Y = (kap/2 lam_Y)^{1/3} X)
  lambda_eff = lam_X - (3/2) kap^{4/3} (2 lam_Y)^{-1/3} > 0  <=>  kap < kap_max = (lam_X/1.1906)^{3/4} lam_Y^{1/4}.
Evaluated with the running couplings at mu = M_Y, |Sigma|_inf, |X|_inf, and beyond; also Young's sufficient condition
and the numerical minimum of V4 on the unit sphere.
Usage: python3 -B gate4_running.py <out.json>
"""
import itertools, json, sys, time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spotmodel as S

P = S.params(); H = S.H_SPOT
NAMES = ['lam_X', 'lam_Y', 'lam_S', 'kap', 'l_XY', 'l_XS', 'l_YS', 'gX2', 'l_pY', 'l_pS', 'lam_phi', 'hD2', 'lam_D', 'l_DX', 'l_DY', 'l_DS']
START = dict(lam_X=1e-3, lam_Y=1.0, lam_S=1.0, kap=P['kap'], l_XY=0., l_XS=0., l_YS=0., gX2=0.01, l_pY=0., l_pS=0.,
             lam_phi=0.03, hD2=1e-4, lam_D=1e-3, l_DX=0., l_DY=0., l_DS=0.)


def basis():
    f = sp.symbols('x1 x2 y1 y2 s1 s2 p d', real=True)
    X = (f[0] + sp.I * f[1]) / sp.sqrt(2); Y = (f[2] + sp.I * f[3]) / sp.sqrt(2); Sg = (f[4] + sp.I * f[5]) / sp.sqrt(2)
    aX, aY, aS = [sp.expand(z * sp.conjugate(z)) for z in (X, Y, Sg)]
    p, d = f[6], f[7]
    mono = dict(lam_X=aX**2, lam_Y=aY**2, lam_S=aS**2, kap=2 * sp.re(sp.expand(sp.conjugate(Y) * X**3)),
                l_XY=aX * aY, l_XS=aX * aS, l_YS=aY * aS, gX2=p**2 * aX, l_pY=p**2 * aY, l_pS=p**2 * aS,
                lam_phi=p**4 / 4, hD2=p**2 * d**2 / 2, lam_D=d**4 / 4, l_DX=d**2 * aX, l_DY=d**2 * aY, l_DS=d**2 * aS)
    T = np.zeros((len(NAMES), 8, 8, 8, 8))
    for k, n in enumerate(NAMES):
        e = sp.expand(mono[n])
        for idx in itertools.combinations_with_replacement(range(8), 4):
            val = float(sp.diff(e, *[f[i] for i in idx]))
            if val != 0.0:
                for perm in set(itertools.permutations(idx)):
                    T[(k,) + perm] = val
    return T


T = basis()
Tflat = T.reshape(len(NAMES), -1).T


def beta(g, resid=False):
    L = np.tensordot(g, T, axes=1)
    B = np.einsum('abef,efcd->abcd', L, L)
    B = (B + B.transpose(0, 2, 1, 3) + B.transpose(0, 3, 2, 1)) / (16 * np.pi**2)
    sol, *_ = np.linalg.lstsq(Tflat, B.ravel(), rcond=None)
    if resid:
        return sol, float(np.linalg.norm(Tflat @ sol - B.ravel()) / np.linalg.norm(B))
    return sol


def V4min(g, n=80, seed=1):
    L = np.tensordot(g, T, axes=1); rng = np.random.default_rng(seed)
    f = lambda v: np.einsum('abcd,a,b,c,d->', L, v, v, v, v) / 24 / np.dot(v, v)**2
    return float(min(minimize(f, rng.normal(size=8), method='BFGS').fun for _ in range(n)))


def kap_max(lx, ly):
    return (lx / (1.5 * 2**(-1 / 3)))**0.75 * ly**0.25


def evaluate(g):
    c = dict(zip(NAMES, g))
    lx, ly, k = c['lam_X'], c['lam_Y'], c['kap']
    lam_eff = lx - 1.5 * k**(4 / 3) * (2 * ly)**(-1 / 3) + c['l_XY'] * (k / (2 * ly))**(2 / 3)   # portal adds +l_XY |X|^2|Y|^2 along the direction
    return dict(couplings=c, kappa_max=kap_max(lx, ly), kappa_over_max=k / kap_max(lx, ly), lambda_eff_aligned=lam_eff,
                young_ratio=ly * lx**3 / 27 / k**4)


def main(out):
    t0 = time.time()
    res = dict(method=__doc__)
    g = np.zeros(len(NAMES)); g[0] = 0.1
    b, r = beta(g, True)
    res['control_lone_complex'] = dict(beta_lamX=b[0], expected_20lam2_over_16pi2=20 * 0.01 / (16 * np.pi**2), resid=r)
    g0 = np.array([START[n] for n in NAMES])
    b0, r0 = beta(g0, True)
    res['beta_at_M_Y'] = dict(zip(NAMES, b0.tolist())); res['projection_residual'] = r0
    print('control', res['control_lone_complex'], 'resid', r0, flush=True)
    print('beta0', {k: '%.3e' % v for k, v in res['beta_at_M_Y'].items()}, flush=True)
    mu0 = P['M'] * H
    scales = dict(M_Y=mu0, absSigma_inf=28.31 * H, absX_inf=80.2 * H, three_absX=3 * 80.2 * H, mu_1e14=1e14, mu_1e15=1e15, mu_1e16=1e16,
                  H_inf=H, mu_1e11=1e11, mu_1e9=1e9)
    ev = lambda t, y: np.max(np.abs(y)) - 4 * np.pi
    ev.terminal = True
    up = solve_ivp(lambda t, y: beta(y), (0, np.log(1e17 / mu0)), g0, method='RK45', rtol=1e-9, atol=1e-14, dense_output=True, events=ev)
    dn = solve_ivp(lambda t, y: beta(y), (0, np.log(1e9 / mu0)), g0, method='RK45', rtol=1e-9, atol=1e-14, dense_output=True)
    pole = float(up.t_events[0][0]) if len(up.t_events[0]) else None
    res['landau_pole_GeV'] = mu0 * np.exp(pole) if pole else None
    rows = {}
    for name, s in scales.items():
        t = np.log(s / mu0)
        if t >= 0:
            if pole is not None and t > pole:
                continue
            gg = up.sol(t)
        else:
            gg = dn.sol(t)
        e = evaluate(gg); e['mu_GeV'] = s; e['V4_min_unit_sphere'] = V4min(gg, n=40)
        rows[name] = e
        c = e['couplings']
        print('%-12s mu=%.3e lamX=%.5e lamY=%.4f lamS=%.4f kap=%.5e lXY=%.2e lXS=%.2e lYS=%.2e kap/max=%.4f lam_eff=%.3e young=%.2f V4min=%.2e'
              % (name, s, c['lam_X'], c['lam_Y'], c['lam_S'], c['kap'], c['l_XY'], c['l_XS'], c['l_YS'], e['kappa_over_max'], e['lambda_eff_aligned'], e['young_ratio'], e['V4_min_unit_sphere']), flush=True)
    res['flow_from_M_Y'] = rows
    # alternative matching: couplings stated at |X|_inf, run down to M_Y
    t_inf = np.log(80.2)
    dn2 = solve_ivp(lambda t, y: beta(y), (0, -t_inf), g0, method='RK45', rtol=1e-9, atol=1e-14)
    e2 = evaluate(dn2.y[:, -1]); e2['note'] = 'couplings stated at |X|_inf = 1.07e13 GeV and run down to M_Y'
    res['stated_at_Xinf_run_down_to_MY'] = e2
    print('stated at X_inf, at M_Y: kap/max = %.4f lam_eff = %.3e' % (e2['kappa_over_max'], e2['lambda_eff_aligned']), flush=True)
    res['gauge_contribution_bound'] = dict(note='U(1) terms in the scalar betas are O(q^2 e^2 lambda/16pi^2) and O(q^4 e^4/16pi^2); for e <= 1e-2 and q <= 6 this is <= 2e-3 relative on lam_S over one decade, negligible for lam_X, lam_Y, kappa')
    res['landau_pole_lamY_analytic_GeV'] = mu0 * np.exp((1 / 1.0 - 1 / (4 * np.pi)) * 16 * np.pi**2 / 20)
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print('pole', res['landau_pole_GeV'], 'analytic', res['landau_pole_lamY_analytic_GeV'], 'done', time.time() - t0)


if __name__ == '__main__':
    main(sys.argv[1])
