#!/usr/bin/env python3
"""One-loop running of the dimensionless scalar couplings and a Coleman-Weinberg estimate at the
inflationary point (own code; estimate only).

Method: write the quartic part of the potential in 8 real fields (X, Y, Sigma as 2 real components each,
phi, D), extract the constant tensor lambda_abcd = d^4 V4, and use the general one-loop result
  16 pi^2 beta_abcd = lambda_abef lambda_efcd + lambda_acef lambda_efbd + lambda_adef lambda_efbc
(pure scalar sector; the U(1) gauge terms with e=1e-7 are O(1e-13) and dropped; SM fields are not coupled).
beta tensors are projected on a closed basis of gauge-invariant, Z2(phi), real-coupling operators,
including the portals that 7737487 sets to zero at tree level. Control: beta(lambda_X) = 20 lambda_X^2/16pi^2
for a lone complex field.
Then: (a) run from mu0 = M_Y (assumed matching scale; stated point has none) and from mu0 = H_inf;
(b) check Young's sufficient condition and the numerical positivity of V4 on the unit sphere along the flow;
(c) RG-improved tree point and a flat-space CW one-loop correction at the inflationary point (H units).
Outputs results/running.json.
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import build, GEV, COUP

OUT = Path(__file__).resolve().parents[1] / 'results'
NAMES = ['lam_X', 'lam_Y', 'lam_S', 'kap', 'l_XY', 'l_XS', 'l_YS', 'gX2', 'l_pY', 'l_pS', 'lam_phi',
         'hD2', 'lam_D', 'l_DX', 'l_DY', 'l_DS']
START = dict(lam_X=1e-3, lam_Y=1.0, lam_S=1.0, kap=1e-3, l_XY=0., l_XS=0., l_YS=0., gX2=0.01, l_pY=0., l_pS=0.,
             lam_phi=0.03, hD2=1e-4, lam_D=1e-3, l_DX=0., l_DY=0., l_DS=0.)


def basis():
    f = sp.symbols('x1 x2 y1 y2 s1 s2 p d', real=True)
    X = (f[0] + sp.I * f[1]) / sp.sqrt(2); Y = (f[2] + sp.I * f[3]) / sp.sqrt(2); S = (f[4] + sp.I * f[5]) / sp.sqrt(2)
    aX, aY, aS = [sp.expand(z * sp.conjugate(z)) for z in (X, Y, S)]
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
    return T, f, mono


T, FS, MONO = basis()
Tflat = T.reshape(len(NAMES), -1).T


def lam_tensor(g):
    return np.tensordot(g, T, axes=1)


def beta_vec(g, return_resid=False):
    L = lam_tensor(g)
    B = np.einsum('abef,efcd->abcd', L, L)
    B = B + B.transpose(0, 2, 1, 3) + B.transpose(0, 3, 2, 1)
    B = B / (16 * np.pi**2)
    sol, res, *_ = np.linalg.lstsq(Tflat, B.ravel(), rcond=None)
    if return_resid:
        return sol, float(np.linalg.norm(Tflat @ sol - B.ravel()) / max(np.linalg.norm(B), 1e-300))
    return sol


def V4_min_on_sphere(g, n=300, seed=3):
    L = lam_tensor(g)
    rng = np.random.default_rng(seed)
    f = lambda v: np.einsum('abcd,a,b,c,d->', L, v, v, v, v) / 24 / np.dot(v, v)**2
    best = np.inf
    for _ in range(n):
        r = minimize(f, rng.normal(size=8), method='BFGS')
        best = min(best, r.fun)
    return float(best)


def run_flow(g0, t_end, t_eval):
    ev = lambda t, g: np.max(np.abs(g)) - 4 * np.pi
    ev.terminal = True
    sol = solve_ivp(lambda t, g: beta_vec(g), (0, t_end), g0, method='RK45', rtol=1e-9, atol=1e-14,
                    t_eval=[t for t in t_eval if t <= t_end], events=ev)
    pole = float(sol.t_events[0][0]) if len(sol.t_events[0]) else None
    return sol, pole


def cw_point(coup, extra, mu_R, c_H=1.0):
    """Tree + flat-space CW one-loop (H units). Returns minimum moduli and angular masses."""
    m = build(c_H=c_H, coup=coup, extra_quartic=extra)
    gX2 = 0.01

    def V1(z):
        ev = np.linalg.eigvalsh(m['hess'](z))
        mphi2 = gX2 * (z[0]**2 + z[1]**2)          # 2 g^2 |X|^2 with |X|^2=(x1^2+x2^2)/2
        allm = np.append(ev, mphi2)
        out = 0.0
        for m2 in allm:
            if abs(m2) > 1e-14:
                out += m2**2 * (np.log(abs(m2) / mu_R**2) - 1.5)
        return out / (64 * np.pi**2)
    tot = lambda z: m['V'](z) + V1(z)
    z0 = np.array([np.sqrt(2) * 23.8, 0, -np.sqrt(2) * 1.89, 0, -np.sqrt(2) * 0.267, 0])
    # minimise over the real (aligned) slice; CP keeps the minimum there
    fr = lambda r: tot(np.array([r[0], 0, r[1], 0, r[2], 0]))
    rt = minimize(lambda r: m['V'](np.array([r[0], 0, r[1], 0, r[2], 0])), z0[[0, 2, 4]], method='Nelder-Mead',
                  options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000))
    r1 = minimize(fr, rt.x, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000))
    out = {}
    for lab, r, fun in [('tree', rt.x, lambda z: m['V'](z)), ('tree+CW', r1.x, tot)]:
        z = np.array([r[0], 0, r[1], 0, r[2], 0]); eps = 1e-3
        Hm = np.zeros((6, 6))
        for i in range(6):
            for j in range(i, 6):
                ei = np.zeros(6); ej = np.zeros(6); ei[i] = eps; ej[j] = eps
                Hm[i, j] = Hm[j, i] = (fun(z + ei + ej) - fun(z + ei - ej) - fun(z - ei + ej) + fun(z - ei - ej)) / (4 * eps**2)
        ang = np.linalg.eigvalsh(Hm[np.ix_([1, 3, 5], [1, 3, 5])])
        rad = np.linalg.eigvalsh(Hm[np.ix_([0, 2, 4], [0, 2, 4])])
        out[lab] = dict(moduli_over_H=(np.abs(r) / np.sqrt(2)).tolist(), angular_mass2=ang.tolist(), radial_mass2=rad.tolist(),
                        V=float(fun(z)))
    out['V1_at_tree_point'] = float(V1(np.array([rt.x[0], 0, rt.x[1], 0, rt.x[2], 0])))
    return out


def main():
    res = {}
    # control
    g = np.zeros(len(NAMES)); g[0] = 0.1
    b, resid = beta_vec(g, True)
    res['control_lone_complex'] = dict(beta_lamX=b[0], expected=20 * 0.01 / (16 * np.pi**2), resid=resid)
    g0 = np.array([START[n] for n in NAMES])
    b0, resid0 = beta_vec(g0, True)
    res['beta_at_start'] = dict(zip(NAMES, b0.tolist())); res['projection_residual'] = resid0
    print('control', res['control_lone_complex'], 'resid', resid0)
    print('beta0', {k: f'{v:.3e}' for k, v in res['beta_at_start'].items()})
    scales = dict(M_Y=1e12, H_inf=4.7e13, m_Y_point=4.65 * 4.7e13, absX_point=1.12e15, T_R=5.8e15, M_GUT=1e16)
    flows = {}
    for lab, mu0 in [('mu0=M_Y', 1e12), ('mu0=H_inf', 4.7e13)]:
        t_eval = sorted(set([np.log(s / mu0) for s in scales.values() if s >= mu0] + list(np.linspace(0, 12, 121))))
        sol, pole = run_flow(g0, 12.0, t_eval)
        rows = {}
        for k, s in scales.items():
            if s < mu0:
                continue
            t = np.log(s / mu0)
            j = np.where(np.isclose(sol.t, t))[0]
            if not len(j):
                continue
            gg = sol.y[:, j[0]]
            rows[k] = dict(mu_GeV=s, couplings=dict(zip(NAMES, gg.tolist())),
                           young_ratio=float(gg[1] * gg[0]**3 / 27 / gg[3]**4),
                           V4_min_unit_sphere=V4_min_on_sphere(gg, n=60))
        flows[lab] = dict(rows=rows, nonperturbative_t=pole,
                          nonperturbative_scale_GeV=(mu0 * np.exp(pole) if pole else None))
        print(lab, 'max|g|=4pi at', flows[lab]['nonperturbative_scale_GeV'])
        for k, r in rows.items():
            c = r['couplings']
            print(f"  {k:10s} mu={r['mu_GeV']:.2e} lamX={c['lam_X']:.5e} lamY={c['lam_Y']:.4f} lamS={c['lam_S']:.4f} kap={c['kap']:.5e}"
                  f" lXY={c['l_XY']:.2e} lXS={c['l_XS']:.2e} lYS={c['l_YS']:.2e} young={r['young_ratio']:.3f} V4min={r['V4_min_unit_sphere']:.3e}")
    # running down from H_inf: what M_Y-scale values does 'lambda=1 at H_inf' mean
    sol = solve_ivp(lambda t, g: beta_vec(g), (0, np.log(1e12 / 4.7e13)), g0, rtol=1e-9, atol=1e-14)
    flows['mu0=H_inf']['values_at_M_Y'] = dict(zip(NAMES, sol.y[:, -1].tolist()))
    res['flows'] = flows
    # RG-improved tree and CW at the inflationary point
    cw = {}
    gm = flows['mu0=M_Y']['rows']['m_Y_point']['couplings']
    coup_run = dict(lam_X=gm['lam_X'], lam_Y=gm['lam_Y'], lam_S=gm['lam_S'], kap=gm['kap'])
    extra_run = dict(XY=gm['l_XY'], XS=gm['l_XS'], YS=gm['l_YS'])
    cw['stated_couplings_muR=4.65H'] = cw_point(dict(COUP), None, 4.65)
    cw['stated_couplings_muR=M_Y'] = cw_point(dict(COUP), None, GEV['M'] / GEV['Hinf'])
    cw['run_from_M_Y_to_4.65H'] = cw_point(coup_run, extra_run, 4.65)
    res['coleman_weinberg'] = cw
    for k, v in cw.items():
        print(k, {kk: (np.round(vv['moduli_over_H'], 4).tolist(), np.round(np.sqrt(np.clip(vv['angular_mass2'], 0, None)), 4).tolist())
                  for kk, vv in v.items() if isinstance(vv, dict)}, 'V1', v['V1_at_tree_point'])
    (OUT / 'running.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
