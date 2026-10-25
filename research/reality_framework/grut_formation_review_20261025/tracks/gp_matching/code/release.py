#!/usr/bin/env python3
"""Homogeneous three-field release from the inflationary point: transfer of the light phase into
late X charge (own code).

Assumptions (supplied, as in the review's W1 branch): instant reheating, radiation domination
H = H_i/a^2 from a=1, the Hubble mass -c_H H^2|X|^2 switched off at a=1, zero initial velocities,
temporal gauge A0=0 with zero initial Gauss charge (e then drops out of homogeneous dynamics).
Conformal variables u = a z, tau with a = tau (H_i = 1):
  u_X'' = -[2 lam_X |u_X|^2 u_X + 3 kap u_Y conj(u_X)^2]
  u_Y'' = -[2 lam_Y |u_Y|^2 u_Y + a^2 M^2 u_Y + kap u_X^3 + 2 mu a u_S conj(u_Y)]
  u_S'' = -[2 lam_S (|u_S|^2 - a^2 h^2/2) u_S + mu a u_Y^2]
Comoving charges N_i = -2 Im(conj(u_i) u_i'); Gauss N_X + 3N_Y + 6N_S = 0.
Late X charge after Y -> 3X decays: N_eff = N_X + 3 N_Y = -6 N_S, sourced only by the psi2 torque.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import GEV, COUP

OUT = Path(__file__).resolve().parents[1] / 'results'
MRED = 2.435e18
YQ = 4.8294009076988175e-9   # required net yield [Q, affleck_dine]


def params(Hinf=GEV['Hinf']):
    return dict(M=GEV['M'] / Hinf, h=GEV['h'] / Hinf, mu=GEV['mu'] / Hinf, **COUP)


def rhs_factory(p):
    lX, lY, lS, k, M, h, mu = p['lam_X'], p['lam_Y'], p['lam_S'], p['kap'], p['M'], p['h'], p['mu']

    def rhs(t, s):
        a = t
        uX, uY, uS = s[0] + 1j * s[1], s[2] + 1j * s[3], s[4] + 1j * s[5]
        fX = 2 * lX * abs(uX)**2 * uX + 3 * k * uY * np.conj(uX)**2
        fY = 2 * lY * abs(uY)**2 * uY + a * a * M * M * uY + k * uX**3 + 2 * mu * a * uS * np.conj(uY)
        fS = 2 * lS * (abs(uS)**2 - a * a * h * h / 2) * uS + mu * a * uY**2
        return [s[6], s[7], s[8], s[9], s[10], s[11],
                -fX.real, -fX.imag, -fY.real, -fY.imag, -fS.real, -fS.imag]
    return rhs


def charges(s):
    u = [s[0] + 1j * s[1], s[2] + 1j * s[3], s[4] + 1j * s[5]]
    v = [s[6] + 1j * s[7], s[8] + 1j * s[9], s[10] + 1j * s[11]]
    return np.array([-2 * np.imag(np.conj(ui) * vi) for ui, vi in zip(u, v)])


def x_energy(s, p):
    uX, vX = s[0] + 1j * s[1], s[6] + 1j * s[7]
    return abs(vX)**2 + p['lam_X'] * abs(uX)**4


def initial(point, delta0=0.0, radial_kick=None):
    R = np.array(point['moduli_over_H'])
    z0 = np.array([R[0], -R[1], -R[2]], complex)
    th = np.array(point['light_mode']['dtheta_per_canonical']) / point['light_mode']['ddelta']
    z = z0 * np.exp(-1j * th * delta0)
    if radial_kick is not None:   # canonical radial amplitude along the given real-component vector
        z = z + np.asarray(radial_kick) / np.sqrt(2)   # canonical real-component displacement
    u = z.copy(); v = z.copy()     # a=1, a'=1, zdot=0  -> u'=z
    return np.ravel([[w.real, w.imag] for w in u] + [[w.real, w.imag] for w in v])


def run(point, p, delta0, a_end, rtol=1e-10, radial_kick=None, n_eval=4000):
    s0 = initial(point, delta0, radial_kick)
    a_eval = np.unique(np.concatenate([np.linspace(1.0, a_end, n_eval), [2., 5., 10., 20., 35., 70., 150., 300., 600., 1000.]]))
    a_eval = a_eval[a_eval <= a_end]
    sol = solve_ivp(rhs_factory(p), (1.0, a_end), s0, method='DOP853', rtol=rtol, atol=1e-12, t_eval=a_eval)
    hist = []
    for j, a in enumerate(sol.t):
        s = sol.y[:, j]; N = charges(s)
        uX, uY, uS = s[0] + 1j * s[1], s[2] + 1j * s[3], s[4] + 1j * s[5]
        vX, vY, vS = s[6] + 1j * s[7], s[8] + 1j * s[9], s[10] + 1j * s[11]
        kX, kY, kS = abs(vX - uX / a)**2, abs(vY - uY / a)**2, abs(vS - uS / a)**2   # a^4 |zdot|^2
        EX = kX + p['lam_X'] * abs(uX)**4
        EY = kY + a * a * p['M']**2 * abs(uY)**2 + p['lam_Y'] * abs(uY)**4
        ES = kS + p['lam_S'] * (abs(uS)**2 - a * a * p['h']**2 / 2)**2
        nmax = 2 * np.sqrt(2 * p['lam_X']) * (EX / (3 * p['lam_X']))**0.75
        hist.append([a, N[0], N[1], N[2], N[0] + 3 * N[1], N[0] + 3 * N[1] + 6 * N[2], N[0] / nmax, EX, EY, ES])
    return np.array(hist), sol.nfev


def summarize(hist, windows=((300, 600), (600, 1000), (1000, 3000))):
    keys = ['a', 'N_X', 'N_Y', 'N_S', 'N_eff', 'gauss', 'eps_X', 'E_X', 'E_Y', 'E_S']
    pick = {}
    for a in [2., 5., 10., 20., 35., 70., 150., 300., 600., 1000., 3000.]:
        j = np.where(np.isclose(hist[:, 0], a))[0]
        if len(j):
            pick[str(a)] = dict(zip(keys, hist[j[0]].tolist()))
    win = {}
    for lo, hi in windows:
        sel = (hist[:, 0] >= lo) & (hist[:, 0] <= hi)
        if sel.sum() > 10:
            h = hist[sel]
            win[f'{lo}-{hi}'] = dict(N_eff_mean=h[:, 4].mean(), N_eff_std=h[:, 4].std(), N_X_mean=h[:, 1].mean(),
                                     N_X_std=h[:, 1].std(), eps_X_mean=h[:, 6].mean(),
                                     EY_over_EX_mean=(h[:, 8] / h[:, 7]).mean(), ES_over_EX_mean=(h[:, 9] / h[:, 7]).mean())
    return dict(at=pick, windows=win, max_abs_gauss=float(np.abs(hist[:, 5]).max()))


def entropy_ai3(Hinf, gstar=106.75):
    rho = 3 * MRED**2 * Hinf**2
    T = (30 * rho / (np.pi**2 * gstar))**0.25
    return 2 * np.pi**2 / 45 * gstar * T**3, T


def main():
    point = json.loads((OUT / 'point_full.json').read_text())['baseline']
    p = params()
    s_i, T_R = entropy_ai3(GEV['Hinf'])
    H3 = GEV['Hinf']**3
    out = dict(assumptions=__doc__.split('\n\n')[1], T_R_GeV=T_R, s_i_GeV3=s_i, Y_per_unit_N=H3 / s_i, runs={})
    for d0 in [0.0, 0.05, -0.05, 0.3, np.pi / 4, np.pi / 2, -np.pi / 2, 2.0, 2.8]:
        hist, nfev = run(point, p, d0, 1000.)
        sm = summarize(hist); sm['nfev'] = nfev; sm['delta0'] = d0
        out['runs'][f'{d0:.4f}'] = sm
        w = sm['windows']['600-1000']
        print(f"delta0={d0:+.4f} N_eff@1000={sm['at']['1000.0']['N_eff']:+.4e} win600-1000 mean={w['N_eff_mean']:+.3e} std={w['N_eff_std']:.3e}"
              f" N_X mean={w['N_X_mean']:+.3e} std={w['N_X_std']:.3e} epsX={w['eps_X_mean']:+.3e} EY/EX={w['EY_over_EX_mean']:.3e}"
              f" ES/EX={w['ES_over_EX_mean']:.3e} gauss={sm['max_abs_gauss']:.1e}", flush=True)
    hist, _ = run(point, p, np.pi / 4, 1000., rtol=1e-12)
    ref = out['runs'][f'{np.pi / 4:.4f}']['at']['1000.0']['N_eff']
    out['tolerance_control_rel_change_N_eff_a1000'] = abs(hist[-1, 4] - ref) / abs(ref)
    a = out['runs']['0.0500']['at']['1000.0']['N_eff']; b = out['runs']['-0.0500']['at']['1000.0']['N_eff']
    out['CP_antisymmetry_rel'] = abs(a + b) / abs(a)
    out['aligned_N_eff_max_abs'] = max(abs(v['N_eff']) for v in out['runs']['0.0000']['at'].values())
    hist, nfev = run(point, p, np.pi / 4, 3000., n_eval=12000)
    sm = summarize(hist); sm['nfev'] = nfev
    out['long_run_pi_over_4'] = sm
    print('long', sm['windows'], 'tol', out['tolerance_control_rel_change_N_eff_a1000'], 'CP', out['CP_antisymmetry_rel'],
          'aligned', out['aligned_N_eff_max_abs'])
    (OUT / 'release.json').write_text(json.dumps(out, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
