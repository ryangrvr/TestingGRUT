"""Part 2b: homogeneous X + phi evolution through the phi condensation (independent code).

Conformal radiation-era potential (a''/a = 0):
  U = lambda_phi/4 phi^4 - lambda_phi/2 M(s)^2 phi^2 + g^2/2 phi^2 r^2 + lambda_X/4 r^4,  M(s) = v a(s),
so that late on phi -> M and the X mass is g M.  M(s) = M_start exp(kappa s) with kappa = eps * omega_r
(eps = d ln a/ds per radial angular frequency of the initial orbit), or the RD form M_start (1 + kappa s).
The X orbit starts as the reviewed pure-quartic orbit (tau=20 or a late-eccentricity orbit with the same q);
phi starts at a small seed with zero velocity.  Final measures (late massive regime, window averaged):
  n_-/q = <((p_X^2/2 + g^2 M^2 r^2/2)/(g M) - q)/2> / q     (free massive antiparticle number per unit charge)
  N_phi/q = <(dphi'^2/2 + m_phi^2 dphi^2/2)/m_phi> / q,  dphi = |phi| - M, m_phi^2 = 2 lambda_phi M^2
(number of phi quanta in the coherent phi oscillation per unit charge)."""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *
from part2_floquet import orbit_info, peri_state

GG = np.sqrt(G2)


def phase_state(y4, Tr, ph, nper=4096):
    par = make_par(mode=0, lamramp=0, T=1e-9, mout=0.0, lamX0=LX, g2=0.0, lamP=0.0)
    y = y8(y4[0], y4[1], 0.0, y4[2], y4[3], 0.0, s=1.0)
    n = int(round(ph * nper))
    if n: c_fixed(par, y, Tr / nper, n, 8, store=False)
    return np.array([y[0], y[1], y[3], y[4]])


def diag(Y, par, q):
    """Instantaneous diagnostics for an array of 8-vectors."""
    Y = np.atleast_2d(Y)
    s = Y[:, 6]
    if par.growth == 0:
        M = par.M0 * np.exp(par.kappa * s); Mp = par.kappa * M
    else:
        M = par.M0 * (1 + par.kappa * s); Mp = par.M0 * par.kappa * np.ones_like(s)
    r2 = Y[:, 0] ** 2 + Y[:, 1] ** 2
    pX2 = Y[:, 3] ** 2 + Y[:, 4] ** 2
    m = GG * M
    nminus = ((0.5 * pX2 + 0.5 * m * m * r2) / m - q) / 2 / q
    mphi = np.sqrt(2 * LP) * M
    sg = np.sign(Y[:, 2]); sg[sg == 0] = 1
    dphi = np.abs(Y[:, 2]) - M
    dphid = sg * Y[:, 5] - Mp
    Nphi = (0.5 * dphid ** 2 + 0.5 * mphi ** 2 * dphi ** 2) / mphi / q
    qq = Y[:, 0] * Y[:, 4] - Y[:, 1] * Y[:, 3]
    return dict(s=s, M=M, nminus=nminus, Nphi=Nphi, phi=Y[:, 2], r2=r2, qrel=qq / q - 1)


def window_avg(t, v):
    return np.trapezoid(v, t) / (t[-1] - t[0])


def transition_run(y4, info, eps, seed, M_start_fac=0.5, M_end_fac=40.0, growth=0, h0=0.08, wb0=20.0, wbM=0.5,
                   n_check=200, final_periods=60, keep_trace=True, phi_dot0=0.0, mode=1):
    q = charge(y4)
    Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']
    M_start = M_start_fac * Mperi; M_end = M_end_fac * Mapo
    kappa = eps * info['omega_r']
    if growth == 0:
        s_end = np.log(M_end / M_start) / kappa
        s_checks = np.log(np.geomspace(M_start, M_end, n_check + 1)[1:] / M_start) / kappa
    else:
        s_end = (M_end / M_start - 1) / kappa
        s_checks = (np.geomspace(M_start, M_end, n_check + 1)[1:] / M_start - 1) / kappa
    par = make_par(mode=mode, M0=M_start, kappa=kappa, growth=growth)
    y = y8(y4[0], y4[1], seed * Mperi if mode == 1 else 0.0, y4[2], y4[3], phi_dot0 * seed * Mperi if mode == 1 else 0.0, s=0.0, ps=0.0)
    E0 = c_energy(par, y)
    trace = []
    for sc in s_checks:
        c_run(par, y, sc, h0, wb0, wbM, 8)
        if keep_trace: trace.append(y.copy())
    Ef = c_energy(par, y)
    # work identity: E(s_end) - E(0) = W = -(ps_end - ps_0)
    work_rel = abs(Ef - E0 + y[7]) / abs(Ef)
    # final averaging window: final_periods X periods at the final mass
    Tw = final_periods * 2 * np.pi / (GG * M_end)
    W = c_run(par, y, s_end + Tw, h0, wb0, wbM, 8, dsample=Tw / 4000, nmax=4400)
    d = diag(W, par, q)
    res = dict(eps=eps, seed=seed, growth=growth, M_start=M_start, M_end=M_end, kappa=kappa,
               nminus=window_avg(d['s'], d['nminus']), Nphi=window_avg(d['s'], d['Nphi']),
               nminus_win_std=float(np.std(d['nminus'])), Nphi_win_std=float(np.std(d['Nphi'])),
               qrel=float(np.max(np.abs(d['qrel']))), work_rel=work_rel)
    if keep_trace:
        T = diag(np.array(trace), par, q)
        res['trace'] = dict(M=T['M'], nminus=T['nminus'], Nphi=T['Nphi'], phi_over_M=T['phi'] / T['M'])
    return res


if __name__ == '__main__':
    y4 = tau20_state(); info = orbit_info(y4)
    eps = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-2
    t0 = time.time()
    st = phase_state(y4, info['T_r'], 0.0)
    r = transition_run(st, info, eps, 1e-6)
    print('eps', eps, 'nminus', r['nminus'], 'Nphi', r['Nphi'], 'std', r['nminus_win_std'], r['Nphi_win_std'], 'qrel', r['qrel'], 'work', r['work_rel'], time.time() - t0)
    tr = r['trace']
    for i in range(0, len(tr['M']), 5):
        print(f"M/Mperi {tr['M'][i]/(np.sqrt(G2/LP)*info['rmin']):9.3f}  n- {tr['nminus'][i]:10.5f}  Nphi {tr['Nphi'][i]:10.5f}  phi/M {tr['phi_over_M'][i]: .4f}")
