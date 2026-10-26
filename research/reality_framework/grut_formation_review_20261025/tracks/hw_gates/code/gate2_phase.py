#!/usr/bin/env python3
"""hw_gates gate 2 (Claude, own code): phase holding at the ls_scaling spot point and what sets the phase value.

(a) 6x6 canonical Hessian (gauge mode projected out) at the aligned minimum and at the three supplied states;
    light relative-phase mode (delta = theta_Sigma - 6 theta_X), its mass and decay constant f_eff.
(b) Conditional-minimum potential V_rel(delta) (Y-relaxed path), its minimum and the de Sitter equilibrium
    P ∝ exp(-8 pi^2 V/3 H^4) (Starobinsky-Yokoyama [R]): probability of the band around pi/2 if nothing holds it.
(c) Spread of a HELD phase: exact Bunch-Davies coarse-grained variance of a mode of mass m at the light-mode mass
    (complex nu for m > 3H/2), sigma_delta = sigma_c/f_eff, for modes outside the horizon at charge creation
    (x_c = k/aH = 1, 1/10, 1/40).
(d) Classical drift of delta from pi/2 once the holder is off while inflation continues (1D: f^2 delta'' + 3 f^2 delta'
    + dV_rel/ddelta = 0 in e-folds), mapped to eps through the own phase scan.
(e) Three holding designs that make a supplied state an exact inflationary minimum, with their post-switch-off states:
    D1 mu-rotation  muP = mu(e^{i delta} - 1)          -> Codex preparation (Y real);
    D2 kappa-rotation kapP = kap(e^{i delta/2} - 1)    -> Y, Sigma rotated by delta/2, delta;
    D3 X^6 Sigma*-rotation c7 (dimension 7, does not enter the Y equation) -> the Y-relaxed preparation.
    Each released with instant removal, exponential removal over Dl, H-tracking persistence and permanence (hwgrel).
(f) Isocurvature: delta psi on CMB scales for the held mode (m > 3H/2: ~ e^{-3N/2}) times the transfer d ln N_X/d delta.
Usage: python3 -B gate2_phase.py <out.json>
"""
import json, sys, time
from math import gamma
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spotmodel as S

P = S.params()
D0 = np.pi / 2
# own phase scan at the spot point [C, ls_scaling B_phase; same model]: delta (rad) -> eps
PHASE_SCAN = [(0.2618, -0.0786), (0.7854, -0.1415), (1.3090, -0.1538), (1.5708, -0.1457), (1.8326, -0.1298),
              (2.3562, -0.0802), (2.8798, -0.0172)]


def eps_of_delta(d):
    xs = [0.0] + [a for a, _ in PHASE_SCAN] + [np.pi]
    ys = [0.0] + [b for _, b in PHASE_SCAN] + [0.0]
    return float(np.interp(abs(d) % (2 * np.pi) if abs(d) <= np.pi else 2 * np.pi - abs(d) % (2 * np.pi), xs, ys))


# ---------------------------------------------------------------- (c) exact BD coarse-grained variance
def P_mode(x, m2):
    nu = mp.sqrt(mp.mpf(9) / 4 - m2)
    return float(mp.re(x**3 * mp.exp(-mp.pi * mp.im(nu)) * abs(mp.hankel1(nu, x))**2 / (8 * mp.pi)))


def sigma2_coarse(m2, xc):
    nu2 = 9 / 4 - m2
    x0 = min(1e-3, xc)
    if nu2 > 0:
        nu = np.sqrt(nu2); Cn = 2**(nu - 1.5) * gamma(nu) / gamma(1.5)
        low = Cn**2 / (4 * np.pi**2) * x0**(3 - 2 * nu) / (3 - 2 * nu)
    else:
        low = float(mp.quad(lambda t: P_mode(t, m2) / t, [1e-9, x0]))
    if xc <= x0:
        return low
    return low + float(mp.quad(lambda t: P_mode(t, m2) / t, mp.linspace(x0, xc, 8)))


# ---------------------------------------------------------------- (b) conditional-minimum potential along delta
def V_rel(delta):
    z, info = S.prep('relaxed', delta, P)
    return info['V'], info


# ---------------------------------------------------------------- (e) designs
def design(kind, delta=D0):
    if kind == 'D1':
        z, _ = S.prep('codex', delta, P)
        return z, dict(muP=P['mu'] * (np.exp(1j * delta) - 1))
    if kind == 'D2':
        z, _ = S.prep('kaprot', delta, P)
        return z, dict(kapP=P['kap'] * (np.exp(1j * delta / 2) - 1))
    if kind == 'D3':
        # c7 such that the Y-relaxed state is stationary: solve for complex c7 and polish the state (6D) jointly
        zr, info = S.prep('relaxed', delta, P)
        a6 = S.alpha6(P, 1.0)
        c7_eft = P['mu'] * P['kap']**2 / (P['M']**4 - 4 * (P['mu'] * P['h'] / np.sqrt(2))**2)   # coefficient of X^6 Sigma*
        # sign convention: V ⊃ 2Re(c7_eft X^6 Sigma*) is the Y-eliminated term; rotation by e^{i delta}: c7 -> c7 e^{-i delta}?
        # choose the rotation that cancels the torque at the relaxed state (both signs tried; keep the stationary one)
        best = None
        for sgn in (+1, -1):
            c7 = c7_eft * (np.exp(sgn * 1j * delta) - 1)
            def F(u):
                return S.V(u[0::2] + 1j * u[1::2], P) + 2 * np.real(c7 * (u[0] + 1j * u[1])**6 * np.conj(u[4] + 1j * u[5]))
            u0 = np.ravel([[c.real, c.imag] for c in zr])
            r = minimize(F, u0, method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
            z = r.x[0::2] + 1j * r.x[1::2]
            # fix the gauge: rotate so that arg X = 0
            ph = np.exp(-1j * np.angle(z[0]) * S.Q)
            z = z * ph
            dev = np.max(np.abs(z - zr) / np.abs(zr))
            g = S.gradient_canon(z, P, dict(c7r=c7.real, c7i=c7.imag))
            if best is None or dev < best[2]:
                best = (z, c7, dev, float(np.linalg.norm(g)), sgn)
        z, c7, dev, gn, sgn = best
        return z, dict(c7=c7, c7_eft=c7_eft, sign=sgn, rel_dev_from_relaxed=dev, grad_norm=gn,
                       Lambda7_GeV=float((1 / abs(c7))**(1 / 3) * S.H_SPOT), relaxed_state=[complex(v) for v in zr])
    raise ValueError(kind)


def release_job(args):
    kind, z0, extra, tag = args
    try:
        r = S.run(np.array(z0, complex), P, a_end=60.0, rtol=1e-10, atol=1e-13, nsamp=30, win0=30.0, extra=extra)
    except Exception as exc:
        return dict(design=kind, tag=tag, status=repr(exc))
    w = r['window']; L = S.ledger(P, w)
    return dict(design=kind, tag=tag, eps=w['eps'], NX=w['NX'], EX=w['EX'], NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig'],
                R_Y=L['R_Y'], R_Sigma=L['R_Sigma'], R_cascade=L['R_cascade'], R_central=L['R_central'],
                R_L=S.R_L(w['eps']), yield_over_YQ=S.yield_over_YQ(w['NX']), thS_window=w['thS'], thG_window=w['thG'],
                minS=r['end']['minS_over_vev'], max_gauss=r['end']['max_gauss'], steps=r['end']['steps'])


def main(out):
    t0 = time.time()
    res = dict(method=__doc__)
    # (a) spectra
    xa, ya, sa = S.aligned(P)
    za = np.array([xa, -ya, -sa])
    spa = S.spectrum(za, P)
    res['aligned'] = dict(moduli=[xa, ya, sa], spectrum=spa, F_g_over_H=float(np.sqrt(S.higgs_norm(za))))
    states = {}
    for kind in ('codex', 'relaxed', 'kaprot'):
        z, info = S.prep(kind, D0, P)
        states[kind] = dict(z=[complex(v) for v in z], info=info, spectrum=S.spectrum(z, P))
    res['supplied_states'] = states
    lp = spa['light_phase']
    print('aligned light phase: m = %.4f H, f_eff = %.3f H, ddelta = %.4f' % (lp['mass_over_H'], lp['f_eff'], lp['ddelta']), flush=True)
    # ls_scaling's estimate 18 eta cH for comparison
    eta = S.alpha6(P, 1.0)['alpha6_GeVm2'] * (P['cH'] / (2 * P['lx'])) / P['lx']
    res['ls_scaling_estimate'] = dict(eta=eta, m2_18_eta_cH=18 * eta * P['cH'], m_over_H=float(np.sqrt(18 * eta * P['cH'])),
                                      note='leading |X| = sqrt(cH/2lx) = 70.7 and f^2 = X^2/18 (no Sigma kinetic weight)')
    # (b) V_rel(delta)
    ds = np.linspace(0, np.pi, 25)
    Vr = []
    for d in ds:
        v, info = V_rel(d)
        Vr.append(dict(delta=float(d), V=v, x=info['x'], s=info['s'], Y=info['Y'], torque_X=info['torque_X']))
    V0 = Vr[0]['V']
    B = 0.5 * (Vr[-1]['V'] - V0)
    kSY = 8 * np.pi**2 / 3
    # curvature at delta = 0 from a parabola through the first points, in delta: compare m^2 f^2
    d2 = np.polyfit(ds[:4], [r['V'] for r in Vr[:4]], 2)[0] * 2
    res['V_rel'] = dict(rows=Vr, B_half_barrier=B, d2V_ddelta2_at_0=float(d2),
                        m2_times_f2_from_hessian=lp['m2'] * lp['f_eff']**2,
                        equilibrium_exponent_at_pi2=float(kSY * (V_rel(D0)[0] - V0)),
                        equilibrium_exponent_full_barrier=float(kSY * 2 * B),
                        note='P(delta) ∝ exp(-8 pi^2 V/3 H^4): the exponent at pi/2 is the log of the suppression relative to delta = 0')
    print('V_rel: B = %.1f, d2V = %.1f, m2 f2 = %.1f, exponent at pi/2 = %.0f' % (B, d2, lp['m2'] * lp['f_eff']**2, res['V_rel']['equilibrium_exponent_at_pi2']), flush=True)
    # (c) spread of a held phase at the light-mode mass, and for a scan of holder masses
    m2 = lp['m2']; f = lp['f_eff']
    spread = {}
    for xc in (1.0, 0.1, 1 / 40):
        s2 = sigma2_coarse(m2, xc)
        spread[f'xc_{xc:g}'] = dict(sigma_canonical_over_H=float(np.sqrt(s2)), sigma_delta_rad=float(np.sqrt(s2) / f))
    scan = []
    for m in (0.6, 1.0, 1.5, 2.0, 2.53, 3.0, 5.0):
        scan.append(dict(m_over_H=m, **{f'sigma_delta_xc_{xc:g}': float(np.sqrt(sigma2_coarse(m * m, xc)) / f) for xc in (1.0, 0.1, 1 / 40)}))
    res['held_spread'] = dict(light_mode_m2=m2, f_eff=f, at_light_mass=spread, vs_mass=scan,
                              band_halfwidth_for_eps_ge_0p10_rad=None)
    # band of delta with |eps| >= 0.10 and >= 0.13 from the own scan
    dd = np.linspace(0, np.pi, 2001)
    ee = np.array([eps_of_delta(d) for d in dd])
    for thr in (0.10, 0.13):
        sel = dd[np.abs(ee) >= thr]
        res['held_spread'][f'band_abs_eps_ge_{thr}'] = [float(sel.min()), float(sel.max())] if sel.size else None
    print('spread', spread, res['held_spread'].get('band_abs_eps_ge_0.1'), flush=True)
    # (d) classical drift once the holder is off (inflation continues): 1D in e-folds
    Vd = np.array([r['V'] for r in Vr]); dV = np.gradient(Vd, ds)
    dVf = lambda d: np.interp(abs(d), ds, dV) * np.sign(d)
    def rhs(N, y):
        return [y[1], -3 * y[1] - dVf(y[0]) / f**2]
    sol = solve_ivp(rhs, (0, 3.0), [D0, 0.0], method='DOP853', rtol=1e-10, atol=1e-12, dense_output=True)
    drift = []
    for N in (0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0):
        d = float(sol.sol(N)[0])
        drift.append(dict(efolds_after_holder_off=N, delta=d, eps_from_scan=eps_of_delta(d)))
    res['classical_drift'] = drift
    print('drift', [(r['efolds_after_holder_off'], round(r['delta'], 3), round(r['eps_from_scan'], 3)) for r in drift], flush=True)
    # (e) designs and releases
    designs = {}
    jobs = []
    for kind in ('D1', 'D2', 'D3'):
        z, d = design(kind)
        sp = S.spectrum(z, P, dict(mpr=np.real(d.get('muP', 0)), mpi=np.imag(d.get('muP', 0)), kpr=np.real(d.get('kapP', 0)),
                                    kpi=np.imag(d.get('kapP', 0)), c7r=np.real(d.get('c7', 0)), c7i=np.imag(d.get('c7', 0))))
        gn = float(np.linalg.norm(S.gradient_canon(z, P, dict(mpr=np.real(d.get('muP', 0)), mpi=np.imag(d.get('muP', 0)),
                                                              kpr=np.real(d.get('kapP', 0)), kpi=np.imag(d.get('kapP', 0)),
                                                              c7r=np.real(d.get('c7', 0)), c7i=np.imag(d.get('c7', 0))))))
        designs[kind] = dict(state=[complex(v) for v in z], params={k: (complex(v) if isinstance(v, (complex, np.complexfloating)) else v)
                                                                     for k, v in d.items() if k != 'relaxed_state'},
                             grad_norm_with_holder=gn, min_mass2_with_holder=min(m['m2'] for m in sp['modes'] if not m['gauge']),
                             light_phase_with_holder=sp['light_phase'])
        cost = {}
        if kind == 'D1':
            cost = dict(abs_muP_over_H=abs(d['muP']), abs_muP_over_mu=abs(d['muP']) / P['mu'], abs_muP_GeV=abs(d['muP']) * S.H_SPOT)
        if kind == 'D2':
            cost = dict(abs_kapP=abs(d['kapP']), abs_kapP_over_kap=abs(d['kapP']) / P['kap'])
        if kind == 'D3':
            cost = dict(abs_c7_Hm3=abs(d['c7']), abs_c7_over_c7_eft=abs(d['c7']) / d['c7_eft'], Lambda7_GeV=d['Lambda7_GeV'],
                        Lambda7_over_MY=d['Lambda7_GeV'] / (P['M'] * S.H_SPOT))
        designs[kind]['cost'] = cost
        print(kind, 'grad with holder', gn, 'min m2', designs[kind]['min_mass2_with_holder'], cost, flush=True)
        ex0 = {k: v for k, v in d.items() if k in ('muP', 'kapP', 'c7')}
        jobs.append((kind, designs[kind]['state'], dict(pmode=0), 'instant'))
        for Dl in (0.1, 0.3, 1.0, 3.0):
            jobs.append((kind, designs[kind]['state'], dict(ex0, pmode=1, Dl=Dl), f'exp_Dl{Dl}'))
        jobs.append((kind, designs[kind]['state'], dict(ex0, pmode=2), 'H_tracking'))
        jobs.append((kind, designs[kind]['state'], dict(ex0, pmode=3), 'permanent'))
    # X Hubble mass persistence (Codex assumption of instant removal) for D3 state
    jobs.append(('D3', designs['D3']['state'], dict(pmode=1, Dl=0.1, cX=P['cH']), 'XHubble_exp_Dl0.1'))
    jobs.append(('D3', designs['D3']['state'], dict(pmode=1, Dl=0.3, cX=P['cH']), 'XHubble_exp_Dl0.3'))
    jobs.append(('D3', designs['D3']['state'], dict(pmode=2, cX=P['cH']), 'XHubble_H_tracking'))
    res['designs'] = designs
    with Pool(4) as pool:
        rel = pool.map(release_job, jobs)
    res['releases'] = rel
    for r in rel:
        if 'eps' in r:
            print('%s %-18s eps %.4f R_c %.2f R_L %s NYH %.3f NYL %.3f Nsig %.2e yield %.3f' % (
                r['design'], r['tag'], r['eps'], r['R_central'], None if r['R_L'] is None else round(r['R_L'], 1), r['NYH'], r['NYL'], r['Nsig'], r['yield_over_YQ']), flush=True)
        else:
            print(r, flush=True)
    # (f) isocurvature of the held phase (m = light mode or holder mass), transfer from the own phase scan
    T_local = abs(np.gradient([b for _, b in PHASE_SCAN], [a for a, _ in PHASE_SCAN])[3]) / abs(PHASE_SCAN[3][1])  # d ln|eps|/d delta at pi/2
    iso = []
    for m in (2.53, 1.5, 1.0):
        for N in (50, 60):
            nu2 = 2.25 - m * m
            if nu2 > 0:
                nu = np.sqrt(nu2); amp = 2**(nu - 1.5) * gamma(nu) / gamma(1.5) / (2 * np.pi * f) * np.exp(-(1.5 - nu) * N)
            else:
                amp = np.exp(-1.5 * N) / (2 * np.pi * f)
            iso.append(dict(m_over_H=m, N=N, dpsi_cmb=float(amp), S_over_Smax=float(T_local * amp / 9.1e-6)))
    res['isocurvature'] = dict(transfer_dlnNX_ddelta_at_pi2=float(T_local), rows=iso, S_max=9.1e-6,
                               note='uncorrelated CDM isocurvature S = T dpsi; S_max 9.1e-6 [Q, as in our earlier tracks]; for m > 3H/2 the superhorizon amplitude decays as e^{-3N/2} (order of magnitude)')
    Path(out).write_text(json.dumps(res, indent=1, default=lambda o: (repr(o) if isinstance(o, complex) else float(o))))
    print('done', time.time() - t0)


if __name__ == '__main__':
    main(sys.argv[1])
