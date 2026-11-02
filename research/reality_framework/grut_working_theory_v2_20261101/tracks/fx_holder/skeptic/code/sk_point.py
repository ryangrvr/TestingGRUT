#!/usr/bin/env python3
"""fx_holder_skeptic (Claude; own code) -- static recomputation of the track's items 1, 2, 4 (bounds) with different
methods:
  * g_I, alpha_I, phi_I(0) in 40-digit mpmath from the GeV pins (unit check in GeV);
  * stationary points at u = 10i (lambda = 1, 0.5) by mpmath.findroot on the five real gradient equations (own
    Wirtinger derivatives, no amplitude ansatz), Hessian by mp.diff + mp.eigsy (not sympy);
  * held-phase 'original' preparations at lambda = 1, 0.5 (own solve);
  * 7x7 Hessian with the inflaton (mixing by mp.diff of dV/dphihat);
  * opposite inflaton branch u = 2 mu - 10 i;
  * EXTENDED scan of the prepared branch in s = phi_I/phi_I(0) at fixed H over [-1.5, 4] (the track scanned [0, 2]), with
    a multi-start global check at each s;
  * joint (phi_I, X, Y, Sigma) quartic positivity by direct minimisation on the unit sphere (the track used AM-GM with
    lambda_Y, not lambda_Y,eff);
  * one-loop running structure Tr M^4 as a function of u (numerical), giving the induced phase-blind portal masses at the
    slice per unit ln Q;
  * decay thresholds at lambda = 0.5.
Output: ../results/sk_point.json
"""
import json, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_common as sc

mp.mp.dps = 40
RES = Path(__file__).resolve().parent.parent / 'results'
RES.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(777)
MC = -10.0                       # xi_X R at R = 12 (prescribed background: R(0) = 12 exactly)
TRACK = dict(g3=8.306742069521301e-07, g5=1.3844570115868836e-06, alpha=2.192144074179184, phi0_GeV=1.988169174553435e+18,
             amps1=[87.19556213287578, 5.754247963285591, 28.38718194436162],
             hess1=[13.091408829452787, 14.589249978666992, 538.1466434275541, 1546.6597219731082, 3239.8852511586983],
             amps05=[90.79986663693572, 6.88772550732482, 28.576239631186507],
             hess05=[10.19148131703187, 18.757349084540536, 471.15523484098924, 1536.2773723501484, 1665.7444251505608],
             held05=[76.46750163083311, 2.3398884331919643, 28.302523769265417], held1=[76.39978266076902, 2.311019885471267, 28.293142125868524],
             opp=[125.9363722351776, 16.621691076846435, 29.686215672449453], opp_delta=2.531956736092635, opp_lightest=22.46282095661046,
             eig7_inflaton_shift=-7.965397230691451e-10, mix_norm=0.0003857498986764756,
             lamI_needed_1=4.8e-25, lamI_needed_05=3.8e-24, runaway_1=2.55e6)


def to_f(z):
    return np.array([complex(c) for c in z])


def gauge_fix(z):
    th = np.angle(z[0])
    return z * np.exp(-1j * th * np.array([1, 3, 6]))


def V_r(r, u, p):
    z = np.array([r[0] + 1j * r[1], r[2] + 1j * r[3], r[4] + 1j * r[5]]) / np.sqrt(2)
    return sc.V_c(z, u, p, MC)


def gr_r(r, u, p):
    z = np.array([r[0] + 1j * r[1], r[2] + 1j * r[3], r[4] + 1j * r[5]]) / np.sqrt(2)
    g = sc.grad_c(z, u, p, MC)
    # dV/d(Re z) = 2 Re g, dV/d(Im z) = 2 Im g; r = sqrt2 (Re, Im) -> factor 1/sqrt2
    return np.array([2 * g[0].real, 2 * g[0].imag, 2 * g[1].real, 2 * g[1].imag, 2 * g[2].real, 2 * g[2].imag]) / np.sqrt(2)


def multistart(u, p, n=120, box=250.0):
    found = []
    for _ in range(n):
        amps = np.exp(rng.uniform(np.log(0.5), np.log(box), 3)); ph = rng.uniform(-np.pi, np.pi, 3)
        z0 = amps * np.exp(1j * ph)
        r0 = np.sqrt(2) * np.array([z0[0].real, z0[0].imag, z0[1].real, z0[1].imag, z0[2].real, z0[2].imag])
        m = minimize(V_r, r0, args=(u, p), jac=gr_r, method='BFGS', options=dict(gtol=1e-8, maxiter=10000))
        z = np.array([m.x[0] + 1j * m.x[1], m.x[2] + 1j * m.x[3], m.x[4] + 1j * m.x[5]]) / np.sqrt(2)
        z = gauge_fix(z)
        found.append((float(m.fun), np.abs(z).tolist(), sc.phases(z)[0]))
    found.sort(key=lambda t: t[0])
    distinct = []
    for f in found:
        if not any(abs(f[0] - d[0]) < 1e-5 * max(1.0, abs(d[0])) for d in distinct):
            distinct.append(f)
    return [dict(V=d[0], amps=d[1], delta=d[2]) for d in distinct[:5]]


def stationary_rec(u, p, guess, hess=True):
    z, res = sc.mp_stationary(u, p, MC, guess)
    rec = dict(amps=[float(abs(c)) for c in z], delta=sc.phases(z)[0], delta_y=sc.phases(z)[1], residual=float(res),
               V=float(sc.mp_V([mp.sqrt(2) * c for zz in z for c in (mp.re(zz), mp.im(zz))], mp.mpc(u), p, MC)))
    if hess:
        _, ev, _ = sc.mp_hessian(z, u, p, MC)
        phys = [float(e) for e in ev if abs(e) > 1e-8]
        rec.update(eig=[float(e) for e in ev], physical=phys, lightest=min(phys))
    return z, rec


def main():
    t0 = time.time()
    out = {}
    # -------------------------------------------------------------- 1. coupling in GeV
    H = mp.mpf('134284273491.30815'); MR = mp.mpf('2.435e18'); MUg = mp.mpf('961402368973.181')
    cp = {}
    for b in (3, 5):
        mI = b * H
        phi0 = mp.sqrt(6) * H * MR / mI                       # 3 H^2 M_r^2 = m_I^2 phi^2 / 2 at rest  [D]
        Gg = (10j * H - MUg) / phi0                            # dimensionless: GeV / GeV
        cp[f'b{b}'] = dict(phi0_GeV=float(phi0), phi0_over_Mr=float(phi0 / MR), g_I=float(abs(Gg)), alpha_I=float(mp.arg(Gg)),
                           gphi_over_H=float(abs(Gg) * phi0 / H),
                           g_I_rel_diff_track=float(abs(Gg) / (TRACK['g3'] if b == 3 else TRACK['g5']) - 1))
    cp['alpha_rel_diff_track'] = cp['b3']['alpha_I'] - TRACK['alpha']
    cp['note'] = 'phi_I(0) from the static Friedmann slice of the quadratic toy; the eta^2 hidden-energy correction (~1e-10) is omitted here'
    out['coupling'] = cp
    G3 = complex(10j - sc.MU) / float(mp.sqrt(6) / 3 * MR / H)          # H_i units, b = 3
    # -------------------------------------------------------------- 2. stationary points, held preparations
    st = {}
    zs = {}
    for lam in (1.0, 0.5):
        p = sc.par(lam)
        guess = np.array([87.2, -5.75, -28.4j]) if lam == 1.0 else np.array([90.8, -6.9, -28.6j])
        z, rec = stationary_rec(10j, p, guess)
        zs[lam] = z
        ref_a = TRACK['amps1'] if lam == 1.0 else TRACK['amps05']; ref_h = TRACK['hess1'] if lam == 1.0 else TRACK['hess05']
        rec['amp_rel_diff_track'] = [a / b_ - 1 for a, b_ in zip(rec['amps'], ref_a)]
        rec['hess_rel_diff_track'] = [a / b_ - 1 for a, b_ in zip(rec['physical'], ref_h)]
        rec['multistart'] = multistart(10j, p, n=150)
        # held-phase original
        zh = sc.mp_held(sc.MU, p, MC, [76.4, -2.1, -0.95, 28.3])
        rec['held_original'] = dict(z=[[float(mp.re(c)), float(mp.im(c))] for c in zh], amps=[float(abs(c)) for c in zh],
                                    amp_rel_diff_track=[float(abs(c)) / b_ - 1 for c, b_ in zip(zh, TRACK['held1'] if lam == 1.0 else TRACK['held05'])])
        rec['z'] = [[float(mp.re(c)), float(mp.im(c))] for c in z]
        # 7x7 with the inflaton: canonical phihat (H_i units), tree mass^2 b^2 = 9; mixing = d/dr_i [2 Re(G S* Y^2)]
        r0 = [mp.sqrt(2) * c for zz in z for c in (mp.re(zz), mp.im(zz))]
        Gm = mp.mpc(G3.real, G3.imag)
        def dVdph(rr):
            x = mp.mpc(rr[0], rr[1]) / mp.sqrt(2); y = mp.mpc(rr[2], rr[3]) / mp.sqrt(2); s = mp.mpc(rr[4], rr[5]) / mp.sqrt(2)
            return 2 * mp.re(Gm * mp.conj(s) * y * y)
        mix = []
        for i in range(6):
            def f(a, i=i):
                rr = list(r0); rr[i] += a; return dVdph(rr)
            mix.append(mp.diff(f, 0))
        Hm, ev, _ = sc.mp_hessian(z, 10j, p, MC)
        H7 = mp.matrix(7, 7)
        for i in range(6):
            for j in range(6):
                H7[i, j] = Hm[i, j]
            H7[i, 6] = H7[6, i] = mix[i]
        H7[6, 6] = 9
        E7, _ = mp.eigsy(H7)
        E7 = sorted([E7[k] for k in range(7)])
        k9 = min(range(7), key=lambda k: abs(E7[k] - 9))
        rec['mix_norm'] = float(mp.sqrt(sum(m ** 2 for m in mix)))
        rec['eig7_inflaton_shift'] = float(E7[k9] - 9)
        rest = [E7[k] for k in range(7) if k != k9]
        rec['eig7_max_hidden_shift'] = float(max(abs(a - b_) for a, b_ in zip(rest, ev)))
        rec['hidden_force_on_inflaton_over_tree'] = float(abs(dVdph(r0)) / (9 * float(mp.sqrt(6) / 3 * MR / H)))
        st[f'lam{lam:g}'] = rec
    out['stationary'] = st

    # -------------------------------------------------------------- 3. opposite branch and conjugate
    p = sc.par(1.0)
    uo = 2 * sc.MU - 10j
    zo, reco = stationary_rec(uo, p, np.array([125.9, -16.6 * np.exp(0.3j), 29.7 * np.exp(2.5j)]))
    # refine with a safer guess from multistart if needed
    ms = multistart(uo, p, n=150)
    reco['multistart'] = ms
    reco['rel_diff_track'] = [a / b_ - 1 for a, b_ in zip(reco['amps'], TRACK['opp'])]
    reco['z'] = [[float(mp.re(c)), float(mp.im(c))] for c in zo]
    out['opposite_branch'] = reco

    # -------------------------------------------------------------- 4. extended s-scan at fixed H (R = 12)
    scan = []
    for direction, svals in (('up', np.round(np.arange(1.0, 4.0001, 0.1), 6)), ('down', np.round(np.arange(1.0, -1.5001, -0.1), 6))):
        guess = to_f(zs[1.0])
        for s in svals:
            u = sc.MU + (10j - sc.MU) * s
            try:
                z, rec = stationary_rec(u, p, guess, hess=True)
                guess = to_f(z)
            except Exception as e:  # noqa
                scan.append(dict(s=float(s), branch_lost=str(e)[:80])); break
            do_ms = abs(s * 10 - round(s * 10)) < 1e-9 and (round(s * 10) % 5 == 0)
            if do_ms:
                ms = multistart(u, p, n=80)
                rec['global_V'] = ms[0]['V']; rec['global_amps'] = ms[0]['amps']; rec['global_delta'] = ms[0]['delta']
                rec['is_global'] = bool(abs(rec['V'] - ms[0]['V']) <= 1e-6 * abs(ms[0]['V']))
                rec['n_distinct_minima'] = len(ms)
            rec.update(s=float(s), abs_u=float(abs(u)), direction=direction, forced_light_Y_mass2=float(sc.M ** 2 - 2 * abs(u) * rec['amps'][2]))
            scan.append(rec)
    out['s_scan_fixed_H'] = scan

    # -------------------------------------------------------------- 5. joint quartic positivity (phi_I, X, Y, Sigma)
    jq = {}
    for lam in (1.0, 0.5):
        g = abs(G3)
        lyeff = lam - 27 * sc.KAP ** 4 / (16 * sc.LX ** 3)
        # Q = (lamI/4) phi^4 + lx x^4 + lam y^4 + lam s^4 - 2 kap y x^3 - 2 g phi s y^2 (worst-case phases) is homogeneous of
        # degree 4; on phi = 1 rescale (x, y, s) = g t:  Q = lamI/4 + g^4 f(t),  f = lx tx^4 + lam ty^4 + lam ts^4 - 2 kap ty tx^3 - 2 ts ty^2,
        # so the joint form is positive iff lamI >= -4 g^4 min f (and the phi = 0 slice, the hidden form, is positive).
        def f(t):
            tx, ty, ts = t
            return sc.LX * tx ** 4 + lam * ty ** 4 + lam * ts ** 4 - 2 * sc.KAP * ty * tx ** 3 - 2 * ts * ty ** 2
        best = 0.0
        for _ in range(60):
            m = minimize(f, rng.uniform(0.1, 30, 3), method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
            best = min(best, m.fun)
        hi = -4 * g ** 4 * best
        jq[f'lam{lam:g}'] = dict(lamI_crit_numeric_quarter_norm=float(hi), min_f=float(best), analytic_min_f=float(-1 / (4 * lam * lyeff ** 2)),
                                 analytic_quarter_norm_with_lamYeff=float(g ** 4 / (lam * lyeff ** 2)),
                                 analytic_quarter_norm_with_lamY=float(g ** 4 / (lam * lam ** 2)), lambda_Y_eff=float(lyeff),
                                 track_quoted=TRACK['lamI_needed_1'] if lam == 1.0 else TRACK['lamI_needed_05'],
                                 runaway_Sigma_lamYeff=float(3 * np.sqrt(lyeff / 2) / g), runaway_Sigma_lamY=float(3 * np.sqrt(lam / 2) / g),
                                 one_loop_quartic_g4_over_16pi2=float(g ** 4 / (16 * np.pi ** 2)))
    out['joint_quartic'] = jq

    # -------------------------------------------------------------- 6. running structure: Tr M^4 vs u (numerical Hessians)
    def trm4(z, u, lam):
        p_ = sc.par(lam)
        Hm, _, _ = sc.mp_hessian(tuple(mp.mpc(c) for c in z), u, p_, 0.0)
        return float(sum(Hm[i, j] ** 2 for i in range(6) for j in range(6)))
    run = {}
    for nm, z in (('Sigma_only', (0, 0, 1.0)), ('Y_only', (0, 1.0, 0)), ('zero', (0, 0, 0))):
        vals = [trm4(z, u, 1.0) for u in (0.0, 3.0, 3j, 3 + 3j)]
        run[nm] = vals
    cS = (run['Sigma_only'][1] - run['Sigma_only'][0] - (run['zero'][1] - run['zero'][0])) / 9.0
    cY = (run['Y_only'][1] - run['Y_only'][0] - (run['zero'][1] - run['zero'][0])) / 9.0
    du2 = 100.0 - sc.MU ** 2
    out['running_structure'] = dict(TrM4_coef_of_u2_Sigma2=cS, TrM4_coef_of_u2_Y2=cY,
                                    check_phase_independence=[run['Sigma_only'][2] - run['Sigma_only'][1], run['Y_only'][3] - run['Y_only'][0]],
                                    note='dV/dlnQ = -Tr M^4/(32 pi^2); induced phase-blind terms cport (|u|^2 - mu^2)(|S|^2 + 2|Y|^2) with cport = 1/(4 pi^2) per unit ln Q',
                                    delta_mS2_slice_per_unit_log=float(cS * du2 / (32 * np.pi ** 2)),
                                    delta_mY2_slice_per_unit_log=float(cY * du2 / (32 * np.pi ** 2)),
                                    track_delta_mY2_per_unit_log=7.7, track_delta_mS2_per_unit_log=1.2,
                                    delta_mY2_phi2_part_only=float(cY * (abs(10j - sc.MU) ** 2) / (32 * np.pi ** 2)),
                                    over_ln_MPl_MY=float(np.log(2.435e18 / (sc.M * sc.H_GEV))))

    # -------------------------------------------------------------- 7. decay thresholds
    dec = {}
    for lam in (1.0, 0.5):
        s_ = sc.HV / np.sqrt(2)
        mL = np.sqrt(sc.M ** 2 - 2 * sc.MU * s_); mS = np.sqrt(4 * lam * s_ ** 2)
        dec[f'lam{lam:g}'] = dict(mYL=float(mL), m_sigma=float(mS), thr_YY=float(2 * mL), thr_SYY=float(mS + 2 * mL), thr_sigma_sigma=float(2 * mS),
                                  closed_for_b5=bool(5 < 2 * mL))
    out['decays'] = dec
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / 'sk_point.json', 'w'), indent=1, default=float)
    print(json.dumps(out['coupling'], indent=1))
    for k, v in out['stationary'].items():
        print(k, v['amps'], v['delta'], v['physical'], 'res', v['residual'], 'reldiff', v['amp_rel_diff_track'], v['hess_rel_diff_track'])
        print('   multistart', v['multistart'][:3])
        print('   held', v['held_original'])
        print('   mix', v['mix_norm'], v['eig7_inflaton_shift'], v['eig7_max_hidden_shift'], v['hidden_force_on_inflaton_over_tree'])
    print('opp', {k: v for k, v in out['opposite_branch'].items() if k not in ('multistart', 'eig')}, out['opposite_branch']['multistart'][:3])
    for d in scan:
        if 'branch_lost' in d:
            print(d); continue
        print('s=%5.2f |u|=%6.2f amps=%s delta=%7.4f lightest=%8.3f mL2f=%8.1f %s' % (d['s'], d['abs_u'], np.round(d['amps'], 3), d['delta'], d['lightest'],
              d['forced_light_Y_mass2'], ('global=%s nmin=%d gV=%.1f V=%.1f gamps=%s' % (d['is_global'], d['n_distinct_minima'], d['global_V'], d['V'], np.round(d['global_amps'], 2))) if 'is_global' in d else ''))
    print(json.dumps(out['joint_quartic'], indent=1))
    print(json.dumps(out['running_structure'], indent=1))
    print(json.dumps(out['decays'], indent=1))
    print('seconds', out['seconds'])


if __name__ == '__main__':
    main()
