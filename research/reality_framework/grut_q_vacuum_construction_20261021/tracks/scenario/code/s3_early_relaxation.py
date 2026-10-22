"""Part 3: early relaxation of a displaced q (Planck-scale initial vacuum energy) with mu = mu0, and the decay of
the oscillation energy into radiation.

Units: M = 1 (K0 = M^2 = 1), KV 0806.2805 shape with chi = 1/M^4 (natural): m^2 = 2/3 M^2, q0 = 3 M^2.
KV time tau = t/sqrt(chi |F0|) = t/sqrt(3); h = sqrt(3) H; r = chi rho = rho.

(a) QV1 proper (Gamma = 0): reproduce KV's <rho_V> = (2/3) d_M^2 M^2/t^2 (equivalently tau^2 <r_V> = (2/9) d_M^2,
    tau <h> = 2/3) with KV's boundary data (a, h, f, r_M) = (1, 1/2, 1/2, r_M) at tau = 1.
(b) QV1+Gamma (effective source [Q/E]): the q-equation gets -Gamma_q qdot (Rayleigh-type dissipation, as in KV
    1601.04676), the metric equation keeps the constant mu, radiation receives Q = Gamma_q qdot^2 (Bianchi-consistent,
    s1). Gamma = N_s m^3/(192 pi M^2), N_s = 4 (Higgs doublet, assumed minimally coupled, xi = 0)
    [Q: Gorbunov-Panin 1009.2448 eq (T4+), reduced M_P; Starobinsky 1980; Vilenkin 1985].
Outputs: transferred fraction f_inj = int Q a^3 dt / (rho_V,i a_i^3); comoving radiation; T_rh; dependence on the
initial displacement y_i = f_i - 1 and on Gamma.
Run: python3 -I s3_early_relaxation.py [outdir]   (about 5 min on 4 cores)
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qv1_core import KV_SHAPE, KVMatter, QV1, integrate_full  # noqa: E402

S3 = math.sqrt(3.0)
M_RED_GEV = 2.435e18
GSTAR = 106.75


def kv_initial(q, y_i, h_i, rM_i, rest=False):
    """KV boundary data (h_i given, Ydot from the constraint), or rest=True: q displaced and at rest (V = 0),
    H from the constraint, so that the only q-sector energy at the start is rho_V(y_i)."""
    Y0 = y_i
    if rest:
        V0 = 0.0
        H0 = q.H_from_constraint(Y0, 0.0, rM_i)
        return [Y0, V0, H0, 0.0, 0.0, rM_i, 0.0]
    H0 = h_i / S3
    V0 = (rM_i + q.rhoV(Y0) - 3 * H0 * H0 * (1 + Y0)) / (3 * H0)
    return [Y0, V0, H0, 0.0, 0.0, rM_i, 0.0]


def run_kv_attractor(wM, rM_i, tau_end=3000.0, rtol=1e-11):
    """Gamma = 0: reproduce KV (5.15)-(5.17)."""
    t0w = time.time()
    q = QV1(math.sqrt(2.0 / 3.0), **KV_SHAPE)
    mat = KVMatter(wM)
    y0 = kv_initial(q, -0.5, 0.5, rM_i)
    T = np.linspace(0.6 * tau_end * S3, tau_end * S3, 200001)
    out = integrate_full(q, mat, y0, S3, T, rtol=rtol, atol=1e-13)
    Y, V, H = out[:, 0], out[:, 1], out[:, 2]
    rV = np.array([q.rhoV(v) for v in Y])
    rM = out[:, 5]
    rM_inf = float(np.mean((T ** 2 / 3 * rM)[-50000:]))
    dM = 1.0 if wM > 0 else math.sqrt(max(0.0, 1 - 2.25 * rM_inf))
    C = np.array([q.constraint(o[0], o[1], o[2], o[5] + o[6])[0] for o in out[::2000]])
    Ca4 = C * np.exp(4 * out[::2000, 3])
    return dict(wM=wM, rM_i=rM_i, tau_window=[0.6 * tau_end, tau_end],
                t_mean_H=float(np.mean(T * H)), KV_t_mean_H=2.0 / 3.0,
                t2_mean_rhoV=float(np.mean(T ** 2 * rV)), KV_t2_mean_rhoV=2.0 / 3.0 * dM ** 2,
                rM_inf_tau2=rM_inf, d_M=dM,
                ratio_rhoV_to_crit=float(np.mean(T ** 2 * rV) / (4.0 / 3.0)),
                max_abs_C_over_3H2=float(np.max(np.abs(C)) / np.max(3 * H ** 2)),
                C_a4_spread=float(np.ptp(Ca4)), wall_s=time.time() - t0w)


def run_decay(y_i, Gamma_factor=1.0, Ns=4, h_i=0.5, rM_i=0.05, wM=1.0 / 3.0, n_decay=14.0, rtol=1e-11,
              pts_per_period=48, rest=False):
    t0w = time.time()
    m = math.sqrt(2.0 / 3.0)
    Gamma = Gamma_factor * Ns * m ** 3 / (192 * math.pi)
    q = QV1(m, Gamma=Gamma, **KV_SHAPE)
    mat = KVMatter(wM)
    y0 = kv_initial(q, y_i, h_i, rM_i, rest=rest)
    rhoV_i = q.rhoV(y_i)
    Ekin_i = 0.75 * y0[1] ** 2
    t_start = S3
    t_end = t_start + n_decay / Gamma
    period = 2 * math.pi / m
    n = int((t_end - t_start) / period * pts_per_period)
    T = np.linspace(t_start, t_end, n + 1)[1:]
    out = integrate_full(q, mat, y0, t_start, T, rtol=rtol, atol=1e-15)
    Y, V, H, N = out[:, 0], out[:, 1], out[:, 2], out[:, 3]
    rM, rr = out[:, 5], out[:, 6]
    a = np.exp(N)
    Q = 1.5 * Gamma * V ** 2 / (1 + Y) ** 2
    TT = np.concatenate([[t_start], T])
    Qa3 = np.concatenate([[1.5 * Gamma * y0[1] ** 2 / (1 + y0[0]) ** 2], Q * a ** 3])
    inj = float(np.trapezoid(Qa3, TT))
    rr_a4_end = float(rr[-1] * a[-1] ** 4)
    # oscillation (Einstein-frame, linearised) energy times a^3
    Eosc_a3 = 0.75 * (V ** 2 + m ** 2 * Y ** 2) * a ** 3
    # average over periods
    k = pts_per_period
    nb = len(T) // k
    Eb = Eosc_a3[:nb * k].reshape(nb, k).mean(axis=1)
    tb = T[:nb * k].reshape(nb, k).mean(axis=1)
    # effective decay rate of Eosc a^3 after the first 3/Gamma
    sel = (tb > t_start + 3 / Gamma) & (tb < t_start + 10 / Gamma)
    Gfit = float(-np.polyfit(tb[sel], np.log(Eb[sel]), 1)[0]) if sel.sum() > 5 else None
    # plateau of Eosc a^3 before decay matters, and late-time radiation domination
    i_g = np.searchsorted(T, t_start + 1.0 / Gamma)
    rho_r_tGamma = float(rr[i_g])
    Trh_GeV = (30 * rho_r_tGamma / (math.pi ** 2 * GSTAR)) ** 0.25 * M_RED_GEV
    # equality time of decay radiation and oscillation energy (linearised Einstein-frame energy)
    Eosc = 0.75 * (V ** 2 + m ** 2 * Y ** 2)
    kk = (len(T) // k) * k
    Eavg = Eosc[:kk].reshape(-1, k).mean(axis=1)
    rravg = rr[:kk].reshape(-1, k).mean(axis=1)
    eq = np.where(rravg > Eavg)[0]
    t_eq = float(tb[eq[0]]) if len(eq) else None
    rr_eq = float(rravg[eq[0]]) if len(eq) else None
    Teq_GeV = (30 * rr_eq / (math.pi ** 2 * GSTAR)) ** 0.25 * M_RED_GEV if rr_eq else None
    # Gorbunov-Panin estimate [Q] T_reh ~ 4.5e-2 g*^-1/4 (N_s mu^3/M_P)^1/2  (in units M = 1)
    T_GP = 4.5e-2 * GSTAR ** -0.25 * math.sqrt(Ns * Gamma_factor * m ** 3) * M_RED_GEV
    C = np.array([q.constraint(o[0], o[1], o[2], o[5] + o[6])[0] for o in out[::max(1, len(T) // 400)]])
    a4 = np.exp(4 * out[::max(1, len(T) // 400), 3])
    late = T > t_end - 2 / Gamma
    return dict(y_i=y_i, rest=rest, Gamma=Gamma, Gamma_over_m=Gamma / m, Ns=Ns, Gamma_factor=Gamma_factor,
                h_i=float(y0[2] * S3), H_i=float(y0[2]), H_i_over_m=float(y0[2] / m), V_i=float(y0[1]),
                Ekin_q_i=Ekin_i, rM_i=rM_i, rho_crit_i=float(3 * y0[2] ** 2 * (1 + y0[0]) + 3 * y0[2] * y0[1]),
                wM=wM, rhoV_i=rhoV_i, rhoV_i_over_M4=rhoV_i, t_end_times_Gamma=n_decay,
                frac_final_radiation_from_q=float(rr_a4_end / (rr_a4_end + rM[-1] * a[-1] ** 4)) if wM > 0 else None,
                f_inj=inj / rhoV_i, injected_comoving=inj, rr_a4_end=rr_a4_end, f_com=rr_a4_end / rhoV_i,
                rM_a4_end=float(rM[-1] * a[-1] ** 4), rM_a4_start=float(rM_i),
                Eosc_a3_first_window=float(Eb[0]), Eosc_a3_at_1_over_Gamma=float(np.interp(t_start + 1 / Gamma, tb, Eb)),
                Gamma_fit_Eosc=Gfit, Gamma_fit_over_Gamma=(Gfit / Gamma if Gfit else None),
                late_tH=float(np.mean(T[late] * H[late])), late_w_note='tH -> 1/2 for radiation domination',
                residual_Eosc_over_rr_end=float(Eavg[-1] / rravg[-1]),
                rho_r_at_t_1_over_Gamma=rho_r_tGamma, T_rh_GeV_at_t_1_over_Gamma=Trh_GeV,
                t_equality_decay_radiation=t_eq, T_at_equality_GeV=Teq_GeV, T_reh_GorbunovPanin_GeV=T_GP,
                max_abs_C_over_3H2=float(np.max(np.abs(C)) / np.max(3 * H ** 2)),
                C_a4_max_drift=float(np.max(np.abs(C * a4 - C[0] * a4[0]))), wall_s=time.time() - t0w,
                series=dict(t=[float(x) for x in tb[::max(1, nb // 150)]],
                            Eosc_a3=[float(x) for x in Eb[::max(1, nb // 150)]],
                            rr_a4=[float(x) for x in (rr * a ** 4)[:nb * k].reshape(nb, k).mean(axis=1)[::max(1, nb // 150)]],
                            tH=[float(x) for x in (T * H)[:nb * k].reshape(nb, k).mean(axis=1)[::max(1, nb // 150)]]))


def _job(a):
    kind, kw = a
    if kind == 'attr':
        return kind, kw, run_kv_attractor(**kw)
    return kind, kw, run_decay(**kw)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                                                  'results')
    T0 = time.time()
    jobs = [('attr', dict(wM=1.0 / 3.0, rM_i=1.0 / 20.0)), ('attr', dict(wM=0.0, rM_i=1.0 / 200.0)),
            ('attr', dict(wM=0.0, rM_i=0.5))]
    for y_i in (-0.9, -0.75, -0.5, -0.25, -0.1, -0.03, 0.1, 0.25, 0.5, 1.0):
        jobs.append(('dec', dict(y_i=y_i)))
    for y_i in (-0.9, -0.75, -0.5, -0.25, -0.1, -0.03, 0.03, 0.1, 0.25, 0.5, 1.0):
        jobs.append(('rest', dict(y_i=y_i, rest=True, rM_i=0.0)))
    jobs.append(('rest', dict(y_i=-0.5, rest=True, rM_i=0.0, Gamma_factor=0.1, n_decay=12.0)))
    jobs.append(('rest', dict(y_i=-0.5, rest=True, rM_i=0.0, Gamma_factor=10.0)))
    for gf in (0.1, 10.0):
        jobs.append(('dec', dict(y_i=-0.5, Gamma_factor=gf, n_decay=12.0 if gf < 1 else 14.0)))
    jobs.append(('dec', dict(y_i=-0.5, rM_i=0.0)))           # no initial matter: pure vacuum start
    jobs.append(('dec', dict(y_i=-0.5, rtol=1e-12)))          # tolerance check
    jobs.append(('dec', dict(y_i=-0.5, rtol=1e-10)))
    res = {'attractor_Gamma0': [], 'decay_vs_displacement': [], 'decay_vs_Gamma': [], 'decay_no_initial_matter': [],
           'decay_tolerance': [], 'rest_vs_displacement': [], 'rest_vs_Gamma': []}
    with ProcessPoolExecutor(max_workers=4) as ex:
        for kind, kw, r in ex.map(_job, jobs):
            print(kind, kw, 'wall %.1f' % r['wall_s'], flush=True)
            if kind == 'attr':
                res['attractor_Gamma0'].append(r)
            elif kind == 'rest':
                r.pop('series', None) if kw['y_i'] != -0.5 else None
                (res['rest_vs_Gamma'] if kw.get('Gamma_factor', 1.0) != 1.0 else res['rest_vs_displacement']).append(r)
            elif kw.get('Gamma_factor', 1.0) != 1.0:
                res['decay_vs_Gamma'].append(r)
            elif kw.get('rM_i', 0.05) == 0.0:
                res['decay_no_initial_matter'].append(r)
            elif 'rtol' in kw:
                res['decay_tolerance'].append({k: r[k] for k in ('y_i', 'f_inj', 'f_com', 'Gamma_fit_over_Gamma',
                                                                 'T_rh_GeV_at_t_1_over_Gamma', 'max_abs_C_over_3H2')}
                                              | {'rtol': kw['rtol']})
            else:
                res['decay_vs_displacement'].append(r)
    res['decay_vs_displacement'].sort(key=lambda d: d['y_i'])
    res['rest_vs_displacement'].sort(key=lambda d: d['y_i'])
    for d in res['decay_vs_displacement'][1:]:
        d.pop('series', None)
    res['meta'] = dict(units='M = 1 (reduced Planck mass), K0 = 1, chi = 1/M^4, m^2 = 2/3 M^2; t in 1/M',
                       Gamma='N_s m^3/(192 pi M^2), N_s = 4, effective source [Q/E]',
                       fractions='f_inj = int Q a^3 dt / (rho_V,i a_i^3) (energy handed to radiation, counted at '
                                 'injection); f_com = rho_r a^4 (end) / (rho_V,i a_i^4) (comoving radiation, a_i = 1)',
                       wall_total_s=time.time() - T0)
    with open(os.path.join(outdir, 'early_relaxation.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    for d in res['decay_vs_displacement']:
        print('y_i %+.2f rhoV_i %.3e f_inj %.4f f_com %.3e Gfit/G %s Trh %.3e GeV' % (
            d['y_i'], d['rhoV_i'], d['f_inj'], d['f_com'], d['Gamma_fit_over_Gamma'], d['T_rh_GeV_at_t_1_over_Gamma']))
    for d in res['rest_vs_displacement']:
        print('rest y_i %+.2f rhoV_i %.3e H_i/m %.3f f_inj %.4f f_com %.3e Trh %.3e GeV' % (
            d['y_i'], d['rhoV_i'], d['H_i_over_m'], d['f_inj'], d['f_com'], d['T_rh_GeV_at_t_1_over_Gamma']))
    print('total', time.time() - T0)


if __name__ == '__main__':
    main()
