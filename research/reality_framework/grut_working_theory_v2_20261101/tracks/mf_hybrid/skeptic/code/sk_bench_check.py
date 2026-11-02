"""Internal-consistency audit of the exported bench_HFa.json / bench_HFb.json histories [C].

Checks (independent of the track's code; reads only the JSON):
- N(t) against the integral of H dt; eps_H against -dH/dt / H^2; K_running = a^2 H; a = e^N;
- R against 6 H^2 (2 - eps) (homogeneous part) and against 3 H^2 (f_sig + f_chi) (fluid part: cycle mean);
- energy fractions: sum, f_SM <= f_rad, positivity;
- spectator: f_s in the fluid part against the full s energy (potential + kinetic, s' from the s column), the jump of
  s at the homogeneous -> fluid switch (t = 0.6) and at the exact -> envelope switch (t ~ 60);
- continuity of H, eps and the fractions at t = 0.6; time ordering; file size.
Usage: python3 -B sk_bench_check.py -> ../results/sk_bench_check.json
"""
import os, math
import numpy as np
import sk_common as C


def audit(name):
    d, P = C.load_bench(name)
    h = C.hist(d)
    t, N, a, H, eps, R = h['t'], h['N'], h['a'], h['H'], h['eps_H'], h['R']
    out = dict(rows=int(len(t)), t_range=[float(t[0]), float(t[-1])],
               size_MB=os.path.getsize(os.path.join(C.BENCH, 'bench_%s.json' % name)) / 1e6,
               t_monotonic=bool(np.all(np.diff(t) > 0)))
    # N vs int H dt (trapezoid, from the first row)
    Nint = N[0] + np.concatenate([[0], np.cumsum(0.5 * (H[1:] + H[:-1]) * np.diff(t))])
    out['max_abs_N_minus_intHdt'] = float(np.max(np.abs(Nint - N)))
    out['max_rel_a_minus_expN'] = float(np.max(np.abs(a / np.exp(N) - 1)))
    out['max_rel_K_minus_a2H'] = float(np.nanmax(np.abs(h['K_running'] / (a * a * H) - 1)))
    hom = t <= 0.6 + 1e-9
    flu = t > 0.6 + 1e-9
    # eps vs -dH/dt/H^2 (central differences; skip the switch rows)
    epsfd = -np.gradient(H, t) / H ** 2
    sel = (t > 0.0) & (np.abs(t - 0.6) > 0.02) & (t < 1e4)
    out['eps_fd_vs_column_median_abs'] = float(np.median(np.abs(epsfd[sel] - eps[sel])))
    out['eps_fd_vs_column_p95_abs'] = float(np.percentile(np.abs(epsfd[sel] - eps[sel]), 95))
    Rh = 6 * H ** 2 * (2 - eps)
    out['R_hom_max_abs_dev'] = float(np.max(np.abs(R[hom] - Rh[hom])))
    Rf = 3 * H ** 2 * (h['f_sig'] + h['f_chi'])
    out['R_fluid_max_rel_dev_from_3H2_fm'] = float(np.nanmax(np.abs(R[flu] / Rf[flu] - 1)))
    fsum = h['f_sig'] + h['f_chi'] + h['f_rad']
    out['fraction_sum_range_hom'] = [float(np.min(fsum[hom & (t > 0)])), float(np.max(fsum[hom & (t > 0)]))]
    out['fraction_sum_range_fluid'] = [float(np.min(fsum[flu])), float(np.max(fsum[flu]))]
    out['f_SM_le_f_rad'] = bool(np.all(h['f_SM'][flu] <= h['f_rad'][flu] + 1e-9))
    # continuity at the switch
    j = int(np.argmax(flu))
    out['switch'] = dict(t=[float(t[j - 1]), float(t[j])], H=[float(H[j - 1]), float(H[j])], eps=[float(eps[j - 1]), float(eps[j])],
                         R=[float(R[j - 1]), float(R[j])], f_chi=[float(h['f_chi'][j - 1]), float(h['f_chi'][j])],
                         s=[float(h['s_over_sstar'][j - 1]), float(h['s_over_sstar'][j])],
                         f_s=[float(h['f_s'][j - 1]), float(h['f_s'][j])])
    # spectator energy in the exact fluid part (0.6 < t < 60): the column holds 1/2 m_s^2 s^2 only
    s, ms = h['s_over_sstar'], h['m_s']
    ex = flu & (h['s_is_envelope'] < 0.5)
    tt, ss, mm = t[ex], s[ex], ms[ex]
    sd = np.gradient(ss, tt)
    rho = 3 * C.mP ** 2 * H[ex] ** 2
    fs_full = 0.5 * P['s_star'] ** 2 * (sd ** 2 + mm ** 2 * ss ** 2) / rho
    fs_col = h['f_s'][ex]
    win = (tt > 5) & (tt < 55)
    out['f_s_exact_part'] = dict(column_min_max=[float(fs_col[win].min()), float(fs_col[win].max())],
                                 full_energy_min_max=[float(fs_full[win].min()), float(fs_full[win].max())],
                                 note='column = potential energy only (oscillates to ~0); full = potential + kinetic')
    je = int(np.argmax(h['s_is_envelope'] > 0.5))
    amp_before = float(np.max(np.abs(s[max(je - 40, 0):je])))
    out['s_envelope_switch'] = dict(t=float(t[je]), envelope=float(s[je]), max_abs_s_before=amp_before)
    # spectator energy fraction at its decay (Gamma_s = 1.26e-6): from the column's late slope f_s / a
    late = t > 1e4
    slope = float(np.median(h['f_s'][late] / a[late]))
    K = d['reheating']['K']
    Gs = 1.2629e-6
    a_dec = math.sqrt(2 * K / Gs)          # Gamma_s t = 1 in radiation, a^2 = 2 K t
    a_dec_H = math.sqrt(K / Gs)            # Gamma_s = H
    out['spectator_fraction_at_decay'] = dict(f_s_over_a_late=slope, a_dec_t_eq_1_over_Gamma=a_dec,
                                              f_s_at_a_dec=slope * a_dec, f_s_at_Gamma_eq_H=slope * a_dec_H,
                                              track_value=d['reheating'].get('f_s_at_decay', 0.0243))
    return out


if __name__ == '__main__':
    res = dict(label='C (skeptic audit of the exported JSON)', note=__doc__)
    for nm in ('HFa', 'HFb'):
        res[nm] = audit(nm)
        print(nm, res[nm])
    C.dump('sk_bench_check.json', res)
