#!/usr/bin/env python3
"""st_response track (Claude, own code): compare our independent tangent runs (runs.py) with Codex's archived
quiet3.npz / quiet5.npz (b73de64, read-only), plus own numerical controls and an rp C-engine background check.
Usage: python3 -B compare.py  -> ../results/st_response.json
"""
import json, os, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import st_tangent as S
import rp

SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
RUNS = SP / 'review/st_response/runs'
RES = HERE.parent / 'results'
SEEDS = [0.03, 1e-6, 1e-9]
LEVELS = [0.1, 1.0]


def load(job):
    return dict(np.load(RUNS / (job + '.npz')))


def crossing(times, a, values, level):
    """First sampled |value| >= level (Codex's definition, one column)."""
    w = np.flatnonzero(np.abs(values) >= level)
    if not len(w):
        return None
    i = int(w[0])
    return dict(tau_bracket=[float(times[max(0, i - 1)]), float(times[i])], a_bracket=[float(a[max(0, i - 1)]), float(a[i])],
                sampled_value=float(np.abs(values[i])))


def codex_bg_obs(bg):
    """Codex archived background columns (v layout of its prepared(): 2:8 z, 8:14 velocities, 15 ln a, 16 H) ->
    a, H, N_X = chi0 chi1' - chi1 chi0' (our own evaluation of the stated observable on the archived numbers)."""
    a = np.exp(bg[:, 15]); H = bg[:, 16]
    x = bg[:, 2] + 1j * bg[:, 3]; vx = bg[:, 8] + 1j * bg[:, 9]
    chi = np.sqrt(2) * a * x; chip = np.sqrt(2) * a ** 2 * (vx + H * x)
    return a, H, (chi.real * chip.imag - chi.imag * chip.real)


def window_mean(times, y, lo=64.0, hi=80.0):
    w = (times >= lo) & (times <= hi)
    return float(np.trapezoid(y[w], times[w]) / (times[w][-1] - times[w][0]))


def diff_metrics(t, mine, ref):
    d = np.abs(mine - ref)
    env = np.maximum.accumulate(np.abs(ref))
    ok = env > 0                                   # the charge response starts at exactly zero
    big = ok & (np.abs(ref) > 1e-3 * env)
    r = np.where(ok, d / np.where(ok, env, 1.0), 0.0)
    i = int(np.argmax(r))
    return dict(max_abs_over_1plus=float(np.max(d / (1 + np.abs(ref)))),
                max_rel_to_running_envelope=float(np.max(r)), at_tau=float(t[i]),
                max_pointwise_rel_where_not_near_zero=float(np.max(d[big] / np.abs(ref[big]))),
                max_rel_to_global_max=float(np.max(d) / np.max(np.abs(ref))),
                final_rel=float(d[-1] / np.abs(ref[-1])) if ref[-1] != 0 else None)


def brackets(times, a, amp, charge):
    out = []
    for s in SEEDS:
        row = dict(seed=s)
        for L in LEVELS:
            row['amplitude_%g' % L] = crossing(times, a, s * amp, L)
        row['local_charge_1'] = crossing(times, a, s * charge, 1.0)
        out.append(row)
    return out


def growth(times, env, t1, t2):
    i1, i2 = np.searchsorted(times, t1), np.searchsorted(times, t2)
    return float(np.log(env[i2] / env[i1]) / (times[i2] - times[i1]))


def compare_history(tag, codexfile, job, pairs):
    C = np.load(CODEX / codexfile)
    tC, ksC = C['times'], C['ks']
    M = load(job)
    t = M['times']; n = len(t)
    assert np.allclose(t, tC[:n], atol=1e-12)
    aC, HC, NXC = codex_bg_obs(C['bg'][:n])
    out = dict(job=job, codex_file=codexfile, tau_end=float(t[-1]))
    out['background'] = dict(
        max_rel_a=float(np.max(np.abs(M['a'] / aC - 1))), max_rel_H=float(np.max(np.abs(M['H'] / HC - 1))),
        max_rel_NX_to_scale=float(np.max(np.abs(M['NX'] - NXC)) / np.max(np.abs(NXC))),
        max_abs_X_component=float(np.max(np.abs(M['Y'][:, 3:5] - C['bg'][:n, 2:4]))),
        own_NX_norm_64_80=window_mean(t, M['NX']), codex_NX_norm_64_80_own_eval=window_mean(t, NXC),
        own_final_a=float(M['a'][-1]), codex_final_a=float(aC[-1]),
        own_late_C=float(M['a'][-1] ** 2 * M['H'][-1]), codex_late_C_own_eval=float(aC[-1] ** 2 * HC[-1]),
        own_final_NX=float(M['NX'][-1]), own_final_NX_change=float(M['NX'][-1] / window_mean(t, M['NX']) - 1))
    norm = abs(window_mean(t, M['NX']))
    out['modes'] = {}
    for jk, kval in pairs:
        ic = int(np.argmin(np.abs(ksC - kval))); assert abs(ksC[ic] - kval) < 1e-12
        jm = int(np.argmin(np.abs(M['ks'] - kval))); assert abs(M['ks'][jm] - kval) < 1e-12
        amp, ampC = M['amp'][:, jm], C['amplitude_per_radian'][:n, ic]
        ch, chC = M['dN'][:, jm] / norm, C['charge_fraction_per_radian'][:n, ic]
        envm, envc = np.maximum.accumulate(np.abs(amp)), np.maximum.accumulate(np.abs(ampC))
        row = dict(k=kval,
                   amplitude=diff_metrics(t, amp, ampC), charge_fraction=diff_metrics(t, ch, chC),
                   charge_sign_check=dict(dot_sign=float(np.sign(np.dot(ch, chC))),
                                          max_rel_if_sign_flipped=float(np.max(np.abs(-ch - chC)) / np.max(np.abs(chC)))),
                   own_brackets=brackets(t, M['a'], amp, ch), codex_column_brackets=brackets(t, aC, ampC, chC),
                   own_max_abs_amp=float(envm[-1]), codex_max_abs_amp=float(envc[-1]),
                   own_max_abs_amp_tau_le_240=float(envm[np.searchsorted(t, 240.0) - 1]) if t[-1] >= 240 else None)
        if t[-1] >= 700:
            row['envelope_growth_rate_200_700'] = dict(own=growth(t, envm, 200, 700), codex=growth(t, envc, 200, 700))
        out['modes']['k=%g' % kval] = row
    return out


def self_compare(jobA, jobB, kval, label):
    A, B = load(jobA), load(jobB)
    n = min(len(A['times']), len(B['times']))
    assert np.allclose(A['times'][:n], B['times'][:n])
    ia = int(np.argmin(np.abs(A['ks'] - kval))); ib = int(np.argmin(np.abs(B['ks'] - kval)))
    normA = abs(window_mean(A['times'], A['NX']))
    normB = abs(window_mean(B['times'], B['NX']))
    return dict(label=label, tau_end=float(A['times'][n - 1]),
                amplitude=diff_metrics(A['times'][:n], A['amp'][:n, ia], B['amp'][:n, ib]),
                charge_fraction=diff_metrics(A['times'][:n], A['dN'][:n, ia] / normA, B['dN'][:n, ib] / normB),
                background_a=float(np.max(np.abs(A['a'][:n] / B['a'][:n] - 1))))


def gauge_effect():
    """Size of the e-dependent term and whether our e-effect equals Codex's (difference of differences)."""
    out = {}
    e0, e7, e6 = load('q3_k002_e0'), load('q3_k3'), load('q3_gauge6')
    Q3, G3 = np.load(CODEX / 'quiet3.npz'), np.load(CODEX / 'gauge3.npz')
    t = e7['times']
    for kval, i7, i6, iq, ig in [(0.02, 0, 0, 0, 0), (2.13, 1, 1, 9, 1), (5.8, 2, 2, 15, 2)]:
        own = e6['amp'][:, i6] - e7['amp'][:, i7]
        cod = G3['amplitude_per_radian'][:, ig] - Q3['amplitude_per_radian'][:, iq]
        env = np.maximum.accumulate(np.abs(Q3['amplitude_per_radian'][:, iq]))
        row = dict(k=kval, own_effect_e6_minus_e7_rel_to_envelope=float(np.max(np.abs(own) / env)),
                   codex_effect_rel_to_envelope=float(np.max(np.abs(cod) / env)),
                   effect_difference_rel_to_effect=float(np.max(np.abs(own - cod)) / np.max(np.abs(cod))))
        if kval == 0.02:
            eff0 = e7['amp'][:, 0] - e0['amp'][:, 0]
            row['own_effect_e7_minus_e0_rel_to_envelope'] = float(np.max(np.abs(eff0) / env))
            row['own_e6_vs_e0_final_ratio'] = float(e6['amp'][-1, 0] / e0['amp'][-1, 0])
        out['k=%g' % kval] = row
    return out


def old_archive():
    """Old 2cbc7f1 PHASE_PROFILES (float32) at k = 2.122989207841491 through tau 80: amplitude and the charge sign."""
    old = np.load(CODEX.parent / 'grut_reheating_modes_20261026' / 'PHASE_PROFILES.npz')
    ks = old['b3_g0.5_ks']; kk = 2.122989207841491
    i = int(np.argmin(np.abs(ks - kk)))
    arc = old['b3_g0.5'][:, i, :].astype(float)
    M = load('q3_old'); main = load('q3_k213')
    norm = abs(window_mean(main['times'], main['NX']))
    ch = M['dN'][:, 0] / norm
    fit = float(np.dot(ch, arc[:, 1]) / np.dot(ch, ch))       # least-squares ratio archive / ours
    return dict(k_archive=float(ks[i]), amplitude=diff_metrics(M['times'], M['amp'][:, 0], arc[:, 0]),
                charge_same_sign=diff_metrics(M['times'], ch, arc[:, 1]),
                charge_opposite_sign=diff_metrics(M['times'], -ch, arc[:, 1]),
                least_squares_ratio_archive_over_ours=fit,
                note='archive is float32; our charge uses q0p1-q1p0 normalised by our 64-80 window mean; the old row used |late_NX|')


def rp_engine_check():
    """Our Python transcription vs the rp C engine (rpreheat.c) at tau = 0.5 n, n = 160..1600; and the yr_b3_g0.5 yield."""
    out = {}
    res, Sx = rp.run('st_bgcheck_b3', mode=0, b=3.0, gam=0.5, end=800.0, nsamp=1600, save=False)
    M = load('q3_k213')
    t = M['times']
    cols = rp.COLS
    out['rp_C_vs_python'] = {}
    for name, mine, sgn in [('a', M['a'], 1), ('H', M['H'], 1), ('NX', M['NX'], -1), ('absX', M['absX'], 1)]:
        ref = np.interp(t[t >= 80], Sx[:, 0], Sx[:, cols.index(name)])
        # samples coincide (0.5 grid from 80): check exact time match
        idx = np.searchsorted(Sx[:, 0], t[t >= 80])
        assert np.allclose(Sx[idx, 0], t[t >= 80], atol=1e-9)
        ref = Sx[idx, cols.index(name)]
        out['rp_C_vs_python'][name] = float(np.max(np.abs(sgn * mine[t >= 80] - ref)) / np.max(np.abs(ref)))
    out['rp_C_checks'] = res['checks']
    r50, _ = rp.run('st_yr_b3_g0.5', mode=0, b=3.0, gam=0.5, end=80.0, save=False)   # runs.py 'extra': end = max(25/g, 80)
    out['yield_yr_b3_g0.5_rerun'] = r50['yield_over_YQ']
    out['yield_yr_b3_g0.5_committed'] = 1.0454206
    return out


def main():
    R = dict(method=__doc__, model_doc=S.__doc__)
    # seed
    C3 = np.load(CODEX / 'quiet3.npz')
    _, q0, info = S.initial(3.0)
    seedC = C3['selected_vectors'][0, 1]
    R['seed'] = dict(own_q0=q0.tolist(), codex_q0=seedC[:6].tolist(), codex_p0_equals_q0=bool(np.array_equal(seedC[:6], seedC[6:])),
                     max_abs_diff=float(np.max(np.abs(q0 - seedC[:6]))), max_rel_to_norm=float(np.max(np.abs(q0 - seedC[:6])) / np.linalg.norm(seedC[:6])),
                     own_A0=float(q0[0] / (np.sqrt(2) * info['z0'][0].real)), prep_residuals=[info['prep']['residual'], info['prep_plus_residual'], info['prep_minus_residual']])
    R['selftest'] = S.selftest()
    R['quiet3_k213'] = compare_history('quiet3', 'quiet3.npz', 'q3_k213', [(0, 2.13)])
    R['quiet3_three_k'] = compare_history('quiet3', 'quiet3.npz', 'q3_k3', [(0, 0.02), (1, 2.13), (2, 5.8)])
    R['quiet5_k2'] = compare_history('quiet5', 'quiet5.npz', 'q5_k2', [(0, 2.0), (1, 2.13)])
    # archive's own max-over-k brackets for context (Codex numbers, recomputed from the npz)
    for name, f in [('quiet3', 'quiet3.npz'), ('quiet5', 'quiet5.npz')]:
        C = np.load(CODEX / f); a = np.exp(C['bg'][:, 15]); ks = C['ks']
        mx = lambda v: np.max(np.abs(v), axis=1)
        rows = []
        for s in SEEDS:
            r = dict(seed=s)
            for L in LEVELS:
                c = crossing(C['times'], a, s * mx(C['amplitude_per_radian']), L)
                if c:
                    i = int(np.searchsorted(C['times'], c['tau_bracket'][1]))
                    c['k'] = float(ks[int(np.argmax(np.abs(C['amplitude_per_radian'][i])))])
                r['amplitude_%g' % L] = c
            r['local_charge_1'] = crossing(C['times'], a, s * mx(C['charge_fraction_per_radian']), 1.0)
            rows.append(r)
        R['codex_archive_max_over_k_' + name] = rows
    R['controls'] = dict(
        tight_vs_main=self_compare('q3_k213_tight', 'q3_k213', 2.13, 'rtol 1e-12 vs 1e-10 (DOP853), tau<=240'),
        rk45_vs_main=self_compare('q3_k213_rk45', 'q3_k213', 2.13, 'RK45 rtol 1e-11 vs DOP853 1e-10, tau<=160'),
        threek_vs_onek=self_compare('q3_k3', 'q3_k213', 2.13, 'k=2.13 inside the 3-k run vs alone (different step sequence), tau<=800'),
        e0_vs_e1e7=self_compare('q3_k213_e0', 'q3_k213', 2.13, 'SENSITIVITY: e=0 vs e=1e-7, tau<=240 (not a numerical error)'))
    R['gauge3_e1e-6'] = compare_history('gauge3', 'gauge3.npz', 'q3_gauge6', [(0, 0.02), (1, 2.13), (2, 5.8)])
    R['quiet5_k002'] = compare_history('quiet5', 'quiet5.npz', 'q5_k002', [(0, 0.02)])
    R['gauge_effect'] = gauge_effect()
    R['old_phase_archive'] = old_archive()
    R['rp_engine'] = rp_engine_check()
    RES.mkdir(exist_ok=True)
    (RES / 'st_response.json').write_text(json.dumps(R, indent=1, default=float))
    print(json.dumps({k: v for k, v in R.items() if k not in ('method', 'model_doc')}, indent=1, default=float))


if __name__ == '__main__':
    main()
