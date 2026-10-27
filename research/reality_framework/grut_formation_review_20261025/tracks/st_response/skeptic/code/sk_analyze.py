#!/usr/bin/env python3
"""st_response skeptic (Claude, own code): analyse sk_runs.py output against Codex's archive (read-only) and against
st_response's stated numbers.  Writes ../results/sk_response.json.

Everything here is recomputed from our own sk_tangent.py histories (third formulation, cosmic-time physical fields,
dual-number Hessian-vector products, simplified geometry, e = 0).  Codex's archive is used only as the comparison
target and for its sample grid / bracket convention (first sampled max-over-k crossing, bracket [t_{i-1}, t_i]).
"""
import json, sys
from pathlib import Path
import numpy as np

SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
RUNS = SP / 'review/st_response_skeptic/runs'
RES = Path(__file__).resolve().parents[1] / 'results'
SEEDS = (0.03, 1e-6, 1e-9)
KOLD = 2.122989207841491


def load(job):
    files = sorted((RUNS / job).glob('chunk_*.npz'))
    parts = [dict(np.load(f)) for f in files]
    out = {}
    for key in ('t', 'a', 'H', 'NX', 'absX', 'A', 'dN'):
        arrs = [p[key] if i == 0 else p[key][1:] for i, p in enumerate(parts)]   # drop duplicated chunk starts
        out[key] = np.concatenate(arrs, axis=0)
    out['ks'] = parts[0]['ks']
    assert np.all(np.diff(out['t']) > 0)
    return out


def at(t, tq):
    j = np.clip(np.searchsorted(t, tq), 1, len(t) - 1)
    idx = np.where(np.abs(t[j - 1] - tq) <= np.abs(t[j] - tq), j - 1, j)
    assert np.max(np.abs(t[idx] - tq)) < 1e-9
    return idx


def window_mean(t, y, lo=64.0, hi=80.0):
    w = (t >= lo) & (t <= hi)
    return float(np.trapezoid(y[w], t[w]) / (t[w][-1] - t[w][0])), float(t[w][0]), float(t[w][-1]), int(w.sum())


def first_cross(t, vals, level):
    """vals (nt, nk): first index where max_k |vals| >= level; returns tau, k-index."""
    m = np.max(np.abs(vals), axis=1); w = np.flatnonzero(m >= level)
    if not len(w):
        return None
    i = int(w[0]); return i, int(np.argmax(np.abs(vals[i])))


def codex_bracket(t, a, vals, ks, level):
    r = first_cross(t, vals, level)
    if r is None:
        return None
    i, j = r
    return dict(k=float(ks[j]), tau_bracket=[float(t[max(0, i - 1)]), float(t[i])], a_upper=float(a[i]),
                sampled_value=float(np.max(np.abs(vals[i]))))


def env_rel(mine, ref):
    env = np.maximum.accumulate(np.max(np.abs(ref), axis=1))
    return float(np.max(np.abs(mine - ref) / env[:, None]))


def history(job, codexfile, gauge=None):
    M = load(job); C = np.load(CODEX / codexfile)
    tC, ks = C['times'], C['ks']; nk = len(ks)
    assert np.allclose(M['ks'][:nk], ks)
    i = at(M['t'], tC)
    out = dict(job=job, codex=codexfile)
    aC = np.exp(C['bg'][:, 15]); HC = C['bg'][:, 16]
    out['background'] = dict(max_rel_a=float(np.max(np.abs(M['a'][i] / aC - 1))),
                             max_rel_H=float(np.max(np.abs(M['H'][i] / HC - 1))),
                             final_a=float(M['a'][-1]), late_C=float(M['a'][-1] ** 2 * M['H'][-1]))
    # Normalisation exactly as Codex (its sample grid) and on our 0.002 grid
    nC = window_mean(tC, M['NX'][i]); nF = window_mean(M['t'], M['NX'])
    out['NX_norm_codex_grid'] = dict(mean=nC[0], spans=[nC[1], nC[2]], samples=nC[3])
    out['NX_norm_fine_grid'] = dict(mean=nF[0], spans=[nF[1], nF[2]], samples=nF[3], rel_shift=nF[0] / nC[0] - 1)
    out['NX_final_change_from_norm'] = float(M['NX'][-1] / nC[0] - 1)
    norm = abs(nC[0])
    A = M['A'][:, :nk]; Q = M['dN'][:, :nk] / norm
    # Agreement with the archive on its own grid
    out['agreement'] = dict(amplitude_env_rel=env_rel(A[i], C['amplitude_per_radian']),
                            charge_env_rel=env_rel(Q[i], C['charge_fraction_per_radian']),
                            per_k_amplitude_env_rel={('%g' % k): env_rel(A[i][:, [j]], C['amplitude_per_radian'][:, [j]]) for j, k in enumerate(ks)},
                            max_abs_amplitude_ours=float(np.max(np.abs(A))), max_abs_amplitude_codex=float(np.max(np.abs(C['amplitude_per_radian']))))
    # Codex-convention brackets from our values on Codex's grid, and the true first continuous crossings (0.002 grid)
    rows = []
    for s in SEEDS:
        r = dict(seed=s)
        for name, vals, lev in (('amplitude_0.1', s * A, 0.1), ('amplitude_1', s * A, 1.0), ('charge_1', s * Q, 1.0)):
            br = codex_bracket(tC, aC, vals[i], ks, lev)
            fc = first_cross(M['t'], vals, lev)
            fine = None
            if fc is not None:
                ii, jj = fc
                fine = dict(tau=float(M['t'][ii]), k=float(ks[jj]), a=float(M['a'][ii]))
                # all k that cross within the Codex bracket or earlier
                per_k = {}
                for j in range(nk):
                    w = np.flatnonzero(np.abs(vals[:, j]) >= lev)
                    per_k['%g' % ks[j]] = float(M['t'][w[0]]) if len(w) else None
                fine['per_k_first'] = per_k
                if br is not None:
                    fine['inside_codex_bracket'] = bool(br['tau_bracket'][0] < fine['tau'] <= br['tau_bracket'][1] + 1e-9)
                    fine['before_codex_lower'] = bool(fine['tau'] <= br['tau_bracket'][0])
            r[name] = dict(codex_convention_bracket_from_ours=br, first_continuous=fine,
                           max_over_tau=float(np.max(np.abs(vals))))
        rows.append(r)
    out['crossings'] = rows
    # Envelope growth at k = 2.13 between tau 200 and 700 (log ratio of running envelopes)
    j = int(np.argmin(np.abs(ks - 2.13)))
    env = np.maximum.accumulate(np.abs(A[:, j])); t = M['t']
    i1, i2 = np.searchsorted(t, 200.0), np.searchsorted(t, 700.0)
    envC = np.maximum.accumulate(np.abs(C['amplitude_per_radian'][:, j]))
    c1, c2 = np.searchsorted(tC, 200.0), np.searchsorted(tC, 700.0)
    out['growth_k213_200_700'] = dict(ours_fine=float(np.log(env[i2] / env[i1]) / (t[i2] - t[i1])),
                                     ours_codex_grid=float(np.log(np.maximum.accumulate(np.abs(A[i][:, j]))[c2] / np.maximum.accumulate(np.abs(A[i][:, j]))[c1]) / (tC[c2] - tC[c1])),
                                     codex=float(np.log(envC[c2] / envC[c1]) / (tC[c2] - tC[c1])))
    # Undersampling diagnostics: how much the archived grid under-reads the running maximum before tau 80 and after
    early = t <= 80
    out['undersampling'] = dict(
        max_abs_A_k213_fine_before80=float(np.max(np.abs(A[early, j]))),
        max_abs_A_k213_codexgrid_before80=float(np.max(np.abs(A[i][tC <= 80, j]))),
        max_over_k_fine_before80=float(np.max(np.abs(A[early]))),
        max_over_k_codexgrid_before80=float(np.max(np.abs(A[i][tC <= 80]))))
    return out, M, i, norm


def old_archive(M):
    old = np.load(CODEX.parent / 'grut_reheating_modes_20261026' / 'PHASE_PROFILES.npz')
    ks = old['b3_g0.5_ks']; idx = int(np.argmin(np.abs(ks - KOLD))); arc = old['b3_g0.5'][:, idx, :].astype(float)
    to = old['times']; i = at(M['t'], to)
    j = int(np.argmin(np.abs(M['ks'] - KOLD)))
    amp = M['A'][i, j]; dN = M['dN'][i, j]
    # Codex verify_response.py convention: mean over the old grid with times >= 64
    w = to >= 64; mean = np.trapezoid(M['NX'][i][w], to[w]) / (to[w][-1] - to[w][0])
    ch = dN / abs(mean)
    fit = float(np.dot(ch, arc[:, 1]) / np.dot(ch, ch))
    scale = lambda d, r: float(np.max(np.abs(d) / (1 + np.abs(r))))
    return dict(k_archive=float(ks[idx]), amplitude_scaled=scale(amp - arc[:, 0], arc[:, 0]),
                charge_same_sign_scaled=scale(ch - arc[:, 1], arc[:, 1]), charge_opposite_sign_scaled=scale(ch + arc[:, 1], arc[:, 1]),
                least_squares_ratio_archive_over_ours=fit, NX_mean_sign=float(np.sign(mean)))


def main():
    out = dict(method='sk_tangent.py: cosmic-time physical-field tangent, dual-number Hessian-vector products, simplified '
                      'geometry (f = 1), e = 0; 0.002 tau grid plus Codex sample times; all 18 momenta')
    q3, M3, i3, n3 = history('q3', 'quiet3.npz')
    q5, M5, i5, n5 = history('q5', 'quiet5.npz')
    q3['old_phase_archive'] = old_archive(M3)
    # gauge3 (e = 1e-6) vs our e = 0: size of the gauge effect in Codex's archive vs the e = 0 / e = 1e-7 difference
    G = np.load(CODEX / 'gauge3.npz'); Cq = np.load(CODEX / 'quiet3.npz')
    sel = [int(np.argmin(np.abs(Cq['ks'] - k))) for k in G['ks']]
    env = np.maximum.accumulate(np.abs(Cq['amplitude_per_radian'][:, sel]), axis=0)
    q3['codex_gauge_effect_1e-6_vs_1e-7'] = {('%g' % k): float(np.max(np.abs(G['amplitude_per_radian'][:, n] - Cq['amplitude_per_radian'][:, sel[n]]) / env[:, n]))
                                             for n, k in enumerate(G['ks'])}
    A3 = M3['A'][i3][:, sel]
    q3['ours_e0_vs_codex_e1e-7_selected'] = {('%g' % k): float(np.max(np.abs(A3[:, n] - Cq['amplitude_per_radian'][:, sel[n]]) / env[:, n]))
                                             for n, k in enumerate(G['ks'])}
    out['quiet3'] = q3; out['quiet5'] = q5
    RES.mkdir(parents=True, exist_ok=True)
    (RES / 'sk_response.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
