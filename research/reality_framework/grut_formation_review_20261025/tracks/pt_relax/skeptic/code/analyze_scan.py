"""pt_relax skeptic (Claude, own code): cutoff dependence of the classical A_- drift.

Reads $RUNDIR/L12_N<N>/part*.json written by scan.py (N = 32, 40, 48, 64) and, as a labelled cross-check,
pt_relax's results/long_run.json (N32 to s = 1004, N48 to s = 254).  Writes results/cutoff_scan.json.
"""
import os, json, glob
import numpy as np

RUN = os.environ['RUNDIR']
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results')
PTR = os.path.join(HERE, '..', '..', 'results', 'long_run.json')


def fd4_max(N, L=12.0):
    dx = L / N
    return float(np.sqrt(3 * 16 / 3 / dx ** 2))          # FD4 symbol maximum, all three axes at pi


def load(N):
    wins, rows, occs = [], [], []
    for f in sorted(glob.glob(os.path.join(RUN, f'L12_N{N}', 'part*.json'))):
        d = json.load(open(f))
        for w in d['wins']:
            if not wins or w['window'][1] > wins[-1]['window'][1] + 1e-9: wins.append(w)
        rows += [r for r in d['rows'] if not rows or r['s'] > rows[-1]['s'] + 1e-9]
        occs += [o for o in d['occs'] if not occs or o['s'] > occs[-1]['s'] + 1e-9]
    return wins, rows, occs


def series(wins, key='all', wn='hann'):
    W = [w for w in wins if w['window'][0] >= 24 - 1e-9]
    s = np.array([np.mean(w['window']) for w in W]); A = np.array([w[wn][key] for w in W])
    return s, A


def summary(s, A):
    out = {}
    for a, b in [(29, 99), (99, 199), (199, 249), (29, 249), (249, 349), (349, 449), (449, 499), (249, 499), (29, 499)]:
        m = (s >= a) & (s <= b)
        if m.sum() >= 3 and s[m].max() >= b - 1e-6:
            out[f'dlnA_ds_{a}_{b}'] = float(np.polyfit(s[m], np.log(A[m]), 1)[0])
            out[f'dA_ds_{a}_{b}'] = float(np.polyfit(s[m], A[m], 1)[0])
    for t in [99, 149, 199, 249, 299, 349, 399, 449, 499]:
        i = np.flatnonzero(abs(s - t) < 1e-6)
        if len(i): out[f'drop_29_to_{t}'] = float(A[0] - A[i[0]]); out[f'A_{t}'] = float(A[i[0]])
    out['A_29'] = float(A[0])
    return out


def main():
    res = dict(runs={}, ptrelax_crosscheck={})
    for N in [32, 40, 48, 64]:
        wins, rows, occs = load(N)
        if not wins: continue
        s, A = series(wins); s8, A8 = series(wins, 'tot8'); sb, Ab = series(wins, 'all', 'bh')
        r = dict(N=N, dx=12.0 / N, omega_lat_max=fd4_max(N), s_end=float(rows[-1]['s']),
                 windows=[dict(s=w['window'], hann=w['hann']['all'], hann_tot8=w['hann']['tot8'], coh=w['hann']['coh'],
                               band=w['hann']['band'], uv=w['hann']['uv'], bh=w['bh']['all'], charge=w['hann']['charge'],
                               uv_quad_frac=w['uv_quad_frac'], phi_frac=w['phi_frac'], dE=w['dE']) for w in wins],
                 all_modes=summary(s, A), tot8=summary(s8, A8),
                 bh_minus_hann_max=float(np.max(np.abs(Ab - A))) if len(Ab) == len(A) else None,
                 max_dE=float(max(abs(w['dE']) for w in wins)), max_charge_err=float(max(abs(w['hann']['charge'] - 1) for w in wins)))
        # cascade front: largest shell with occupation >= 1 and minimum occupation, at s ~ 100, 200, 250
        fr = {}
        for o in occs:
            if any(abs(o['s'] - t) < 1e-6 for t in [34, 100, 150, 200, 250, 300, 400, 500]):
                ox = np.array(o['X']); op = np.array(o['phi'])
                fr[f"s{o['s']:.0f}"] = dict(occX_min=float(ox[1:].min()), occphi_min=float(op[1:].min()),
                                           occX_last=float(ox[-1]), occX_shell_of_min=int(np.argmin(ox[1:]) + 1))
        r['occupations'] = fr
        res['runs'][N] = r
        print(f"N{N} dx {12/N:.4f} wmax {fd4_max(N):.2f} s_end {r['s_end']:.0f} maxdE {r['max_dE']:.1e} chg {r['max_charge_err']:.1e} BH-Hann {r['bh_minus_hann_max']}")
        print('   all ', {k: round(v, 6) for k, v in r['all_modes'].items()})
        print('   tot8', {k: round(v, 6) for k, v in r['tot8'].items()})
        print('   occ ', fr)
    # pt_relax windows (reproduction data from that track) at the same windows, for N32/N48
    if os.path.exists(PTR):
        d = json.load(open(PTR))
        for tag, N in [('IR12_s1_w1', 32), ('IR12_fine_s1_w1', 48)]:
            W = [w for w in d[tag]['windows'] if w['s'][0] >= 24 - 1e-9 and w['s'][1] <= 254 + 1e-9]
            s = np.array([np.mean(w['s']) for w in W]); A = np.array([w['total_all_modes'] for w in W])
            res['ptrelax_crosscheck'][N] = summary(s, A)
            mine = res['runs'].get(N)
            if mine:
                diffs = []
                for w in mine['windows']:
                    j = [x for x in W if abs(x['s'][0] - w['s'][0]) < 1e-6]
                    if j: diffs.append(abs(j[0]['total_all_modes'] - w['hann']))
                res['ptrelax_crosscheck'][N]['max_abs_diff_vs_mine'] = float(max(diffs)) if diffs else None
            print('pt_relax', N, {k: round(v, 6) for k, v in res['ptrelax_crosscheck'][N].items() if v is not None})
    # continuum estimate of the drop(29 -> t) from N = 40, 48, 64 only (N32's late acceleration is its cascade front
    # reaching the lattice edge, omega_max = 10.7): geometric (Aitken) extrapolation of the three-point sequence,
    # and the quantum-unreliable part: energy and A_- in shells with classical occupation < 1 (N64).
    fits = {}
    for t in [99, 149, 199, 249]:
        y = [res['runs'][N]['all_modes'].get(f'drop_29_to_{t}') for N in (40, 48, 64) if N in res['runs']]
        if len(y) == 3 and None not in y:
            d1, d2 = y[1] - y[0], y[2] - y[1]
            r = d2 / d1 if d1 else 0.0
            aitken = y[2] - d2 * r / (1 - r) if 0 < r < 1 else y[2]
            fits[f'drop_29_to_{t}'] = dict(N40=y[0], N48=y[1], N64=y[2], ratio=r, aitken_continuum=aitken)
            print('drop', t, [round(v, 5) for v in y], 'ratio', round(r, 3), 'continuum ~', round(aitken, 5))
    if 64 in res['runs']:
        wins, rows, occs = load(64)
        o = [x for x in occs if abs(x['s'] - 250) < 1e-6][0]
        ox, op = np.array(o['X']), np.array(o['phi']); EXs, EPs = np.array(o['E_X']), np.array(o['E_phi'])
        Et = [r for r in rows if abs(r['s'] - 250) < 1e-6][0]['rho']
        lowX = ox < 1; lowP = op < 1; lowX[0] = lowP[0] = False
        sm = np.array([w for w in wins if abs(w['window'][1] - 254) < 1e-6][0]['hann']['shell_minus'])
        res['N64_quantum_unreliable_s250'] = dict(
            energy_frac_X_shells_occ_lt1=float(EXs[lowX].sum() / Et), energy_frac_phi_shells_occ_lt1=float(EPs[lowP].sum() / Et),
            Aminus_in_X_shells_occ_lt1=float(sm[:len(lowX)][lowX].sum()), first_X_shell_occ_lt1=int(np.argmax(lowX)) if lowX.any() else None)
        print('N64 s250 occ<1:', res['N64_quantum_unreliable_s250'])
    late = {}
    for t in [299, 349, 399, 449, 499]:
        y = {N: res['runs'][N]['all_modes'].get(f'drop_29_to_{t}') for N in (48, 64) if N in res['runs']}
        if len(y) == 2 and None not in y.values():
            late[f'drop_29_to_{t}'] = dict(N48=y[48], N64=y[64], N64_minus_N48=y[64] - y[48])
            print('late drop', t, {k: round(v, 5) for k, v in late[f'drop_29_to_{t}'].items()})
    res['late_N48_vs_N64'] = late
    # same-window split at s = 494-504 and energy bookkeeping on N48/N64
    sp = {}
    for N in (48, 64):
        if N not in res['runs']: continue
        W = [w for w in res['runs'][N]['windows'] if abs(w['s'][1] - 504) < 1e-6]
        if W:
            w = W[0]
            sp[N] = dict(all=w['hann'], tot8=w['hann_tot8'], coh=w['coh'], band=w['band'], uv=w['uv'],
                         phi_frac=w['phi_frac'], uv_quad_frac=w['uv_quad_frac'], dE=w['dE'])
            print('split 494-504 N%d' % N, {k: round(v, 5) for k, v in sp[N].items()})
    res['split_494_504'] = sp
    res['continuum_extrapolation'] = fits
    json.dump(res, open(os.path.join(OUT, 'cutoff_scan.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
