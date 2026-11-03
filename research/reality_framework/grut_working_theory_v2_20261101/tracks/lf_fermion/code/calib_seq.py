"""Sequential (instantaneous) calibration of the capacity multiplier kappa_j(t) of the local closure: window by window
(0.05 H^-1), with the history before the window fixed, kappa is chosen so that the line-mean closure energy at the
window end equals the line fermion energy (lfline, non-back-reacting, same recorded background).  Model choice.
Claude, lf_fermion (issue #13 context). Usage: python3 -B calib_seq.py RUN -> ../results/calib_seq_<RUN>.json"""
import json, math, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
import lines_io, lfcommon as LC
import closure_cmp as CC

SPI = {'N1': 0, 'N2': 1, 'N3': 2, 'N4': 3}


def prep(rd, D, lines, t_end, homog=False):
    P = CC.par(rd)
    u0, du, ch = CC.tables(rd)
    ug = u0 + du * np.arange(len(ch[0]['G']))
    t, H = D['t'], D['H']
    nt = int(np.searchsorted(t, t_end)) + 1
    F = D['F'][:nt][:, lines]
    if homog:
        F = np.repeat(F.mean(axis=2, keepdims=True), F.shape[2], axis=2)
    s, sd, x, xd = F[..., 0], F[..., 1], F[..., 2], F[..., 3]
    d = s - P['sct']
    V = 0.25 * P['Lv2'] * (x * x - 1) ** 2 + 0.5 * P['Gv2'] * s * s * x * x - P['cbv'] * x + P['c1v'] * d - 0.5 * P['m2c'] * d * d + P['Vc']
    M2 = P['Gv2'] * s * s - P['Lv2']
    with np.errstate(divide='ignore', invalid='ignore'):
        V = V + np.where((s > P['sct']) & (M2 != 0), M2 * M2 * (np.log(np.abs(M2) / P['mu2']) - 1.5) / (64 * math.pi ** 2 * P['vH'] ** 2), 0.0)
    Et = 0.5 * (sd * sd + xd * xd) + V
    Es = np.maximum(0.5 * sd * sd + 0.5 * P['Gv2'] * x * x * (s - P['s0t']) ** 2, 0.0)
    Ex = np.maximum(Et - Es, 0.0)
    out = []
    for c in ch:
        E = Ex if c['w'] else Es
        vel2 = xd * xd if c['w'] else sd * sd
        uE = np.log(np.maximum(E, 1e-300))
        out.append(dict(G=np.where(E > 1e-300, np.interp(uE, ug, c['G']), 0.0), cap1=np.interp(uE, ug, c['F']) * E,
                        vel2=vel2, Gp=c['Gp'], on=E > 1e-300))
    return t[:nt], H[:nt], out


def advance(ch, t, H, n0, n1, kap, B, TE):
    for n in range(n0, n1):
        dt = t[n + 1] - t[n]
        cap = kap * ch['cap1'][n]
        bl = np.where(cap > 0, np.maximum(0.0, 1 - B / np.where(cap > 0, cap, 1)), 0.0)
        dE = np.where(ch['on'][n], np.maximum(ch['G'][n] * bl, ch['Gp']), 0.0) * ch['vel2'][n] * dt
        B = B * math.exp(-6 * H[n] * dt) + dE
        TE = TE * math.exp(-4 * H[n] * dt) + dE
    return B, TE


def seq_fit(t, H, ch, tq, q, tw):
    B = np.zeros(ch['G'].shape[1:])
    TE = np.zeros_like(B)
    n_prev = 0
    kaps, tes = [], []
    for i, tb in enumerate(tw):
        n1 = min(int(np.searchsorted(t, tb)), len(t) - 1)
        target = float(np.interp(t[n1], tq, q))
        def f(lk):
            _, te = advance(ch, t, H, n_prev, n1, math.exp(lk), B, TE)
            return te.mean() - target
        lo, hi = math.log(0.1), math.log(1e4)
        flo, fhi = f(lo), f(hi)
        if flo >= 0:
            lk = lo
        elif fhi <= 0:
            lk = hi
        else:
            for _ in range(22):
                mid = 0.5 * (lo + hi)
                if f(mid) > 0:
                    hi = mid
                else:
                    lo = mid
            lk = 0.5 * (lo + hi)
        B, TE = advance(ch, t, H, n_prev, n1, math.exp(lk), B, TE)
        kaps.append(math.exp(lk))
        tes.append(float(TE.mean()))
        n_prev = n1
    return kaps, tes


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else 'REC_P0'
    rd = LC.RUNS / run
    D = lines_io.read_lines(rd)
    cmpd = json.loads((LC.TRACK / 'results' / ('cmp_%s.json' % run)).read_text())
    tg = np.array(cmpd['t'])
    tw = [x for x in tg if x >= 0.2]
    res = dict(run=run, label='C', note=__doc__, t_window_end=tw, species={})
    for nm, d in cmpd['species'].items():
        j = SPI[nm]
        t, H, ch = prep(rd, D, d['lines'], 1.1)
        k_in, te_in = seq_fit(t, H, ch[j], tg, np.array(d['Q_mean']), tw)
        # +-1 sem bands
        k_lo, _ = seq_fit(t, H, ch[j], tg, np.array(d['Q_mean']) - np.array(d['Q_sem']), tw)
        k_hi, _ = seq_fit(t, H, ch[j], tg, np.array(d['Q_mean']) + np.array(d['Q_sem']), tw)
        row = dict(kappa_inhom=k_in, kappa_inhom_minus_sem=k_lo, kappa_inhom_plus_sem=k_hi, closure_energy=te_in)
        if nm in cmpd['homog']:
            h = cmpd['homog'][nm]
            th, Hh, chh = prep(rd, D, [h['line']], 1.1, homog=True)
            k_h, _ = seq_fit(th, Hh, chh[j], tg, np.array(h['Q_hom']), tw)
            row['kappa_homog'] = k_h
        res['species'][nm] = row
        print(nm, 'inhom:', ' '.join('%.2f:%.2g' % (a, b) for a, b in zip(tw, k_in)))
        if 'kappa_homog' in row:
            print('   homog:', ' '.join('%.2f:%.2g' % (a, b) for a, b in zip(tw, row['kappa_homog'])))
        print('   -sem :', ' '.join('%.2g' % b for b in k_lo))
        print('   +sem :', ' '.join('%.2g' % b for b in k_hi))
    Path(__file__).resolve().parents[1].joinpath('results', 'calib_seq_%s.json' % run).write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
