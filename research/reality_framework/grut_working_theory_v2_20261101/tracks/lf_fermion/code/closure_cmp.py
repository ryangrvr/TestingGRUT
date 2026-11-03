"""Non-back-reacting comparison: fermion production on recorded lattice lines (lfline) against the closure family
(lf3d passive tracers, and a Python copy of the same closure on the line-mean field for the homogeneous control).
Claude, lf_fermion (issue #13 context). Usage: python3 -B closure_cmp.py RUN -> ../results/cmp_<RUN>.json"""
import json, math, sys, glob
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
import lines_io, lfcommon as LC

SPI = {'N1': 0, 'N2': 1, 'N3': 2, 'N4': 3}
LAM = [0.0, 0.25, 0.5, 0.75, 1.0]


def par(rd):
    P = {}
    for l in (rd / 'run.par').read_text().splitlines():
        sp = l.split()
        if len(sp) == 2:
            try:
                P[sp[0]] = float(sp[1])
            except ValueError:
                pass
    return P


def tables(rd):
    lines = (rd / 'tables.txt').read_text().split('\n')
    n, u0, du = lines[0].split()
    n, u0, du = int(n), float(u0), float(du)
    ch, pos = [], 1
    for j in range(4):
        w, Gp, Gr, M = lines[pos].split()
        arr = np.array([[float(x) for x in lines[pos + 1 + k].split()] for k in range(n)])
        ch.append(dict(w=int(w), Gp=float(Gp), M=float(M), G=arr[:, 0], F=arr[:, 1]))
        pos += 1 + n
    return u0, du, ch


def Vpot(P, s, x):
    d = s - P['sct']
    V = 0.25 * P['Lv2'] * (x * x - 1) ** 2 + 0.5 * P['Gv2'] * s * s * x * x - P['cbv'] * x + P['c1v'] * d - 0.5 * P['m2c'] * d * d + P['Vc']
    M2 = P['Gv2'] * s * s - P['Lv2']
    if s > P['sct'] and M2 != 0:
        V += M2 * M2 * (math.log(abs(M2) / P['mu2']) - 1.5) / (64 * math.pi ** 2 * P['vH'] ** 2)
    return V


def homog_closure(rd, D, line):
    """the lf3d closure (P0 = GM for a homogeneous field) on the line-mean field, passive, no decays (as the tracers)."""
    P = par(rd)
    u0, du, ch = tables(rd)
    t, a, H = D['t'], D['a'], D['H']
    F = D['F'][:, line].mean(axis=1)            # sigma, sigma_dot, chi, chi_dot (line mean)
    B = np.zeros(4)
    TE = np.zeros(4)
    out = np.zeros((len(t), 4))
    for n in range(len(t) - 1):
        s, sd, x, xd = F[n]
        Et = 0.5 * (sd * sd + xd * xd) + Vpot(P, s, x)
        Es = max(0.5 * sd * sd + 0.5 * P['Gv2'] * x * x * (s - P['s0t']) ** 2, 0.0)
        Ex = max(Et - Es, 0.0)
        dt = t[n + 1] - t[n]
        for j, c in enumerate(ch):
            E = Ex if c['w'] else Es
            vel2 = xd * xd if c['w'] else sd * sd
            if E <= 1e-300:
                continue
            u = math.log(E)
            G = float(np.interp(u, u0 + du * np.arange(len(c['G'])), c['G']))
            cap = float(np.interp(u, u0 + du * np.arange(len(c['F'])), c['F'])) * E
            bl = max(0.0, 1 - B[j] / cap) if cap > 0 else 0.0
            dE = max(G * bl, c['Gp']) * vel2 * dt
            B[j] = B[j] * math.exp(-6 * H[n] * dt) + dE
            TE[j] = TE[j] * math.exp(-4 * H[n] * dt) + dE
        out[n + 1] = TE
    return t, out


def load_q(run, tag, sp, v):
    res = {}
    for f in sorted(glob.glob(str(LC.SP / 'lfq' / run / ('%s_l*_%s.json' % (tag, sp))))):
        d = json.loads(Path(f).read_text())
        m = d['meas']
        t = np.array([r['t'] for r in m])
        for nm in [k for k in m[0] if k in SPI]:
            for meas in ('e_proj_cut', 'e_chir_cut'):
                e = np.array([r[nm][meas] for r in m])
                b = e[t < 0.15].mean()
                res.setdefault((nm, meas), {})[d['line']] = (t, (e - b) / v ** 2, b / v ** 2)
    return res


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else 'REC_P0'
    rd = LC.RUNS / run
    v = json.loads((rd / 'info.json').read_text())['P']['v']
    D = lines_io.read_lines(rd)
    T = lines_io.read_trlines(rd, D['N'], D['nl'])
    tg = np.round(np.arange(0.05, 1.1001, 0.05), 3)
    out = dict(run=run, label='C', note=__doc__, t=list(tg), species={}, homog={})
    for sp in ('chi', 'sig'):
        Q = load_q(run, 'base', sp, v)
        Qh = load_q(run, 'hom', sp, v)
        names = sorted({k[0] for k in Q})
        for nm in names:
            j = SPI[nm]
            lines = sorted(Q[(nm, 'e_proj_cut')].keys())
            q_by_line, ch_by_line, tr_by_line = [], [], []
            for l in lines:
                tq, eq, _ = Q[(nm, 'e_proj_cut')][l]
                q_by_line.append(np.interp(tg, tq, eq))
                tq2, ec, _ = Q[(nm, 'e_chir_cut')][l]
                ch_by_line.append(np.interp(tg, tq2, ec))
                tr = T['TE'][:, :, j, l, :].mean(-1)       # (nt, NL)
                tr_by_line.append(np.array([np.interp(tg, T['t'], tr[:, i]) for i in range(len(LAM))]).T)
            q = np.array(q_by_line); chq = np.array(ch_by_line); tr = np.array(tr_by_line)   # (nl, nt), (nl, nt, NL)
            qm, qs = q.mean(0), q.std(0) / math.sqrt(len(lines))
            trm = tr.mean(0)
            # output-level interpolation parameter: lambda_out = (E_P0 - E_Q)/(E_P0 - E_GM)
            lam_out = (trm[:, 0] - qm) / np.where(np.abs(trm[:, 0] - trm[:, 4]) > 1e-12, trm[:, 0] - trm[:, 4], np.nan)
            out['species'][nm] = dict(lines=lines, Q_mean=list(qm), Q_sem=list(qs), Q_chir_mean=list(chq.mean(0)),
                                      tracer_mean={str(L): list(trm[:, i]) for i, L in enumerate(LAM)},
                                      ratio_Q_P0=list(qm / trm[:, 0]), ratio_Q_GM=list(qm / trm[:, 4]),
                                      lambda_out=list(lam_out), Q_by_line=[list(x) for x in q],
                                      P0_by_line=[list(x) for x in tr[:, :, 0]])
            if (nm, 'e_proj_cut') in Qh:
                for l, (tq, eq, base) in Qh[(nm, 'e_proj_cut')].items():
                    th, hc = homog_closure(rd, D, l)
                    out['homog'][nm] = dict(line=l, Q_hom=list(np.interp(tg, tq, eq)), closure_hom=list(np.interp(tg, th, hc[:, j])),
                                            Q_inhom_same_line=list(np.interp(tg, *Q[(nm, 'e_proj_cut')][l][:2])) if l in Q[(nm, 'e_proj_cut')] else None)
    Path(__file__).resolve().parents[1].joinpath('results', 'cmp_%s.json' % run).write_text(json.dumps(out, indent=1, default=float))
    for nm, d in out['species'].items():
        print(nm, 'lines', d['lines'])
        for i, tt in enumerate(tg):
            if abs(tt * 20 - round(tt * 20)) < 1e-6 and round(tt * 100) % 10 == 0 or tt in (0.25, 0.35, 0.45):
                h = out['homog'].get(nm)
                hs = '' if h is None else ' | hom: Q %.3g clos %.3g' % (h['Q_hom'][i], h['closure_hom'][i])
                print('  t %.2f Q %.4g +- %.2g (chir %.4g)  P0 %.4g  l.5 %.4g  GM %.4g  Q/P0 %.2f Q/GM %.2f lam_out %.2f%s' % (
                    tt, d['Q_mean'][i], d['Q_sem'][i], d['Q_chir_mean'][i], d['tracer_mean']['0.0'][i], d['tracer_mean']['0.5'][i],
                    d['tracer_mean']['1.0'][i], d['ratio_Q_P0'][i], d['ratio_Q_GM'][i], d['lambda_out'][i], hs))


if __name__ == '__main__':
    main()
