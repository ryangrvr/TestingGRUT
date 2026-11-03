"""Regression R1 (homogeneous, constant amplitude, no expansion; units m = 1): the line solver lf1d against
mf_hybrid's exact single half-period production (results/axial_tables.json) and the mf_hybrid skeptic's multi-period
saturation F_tot(n) (skeptic/results/sk_axial.json, sk_axial_Fsat_vs_A.json).  Read-only reuse of those JSON files.
Usage: python3 -B regress_homog.py [dt] -> ../results/regress_homog.json (Claude, lf_fermion; issue #13 context)."""
import json, math, sys, time
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
import lf1d

TR = Path(__file__).resolve().parents[2]
MFH = TR / 'mf_hybrid'
tab = json.load(open(MFH / 'results' / 'axial_tables.json'))
sk = json.load(open(MFH / 'skeptic' / 'results' / 'sk_axial.json'))['multi_period']
skf = json.load(open(MFH / 'skeptic' / 'results' / 'sk_axial_Fsat_vs_A.json'))


def eps_table(Mh, A):
    key = min(tab, key=lambda s: abs(float(s) - Mh))
    return float(np.exp(np.interp(math.log(A), np.log(tab[key]['A']), np.log(tab[key]['eps']))))


def run(A, Mhs, fm, nhalf, dt, dk=0.04):
    kmax = 1.3 * A + 3 * max(Mhs) + 4
    L = 2 * math.pi / dk
    N = 1 << int(math.ceil(math.log2(2 * kmax / dk)))
    ln = lf1d.Line(N, L, Mhs, n_r=1, theta0=2 * A, phases='ones')
    nst = int(round(math.pi / dt))
    h = math.pi / nst
    S = 2 * A * fm
    rows = []
    t = 0.0
    for j in range(nhalf):
        for i in range(nst):
            ln.step(t, h, 1.0, 2 * A * math.cos(t + 0.5 * h))
            t += h
        th = 2 * A * math.cos(t)
        n = ln.occupation(1.0, th)
        e = ln.eps3d(n, 1.0)
        rows.append(dict(nhalf=j + 1, F_tot=list(e / (0.5 * S * S)), eps=list(e),
                         nmax=list(np.max(n, axis=1)), nbar_band=list([float(np.mean(n[q][np.abs(ln.k) < A])) for q in range(len(Mhs))])))
    return dict(A=A, Mh=Mhs, fm=fm, N=N, L=L, dt=h, rows=rows)


if __name__ == '__main__':
    dt = float(sys.argv[1]) if len(sys.argv) > 1 else 0.004
    out = dict(label='C', note=__doc__, dt=dt, cases=[])
    t0 = time.time()
    # chi species at the slice amplitude A = 42.3 (16 half periods) and at A = 10, 20 (saturation points)
    for A, nh in ((42.3, 16), (20.0, 16), (10.0, 16)):
        r = run(A, [1.6, 0.8, 0.48], 8.0, nh, dt)
        cmp = []
        for q, (nm, Mh) in enumerate((('N2', 1.6), ('N3', 0.8), ('N4', 0.48))):
            F1 = r['rows'][0]['eps'][q]
            e_tab = eps_table(Mh, A)
            Fs_mine = float(np.mean([row['F_tot'][q] for row in r['rows'][6:]]))
            ref = None
            if abs(A - 42.3) < 1e-6:
                rr = sk['chi->' + nm]['rows']
                ref = dict(F_tot=[x['F_tot'] for x in rr], Fsat=float(np.mean([x['F_tot'] for x in rr[6:]])))
            else:
                kk = '%s_A%g' % (nm, A)
                if kk in skf:
                    ref = dict(Fsat=skf[kk]['Fsat'], F1=skf[kk]['F1'])
            cmp.append(dict(species=nm, Mh=Mh, eps1_line=F1, eps1_table=e_tab, ratio_eps1=F1 / e_tab,
                            F_tot_line=[row['F_tot'][q] for row in r['rows']], Fsat_line=Fs_mine, ref=ref))
            print('A %.1f %s eps1 line %.5g table %.5g ratio %.4f | Fsat line %.4f ref %s' % (
                A, nm, F1, e_tab, F1 / e_tab, Fs_mine, None if ref is None else '%.4f' % ref['Fsat']), flush=True)
        r['compare'] = cmp
        out['cases'].append(r)
    # sigma -> N1 at its slice amplitude 12.34 (f/m = 2)
    r = run(12.34, [250. / 600.], 2.0, 16, dt)
    rr = sk['sigma->N1']['rows']
    F1 = r['rows'][0]['eps'][0]
    r['compare'] = [dict(species='N1', eps1_line=F1, eps1_table=eps_table(250 / 600., 12.34),
                         F_tot_line=[row['F_tot'][0] for row in r['rows']],
                         Fsat_line=float(np.mean([row['F_tot'][0] for row in r['rows'][6:]])),
                         ref=dict(F_tot=[x['F_tot'] for x in rr], Fsat=float(np.mean([x['F_tot'] for x in rr[6:]]))))]
    c = r['compare'][0]
    print('N1 A 12.34 eps1 ratio %.4f Fsat line %.4f ref %.4f' % (c['eps1_line'] / c['eps1_table'], c['Fsat_line'], c['ref']['Fsat']))
    out['cases'].append(r)
    out['wall_s'] = time.time() - t0
    Path(__file__).resolve().parents[1].joinpath('results', 'regress_homog%s.json' % ('' if dt == 0.004 else '_dt%g' % dt)).write_text(json.dumps(out, indent=1))
