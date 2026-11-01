"""Skeptic part E: aggregate every transition run (own DOP853 code, plus the track-code cross-run on the skeptic's random
design) into results/skC_transition.json.  Adds: per-design statistics, pooled slow-limit values, seed-trend fits,
eps trend, energy-spread trend (ergodic test), uniform-on-[0, n_max] test of the n_- distribution."""
import sys, os, json, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *
from scipy import stats as sst

RUNS = os.path.join(SCR, 'runs'); TCODE = os.path.join(SCR, 'trackcode')
R = np.sqrt(2 * LP) / GG
J0 = dict(tau20=1.7637252793969966, late=0.48049184685032387)


def load(pat, d=RUNS):
    out = []
    for fn in sorted(glob.glob(os.path.join(d, pat))):
        r = json.load(open(fn)); r['file'] = os.path.basename(fn); out.append(r)
    return out


def design(r):
    sd, ph = r['seed'], r['phase']
    on16 = abs(ph * 16 - round(ph * 16)) < 1e-9
    if on16 and any(abs(sd / s - 1) < 1e-9 for s in [1e-2, 1e-6, 1e-10]): return 'grid16'
    if on16: return f'seed1e-{-np.log10(sd):.0f}'
    return 'rand'


def st(a):
    a = np.asarray(a, float)
    return dict(n=len(a), mean=a.mean(), sem=a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0,
                std=a.std(ddof=1) if len(a) > 1 else 0.0, min=a.min(), max=a.max())


def block(rows, orb):
    n = np.array([r['nminus'] for r in rows]); N = np.array([r['Nphi'] for r in rows])
    E = 2 * n + R * N
    esc = [r['M_escape_over_Mperi'] for r in rows if r.get('M_escape_over_Mperi') is not None]
    return dict(nminus=st(n), Nphi=st(N), sum=st(n + N), E_final_over_mXq=st(E), nminus_over_J0=float(n.mean() / J0[orb]),
                Nphi_over_J0=float(N.mean() / J0[orb]), M_escape_over_Mperi=st(esc) if esc else None,
                max_qrel=max(r['qrel'] for r in rows), max_work_rel=max(r['work_rel'] for r in rows))


out = dict(method='own adaptive DOP853 C code (sk_dop.c), rtol 1e-11; per-run files in $SP/review/pt_matching_skeptic/runs',
           groups={}, pooled={}, seed_trend={}, eps_trend={}, uniform_test={}, track_code_cross_run={}, bias={}, rd_growth={})
allrows = load('*.json')
for r in allrows:
    r['design'] = design(r)
    r['mend'] = r['file'].split('_mend')[1].replace('.json', '') if '_mend' in r['file'] else '40'
keys = sorted(set((r['orbit'], r['growth'], r['beta'], r['eps'], r['design'], r['mend']) for r in allrows))
for orb, gr, be, eps, des, me in keys:
    rows = [r for r in allrows if (r['orbit'], r['growth'], r['beta'], r['eps'], r['design'], r['mend']) == (orb, gr, be, eps, des, me)]
    out['groups'][f'{orb}|growth{gr}|beta{be:g}|eps{eps:.0e}|{des}|mend{me}'] = block(rows, orb)

base = lambda orb: [r for r in allrows if r['orbit'] == orb and r['growth'] == 0 and r['beta'] == 0]
for orb in ['tau20', 'late']:
    rows = [r for r in base(orb) if r['eps'] <= 1e-3 and r['design'] in ('rand', 'grid16')]
    out['pooled'][orb + '_eps_le_1e-3_rand_and_grid'] = block(rows, orb)
    # seed trend inside the random design (log-uniform seeds 1e-12..1e-2) and across the 16-phase tiny-seed sets
    rr = [r for r in rows if r['design'] == 'rand']
    if rr:
        x = np.log10([r['seed'] for r in rr]); y = np.array([r['nminus'] for r in rr])
        lr = sst.linregress(x, y)
        out['seed_trend'][orb] = dict(slope_per_decade=lr.slope, sem=lr.stderr, n=len(rr))
    # eps trend over all eps <= 1e-3
    x = np.log10([r['eps'] for r in rows]); y = np.array([r['nminus'] for r in rows])
    if len(set(x)) > 1:
        lr = sst.linregress(x, y); out['eps_trend'][orb] = dict(slope_per_decade=lr.slope, sem=lr.stderr, n=len(rows))
    # uniform-on-[0, a] test of the n_- distribution (microcanonical exit picture), a = max over the pool
    y = np.array([r['nminus'] for r in rows]); a = y.max()
    ks = sst.kstest(y / a, 'uniform')
    out['uniform_test'][orb] = dict(mean_over_max=float(y.mean() / a), std_over_max=float(y.std() / a), ks_stat=ks.statistic,
                                    ks_p=ks.pvalue, uniform_expect=dict(mean_over_max=0.5, std_over_max=1 / np.sqrt(12)),
                                    corr_nminus_Nphi=float(np.corrcoef(y, [r['Nphi'] for r in rows])[0, 1]))
# eps-resolved table incl. small eps (MEND 10)
for orb in ['tau20', 'late']:
    tab = []
    for eps in sorted(set(r['eps'] for r in base(orb)), reverse=True):
        rows = [r for r in base(orb) if r['eps'] == eps and not r['design'].startswith('seed')]
        b = block(rows, orb); tab.append(dict(eps=eps, n=b['nminus']['n'], nminus=b['nminus']['mean'], sem=b['nminus']['sem'],
                                             Nphi=b['Nphi']['mean'], E_std=b['E_final_over_mXq']['std'],
                                             M_escape=b['M_escape_over_Mperi']['mean'] if b['M_escape_over_Mperi'] else None))
    out['eps_trend'][orb + '_table'] = tab
# tiny seeds
for s in ['seed1e-20', 'seed1e-30', 'seed1e-60']:
    rows = [r for r in base('tau20') if r['design'] == s]
    if rows: out['seed_trend']['tau20_eps1e-3_' + s] = block(rows, 'tau20')
# bias and RD growth
for be in sorted(set(r['beta'] for r in allrows if r['beta'] != 0)):
    for eps in [1e-3, 1e-4, 1e-5]:
        rows = [r for r in allrows if r['beta'] == be and r['eps'] == eps]
        if rows: out['bias'][f'beta{be:g}_eps{eps:.0e}'] = block(rows, 'tau20')
for me in ['5', '40']:
    rows = [r for r in allrows if r['growth'] == 1 and r['mend'] == me]
    if rows: out['rd_growth'][f'eps2e-3_mend{me}'] = block(rows, 'tau20')
# track-code cross-run (pt_matching's own integrator run read-only on the skeptic random design)
for orb in ['late', 'tau20']:
    rows = load(f'{orb}_eps1e-03_*.json', TCODE)
    if rows:
        n = np.array([r['nminus'] for r in rows]); N = np.array([r['Nphi'] for r in rows])
        out['track_code_cross_run'][orb + '_eps1e-3_rand'] = dict(nminus=st(n), Nphi=st(N))
# the track's own published per-run values (grid design) for comparison
for orb, eps in [('tau20', 1e-3), ('tau20', 1e-4), ('late', 1e-3), ('late', 1e-4)]:
    rows = [json.load(open(f)) for f in glob.glob(f'{TRACK_SCAN}/{orb}_g0_eps{eps:.0e}_ph*_h0.08.json')]
    out['track_code_cross_run'][f'track_published_{orb}_eps{eps:.0e}_grid16'] = dict(nminus=st([r['nminus'] for r in rows]),
                                                                                        Nphi=st([r['Nphi'] for r in rows]))
late_both = [r['nminus'] for r in base('late') if r['eps'] == 1e-3 and r['design'] == 'rand'] + \
            [r['nminus'] for r in load('late_eps1e-03_*.json', TCODE)]
out['late_eps1e-3_rand_both_codes'] = st(late_both)
late_bothN = [r['Nphi'] for r in base('late') if r['eps'] == 1e-3 and r['design'] == 'rand'] + \
             [r['Nphi'] for r in load('late_eps1e-03_*.json', TCODE)]
out['late_eps1e-3_rand_both_codes_Nphi'] = st(late_bothN)
save_json('skC_transition.json', out)
for k, v in out['pooled'].items():
    print(k, {kk: (round(vv['mean'], 4), round(vv['sem'], 4)) if isinstance(vv, dict) and 'mean' in vv else vv for kk, vv in v.items() if kk in ('nminus', 'Nphi', 'sum', 'nminus_over_J0', 'Nphi_over_J0')})
print('seed', out['seed_trend'].get('tau20'), out['seed_trend'].get('late'))
print('eps', out['eps_trend'].get('tau20'), out['eps_trend'].get('late'))
for t in out['eps_trend']['tau20_table']: print('  ', t)
for t in out['eps_trend']['late_table']: print('  ', t)
print('uniform', out['uniform_test'])
print('late both codes', out['late_eps1e-3_rand_both_codes'], out['late_eps1e-3_rand_both_codes_Nphi'])
for k, v in out['bias'].items(): print('bias', k, round(v['nminus']['mean'], 4), round(v['nminus']['sem'], 4), round(v['Nphi']['mean'], 4))
for k, v in out['rd_growth'].items(): print('rd', k, round(v['nminus']['mean'], 4), round(v['nminus']['sem'], 4), v['nminus']['n'])
