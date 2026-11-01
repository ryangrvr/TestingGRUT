"""Aggregate skeptic transition runs; compare with the track's per-run files (run-by-run where non-chaotic)."""
import sys, os, json, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *

OUT = os.path.join(SCR, 'runs')
R = np.sqrt(2 * LP) / GG   # m_phi / m_X


def load(pat):
    rows = []
    for fn in sorted(glob.glob(os.path.join(OUT, pat))):
        r = json.load(open(fn)); r['file'] = os.path.basename(fn); rows.append(r)
    return rows


def stats(a):
    a = np.asarray(a, float)
    if len(a) == 0: return None
    return dict(n=len(a), mean=float(a.mean()), sem=float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0,
                std=float(a.std(ddof=1)) if len(a) > 1 else 0.0, min=float(a.min()), max=float(a.max()),
                median=float(np.median(a)))


def summary(rows):
    n = np.array([r['nminus'] for r in rows]); N = np.array([r['Nphi'] for r in rows])
    esc = [r['M_escape_over_Mperi'] for r in rows if r['M_escape_over_Mperi'] is not None]
    w5 = np.array([r['windows']['5']['nminus'] for r in rows])
    return dict(nminus=stats(n), Nphi=stats(N), sum=stats(n + N), E_final_per_mXq=stats(2 * n + R * N),
                M_escape_over_Mperi=stats(esc), max_qrel=max(r['qrel'] for r in rows), max_work_rel=max(r['work_rel'] for r in rows),
                max_abs_n40_minus_n5=float(np.max(np.abs(n - w5))),
                max_abs_actualmass_minus_gM=float(max(abs(r['windows']['40']['nminus_actualmass'] - r['windows']['40']['nminus']) for r in rows)))


def track_rows(orb, eps, h='0.08'):
    rows = []
    for fn in glob.glob(f'{TRACK_SCAN}/{orb}_g0_eps{eps:.0e}_ph*_h{h}.json'):
        r = json.load(open(fn)); rows.append(r)
    return rows


def compare_runbyrun(orb, eps):
    mine = load(f'{orb}_g0_b0_eps{eps:.0e}_ph*_rt1e-11.json')
    trk = {(round(r['phase'], 4), f"{r['seed']:.0e}"): r for r in track_rows(orb, eps)}
    d = []
    for r in mine:
        key = (round(r['phase'], 4), f"{r['seed']:.0e}")
        if key in trk:
            d.append((r['nminus'] - trk[key]['nminus'], r['Nphi'] - trk[key]['Nphi']))
    d = np.array(d)
    return dict(n=len(d), max_abs_dn=float(np.max(np.abs(d[:, 0]))), max_abs_dN=float(np.max(np.abs(d[:, 1]))),
                median_abs_dn=float(np.median(np.abs(d[:, 0])))) if len(d) else None


if __name__ == '__main__':
    out = {}
    groups = sorted(set(r['file'].split('_ph')[0] + '|' + r['file'].split('_rt')[1] for r in load('*.json')))
    for g in groups:
        pre, rt = g.split('|')
        rows = load(f'{pre}_ph*_rt{rt}')
        out[pre + '_rt' + rt.replace('.json', '')] = summary(rows)
        s = out[pre + '_rt' + rt.replace('.json', '')]
        print(f"{pre:32s} n={s['nminus']['n']:4d} n-/q {s['nminus']['mean']:.4f} +- {s['nminus']['sem']:.4f} (sd {s['nminus']['std']:.3f}, max {s['nminus']['max']:.3f}) "
              f"Nphi {s['Nphi']['mean']:.4f} +- {s['Nphi']['sem']:.4f} | E {s['E_final_per_mXq']['mean']:.3f} sd {s['E_final_per_mXq']['std']:.3f} | "
              f"esc {s['M_escape_over_Mperi']['mean'] if s['M_escape_over_Mperi'] else float('nan'):.3f} | q {s['max_qrel']:.1e} W {s['max_work_rel']:.1e} | n40-n5 {s['max_abs_n40_minus_n5']:.3f}")
    for orb, epss in [('tau20', [1e-2, 1e-3, 1e-4, 3e-5]), ('late', [1e-2, 1e-3, 1e-4])]:
        for e in epss:
            trk = track_rows(orb, e)
            if trk:
                c = compare_runbyrun(orb, e)
                ts = dict(nminus=stats([r['nminus'] for r in trk]), Nphi=stats([r['Nphi'] for r in trk]))
                out[f'track_{orb}_{e:.0e}'] = dict(track=ts, runbyrun=c)
                print('track', orb, e, f"n- {ts['nminus']['mean']:.4f} +- {ts['nminus']['sem']:.4f}  Nphi {ts['Nphi']['mean']:.4f}", 'run-by-run', c)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        save_json('skC_transition.json', out)
