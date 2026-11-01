"""Compare the skeptic's runs with the track's per-run files (same orbit, eps, phase, M0 = 0.5 Mperi, epst = 0.002,
exponential growth, M_end = 40 Mapo).  Prints and saves a run-by-run table."""
import glob, json, os, numpy as np
from sk_common import *

rows = []
for fn in sorted(glob.glob(os.path.join(RUNS, '*_g0_et0.002_me40_ms0.5_n32_m8.json'))):
    d = json.load(open(fn)); j = d['job']
    pat = os.path.join(TRACK_RUNS, f"{j['orbit']}_m3_et0.002_b0_phi0qs_eps{j['eps']:.0e}_ph{j['phase']:.6f}_rt1e-11.json")
    tr = json.load(open(pat)) if os.path.exists(pat) else None
    row = dict(orbit=j['orbit'], eps=j['eps'], phase=j['phase'], sk_nminus=d['nminus'], sk_Nphi=d['Nphi'], sk_hext=d['max_abs_hext'],
               sk_qrel=d['max_qrel'], sk_steps=d['steps'], J0=d['J0_per_q'],
               sk_n10=d['windows'].get('10', {}).get('nminus'), sk_n5=d['windows'].get('5', {}).get('nminus'))
    if tr:
        row.update(tr_nminus=tr['nminus'], tr_Nphi=tr['Nphi'], tr_work=tr['work_rel'], tr_qrel=tr['qrel'], tr_steps=tr['steps'],
                   d_nminus=d['nminus'] - tr['nminus'], d_Nphi=d['Nphi'] - tr['Nphi'], tr_n10=tr['windows']['10']['nminus'])
    rows.append(row)
rows.sort(key=lambda r: (r['orbit'], -r['eps'], r['phase']))
print(f"{'orbit':6s} {'eps':8s} {'phase':6s} | {'sk n-/q':10s} {'track n-/q':10s} {'diff':9s} | {'sk Nphi/q':10s} {'track':10s} {'diff':9s} | sk Hext  sk q   | (J0-n)/N sk")
for r in rows:
    if 'tr_nminus' in r:
        rat = (r['J0'] - r['sk_nminus']) / r['sk_Nphi'] if r['sk_Nphi'] > 1e-12 else float('nan')
        print(f"{r['orbit']:6s} {r['eps']:8.0e} {r['phase']:6.3f} | {r['sk_nminus']:10.5f} {r['tr_nminus']:10.5f} {r['d_nminus']:+9.1e} | {r['sk_Nphi']:10.5f} {r['tr_Nphi']:10.5f} {r['d_Nphi']:+9.1e} | {r['sk_hext']:7.0e} {r['sk_qrel']:7.0e} | {rat:7.4f}")
    else:
        print(f"{r['orbit']:6s} {r['eps']:8.0e} {r['phase']:6.3f} | {r['sk_nminus']:10.5f} {'(no track run)':10s} {'':9s} | {r['sk_Nphi']:10.5f}")
save_json(os.path.join(RESULTS, 'sk_reproduction.json'), rows)
