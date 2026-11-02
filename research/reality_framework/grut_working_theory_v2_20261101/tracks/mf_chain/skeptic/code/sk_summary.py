#!/usr/bin/env python3
"""mf_chain SKEPTIC: collect every number quoted in ../../SKEPTIC.md into ../results/sk_summary.json.
Re-runs the cheap independent checks (release regression on two histories, ~45 s; vector modes at 5 k, ~30 s) and reads the
full-chain outputs of sk_chainrun.py / sk_c0chain.py / sk_idealc0.py from the skeptic's scratch (small summaries copied).
Usage: python3 -B sk_summary.py [--fast]   (--fast skips the re-runs and keeps the previous values)
"""
import json, math, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
TRK = HERE.parents[2]
MFC = TRK / 'mf_chain' / 'results' / 'engine'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain_skeptic')
SKE = SP / 'engine'
KEYS = ('K', 'N_X', 'E_X', 'eps', 'E_over_Ec', 'N_YL', 'N_YH', 'N_sigma', 'S_per_Q_e_inf', 'R_Y', 'e_crit_bulk', 'e_crit_Q2e5')


def chain_summary(path):
    o = json.load(open(path))
    s = {k: o[k] for k in KEYS}
    for lab in ('e1e-6', 'e3e-6'):
        s[lab] = {k: o[lab][k] for k in ('R_vec', 'R_Y', 'R_casc', 'W', 'S_per_Q', 'branch', 'M', 'H_i_req_GeV', 'r')}
    s['c0'] = o['spec'].get('c0'); s['bg'] = o['spec'].get('bg', '').split('/')[-1]
    return s


def main(fast=False):
    out_p = RES / 'sk_summary.json'
    old = json.load(open(out_p)) if out_p.exists() else {}
    out = dict(note='mf_chain skeptic (Claude). Every history, holder amplitude, map and ledger is a model choice; nothing '
                    'validates or excludes anything; data comparisons are retrospective.')
    # 1. independent release vs engine
    if fast and 'release_regression' in old:
        out['release_regression'] = old['release_regression']
    else:
        import sk_release as sr
        rr = {}
        for lab, bench, eng in (('HFa exported', TRK / 'mf_hybrid' / 'results' / 'bench_HFa.json', MFC / 'release' / 'release_HFa_raw.json'),
                                ('HFa Pauli P0 sc 2.40e11', TRK / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_P0_it2.json',
                                 MFC / 'release' / 'release_itHFa_P0_it2.json')):
            o, _ = sr.run(str(bench), T=100.0, tag=lab)
            e = json.load(open(eng))['dense']['late_mean']
            m = o['late_mean']
            rr[lab] = dict(mine=dict(NX=m['NX'], EX=m['EX'], eps=m['eps'], NYL=m['NYL'], NYH=m['NYH'], Nsig=m['Nsig']),
                           engine=dict(NX=e['NX'], EX=e['EX'], eps=e['eps'], NYL=e['NYL'], NYH=e['NYH'], Nsig=e['Nsig']),
                           rel_diff={k: m[k] / e[k] - 1 for k in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig')},
                           gauss=o['gauss_max'], prep_amps=o['prep']['amps'])
        out['release_regression'] = rr
    # 2. independent vector modes
    if fast and 'vector_check' in old:
        out['vector_check'] = old['vector_check']
    else:
        import sk_vec as sv
        bench = TRK / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_P0_it2.json'
        B = sv.background_F(str(bench))
        v = json.load(open(MFC / 'vectors' / 'vectors_itHFa_P0_it2.json'))
        sp = [x for x in v['by_e'] if abs(x['e'] - 1e-6) < 1e-18][0]['spectrum']
        rows = []
        for k in (0.01, 0.0631, 0.15849, 0.39811, 1.0):
            m = sv.mode(k, 1e-6, B)
            tr = min(sp, key=lambda r: abs(math.log(r['k'] / k)))
            rows.append(dict(k=k, n_mine=m.get('n2'), n_track=tr['n_k'], ratio=m.get('n2') / tr['n_k'], wronskian=m['wronskian']))
        out['vector_check'] = dict(history='HFa Pauli P0 sc 2.40e11', e=1e-6, rows=rows)
    # 3. maps, ledger, leak, scans (files written by the other scripts)
    for nm in ('sk_pmap', 'sk_ledger', 'sk_leak', 'sk_scan_HFa_raw', 'sk_scan_P0sc', 'sk_scan_P0H110', 'sk_scan_P0H084', 'sk_scan_HFb_raw'):
        p = RES / (nm + '.json')
        if p.exists():
            d = json.load(open(p))
            if nm.startswith('sk_scan'):
                d = {k: dict(eps=abs(x['late_mean']['eps']), N_YL=x['late_mean']['NYL'], **{kk: x['ledger'][kk] for kk in
                     ('E_Ec', 'R_Y', 'SQ_1e6', 'SQ_3e6', 'SQ_inf', 'e_crit')}) for k, x in d['variants'].items()}
            elif nm == 'sk_ledger':
                d = d['cases']
            elif nm == 'sk_leak':
                d = {k: dict(S_per_Q_noleak=x['S_per_Q_noleak'], **{c: dict(R_leak=r['R_leak'], S_per_Q=r['S_per_Q_with_leak'])
                                                                     for c, r in x['rows'].items()}) for k, x in d['cases'].items()}
            out[nm] = d
    for nm in ('sk_scan_P0H084_imag', 'sk_scan_P0H084_phase'):
        p = RES / (nm + '.json')
        if p.exists():
            out[nm] = json.load(open(p))
    # 4. full-chain runs in scratch
    fc = {}
    for p in sorted((SKE / 'chain').glob('chain_*.json')):
        tag = p.stem[len('chain_'):]
        try:
            fc[tag] = chain_summary(p)
        except Exception as ex:
            fc[tag] = dict(error=repr(ex))
    out['full_chain_runs'] = fc
    for p in sorted((SKE / 'chain').glob('iso_*.json')):
        out.setdefault('iso_runs', {})[p.stem] = json.load(open(p))
    # 5. isocurvature at the strong-holder self-consistent point (c0 = 1.386, 1.414 around 1.40, H_i = 8.40e10 GeV)
    try:
        a = json.load(open(SKE / 'chain' / 'chain_skP0H084_C1.386.json')); b = json.load(open(SKE / 'chain' / 'chain_skP0H084_C1.414.json'))
        dl = math.log(1.414 / 1.386); iso = {}
        for lab in ('e1e-6', 'e3e-6'):
            lo = lambda o: math.log(o[lab]['M'] * abs(o['N_X']) * o['K'] ** -1.5)
            D = (lo(b) - lo(a)) / dl
            H = 8.40e10; x = (H / (2 * math.pi * 7.305e15)) ** 2 * D * D / 2.1e-9
            iso[lab] = dict(D=D, beta_iso_at_sstar_7p3e15=x / (1 + x), s_star_min_GeV=abs(D) * H / (2 * math.pi * math.sqrt(0.038 * 2.1e-9 / 0.962)))
        out['iso_strong_holder'] = iso
    except FileNotFoundError:
        pass
    # 6. beta_iso at the e = 3e-6 unblocked fixed point (track map), H_i of that history
    try:
        b = json.load(open(SP / 'hist' / 'bench_itHFa_track_e3_it2.json'))
        H = b['H_i_GeV']; D = json.load(open(SKE / 'chain' / 'iso_e3U.json'))['e3e-6']['D']
        x = (H / (2 * math.pi * 7.305e15)) ** 2 * D * D / 2.1e-9
        out['iso_unblocked_3e6_fixed_point'] = dict(H_i_GeV=H, D_trackmap=D, beta_iso=x / (1 + x),
                                                    s_star_min_GeV=abs(D) * H / (2 * math.pi * math.sqrt(0.038 * 2.1e-9 / 0.962)))
    except FileNotFoundError:
        pass
    cmb = SP / 'cmb_check.json'
    if cmb.exists():
        out['cmb_capacity'] = json.load(open(cmb))
    json.dump(out, open(out_p, 'w'), indent=1, default=float)
    return out


if __name__ == '__main__':
    o = main('--fast' in sys.argv)
    print(json.dumps(o.get('release_regression', {}), indent=1)[:1500])
    print(json.dumps(o.get('vector_check', {}), indent=1)[:1200])
    print(json.dumps(o.get('iso_strong_holder', {}), indent=1))
    print(json.dumps(o.get('iso_unblocked_3e6_fixed_point', {}), indent=1))
    print(sorted(o['full_chain_runs']))
