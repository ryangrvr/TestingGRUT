#!/usr/bin/env python3
"""gp_phase: collect the result files into results/summary.json and print the key numbers."""
import json
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1] / 'results'


def load(name):
    p = R / name
    return json.loads(p.read_text()) if p.exists() else None


def main():
    out = {}
    rs = load('real_subspace.json')
    if rs:
        out['a_real_subspace'] = dict(exact_zero_max_abs_charge=max(z['max_abs_a3n'] for z in rs['exact_zero_runs']),
                                      codex_pi24_max_rel_diff=max(c['max_rel_diff'] for c in rs['codex_pi24_comparison']),
                                      cp_odd_growth_rate=rs.get('cp_odd_envelope_growth_rate_per_unit_conformal_time'),
                                      eps_at_tau300_from_seed_1e15=rs['cp_odd_growth']['1e-15']['epsX'][-1])
    ops = load('operators.json')
    if ops:
        out['b_phase_count'] = {k: {r['dim_max']: r['physical_phases'] for r in v} for k, v in ops.items() if k.startswith('Z2')}
    ps = load('phase_scan.json')
    if ps:
        cfgs = [c for c in ps['configs'] if 'eps_mean' in c]
        allabs = max(c['eps_abs_max'] for c in cfgs)
        nruns = sum(c['n_ok'] for c in cfgs)
        nmeet = sum(round(c['frac_meeting_requirement'] * c['n_ok']) for c in cfgs)
        sig = [(abs(c['eps_mean']) / c['eps_sem'], c['op'], c['r'], c['delta'], c['eps_mean'], c['eps_sem']) for c in cfgs if c['op']]
        sig.sort(reverse=True)
        base = next(c for c in cfgs if c['op'] == 0)
        # CP-odd combination per (op, r): mean(+pi/2) - mean(-pi/2)
        cpodd = []
        for c in cfgs:
            if c['op'] and abs(c['delta'] - np.pi / 2) < 1e-9:
                m = next((d for d in cfgs if d['op'] == c['op'] and d['r'] == c['r'] and abs(d['delta'] + np.pi / 2) < 1e-9), None)
                if m:
                    cpodd.append(dict(op=c['op'], r=c['r'], half_difference=(c['eps_mean'] - m['eps_mean']) / 2,
                                      sem=0.5 * np.hypot(c['eps_sem'], m['eps_sem'])))
        out['c_full_model'] = dict(n_configs=len(cfgs), n_runs=nruns, max_abs_eps_any_run=allabs, runs_meeting_requirement=nmeet,
                                   baseline=dict(mean=base['eps_mean'], sem=base['eps_sem'], std=base['eps_std'], abs_mean=base['eps_abs_mean'], abs_max=base['eps_abs_max']),
                                   std_range=[min(c['eps_std'] for c in cfgs), max(c['eps_std'] for c in cfgs)],
                                   abs_mean_range=[min(c['eps_abs_mean'] for c in cfgs), max(c['eps_abs_mean'] for c in cfgs)],
                                   mean_range=[min(c['eps_mean'] for c in cfgs), max(c['eps_mean'] for c in cfgs)],
                                   largest_significance=sig[:5], n_configs_beyond_3sem=sum(s[0] > 3 for s in sig), n_cp_configs=len(sig),
                                   cp_odd_half_differences=cpodd,
                                   max_abs_cp_odd=max(abs(x['half_difference']) for x in cpodd) if cpodd else None,
                                   no_local_minimum=[(c['op'], c['r'], c['delta']) for c in ps['configs'] if 'status' in c],
                                   max_gauss=max(c['max_gauss'] for c in cfgs), max_ledger=max(c['max_ledger'] for c in cfgs))
    ct = load('controls.json')
    if ct:
        out['c_controls'] = {k: {kk: dict(mean=vv['eps_mean'], sem=vv['eps_sem'], std=vv['eps_std'], absmax=vv['eps_abs_max']) for kk, vv in v.items()}
                             for k, v in ct.items() if isinstance(v, dict) and k != 'eps_requirement'}
        out['c_DRT'] = [dict(r=d['r'], delta=d['delta'], mean=d['eps_mean'], sem=d['eps_sem'], std=d['eps_std'], A0=d['A_initial'], det=d['deterministic_eps'])
                        for d in ct['DRT_hubble_A_term']]
    sf = load('single_field.json')
    if sf:
        out['c_single_field_envelope'] = dict(fine_best=sf['fine_best'], max_abs_eps_by_eta=sf['max_abs_eps_by_eta'],
                                              controls={k: v.get('eps_final') for k, v in sf['controls'].items()},
                                              n_coarse_runaway=sf['n_coarse_runaway'])
    (R / 'summary.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float)[:6000])


if __name__ == '__main__':
    main()
