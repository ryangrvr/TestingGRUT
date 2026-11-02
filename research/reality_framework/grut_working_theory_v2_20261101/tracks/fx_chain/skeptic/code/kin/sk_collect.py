"""fx_chain_skeptic: collect the skeptic network runs (scratchpad kin_<case>.json) into ../../results/sk_kin.json,
with a decimated condensed-fraction / pair trajectory per run (object survival through the cooling)."""
import json, glob, os
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/fx_chain_skeptic'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'results', 'sk_kin.json')
rows = {}
for f in sorted(glob.glob(SP + '/kin_*.json')):
    d = json.load(open(f)); n = d.get('name', os.path.basename(f)[4:-5])
    if 'error' in d:
        rows[n] = d; continue
    tr = d['trajectories']
    c = [r['condensed_fraction'] for r in tr]
    zero = [r['T_over_m'] for r in tr if r['condensed_fraction'] <= 1e-12]
    i = min(range(len(c)), key=lambda j: c[j])
    rows[n] = dict(R=d['R'], formation=d['formation'], x0=d['x0'], a0_over_aon=None, start_T_over_m=d['initial']['t'],
                   start_cond=d['initial']['cond_fraction'], rho_over_mX4=d['initial']['rho_over_mX4'], V0_GeV4=d['V0_GeV4'],
                   YQ=d['YQ'], capture_scale=d['capture_scale'], cubic=d['cubic_applied'], stop=d['stop'],
                   T_stop=d['final']['T_over_m'], p=d['p'], test_term_2p13356=d['frozen_test_term'], eta_min_Q2e5_calibration=d['eta_min'],
                   cond_min=c[i], T_at_cond_min=tr[i]['T_over_m'], cond_end=c[-1],
                   cond_zero_window_T=[min(zero), max(zero)] if zero else None,
                   checks=d['checks'],
                   traj=[dict(T=r['T_over_m'], cond=r['condensed_fraction'], p=r['pairs_per_net'], Xbar=r['Xbar_gas_per_net'])
                         for r in tr[::max(1, len(tr) // 14)]] + [dict(T=tr[-1]['T_over_m'], cond=c[-1], p=tr[-1]['pairs_per_net'])])
json.dump(rows, open(OUT, 'w'), indent=1, default=float)
for n, v in rows.items():
    if 'p' in v:
        print('%-36s R=%.3f %-8s p=%.4e cond_min=%.4f@T=%.3f cond_end=%.4f zero=%s stop=%s' % (n, v['R'], v['formation'], v['p'], v['cond_min'], v['T_at_cond_min'], v['cond_end'], v['cond_zero_window_T'], v['stop']))
