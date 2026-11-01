"""fx_chain: collect the m1_kin.py case outputs (scratchpad kin_<case>.json) into results/kin_m1.json.
Trajectories are kept (decimated) only for the main cases; the channel budget for the main R = 4.97 case."""
import json, glob, os, math
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'results')
MAIN = {'m1_gas_relaxed_R497', 'm1_gas_relaxed_R570', 'm1_gas_relaxed_R114', 'm1_inst_relaxed_R497', 'repro_gas_R497'}
rows = {}
for f in sorted(glob.glob(SP + '/wt2/fx_chain/kin_*.json')):
    d = json.load(open(f)); n = d.get('name', os.path.basename(f)[4:-5])
    if 'error' in d:
        rows[n] = d; continue
    keep = {k: d[k] for k in ('R', 'formation', 'x0', 'a0_over_aref', 'initial', 'final', 'stop', 'checks', 'p', 'f_old_convention',
                              'frozen_test_term', 'eta_min', 'YQ', 'V0_GeV4', 'omega_cap', 'e_cap', 'capture_scale',
                              'cubic_applied', 'rtol', 'case', 'seconds') if k in d}
    if n in MAIN:
        keep['trajectory'] = d['trajectories'][::3]
    if 'channel_budget' in d:
        keep['channel_budget'] = d['channel_budget'][::2]
    rows[n] = keep
# R scan (M1, gas-only, relaxed start) and the R at which 2.13356 p = 0.01
sc = sorted((v['R'], v['p']) for k, v in rows.items() if k.startswith('m1_gas_relaxed_Rx') and 'p' in v)
Rc = None
for (r1, p1), (r2, p2) in zip(sc[:-1], sc[1:]):
    if 2.13356 * p1 < 0.01 <= 2.13356 * p2:
        Rc = r1 + (r2 - r1) * (math.log(0.01 / 2.13356) - math.log(p1)) / (math.log(p2) - math.log(p1))
# e at which R_tot(e) = Rc: R_tot = R_Y + (R_vec(1e-6) + R_casc(1e-6)) (1e-6/e)
HO = json.load(open(os.path.join(OUT, 'heat_onset.json')))['R']
e_band = 1e-6 * (HO['vec_e1e6_band'] + HO['casc_e1e6']) / (Rc - HO['Y_L']) if Rc else None
e_broad = 1e-6 * (HO['vec_e1e6_broad'] + HO['casc_e1e6']) / (Rc - HO['Y_L']) if Rc else None
json.dump(dict(cases=rows, R_scan=sc, R_crit_eta1=Rc, e_min_band=e_band, e_min_broad=e_broad,
               note='p = P/n_Q after annihilation (pairs per net charge); f = 2.02855 p; frozen test term 2.13356 p; '
                    'eta_min = 1 - (0.01 - 2.13356 p)/0.0667805. Conditional moment recipes (REVIEW_KINETIC).'),
          open(os.path.join(OUT, 'kin_m1.json'), 'w'), indent=1, default=float)
print('R scan', sc); print('R_crit(eta=1) = %s, e_min band %s broad %s' % (Rc, e_band, e_broad))
for k, v in rows.items():
    if 'p' in v: print('%-34s R=%7.4f %-8s x0=%+.3f p=%.4e test=%.4e eta>=%.4f stop=%s' % (k, v['R'], v['formation'], v['x0'], v['p'], v['frozen_test_term'], v['eta_min'], v['stop']))
