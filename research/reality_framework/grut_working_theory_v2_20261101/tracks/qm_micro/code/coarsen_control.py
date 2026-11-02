"""qm_micro item 3 numerical controls (Claude, 2 Nov 2026): the S/Q = 3.0, Q0 = 4e11 run with (i) an initial width
sigma_lnQ = 1.5 and (ii) steps 3x finer and 2x finer re-resolution (coarsen_fine.py, generated from coarsen.py by
string replacement of the three control constants)."""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coarsen as c
import coarsen_fine as cf
SCR = os.environ.get('QM_SCR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/qm_micro')


def summ(r):
    f, z = r['freeze'], r['final']
    return dict(T_f=f['T'], f_v_freeze=f['f_v'], Q_freeze=f['Q_charge_weighted'], f_v_final=z['f_v'], Q_final=z['Q_charge_weighted'],
                f_v_b=r['variant_b']['f_v_final'], sigma=z['sigma_lnQ_charge'], steps=r['steps'], refines=r['refines'])


out = {}
H = c.History(os.path.join(SCR, 'hist_3.0.json')); Tsp, _ = c.T_spinodal(3.0)
out['baseline'] = summ(c.run(H, 4e11, T_start=Tsp, verbose=False)); print('baseline', out['baseline'], flush=True)
out['sig0_1.5'] = summ(c.run(H, 4e11, sig0=1.5, T_start=Tsp, verbose=False)); print('sig0', out['sig0_1.5'], flush=True)
H2 = cf.History(os.path.join(SCR, 'hist_3.0.json'))
out['fine'] = summ(cf.run(H2, 4e11, T_start=Tsp, verbose=False, nmax=600)); print('fine', out['fine'], flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
