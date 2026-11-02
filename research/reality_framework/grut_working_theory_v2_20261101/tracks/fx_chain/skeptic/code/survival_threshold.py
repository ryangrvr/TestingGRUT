"""fx_chain_skeptic: total heat R (and e) below which objects of a given omega survive the whole cooling, from
s/n_Q(R) of the pre-onset Hartree plasma (hartree_indep.solve, phi in equilibrium) and sigma_min of survival.py.
e from R_tot(e) = R_Y + R_vec(1e-6)(1e-6/e) + R_casc(1e-6)(1e-6/e)^(1/2) (kin_source: cascade vectors have a fixed
number N_YH, so R_casc ~ m_A a_decay ~ e^(-1/2)).  [C]"""
import json, math, sys, os
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hartree_indep as HI

surv = json.load(open(sys.argv[1]))
Rg = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.9665]
sn = []
g = None
for R in Rg:
    st, g = HI.solve(R, True, 0.0, g)
    sn.append(st['S'] / HI.NX)
    print('R=%.3f s/n_Q=%.5f T=%.4f' % (R, sn[-1], st['T']), flush=True)
lR = np.log1p(Rg); ls = np.log(sn)
YL, VB, VBR, CAS = 0.8055701153147813, 4.053801019129129, 4.790544357727859, 0.10716449056156367
out = dict(R_grid=Rg, s_per_net=sn, thresholds={})
for name, c in surv['cases'].items():
    smin = c['sigma_min']
    if smin <= sn[0]:
        out['thresholds'][name] = dict(sigma_min=smin, R_crit=None, note='does not survive even at R = 0'); continue
    Rc = math.expm1(float(np.interp(math.log(smin), ls, lR)))
    eb = 1e-6 / brentq(lambda x: YL + VB * x + CAS * math.sqrt(x) - Rc, 1e-6, 100) if Rc > YL else None
    er = 1e-6 / brentq(lambda x: YL + VBR * x + CAS * math.sqrt(x) - Rc, 1e-6, 100) if Rc > YL else None
    out['thresholds'][name] = dict(omega=c['omega'], sigma_min=smin, T_at_min=c['T_at_min'], R_crit=Rc, e_min_band=eb, e_min_broad=er)
    print(name, out['thresholds'][name], flush=True)
json.dump(out, open(sys.argv[2], 'w'), indent=1)
