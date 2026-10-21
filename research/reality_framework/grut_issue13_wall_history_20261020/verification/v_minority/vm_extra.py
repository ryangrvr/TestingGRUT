"""v_minority: extra spot checks of closed forms (verifier's own code; uses vm_cosmo, vm_fp only)."""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vm_cosmo as C
import vm_fp as F
OUT = os.path.dirname(os.path.abspath(__file__))
target = 0.12 * C.rho_c_over_h2 / (2 * math.pi ** 2 / 45 * 3.909 * C.T0 ** 3)
R = {}
# count closed forms: N(T_re < T_X) ≈ (Lambda1/2) e^{3 n'_X}/3 < 1  <=> Lambda1 < 6 e^{-3 n'_X}
for lab, TX in [('Tdom_iso', C.T_dom_iso()), ('1MeV', 1e-3), ('Trec', 0.256e-9)]:
    npx = C.nprime(TX)
    Lmax = 6 * math.exp(-3 * npx)
    Hx = brentq(lambda h: math.log(F.lambda1_mfpt(h)[0] / Lmax), 120, 1500, xtol=0.01)
    R['count_closed_' + lab] = dict(nprime=npx, Lambda1_max=Lmax, H=Hx)
# relic at 241 GeV relaxed N*=50 (Omega h^2 = 0.12 * rho_s/bound)
G, Rr, rel, fr, ng = F.rates_for(241.0)
r = F.relic_rho_over_s(241.0, 50, rel)[0]
R['relic_241_relaxed_N50_log10_Omegah2'] = math.log10(0.12 * r / target)
R['Lambda1_241'] = Rr['Lambda1']
# closed-form relic threshold Lambda1 <= 2 * bound / (C_E sigma H/s at T_dom)
CE = 3 * 7.7866 / (4 * math.pi)
Td = C.T_dom_iso()
L = 2 * target / (CE * C.sigma * C.H_of_T(Td) / C.s_of_T(Td))
R['closed_relic_Lambda1_max'] = L
R['closed_relic_H'] = brentq(lambda h: math.log(F.lambda1_mfpt(h)[0] / L), 250, 600, xtol=0.01)
with open(os.path.join(OUT, 'results_extra.json'), 'w') as f:
    json.dump(R, f, indent=1)
print(json.dumps(R, indent=1))
