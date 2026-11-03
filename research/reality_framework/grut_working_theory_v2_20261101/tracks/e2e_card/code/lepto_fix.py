"""e2e_card: corrected rerun of e2e_sm's leptogenesis (e2e_sm/code/sm_lepto.py imported read-only, unchanged).

Finding [C]: e2e_sm reads the 'non-thermal' N_j population X_NT from the history's f_Nj columns. In the e2e_kinetic
export f_N4 equals the thermal-equilibrium N4 energy density to < 2% at t <= 300 and t >= 8000 H_i^-1
(and is x1.2-3.6 above it at t = 1e3-7e3). sm_lepto adds eps*Gamma*X_NT with no inverse-decay counterpart, so the
equilibrium part (which yields no asymmetry) is counted as a source. Its late, weakly washed-out tail (t >~ 1e4) is
what made eta_B/delta = 1.7e-4.

Fix (model choice [A]): X_NT,j -> the excess over equilibrium, max(rho_Nj - rho_Nj^eq(T_thermal), 0)/E_j.
The thermal N_j are already evolved separately by sm_lepto (X_th, relaxing to X_eq by inverse decays).
The patched history is written to scratch and passed to sm_lepto.run unchanged.
"""
import json, math, os, sys, copy
import numpy as np
from scipy.special import kve

TR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(TR, 'e2e_sm', 'code'))
import sm_lepto as L  # noqa: E402  (read-only reuse)

SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_card'
MJ = [250.0, 400.0, 200.0, 120.0]


def rho_eq(M, T):
    if T <= 0:
        return 0.0
    z = min(M / T, 600.0)
    n = (2.0 / (2 * math.pi**2)) * M * M * T * kve(2, z) * math.exp(-z)
    return n * (M * kve(1, z) / kve(2, z) + 3 * T)


def patch(path, H_GeV, tag):
    d = json.load(open(path))
    h = d['history']
    MP = 2.435e18 / H_GeV
    H = np.array(h['H']); rho = 3 * H**2 * MP**2
    T = [x if x is not None else 0.0 for x in h['T_thermal']]
    diag = {}
    for j in range(4):
        f = np.array(h['f_N%d' % (j + 1)], float)
        feq = np.array([rho_eq(MJ[j], T[i]) for i in range(len(T))]) / rho
        h['f_N%d' % (j + 1)] = np.maximum(f - feq, 0.0).tolist()
        diag['N%d' % (j + 1)] = dict(max_ratio_to_eq=float(np.max(np.where(feq > 0, f / np.maximum(feq, 1e-300), 0))))
    out = os.path.join(SCR, 'excess_' + tag + '.json')
    json.dump(d, open(out, 'w'))
    return out, diag


SEE = os.path.join(TR, 'e2e_sm', 'results', 'seesaw_H1.3e11.json')
KIN = os.path.join(TR, 'e2e_kinetic', 'results')
CH = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_chain/kin'
cases = [
    ('kin_HFa_1.3e11', os.path.join(KIN, 'fluid_HFa.json'), 1.3e11),
    ('final_GM_1.3e11', os.path.join(CH, 'fluid_hoHFa_GM.json'), 1.3e11),
    ('final_GM_8.4e10', os.path.join(CH, 'fluid_hoH84_GM.json'), 8.4e10),
    ('final_P0_1.3e11', os.path.join(CH, 'fluid_hoHFa_P0.json'), 1.3e11),
    ('final_P0_8.4e10', os.path.join(CH, 'fluid_hoH84_P0.json'), 8.4e10),
    ('final_GMg03_1.3e11', os.path.join(CH, 'fluid_hoHFa_GM_g03.json'), 1.3e11),
]
if __name__ == '__main__':
    res = {}
    for name, p, Hg in cases:
        if not os.path.exists(p):
            continue
        row = {}
        for mode in ('as_e2e_sm', 'excess'):
            hp = p
            if mode == 'excess':
                hp, diag = patch(p, Hg, name)
                row['eq_ratio_max'] = diag
            for ws in (1.0, 0.5, 2.0):
                r = L.run(hp, SEE, 'S1_mt050', Hg, washout_scale=ws)
                row['%s_w%g' % (mode, ws)] = dict(etaB_per_delta=r['etaB_per_unit_delta'],
                                                 thermal_only=r['etaB_thermal_only_per_unit_delta'],
                                                 delta_needed=r['delta_needed'])
        res[name] = row
        print(name, json.dumps({k: (v['etaB_per_delta'], v['delta_needed']) for k, v in row.items() if k != 'eq_ratio_max'}))
    json.dump(res, open(os.path.join(SCR, 'lepto_fix.json'), 'w'), indent=1)
