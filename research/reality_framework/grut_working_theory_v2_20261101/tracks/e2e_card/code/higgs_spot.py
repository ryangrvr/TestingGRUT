"""e2e_card: spot checks of e2e_higgs (own short code; reads e2e_higgs's cached flat lambda_eff tables read-only).

1. Hawking-Moss estimate of xi_min during inflation at constant H, own criterion:
   V = 6 xi H^2 h^2 + lambda_eff(h) h^4 / 4 (xi constant, R = 12 H^2), barrier V_b;
   safe if 8 pi^2 V_b / (3 H^4) > 3 N0 (no HM jump in any of e^{3 N0} Hubble patches) [E].
2. M_t,crit (gravitational only) at the final M2 scales, quadratic in ln H_i through e2e_higgs's 5e10/1.3e11/2.4e11
   values [E, interpolation].
3. N-gas Higgs medium mass at t = 0.2-3 on the e2e_kinetic history, own estimate m^2 = sum_j y_j^2 rho_Nj / E_j^2
   with e2e_sm S1 y_j^2 and E_j = max(3/4 A_j m_parent, M_j, m_parent/2).
"""
import json, math, os
import numpy as np

TR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SMC = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_higgs/smcache'
out = {}


def lam_tab(Mt, Mh=125.2, als=0.118):
    d = json.load(open(os.path.join(SMC, 'sm_%.4f_%.3f_%.4f.json' % (Mt, Mh, als))))
    lh, le = np.log(d['hs']), np.array(d['le'])
    return lambda h: np.interp(np.log(h), lh, le)


def xi_min(H, Mt, N0=57.2):
    lam = lam_tab(Mt)
    hs = np.exp(np.linspace(math.log(1e8), math.log(1e18), 4000))

    def ok(xi):
        V = 6 * xi * H**2 * hs**2 + lam(hs) * hs**4 / 4
        i = int(np.argmax(V))
        if i == len(hs) - 1:
            return True
        return 8 * math.pi**2 * V[i] / (3 * H**4) > 3 * N0
    lo, hi = 1e-4, 1.0
    for _ in range(50):
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return hi


out['xi_min_HM_own'] = {}
for Hg in (1.3e11, 1.123e11, 9.89e10):
    out['xi_min_HM_own']['%.3g' % Hg] = {str(Mt): round(xi_min(Hg, Mt), 4) for Mt in (172.57, 172.86)}
out['xi_min_e2e_higgs_const_tree_1.3e11'] = {'172.57': 0.0487, '172.86': 0.0531}

s = json.load(open(os.path.join(TR, 'e2e_higgs', 'results', 'summary.json')))
pg = s['post_grav']
lnH = np.log([5e10, 1.3e11, 2.4e11])
out['Mt_crit_grav_at_M2'] = {}
for hist in ('KIN', 'LAT_P0', 'LAT_GM', 'LATF_P0'):
    y = [pg['%s@%s' % (hist, k)]['Mt_crit'] for k in ('5e+10', '1.3e+11', '2.4e+11')]
    c = np.polyfit(lnH, y, 2)
    out['Mt_crit_grav_at_M2'][hist] = {'M2X_1.123e11': round(float(np.polyval(c, math.log(1.123e11))), 3),
                                       'M2Q_9.89e10': round(float(np.polyval(c, math.log(9.89e10))), 3)}
for als in ('0.1171', '0.1189'):
    y = [pg['KIN@%s' % k]['alpha_s'][als]['Mt_crit'] for k in ('5e+10', '1.3e+11', '2.4e+11')]
    c = np.polyfit(lnH, y, 2)
    out['Mt_crit_grav_at_M2']['KIN_alphas_' + als] = {'M2X_1.123e11': round(float(np.polyval(c, math.log(1.123e11))), 3),
                                                       'M2Q_9.89e10': round(float(np.polyval(c, math.log(9.89e10))), 3)}

# N-gas medium mass on the kinetic history
d = json.load(open(os.path.join(TR, 'e2e_kinetic', 'results', 'fluid_HFa.json')))
h = d['history']
card = json.load(open(os.path.join(TR, 'e2e_sm', 'results', 'sm_card.json')))
y2 = np.array(card['singlets']['y2']['value'])
M = np.array([250.0, 400.0, 200.0, 120.0]); mpar = np.array([599.99, 249.995, 249.995, 249.995])
MP = 2.435e18 / 1.3e11
t = np.array(h['t']); H = np.array(h['H'])
med = {}
for tt in (0.2, 0.3, 0.5, 1.0, 3.0):
    i = int(np.searchsorted(t, tt))
    rho = 3 * H[i]**2 * MP**2
    m2 = 0.0
    for j in range(4):
        A = h['A_N%d' % (j + 1)][i] or 0.0
        E = max(0.75 * A * mpar[j], M[j], mpar[j] / 2)
        m2 += y2[j] * h['f_N%d' % (j + 1)][i] * rho / E**2
    med['%g' % tt] = dict(mH2_Ngas_over_Hi2=float(m2), kinetic_medium_column=h['m_H2_medium'][i])
out['N_gas_medium_mass'] = med
json.dump(out, open(os.path.join(TR, 'e2e_card', 'results', 'higgs_spot.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
