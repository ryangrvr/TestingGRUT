"""e2e_card: independent spot checks of e2e_sm's neutrino sector (own code; inputs from e2e_sm/results/sm_card.json).

1. light masses and sum (NuFIT 6.0 NO splittings as quoted by e2e_sm), CP-conserving sum rule;
2. Casas-Ibarra: m~_j from R0, light-mass reconstruction from R(delta) = R0 exp(delta A), R R^T = 1;
3. rest widths Gamma_j = (y^2)_jj M_j / 8 pi, y^2 = m~ M / v^2;
4. CP asymmetries eps_j(delta) from the Covi-Roulet-Vissani unflavoured formula (own implementation),
   eps_j = -(1/8 pi) sum_k Im[(h^dag h)_jk^2]/(h^dag h)_jj * f(M_k^2/M_j^2),
   f(x) = sqrt(x) [1 - (1 + x) ln((1 + x)/x) + 1/(1 - x)]  (vertex + self-energy, Majorana);
5. radial-mode input theta-dot/f of sigma from the mf_hybrid bench history (max |sigma-dot|/f1).
"""
import json, math, os
import numpy as np
from scipy.linalg import expm

TR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
card = json.load(open(os.path.join(TR, 'e2e_sm', 'results', 'sm_card.json')))
VEW = 174.1
Hi = 1.3e11
out = {}
dm21, dm31 = 7.49e-5, 2.513e-3
m = np.array([0.0, math.sqrt(dm21), math.sqrt(dm31)])
out['masses_eV'] = m.tolist(); out['sum_eV'] = float(m.sum())
sing = card['singlets']
Mh = np.array(sing['M_over_Hi']['value']); MG = Mh * Hi
eta = np.array(sing['CP_parities_at_delta0']['value']); mt = np.array(sing['mtilde_eV']['value'])
out['sum_rule_sum_eta_mt'] = float(np.dot(eta, mt))
R0 = np.array(sing['CI_R0_re']['value']) + 1j * np.array(sing['CI_R0_im']['value'])
A_real = np.array(sing['CI_A_odd']['value'])
# e2e_sm builds A = where(eta_j eta_k > 0, 1j, 1) * A_ODD (sm_seesaw.family); the card stores only A_ODD
A = np.where(np.outer(eta, eta) > 0, 1j, 1.0) * A_real
out['R0R0T_minus_1_max'] = float(np.max(np.abs(R0 @ R0.T - np.eye(3))))
out['mt_from_R0'] = (np.abs(R0)**2).T.dot(m).tolist()
y2 = mt * MG / VEW**2 * 1e-9
out['Gamma_over_Hi'] = (y2 * Mh / (8 * math.pi)).tolist()


def f_crv(x):
    return math.sqrt(x) * (1 - (1 + x) * math.log((1 + x) / x) + 1 / (1 - x))


def eps_of(delta):
    R = R0 @ expm(delta * A)
    # Yukawa (3 x 4) up to the PMNS rotation (drops out unflavoured): Y v = i sqrt(m) R sqrt(M)  (eV-GeV mixed: use ratio)
    Yv = 1j * np.diag(np.sqrt(m * 1e-9)) @ R @ np.diag(np.sqrt(MG))
    h = Yv / VEW
    hh = h.conj().T @ h
    eps = []
    for j in range(4):
        s = 0.0
        for k in range(4):
            if k == j:
                continue
            s += (hh[j, k]**2).imag * f_crv(MG[k]**2 / MG[j]**2)
        eps.append(-s / (8 * math.pi * hh[j, j].real))
    mlight = np.linalg.svd(Yv @ np.diag(1 / MG) @ Yv.T, compute_uv=False)  # GeV -> light masses
    return np.array(eps), np.sort(mlight * 1e9), (np.abs(R)**2).T.dot(m), float(np.max(np.abs(R @ R.T - np.eye(3))))


for d in (1e-4, 1e-3, 1e-2):
    e, ml, mtd, err = eps_of(d)
    out['delta_%g' % d] = dict(direction='e2e_sm complex A', eps_over_delta=(e / d).tolist(), light_masses_eV=ml.tolist(), mt=mtd.tolist(), RRT_err=err)
Asave = A; A = A_real
e, ml, mtd, err = eps_of(1e-3)
out['real_A_direction_delta_1e-3'] = dict(eps_over_delta=(e / 1e-3).tolist(), mt=mtd.tolist())
A = Asave
out['e2e_sm_eps_per_delta'] = card['CP_phase']['eps_per_unit_delta']['value']
# Davidson-Ibarra bound for the lightest singlet, for orientation (hierarchical-limit formula [Q, hep-ph/0202239])
out['DI_bound_N4'] = float(3 * MG[3] * 0.0501e-9 / (16 * math.pi * VEW**2))

# sigma-dot / f1 from the bench
b = json.load(open(os.path.join(TR, 'mf_hybrid', 'results', 'bench_HFa.json')))
hist = b['history']; cols = b['history_columns']
if isinstance(hist, dict):
    t = np.array(hist['t']); sg = np.array([x if x is not None else np.nan for x in hist['sigma_MP']])
else:
    names = list(cols.keys()) if isinstance(cols, dict) else cols
    arr = np.array([[x if x is not None else np.nan for x in r] for r in hist], float)
    t = arr[:, names.index('t')]; sg = arr[:, names.index('sigma_MP')]; ch = arr[:, names.index('chi_MP')]
MPH = 2.435e18 / Hi
ok = np.isfinite(sg) & (t > 0) & (t < 0.6)
sd = np.gradient(sg[ok] * MPH, t[ok])
out['sigma_dot_over_f1_max_t>0'] = float(np.nanmax(np.abs(sd)) / 1200.0)
out['sigma_dot_over_f1_max_t>0.19'] = float(np.nanmax(np.abs(sd[t[ok] > 0.19])) / 1200.0)
out['bench_dt_min'] = float(np.min(np.diff(t[ok])))
json.dump(out, open(os.path.join(TR, 'e2e_card', 'results', 'sm_spot.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
