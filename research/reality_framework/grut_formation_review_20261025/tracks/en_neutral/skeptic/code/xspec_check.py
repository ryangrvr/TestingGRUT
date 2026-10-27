"""en_neutral_skeptic: own X-fluctuation spectrum (where S = <|delta chi_X|^2> sits) + period-mean S/<|chi_X|^2> two ways
(mean of ratio vs ratio of means) from the N48 trace, + max U_phi/rho0. OWN code. usage: python3 -I xspec_check.py <rundir> <out>"""
import sys, json, struct
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K
rd, out = Path(sys.argv[1]), Path(sys.argv[2])
tr = K.read_trace(rd / 'state.tsv'); rho0 = float(tr['rho'][0])
a0 = float(tr['a'][0]); C = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
N, L = 48, 12.0; kl = np.sqrt(K.lattice_k2_grid(N, L))
res = dict(label='OWN', xspec={})
for path in sorted(rd.glob('state_s*.bin')):
    with open(path, 'rb') as f:
        n_, L_, s_, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f4').reshape(6, n_, n_, n_).astype(np.float64)
    tau = 20 + a0 * s_ + 0.5 * C * s_ * s_
    if min(abs(tau - np.array([60, 70, 80, 100, 150, 200, 400, 800]))) > 1e-6: continue
    pw = sum(abs(np.fft.fftn(a[c]) / N ** 3) ** 2 for c in (0, 1)); pw[0, 0, 0] = 0; tot = pw.sum()
    res['xspec'][f'{tau:.0f}'] = dict(S=float(tot), lt_1p8=float(pw[kl < 1.8].sum() / tot), b1p8_4=float(pw[(kl >= 1.8) & (kl < 4)].sum() / tot),
                                      b4_8=float(pw[(kl >= 4) & (kl < 8)].sum() / tot), gt_8=float(pw[kl >= 8].sum() / tot))
s, tau = tr['eta'], tr['tau']; rb2 = tr['q0'] ** 2 + tr['q1'] ** 2; X2 = tr['X2']; S = X2 - rb2
per = np.flatnonzero((rb2[1:-1] < rb2[:-2]) & (rb2[1:-1] <= rb2[2:])) + 1
rows = []
for a_, b_ in zip(per[:-1], per[1:]):
    sl = slice(a_, b_ + 1)
    rows.append(dict(tau_start=float(tau[a_]), tau_end=float(tau[b_]),
                     mean_of_ratio=float(np.trapezoid(S[sl] / X2[sl], s[sl]) / (s[b_] - s[a_])),
                     ratio_of_means=float(np.trapezoid(S[sl], s[sl]) / np.trapezoid(X2[sl], s[sl])),
                     g2S_mean=float(K.G2 * np.trapezoid(S[sl], s[sl]) / (s[b_] - s[a_]))))
res['period_S_over_X2'] = rows
res['max_Uphi_over_rho0'] = float(np.max(tr['Uphi']) / rho0)
res['max_Umix_over_Ephi_plus_mix_late'] = float(np.mean((tr['Umix'] / (tr['Kphi'] + tr['Gphi'] + tr['Uphi'] + tr['Umix']))[tau >= 700]))
out.write_text(json.dumps(res, indent=1) + '\n')
for k, v in res['xspec'].items(): print(k, {a: round(b, 4) for a, b in v.items()})
for r in rows[:9]: print(r)
print(rows[-1]); print('max Uphi/rho0', res['max_Uphi_over_rho0'], 'late Umix share', res['max_Umix_over_Ephi_plus_mix_late'])
