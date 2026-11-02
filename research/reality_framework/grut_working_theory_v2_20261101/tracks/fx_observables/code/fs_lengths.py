"""Comoving rms free-streaming length after momentum decoupling, lambda = sqrt(3) sigma_p int_{a_m}^{1} dtau/a, for the
gas scenarios of vapour_transfer.json (own; background of vapour_transfer.py).  Writes results/free_streaming_lengths.json."""
import sys, json, math, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from scipy.integrate import quad
import vapour_transfer as vt
RES = os.path.join(HERE, '..', 'results')


def S(a0, a1=1.0):
    return quad(lambda la: 1 / (math.exp(la)**2 * vt.H_Mpc(math.exp(la))), math.log(a0), math.log(a1), limit=400)[0]


d = json.load(open(os.path.join(RES, 'vapour_transfer.json')))
A = d['A']
out = {'A': dict(lambda_fs_rms_pc=math.sqrt(3) * A['sigma_p_today'] * S(A['a_m']) * 1e6, r_s_pc=A['r_s_Mpc'] * 1e6)}
for r in d['B']:
    if r['S_over_Q'] == 1.235 and 'a_m' in r:
        out['Q%.0e_s%g' % (r['Q'], r['sticking'])] = dict(lambda_fs_rms_kpc=math.sqrt(3) * r['sigma_p_today'] * S(r['a_m']) * 1e3,
                                                         r_s_kpc=r['r_s_Mpc'] * 1e3, zeta=r['sticking'] * r['Q']**(-1 / 3))
json.dump(out, open(os.path.join(RES, 'free_streaming_lengths.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
