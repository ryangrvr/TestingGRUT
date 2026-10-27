"""Generate the tau=20 fields of further cases with Codex's prepare.initialize (reproduction) and apply our
whitened-covariance test (own code). Usage: python3 -B whiten_cases.py <packet_dir> <out.json>"""
import sys, json, os
here = os.path.dirname(os.path.abspath(__file__))
pk, out = sys.argv[1:3]
sys.path.insert(0, os.path.join(pk, 'code'))
from prepare import initialize  # Codex code, run in place
sys.path.insert(0, here)
import numpy as np
from whiten import whiten_test
res = {}
for tag, N, L, seed, width in [('IR12_s1_w1', 32, 12., 1, 1), ('IR12_s2_w1', 32, 12., 2, 1), ('IR12_s1_w3', 32, 12., 1, 3),
                               ('IR12_fine_s1_w1', 48, 12., 1, 1), ('IR24_s1_w1', 64, 24., 1, 1), ('IR24_s2_w1', 64, 24., 2, 1)]:
    z, params = initialize(N, L, seed, width)
    d = np.load(os.path.join(pk, f'PREPARATION_L{L:g}.npz'))
    r = whiten_test(z, d, f'factor_width{width}', L, N)
    r['q_RMS_over_mean'] = float(params['q_RMS_over_mean']); r['p_RMS_over_mean'] = float(params['p_RMS_over_mean'])
    r['nQ'] = float(params['nQ']); r['rho'] = float(params['rho'])
    pub = json.load(open(os.path.join(pk, f'RUN_{tag}.json')))['initial'] if os.path.exists(os.path.join(pk, f'RUN_{tag}.json')) else {}
    r['published_q_RMS_over_mean'] = pub.get('q_RMS_over_mean'); r['published_nQ'] = pub.get('nQ')
    res[tag] = r
    print(tag, 'pooled', r['pooled_mean_abs2'], 'p', r['pooled_chi2_two_sided_p'], 'diag std(>=50)', r['diag_std_shells_ge50'], r['expected_std_shells_ge50'],
          'qRMS', r['q_RMS_over_mean'], r['published_q_RMS_over_mean'], 'nQ', r['nQ'], r['published_nQ'], flush=True)
json.dump(res, open(out, 'w'), indent=1)
