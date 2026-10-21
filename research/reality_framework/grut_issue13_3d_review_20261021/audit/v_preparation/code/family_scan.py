"""Independent family scan of Codex's sparse generator (run.py modes()) over seeds 1..150 at L=76.8 and 115.2.
Computes C_uS, C_wS and the leapfrog quasi-linear S0/sigma_f at total rms 0.5; checks same-seed reuse across L.
Usage: python3 -I family_scan.py <codex_raw_dir> <out_json> [nseed]"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import leapfrog_ql as lq
raw, outp = Path(sys.argv[1]), Path(sys.argv[2]); NS = int(sys.argv[3]) if len(sys.argv) > 3 else 150
def gen(L, seed):  # transcribed from Codex run.py modes()
    rng = np.random.default_rng(seed); kv = []
    while len(kv) < 64:
        k = np.rint(rng.normal(0, L / (2 * np.pi), 3)).astype(int)
        if not k.any() or np.linalg.norm(k) * 2 * np.pi / L > 2.5: continue
        f = np.flatnonzero(k)[0]
        if k[f] < 0: k = -k
        r = tuple(int(x) for x in k)
        if r not in kv: kv.append(r)
    co = rng.normal(size=(64, 6))
    for f in range(3):
        co[:, 2 * f:2 * f + 2] /= np.sqrt(np.sum(co[:, 2 * f:2 * f + 2] ** 2) / 2)
    return np.array(kv), co
chk = {s: float(np.max(np.abs(np.column_stack(gen(76.8, s)) - np.loadtxt(raw / f'modes_L76.8_s{s}.txt')))) for s in (101, 202, 303)}
chk['L115.2_s101'] = float(np.max(np.abs(np.column_stack(gen(115.2, 101)) - np.loadtxt(raw / 'modes_L115.2_s101.txt'))))
res = {'replication_max_abs_diff': chk}
per = {}
for L, N in ((76.8, 192), (115.2, 288)):
    rows = []
    for s in range(1, NS + 1):
        n, c = gen(L, s)
        o, _ = lq.evolve(n, c, N, 0.4, 2.0, T=30.0)
        tot = np.sqrt(o['sigma_f'] ** 2 + o['S0'] ** 2); j = int(np.argmax(tot >= 0.5))
        rows.append({'seed': s, 'C_uS': float(np.sum(c[:, 0] * c[:, 4] + c[:, 1] * c[:, 5]) / 2), 'C_wS': float(np.sum(c[:, 2] * c[:, 4] + c[:, 3] * c[:, 5]) / 2),
                     'ratio': float(o['S0'][j] / o['sigma_f'][j]), 'k_min': float(np.min(np.linalg.norm(n, axis=1)) * 2 * np.pi / L), 'coef': c})
    per[L] = rows
    r = np.array([x['ratio'] for x in rows]); cu = np.array([x['C_uS'] for x in rows]); cw = np.array([x['C_wS'] for x in rows])
    res[f'L{L:g}'] = {'sd_C_uS': float(cu.std()), 'sd_C_wS': float(cw.std()), 'frac_abs_ratio_ge_0.34': float(np.mean(np.abs(r) >= 0.34)),
                      'frac_abs_ratio_ge_0.25': float(np.mean(np.abs(r) >= 0.25)), 'frac_abs_ratio_le_0.13': float(np.mean(np.abs(r) <= 0.13)),
                      'mean_signed_ratio': float(r.mean()), 'median_abs_ratio': float(np.median(np.abs(r)))}
    print(L, res[f'L{L:g}'], flush=True)
same = [bool(np.array_equal(a['coef'], b['coef'])) for a, b in zip(per[76.8], per[115.2])]
ra = np.array([x['ratio'] for x in per[76.8]]); rb = np.array([x['ratio'] for x in per[115.2]])
res['same_coefficients_across_L_fraction'] = float(np.mean(same))
res['class_agreement_across_L_(>=0.34)'] = float(np.mean((np.abs(ra) >= 0.34) == (np.abs(rb) >= 0.34)))
res['corr_ratio_L76.8_vs_L115.2'] = float(np.corrcoef(ra, rb)[0, 1])
res['binomial_sd_of_fraction'] = float(np.sqrt(0.14 * 0.86 / NS))
res['expected_sd_C_for_128_gaussian_pairs'] = float(1 / np.sqrt(128))
for x in per[76.8] + per[115.2]: x.pop('coef')
res['rows'] = {str(L): v for L, v in per.items()}
print({k: v for k, v in res.items() if k != 'rows'})
outp.write_text(json.dumps(res, indent=1))
