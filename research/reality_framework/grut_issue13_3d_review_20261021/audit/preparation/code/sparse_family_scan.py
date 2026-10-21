"""Distribution of the sourced zero mode over Codex's sparse 64-mode family (seeds 1..NSEED) at L=76.8 and 115.2.

Usage: python3 -I sparse_family_scan.py <codex_raw_dir> <out_json> [nseed]
Replicates run.py modes() exactly (checked bitwise against Codex's s101/s202/s303 tables), then evaluates the
quasi-linear S0/sigma_f at total rms 0.5 (kappa=2). Class thresholds are EMPIRICAL from Codex runs + this audit's
interventions at kappa=2: S0/sigma_f >= 0.34 -> minority removed by a=2.2 (noW 0.34, noU 0.59, s101 0.91 all removed);
<= 0.13 -> network retained (s303 0.13, s202 0.04); 0.13-0.34 untested.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import quasilinear as ql  # noqa: E402

TG = np.round(np.arange(0, 48.0001, 0.05), 10)


def codex_modes(L, seed):
    rng = np.random.default_rng(seed)
    kvec = []
    while len(kvec) < 64:
        k = np.rint(rng.normal(0, L / (2 * np.pi), 3)).astype(int)
        if not k.any() or np.linalg.norm(k) * 2 * np.pi / L > 2.5:
            continue
        first = np.flatnonzero(k)[0]
        if k[first] < 0:
            k = -k
        row = tuple(int(t) for t in k)
        if row not in kvec:
            kvec.append(row)
    co = rng.normal(size=(64, 6))
    for f in range(3):
        norm = np.sqrt(np.sum(co[:, 2 * f:2 * f + 2] ** 2) / 2)
        co[:, 2 * f:2 * f + 2] /= norm
    return np.array(kvec), co


def ratio(n, c, N, dx, kappa=2.0):
    F = kappa ** (1 / 3)
    k2 = ql.khat2_int(n, N, dx)
    cu, su, cw, sw, cS, sS = c.T
    pS = (cS ** 2 + sS ** 2) / 2
    CuS = (cu * cS + su * sS) / 2
    CwS = (cw * cS + sw * sS) / 2
    g, _ = ql.s_growth(k2, kappa, TG)
    Au, Aw = ql.x_response(k2, kappa, TG)
    J = -np.sum(ql.EPS_X * F * ql.EPS_S * g * (Au * CuS + Aw * CwS), axis=1)
    chi0 = ql.zero_mode(TG, J, kappa)
    a = ql.a_of(TG)
    sig = ql.EPS_S * np.sqrt(np.sum(g ** 2 * pS, axis=1)) / a
    S0 = chi0 / a
    tot = np.sqrt(sig ** 2 + S0 ** 2)
    j = int(np.argmax(tot >= 0.5))
    kmin = float(np.min(np.linalg.norm(n, axis=1)) * 2 * np.pi / (N * dx))
    return float(S0[j] / sig[j]), float(a[j]), kmin, float(CuS.sum()), float(CwS.sum())


def main():
    raw, outp = Path(sys.argv[1]), Path(sys.argv[2])
    nseed = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    check = {}
    for s in (101, 202, 303):
        n, c = codex_modes(76.8, s)
        ref = np.loadtxt(raw / f'modes_L76.8_s{s}.txt')
        check[s] = float(np.max(np.abs(np.column_stack((n, c)) - ref)))
    n, c = codex_modes(115.2, 101)
    check['L115.2_s101'] = float(np.max(np.abs(np.column_stack((n, c)) - np.loadtxt(raw / 'modes_L115.2_s101.txt'))))
    print('replication max abs diff', check, flush=True)
    res = {'replication_check_max_abs_diff': check, 'L76.8': [], 'L115.2': []}
    for L, N in ((76.8, 192), (115.2, 288)):
        for s in range(1, nseed + 1):
            n, c = codex_modes(L, s)
            r, a, kmin, cu, cw = ratio(n, c, N, 0.4)
            res[f'L{L:g}'].append({'seed': s, 'S0_over_sigma_f': r, 'a_rms0.5': a, 'k_min': kmin, 'C_uS': cu, 'C_wS': cw})
        x = np.array([e['S0_over_sigma_f'] for e in res[f'L{L:g}']])
        km = np.array([e['k_min'] for e in res[f'L{L:g}']])
        res[f'summary_L{L:g}'] = {'n': len(x), 'frac_abs_ge_0.34': float(np.mean(np.abs(x) >= 0.34)),
                                  'frac_abs_le_0.13': float(np.mean(np.abs(x) <= 0.13)),
                                  'median_abs': float(np.median(np.abs(x))), 'mean_signed': float(np.mean(x)),
                                  'corr_abs_ratio_vs_kmin': float(np.corrcoef(np.abs(x), km)[0, 1]),
                                  'sd_C_uS': float(np.std([e['C_uS'] for e in res[f'L{L:g}']])),
                                  'sd_C_wS': float(np.std([e['C_wS'] for e in res[f'L{L:g}']]))}
        print(L, res[f'summary_L{L:g}'], flush=True)
    outp.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
