#!/usr/bin/env python3
"""REPRODUCTION comparison (own code): reproduced JSON outputs vs the published packet JSON.

Every numeric leaf is compared (relative difference |a-b|/max(|b|,1e-300)); strings (hashes) are
compared for equality.  Timing fields are ignored.  Also recomputes the headline numbers directly
from the regenerated raw arrays, independently of Codex's summary scripts.
Usage: python3 compare_reproduction.py SCRATCH_REVIEW_DIR OUT.json
  SCRATCH_REVIEW_DIR holds tree/ (reproduced) and ff49e5d/, 9636950/, bb211d1/ (pristine git-archive copies).
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np

IGNORE = {'seconds', 'elapsed_seconds', 'executable_stdout'}
RF = 'research/reality_framework'


def walk(a, b, path, acc):
    if isinstance(b, dict):
        if not isinstance(a, dict):
            acc['structural'].append(path); return
        for k in b:
            if k in IGNORE:
                continue
            if k not in a:
                acc['missing'].append(path + '/' + k); continue
            walk(a[k], b[k], path + '/' + k, acc)
    elif isinstance(b, list):
        if not isinstance(a, list) or len(a) != len(b):
            acc['structural'].append(path); return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f'{path}[{i}]', acc)
    elif isinstance(b, bool) or b is None:
        if a != b:
            acc['string_mismatch'].append((path, a, b))
    elif isinstance(b, (int, float)):
        acc['n'] += 1
        if a == b:
            acc['exact'] += 1
        else:
            d = abs(a - b) / max(abs(b), 1e-300)
            if d > acc['max_rel'][0]:
                acc['max_rel'] = (d, path, a, b)
    elif isinstance(b, str):
        acc['strings'] += 1
        if a != b:
            acc['string_mismatch'].append((path, a[:80], b[:80]))


def compare(rep, pub):
    acc = dict(n=0, exact=0, strings=0, max_rel=(0., '', None, None), string_mismatch=[], missing=[], structural=[])
    walk(json.loads(rep.read_text()), json.loads(pub.read_text()), '', acc)
    acc['string_mismatch'] = acc['string_mismatch'][:20]
    return acc


def main(S, out):
    S = Path(S)
    res = dict(label='REPRODUCTION comparison: reproduced (g++ -O2, local Python) vs published JSON', files={})
    threshold = max(p.stat().st_mtime for sha in ['ff49e5d', '9636950', 'bb211d1'] for p in (S / sha).rglob('*.json')) + 3600
    for sha, pkt in [('ff49e5d', 'grut_causal_receiver_20261024'), ('9636950', 'grut_soft_receiver_20261023'),
                     ('bb211d1', 'grut_core_conversion_20261023')]:
        pubdir = S / sha / RF / pkt; repdir = S / 'tree' / RF / pkt
        for pub in sorted(pubdir.glob('*.json')):
            if pub.name in ('MANIFEST.json', 'ENVIRONMENT.json', 'RAW_DATA_MANIFEST.json'):
                continue
            rep = repdir / pub.name
            key = f'{sha}/{pub.name}'
            # git archive stamps files with commit times; anything rewritten by the reproduction is newer
            regenerated = rep.exists() and rep.stat().st_mtime > threshold
            if not regenerated:
                res['files'][key] = 'NOT REGENERATED'; continue
            c = compare(rep, pub)
            c['byte_identical'] = rep.read_bytes() == pub.read_bytes()
            res['files'][key] = c
    # ---- headline numbers recomputed from regenerated raw arrays (own code) ----
    C = S / 'tree' / RF / 'grut_causal_receiver_20261024'
    head = {}
    def hist(name):
        return np.load(C / 'data' / f'{name}_history.npz')['history']
    ret = {}
    for name, tmax in [('P025', 12000), ('P025_time', 12000), ('P0125', 12000), ('P00625', 12000),
                       ('P025_remote', 24000), ('P025_remote_time', 12000), ('I025', 12000)]:
        try:
            H = hist(name)
        except FileNotFoundError:
            continue
        late = H[:, 0] > .8 * tmax
        e0, q0 = H[0, 1], H[0, 2]
        ret[name] = dict(core60_late_mean_fraction=float(H[late, 5].mean() / q0), window=[float(H[late, 0].min()), float(H[late, 0].max())],
                         max_E_plus_receiver_err=float(np.max(abs((H[:, 1] + H[:, 3]) / e0 - 1))),
                         max_Q_plus_receiver_err=float(np.max(abs((H[:, 2] + H[:, 4]) / q0 - 1))))
        if name == 'I025':
            ret[name]['cut_energy_closure'] = float(np.max(abs((H[:, 8] + H[:, 10]) / H[0, 8] - 1)))
            ret[name]['cut_charge_closure'] = float(np.max(abs((H[:, 9] + H[:, 11]) / H[0, 9] - 1)))
            ret[name]['final_net_Eflux'] = float(H[-1, 10]); ret[name]['final_net_Qflux'] = float(H[-1, 11])
    head['retention_from_raw'] = ret
    if all(n in ret for n in ['P025_time', 'P0125', 'P00625']):
        v = [ret[n]['core60_late_mean_fraction'] for n in ['P025_time', 'P0125', 'P00625']]
        d1, d2 = v[0] - v[1], v[1] - v[2]
        head['grid'] = dict(values=v, spread_max_over_min=max(v) / min(v) - 1, two_finest_relative=abs(v[2] / v[1] - 1),
                            observed_order=float(np.log2(abs(d1 / d2))) if d1 * d2 > 0 else None,
                            richardson_order4_from_two_finest=v[2] + (v[2] - v[1]) / 15)
    if 'P025' in ret and 'P025_time' in ret:
        head['time_step_control_P025_vs_P025_time'] = ret['P025']['core60_late_mean_fraction'] / ret['P025_time']['core60_late_mean_fraction'] - 1
    if 'P025_remote_time' in ret and 'P025_time' in ret:
        head['layer_move_360_to_720'] = ret['P025_remote_time']['core60_late_mean_fraction'] / ret['P025_time']['core60_late_mean_fraction'] - 1
    # bitwise identity I025 vs P025_time, directly on arrays (also physical history columns 0..7)
    try:
        ident = {}
        for key in ['trace', 'final', 'memory']:
            a = np.load(C / 'data' / f'I025_{key}.npz')[key]; b = np.load(C / 'data' / f'P025_time_{key}.npz')[key]
            ident[key] = bool(np.array_equal(a, b)) and a.tobytes() == b.tobytes()
        a = hist('I025')[:, :8]; b = hist('P025_time')
        ident['history_physical_columns_0_7'] = a.tobytes() == np.ascontiguousarray(b).tobytes()
        head['I025_vs_P025_time_bitwise'] = ident
    except FileNotFoundError:
        pass
    # raw-array hashes vs published hashes
    raw = {}
    pubC = S / 'ff49e5d' / RF / 'grut_causal_receiver_20261024'
    for p in sorted(pubC.glob('*_RESULTS.json')):
        r = json.loads(p.read_text())
        if 'raw_numeric_sha256' not in r:
            continue
        for key, h in r['raw_numeric_sha256'].items():
            f = C / 'data' / f"{r['name']}_{key}.npz"
            raw[f"{r['name']}_{key}"] = (hashlib.sha256(np.load(f)[key].tobytes()).hexdigest() == h) if f.exists() else 'missing'
    head['ff49e5d_raw_hash_matches'] = dict(matched=sum(v is True for v in raw.values()), total=len(raw),
                                            mismatched=[k for k, v in raw.items() if v is not True])
    B = S / 'tree' / RF / 'grut_core_conversion_20261023'; pubB = S / 'bb211d1' / RF / 'grut_core_conversion_20261023'
    braw = {}
    for p in sorted(pubB.glob('[HN]*_RESULTS.json')):
        r = json.loads(p.read_text()); f = B / 'data' / f"{r['name']}_trace.npz"
        braw[r['name']] = (hashlib.sha256(np.load(f)['trace'].tobytes()).hexdigest() == r['trace_sha256']) if f.exists() else 'missing'
    head['bb211d1_trace_hash_matches'] = braw
    # bb211d1 nodal retention statistics recomputed from regenerated histories (own code)
    nret = {}
    for p in sorted(B.glob('N*_RESULTS.json')):
        r = json.loads(p.read_text()); f = B / 'data' / f"{r['name']}_history.csv"
        if not f.exists():
            continue
        H = np.loadtxt(f, delimiter=',', skiprows=1); late = H[:, 0] > .8 * r['parameters']['tmax']
        nret[r['name']] = float(H[late, 12].mean() / H[0, 2])
    head['bb211d1_nodal_retention_from_raw'] = nret
    if all(n in nret for n in ['N025', 'N0125', 'N00625']):
        v = [nret[n] for n in ['N025', 'N0125', 'N00625']]
        head['bb211d1_spatial'] = dict(values=v, relative_spread_over_mean=(max(v) - min(v)) / np.mean(v))
    pairs = [('N025', 'N025_time'), ('N0125', 'N0125_time'), ('N0125', 'N0125_abs'), ('N0125', 'N0125_box'), ('N0125_box', 'N0125_far')]
    head['bb211d1_controls'] = {f'{a}->{b}': nret[b] / nret[a] - 1 for a, b in pairs if a in nret and b in nret}
    res['headline_from_raw'] = head
    Path(out).write_text(json.dumps(res, indent=1, default=str) + '\n')
    for k, v in res['files'].items():
        if isinstance(v, str):
            print(k, v)
        else:
            print(k, 'numeric', v['n'], 'exact', v['exact'], 'maxrel %.3g' % v['max_rel'][0], v['max_rel'][1][:60],
                  'strmismatch', len(v['string_mismatch']), 'byte_identical', v['byte_identical'])
    print(json.dumps(head, indent=1, default=str))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
