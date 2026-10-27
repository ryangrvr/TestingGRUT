"""Reproduction comparison (st_rerun): regenerated Codex b73de64 outputs vs the published copies.

Usage: python3 -I compare_outputs.py <published_dir> <regenerated_packet_dir> <out_json>
Walks every JSON leaf (numbers, strings, bools, None) and every NPZ array and reports the
largest relative/absolute differences. Read-only on both directories.
"""
import json, sys, math
from pathlib import Path
import numpy as np

JSONS = ['RESPONSE_CHECKS.json', 'RESPONSE_RESULTS.json', 'RESULTS.json', 'SOURCE_CLOCKS.json',
         'QUARTIC_SCAN.json', 'QUARTIC_CHECKS.json', 'VERIFICATION_RESULTS.json']
NPZS = ['quiet3.npz', 'quiet5.npz', 'gauge3.npz', 'gauge5.npz']


def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(dict(path=path + '/' + k, kind='missing_key', published=k in a, regenerated=k in b))
                continue
            walk(a[k], b[k], path + '/' + k, out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(dict(path=path, kind='length', published=len(a), regenerated=len(b)))
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f'{path}[{i}]', out)
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str):
        if a != b:
            out.append(dict(path=path, kind='value', published=a, regenerated=b))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a != b:
            rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
            out.append(dict(path=path, kind='number', published=a, regenerated=b, abs_diff=abs(a - b), rel_diff=rel))
    else:
        if a != b:
            out.append(dict(path=path, kind='type', published=repr(a)[:80], regenerated=repr(b)[:80]))


def count_leaves(x):
    if isinstance(x, dict):
        return sum(count_leaves(v) for v in x.values())
    if isinstance(x, list):
        return sum(count_leaves(v) for v in x)
    return 1


def main():
    pub, reg, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    report = {}
    for name in JSONS:
        a = json.loads((pub / name).read_text()); b = json.loads((reg / name).read_text())
        diffs = []; walk(a, b, '', diffs)
        nums = [d for d in diffs if d['kind'] == 'number']
        other = [d for d in diffs if d['kind'] != 'number']
        nums.sort(key=lambda d: -d['rel_diff'])
        report[name] = dict(leaves=count_leaves(a), differing_numbers=len(nums), non_numeric_differences=other,
                            max_rel_diff=nums[0]['rel_diff'] if nums else 0.0,
                            largest=nums[:12])
    for name in NPZS:
        a = np.load(pub / name); b = np.load(reg / name); rows = {}
        for k in sorted(set(a.files) | set(b.files)):
            if k not in a.files or k not in b.files:
                rows[k] = 'missing'; continue
            x = np.asarray(a[k], float); y = np.asarray(b[k], float)
            if x.shape != y.shape:
                rows[k] = dict(shape=[x.shape, y.shape]); continue
            d = np.abs(x - y)
            rows[k] = dict(shape=list(x.shape), identical=bool(np.array_equal(x, y)), max_abs=float(d.max()),
                           max_scaled=float(np.max(d / (1 + np.abs(x)))))
        report[name] = rows
    outp.write_text(json.dumps(report, indent=1) + '\n')
    for name in JSONS:
        r = report[name]
        print(name, 'leaves', r['leaves'], 'differing numbers', r['differing_numbers'], 'max rel', r['max_rel_diff'],
              'non-numeric', len(r['non_numeric_differences']))
    for name in NPZS:
        print(name, {k: (v if isinstance(v, str) else (v.get('identical'), v.get('max_scaled'))) for k, v in report[name].items()})


if __name__ == '__main__':
    main()
