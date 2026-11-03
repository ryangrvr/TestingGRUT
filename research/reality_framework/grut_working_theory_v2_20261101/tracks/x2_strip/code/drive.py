import json, os, subprocess, sys
from pathlib import Path
TR = Path('/home/user/TestingGRUT/research/reality_framework/grut_working_theory_v2_20261101/tracks/x2_strip')
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/x2_strip')
node = sys.argv[1]; base = 'x2_' + node
c0s = [1.00, 1.04, 1.08, 1.12, 1.16, 1.20] + [round(1.30 + 0.01 * i, 2) for i in range(16)] + [1.355, 1.365, 1.375, 1.385, 1.395]
keys = ['%.3f:%+.0f' % (c, d) for d in (-4, -2, 0, 2, 4) for c in sorted(set(c0s))]
def have():
    p = TR / 'results' / ('scan_%s.json' % base)
    return set(json.load(open(p)).get('runs', {})) if p.exists() else set()
todo = [k for k in keys if k not in have()]
env = dict(os.environ, LFC_SCRATCH=str(SCR / ('sp_' + node)), PYTHONDONTWRITEBYTECODE='1')
for i in range(0, len(todo), 9):
    chunk = [k.replace(':+', ':') for k in todo[i:i + 9]]
    r = subprocess.run([sys.executable, '-B', str(TR / 'code' / 'x2scan.py'), node, base] + chunk, env=env, capture_output=True, text=True)
    print(node, i, len(todo), 'rc', r.returncode, flush=True)
    if r.returncode != 0:
        print(r.stderr[-2000:], flush=True); break
print(node, 'DONE', len(have()), flush=True)
