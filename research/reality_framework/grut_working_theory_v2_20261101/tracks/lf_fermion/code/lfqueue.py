"""Foreground 2-worker queue for resumable lfline jobs (each job checkpoints; rerun the same command to continue).
Usage: python3 -B lfqueue.py BUDGET_S JOB [JOB ...], JOB = 'RUN:LINE:SP:NR:SEED:TEND[:Nc:Nf:tref:cfl:tag]'.
Claude, lf_fermion (issue #13 context)."""
import subprocess, sys, time, json
from pathlib import Path
sys.dont_write_bytecode = True
import lfcommon as LC

budget = float(sys.argv[1])
jobs = [j.split(':') for j in sys.argv[2:]]
t0 = time.time()


def finished(j):
    run, line, sp, tend = j[0], j[1], j[2], float(j[5])
    tag = j[10] if len(j) > 10 else 'base'
    f = LC.SP / 'lfq' / run / ('%s_l%s_%s.json' % (tag, line, sp))
    if not f.exists():
        return False
    m = json.loads(f.read_text())['meas']
    return bool(m) and m[-1]['t'] >= tend - 2e-3


pending = [j for j in jobs if not finished(j)]
running = []
while (pending or running):
    left = budget - (time.time() - t0)
    while pending and len(running) < 2 and left > 60:
        j = pending.pop(0)
        args = j[:6] + (j[6:] if len(j) > 6 else ['4096', '8192', '0.34', '0.5', 'base'])
        while len(args) < 11:
            args.append(['4096', '8192', '0.34', '0.5', 'base'][len(args) - 6])
        p = subprocess.Popen(['python3', '-B', 'lfline.py'] + args + ['%.0f' % (left - 40)], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, env=dict(__import__('os').environ, OMP_NUM_THREADS='1'))
        running.append((p, j))
    time.sleep(2)
    for p, j in list(running):
        if p.poll() is not None:
            print(p.stdout.read().strip()[-300:], flush=True)
            running.remove((p, j))
            if not finished(j) and budget - (time.time() - t0) > 90:
                pending.insert(0, j)
    if not running and (budget - (time.time() - t0)) <= 60:
        break
print('left pending:', [':'.join(j[:3]) for j in pending + [j for _, j in running]])
