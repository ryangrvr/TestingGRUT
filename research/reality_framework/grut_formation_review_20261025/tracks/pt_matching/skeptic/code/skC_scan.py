"""Skeptic part C: transition scan with the independent adaptive DOP853 code (sk_dop.c).
Usage: python3 skC_scan.py <orbit tau20|late> <eps> <design> [growth] [beta] [rtol]
  design: 'grid16' = track's 16 phases x seeds (1e-2,1e-6,1e-10);
          'rand<N>' = N runs with uniform random phase and log-uniform seed in [1e-12, 1e-2] (rng seed fixed);
          'seed<exp>' = 16 phases at seed 10^-exp;  'grid16s1' = 16 phases at seed 1e-6.
One JSON per run in $SP/review/pt_matching_skeptic/runs (restartable); 2 worker processes."""
import sys, os, json, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *

OUT = os.path.join(SCR, 'runs')
LATE_R0 = 31.494634874218747


def orbit(name):
    y4 = tau20_state()
    if name == 'late':
        q = charge(y4); y4 = np.array([LATE_R0, 0.0, 0.0, q / LATE_R0])
    info = orbit_info(y4, nsamp=1024)
    for k in ('n', 'w', 'P'): info.pop(k)
    return y4, info


def design(d):
    if d == 'grid16':
        return [(ph, sd) for sd in [1e-2, 1e-6, 1e-10] for ph in np.arange(16) / 16]
    if d == 'grid16s1':
        return [(ph, 1e-6) for ph in np.arange(16) / 16]
    if d.startswith('seed'):
        e = float(d[4:]); return [(ph, 10.0 ** (-e)) for ph in np.arange(16) / 16]
    if d.startswith('rand'):
        n = int(d[4:]); rng = np.random.default_rng(20261027)
        return [(float(rng.uniform()), float(10 ** rng.uniform(-12, -2))) for _ in range(n)]
    raise ValueError(d)


def job(a):
    name, eps, ph, sd, growth, beta, rtol = a
    mend = float(os.environ.get('MEND', '40'))
    tag = '' if mend == 40.0 else f'_mend{mend:g}'
    fn = os.path.join(OUT, f'{name}_g{growth}_b{beta:g}_eps{eps:.0e}_ph{ph:.6f}_seed{sd:.3e}_rt{rtol:.0e}{tag}.json')
    if os.path.exists(fn):
        return fn
    y4, info = orbit(name)
    st = phase_state(y4, info['T_r'], ph)
    t0 = time.time()
    r = run_measure(st, info, eps, sd, growth=growth, beta=beta, rtol=rtol, M_end_fac=mend,
                    checkpoints=tuple(c for c in (5.0, 10.0, 20.0) if c < mend))
    r.update(orbit=name, phase=ph, J0_per_q=info['J_per_q'], runtime=time.time() - t0)
    with open(fn, 'w') as f:
        json.dump(r, f, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x.tolist())
    return fn


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    name, eps, des = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    growth = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    beta = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    rtol = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-11
    jobs = [(name, eps, ph, sd, growth, beta, rtol) for ph, sd in design(des)]
    t0 = time.time()
    with Pool(2) as p:
        for _ in p.imap_unordered(job, jobs):
            pass
    print('done', name, eps, des, growth, beta, len(jobs), 'runs', round(time.time() - t0, 1), 's')
