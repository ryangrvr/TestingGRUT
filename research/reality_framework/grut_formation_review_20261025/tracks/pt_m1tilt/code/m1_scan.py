"""Track pt_m1tilt: transition scan under the GRUT-M1 potential (and controls) with the own adaptive DOP853 code.
Usage: python3 m1_scan.py <orbit tau20|late> <eps> <design> [--mode 3] [--epst 0.002] [--beta 0] [--phi0 qs]
                          [--rtol 1e-11] [--mend 40] [--nproc 2] [--pert 0]
  design: 'rand<N>'  = N runs with uniform random phase (rng seed fixed, same phases for every eps/orbit);
          'pair<N>'  = N runs at phases k/N (phase-smoothness check), plus neighbours at k/N + 1e-3;
          'ph<x>'    = one run at phase x.
  --pert d: multiply chi0 by (1 + d) (seed-perturbation pair; d = 0 is the unperturbed run).
One JSON per run in $SP/review/pt_m1tilt/runs (restartable); at most 2 worker processes."""
import sys, os, json, time, argparse, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1_common import *

OUT = os.path.join(SCR, 'runs')


def design(d):
    if d.startswith('rand'):
        n = int(d[4:]); rng = np.random.default_rng(20261101)
        return [float(rng.uniform()) for _ in range(n)]
    if d.startswith('pair'):
        n = int(d[4:]); ph = [k / n for k in range(n)]
        return ph + [p + 1e-3 for p in ph]
    if d.startswith('ph'):
        return [float(d[2:])]
    raise ValueError(d)


def fname(a):
    name, eps, ph, mode, epst, beta, phi0, rtol, mend, pert = a
    tag = '' if mend == 40.0 else f'_mend{mend:g}'
    tag += '' if pert == 0 else f'_pert{pert:.0e}'
    return os.path.join(OUT, f'{name}_m{mode}_et{epst:g}_b{beta:g}_phi0{phi0}_eps{eps:.0e}_ph{ph:.6f}_rt{rtol:.0e}{tag}.json')


def job(a):
    name, eps, ph, mode, epst, beta, phi0, rtol, mend, pert = a
    fn = fname(a)
    if os.path.exists(fn):
        return fn
    y4, info = orbit(name)
    st = phase_state(y4, info['T_r'], ph)
    t0 = time.time()
    p0 = 'qs' if phi0 == 'qs' else float(phi0)
    r = run_measure(st, info, eps, p0, mode=mode, epst=epst, beta=beta, rtol=rtol, M_end_fac=mend,
                    checkpoints=tuple(c for c in (5.0, 10.0, 20.0) if c < mend), phi0_rel=(1 + pert) if pert else None)
    r.update(orbit=name, phase=ph, pert=pert, J0_per_q=info['J_per_q'], runtime=time.time() - t0)
    with open(fn, 'w') as f:
        json.dump(r, f, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x.tolist())
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('orbit'); ap.add_argument('eps', type=float); ap.add_argument('design')
    ap.add_argument('--mode', type=int, default=3); ap.add_argument('--epst', type=float, default=EPST_M1)
    ap.add_argument('--beta', type=float, default=0.0); ap.add_argument('--phi0', default='qs')
    ap.add_argument('--rtol', type=float, default=1e-11); ap.add_argument('--mend', type=float, default=40.0)
    ap.add_argument('--nproc', type=int, default=2); ap.add_argument('--pert', type=float, default=0.0)
    A = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    phi0 = 'qs' if A.phi0 == 'qs' else f'{float(A.phi0):g}'
    jobs = [(A.orbit, A.eps, ph, A.mode, A.epst, A.beta, phi0, A.rtol, A.mend, A.pert) for ph in design(A.design)]
    t0 = time.time()
    with Pool(A.nproc) as p:
        for fn in p.imap_unordered(job, jobs):
            pass
    print('done', A.orbit, A.eps, A.design, 'mode', A.mode, 'epst', A.epst, 'beta', A.beta, 'phi0', phi0, len(jobs), 'runs',
          round(time.time() - t0, 1), 's')
