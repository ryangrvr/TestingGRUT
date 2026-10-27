"""Part 2c: scan of the phi-condensation transition over eps, initial orbital phase and phi seed.
Usage: python3 part2_scan.py <orbit: tau20|late> <eps> <nphase> [i0 i1] ; writes one JSON per run to
$SP/review/pt_matching/scan/ (restartable; existing runs are skipped).  Two worker processes."""
import sys, os, json, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *
from part2_floquet import orbit_info, peri_state
from part2_transition import transition_run, phase_state

OUT = os.path.join(SCRATCH, 'scan')
SEEDS = [1e-2, 1e-6, 1e-10]
LATE_R0 = None


def orbit(name):
    y4 = tau20_state()
    if name == 'late':
        fl = json.load(open(os.path.join(RESULTS, 'floquet.json')))
        y4 = peri_state(charge(y4), fl['late_A_match']['r0'])
    return y4, orbit_info(y4)


def job(args):
    name, eps, ph, seed, growth, h0, mend = args
    tag = '' if mend == 40.0 else f'_mend{mend:g}'
    fn = os.path.join(OUT, f'{name}_g{growth}_eps{eps:.0e}_ph{ph:.4f}_seed{seed:.0e}_h{h0}{tag}.json')
    if os.path.exists(fn):
        return fn
    y4, info = orbit(name)
    st = phase_state(y4, info['T_r'], ph)
    t0 = time.time()
    r = transition_run(st, info, eps, seed, growth=growth, h0=h0, M_end_fac=mend)
    tr = r.pop('trace')
    Mperi = np.sqrt(G2 / LP) * info['rmin']
    big = np.where(np.abs(tr['phi_over_M']) > 0.1)[0]
    r.update(orbit=name, phase=ph, J0_per_q=info['J_per_q'], A0_per_q=info['A_minus_per_q'],
             M_escape_over_Mperi=float(tr['M'][big[0]] / Mperi) if len(big) else None, runtime=time.time() - t0,
             trace=dict(M_over_Mperi=(tr['M'] / Mperi).tolist(), nminus=tr['nminus'].tolist(), Nphi=tr['Nphi'].tolist(), phi_over_M=tr['phi_over_M'].tolist()))
    with open(fn, 'w') as f:
        json.dump(r, f, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x.tolist())
    return fn


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    name, eps, nph = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    growth = int(os.environ.get('GROWTH', '0')); h0 = float(os.environ.get('H0', '0.08'))
    seeds = [float(x) for x in os.environ.get('SEEDS', ','.join(str(s) for s in SEEDS)).split(',')]
    mend = float(os.environ.get('MEND', '40'))
    jobs = [(name, eps, ph, sd, growth, h0, mend) for sd in seeds for ph in np.arange(nph) / nph]
    if len(sys.argv) > 5:
        jobs = jobs[int(sys.argv[4]):int(sys.argv[5])]
    t0 = time.time()
    with Pool(2) as p:
        for fn in p.imap_unordered(job, jobs):
            pass
    print('done', name, eps, len(jobs), 'jobs', time.time() - t0, 's')
