# Runs the pt_matching track's own integrator (read-only import) on the skeptic's random design; output to skeptic scratch.
import sys, os, json, time, numpy as np
from multiprocessing import Pool
TR = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks/pt_matching/code'
SK = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks/pt_matching/skeptic/code'
OUT = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/pt_matching_skeptic/trackcode'
def job(a):
    sys.path.insert(0, TR)
    from pm_common import tau20_state, charge
    from part2_floquet import orbit_info, peri_state
    from part2_transition import transition_run, phase_state
    i, ph, sd, eps, orb = a
    fn = f'{OUT}/{orb}_eps{eps:.0e}_{i:04d}.json'
    if os.path.exists(fn): return
    y4 = tau20_state()
    if orb == 'late': y4 = peri_state(charge(y4), 31.49463487422226)
    info = orbit_info(y4)
    st = phase_state(y4, info['T_r'], ph)
    r = transition_run(st, info, eps, sd, keep_trace=False)
    json.dump(dict(phase=ph, seed=sd, eps=eps, nminus=r['nminus'], Nphi=r['Nphi'], qrel=r['qrel'], work_rel=r['work_rel']), open(fn, 'w'))
if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    orb, eps, n = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    rng = np.random.default_rng(20261027)
    pairs = [(float(rng.uniform()), float(10 ** rng.uniform(-12, -2))) for _ in range(n)]
    t0 = time.time()
    with Pool(2) as p: p.map(job, [(i, ph, sd, eps, orb) for i, (ph, sd) in enumerate(pairs)])
    import glob
    rows = [json.load(open(f)) for f in glob.glob(f'{OUT}/{orb}_eps{eps:.0e}_*.json')]
    a = np.array([r['nminus'] for r in rows]); b = np.array([r['Nphi'] for r in rows])
    print(orb, eps, len(a), 'track-code n- %.4f +- %.4f  Nphi %.4f +- %.4f' % (a.mean(), a.std(ddof=1)/np.sqrt(len(a)), b.mean(), b.std(ddof=1)/np.sqrt(len(b))), time.time()-t0)
