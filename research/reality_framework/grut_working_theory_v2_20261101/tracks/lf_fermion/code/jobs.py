"""lf_fermion lattice jobs (copy-free driver of lf3d through lfcommon). Usage: python3 -B jobs.py TAG [budget_s].
Recording runs (REC_*) carry passive closure tracers and line dumps; LF_* runs use the calibrated mixed closure.
Claude, lf_fermion track (issue #13 context). Every closure, seed and line choice is a model choice."""
import sys
sys.dont_write_bytecode = True
import lfcommon as LC

LINES = [(0, 5, 17), (0, 29, 43), (0, 51, 9), (0, 37, 60), (1, 11, 33), (1, 47, 2), (1, 22, 55), (1, 60, 28),
         (2, 3, 40), (2, 19, 8), (2, 42, 24), (2, 56, 50)]
REC = ['tracers 1', 'line_every 19'] + ['line %d %d %d' % l for l in LINES]
JOBS = {
    'REC_P0': dict(closure='pauli', Bglobal=0, t_end=1.5, extra=REC),
    'REC_GM': dict(closure='pauli', Bglobal=2, t_end=1.5, extra=REC),
    'REC_H84_P0': dict(closure='pauli', Bglobal=0, t_end=1.1, H_GeV=8.4e10, extra=REC),
}


def lf_job(lam, **kw):
    d = dict(closure='pauli', Bglobal=3, t_end=3.0, extra=['lam_pb %.6f' % lam])
    d.update(kw)
    return d


if __name__ == '__main__':
    tag = sys.argv[1]
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 560.0
    if tag.startswith('KF'):
        # KF<var>[_H84][_N96|_L63]: calibrated capacity-multiplier closure (local blocking, kappa_j(t) from
        # results/kappa_REC_P0.json), recorded lines and tracers to t = 3 for the self-consistency check
        import kappa_table as KT
        parts = tag.split('_')
        var = parts[0][2:] or 'mid'
        kw = dict(closure='pauli', Bglobal=0, t_end=3.0, extra=REC + KT.par_lines(KT.table('REC_P0')[var]))
        if 'H84' in parts:
            # fragmentation is later at 8.4e10 (e2e_lattice: f_coh < 0.1 at 0.74): kappa fitted on REC_H84_P0 with
            # the switch window 0.65-0.75 [A]
            T84 = KT.table('REC_H84_P0', (0.65, 0.75))[var]
            kw.update(H_GeV=8.4e10, extra=['tracers 0'] + KT.par_lines(T84))
        if 'nl' in parts and 'H84' not in parts:
            kw['extra'] = ['tracers 0'] + KT.par_lines(KT.table('REC_P0')[var])
        if 's2' in parts:
            kw['seed'] = 2
        if 'N96' in parts:
            kw.update(N=96, L=0.42, extra=['tracers 0'] + KT.par_lines(KT.table('REC_P0')[var]))
        if 'L63' in parts:
            kw.update(N=96, L=0.63, extra=['tracers 0'] + KT.par_lines(KT.table('REC_P0')[var]))
        JOBS[tag] = kw
    if tag.startswith('LF_'):
        # LF_<lam>[_suffix]: calibrated mixed closure at lam
        parts = tag.split('_')
        kw = {}
        if 'N96' in parts:
            kw.update(N=96, L=0.42)
        if 'L63' in parts:
            kw.update(N=96, L=0.63)
        if 'H84' in parts:
            kw.update(H_GeV=8.4e10)
        JOBS[tag] = lf_job(float(parts[1]), **kw)
    job = dict(JOBS[tag])
    rd = LC.RUNS / tag
    if not (rd / 'run.par').exists():
        LC.make_run(tag, gw=0, snaps=[3.0] if job.get('t_end', 3) >= 3 else [], **job)
    done = LC.run(tag, budget=budget, threads=2)
    print('done' if done else 'partial', tag)
