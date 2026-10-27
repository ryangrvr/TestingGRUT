"""Part 2d: controls for the transition result.
 (a) slaved mass: X mass g M(s) imposed, no phi field (mode 2)  -> J_r must be conserved;
 (b) Born-Oppenheimer: phi at its instantaneous minimum (mode 3) -> X moves in a slowly changing central potential;
 (c) step-size halving of the full model (statistics), (d) RD-linear growth instead of exponential.
Independent code."""
import sys, os, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *
from part2_floquet import orbit_info, peri_state
from part2_transition import transition_run, phase_state
import json


def initial_J(st, mode, M_start):
    q = charge(st); r2 = st[0] ** 2 + st[1] ** 2
    m2 = G2 * M_start ** 2 if mode == 2 else 0.0
    E = 0.5 * (st[2] ** 2 + st[3] ** 2) + 0.5 * m2 * r2 + 0.25 * LX * r2 * r2
    return orbit_quadrature(E, q, LX, m2)['J_r'] / q


def job(a):
    orb, mode, eps, ph = a
    y4 = tau20_state()
    if orb == 'late':
        fl = json.load(open(os.path.join(RESULTS, 'floquet.json'))); y4 = peri_state(charge(y4), fl['late_A_match']['r0'])
    info = orbit_info(y4)
    st = phase_state(y4, info['T_r'], ph)
    Mstart = 0.5 * np.sqrt(G2 / LP) * info['rmin']
    J0 = initial_J(st, mode, Mstart)
    r = transition_run(st, info, eps, 0.0, mode=mode, keep_trace=False)
    return dict(orbit=orb, mode=mode, eps=eps, phase=ph, J0=J0, nminus=r['nminus'], dev=r['nminus'] - J0, qrel=r['qrel'], work_rel=r['work_rel'])


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(orb, mode, eps, ph) for orb in ['tau20', 'late'] for mode in [2, 3] for eps in [1e-2, 1e-3, 1e-4] for ph in np.arange(8) / 8]
    with Pool(2) as p:
        rows = p.map(job, jobs)
    summ = []
    for orb in ['tau20', 'late']:
        for mode in [2, 3]:
            for eps in [1e-2, 1e-3, 1e-4]:
                sel = [r for r in rows if r['orbit'] == orb and r['mode'] == mode and r['eps'] == eps]
                d = np.array([r['dev'] for r in sel]); J0 = np.mean([r['J0'] for r in sel])
                summ.append(dict(orbit=orb, control=['', '', 'slaved_mass', 'born_oppenheimer'][mode], eps=eps, J0_mean=J0,
                                 nminus_mean=np.mean([r['nminus'] for r in sel]), max_abs_dev=np.abs(d).max(), mean_dev=d.mean(),
                                 max_qrel=max(r['qrel'] for r in sel), max_work_rel=max(r['work_rel'] for r in sel)))
                print(summ[-1])
    save_json('transition_controls.json', dict(rows=rows, summary=summ))
    print('elapsed', time.time() - t0)
