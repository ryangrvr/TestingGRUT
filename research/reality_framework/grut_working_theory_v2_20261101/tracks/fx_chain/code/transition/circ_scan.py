"""fx_chain item 3 (Claude): does the M1 phi transition preserve a CIRCULAR X condensate?

Uses the copied pt_m1tilt engine (m1_common.py + m1_dop.c, own adaptive DOP853; mode 3 = GRUT-M1 potential
V = lX/4 r^4 + g^2/2 chi^2 r^2 + lP/4 (chi^2 - M^2)^2 - eps_t M (chi - M)^3, conformal chi = a phi, M = a v), the same
drift M = M0 exp(kappa s), kappa = eps omega_r, start M0 = 0.5 M_c (M_c = g r/sqrt(lP)), phi at the instantaneous
tilted minimum, the same end window (60 X periods at 40 M_c) and the same measures n_-/q and N_phi/q.

Orbits: the tau-20 charge q = 8117.8816 on a circular orbit r0 = (q/sqrt(lX))^(1/3) (J_r = 0) and nearly circular
orbits r(0) = r0 (1 + delta), p_theta = q (radial amplitude ~ delta).  omega_r for the drift is the small-oscillation
radial frequency sqrt(6 lX) r0 of the circular orbit (eccentric cases: from the exact radial quadrature).
[C] homogeneous classical X + phi; not a cosmological result.
Usage: python3 circ_scan.py [quick]
"""
import sys, json, time, os
import numpy as np
from multiprocessing import Pool
import m1_common as C

LX = C.LX


def circ_state(delta=0.0):
    q = C.charge(C.tau20_state())
    r0 = (q / np.sqrt(LX)) ** (1 / 3)
    y4 = np.array([r0 * (1 + delta), 0.0, 0.0, q / (r0 * (1 + delta))])
    if delta == 0.0:
        info = dict(q=q, E=C.energy4(y4), J_per_q=0.0, T_r=2 * np.pi / (np.sqrt(6 * LX) * r0), omega_r=np.sqrt(6 * LX) * r0,
                    Omega=np.sqrt(LX) * r0, rmin=r0, rmax=r0, r2mean=r0 * r0)
    else:
        info = C.orbit_info(y4)
    return y4, info


def job(args):
    delta, eps, ph, epst, extra = args
    y4, info = circ_state(delta)
    if ph and delta:
        y4 = C.phase_state(y4, info['T_r'], ph)
    t0 = time.time()
    kw = dict(epst=epst, mode=3, atol=1e-12)
    kw.update(extra)
    r = C.run_measure(y4, info, eps, 'qs', **kw)
    r.pop('trace_M_chiOverM', None)
    return dict(delta=delta, eps=eps, phase=ph, epst=epst, extra=extra, J_r_over_q_initial=info['J_per_q'],
                omega_r=info['omega_r'], rmin=info['rmin'], rmax=info['rmax'], nminus=r['nminus'], Nphi=r['Nphi'],
                windows=r['windows'], work_rel=r['work_rel'], qrel=r['qrel'], steps=r['steps'], sec=time.time() - t0,
                M_half=r.get('M_half_over_Mperi'), M_settle=r.get('M_settle_over_Mperi'), chi0_over_Mperi=r['chi0_over_Mperi'])


def main():
    grp = sys.argv[1] if len(sys.argv) > 1 else 'A'
    E = C.EPST_M1
    jobs = []
    if grp == 'quick':
        jobs = [(0.0, 1e-3, 0.0, E, {}), (0.0, 1e-4, 0.0, E, {})]
    elif grp == 'A':                                                     # exactly circular, all rates + controls
        for eps in (1e-6, 3e-6, 1e-3, 1e-4, 3e-5, 1e-5):
            jobs.append((0.0, eps, 0.0, E, {}))
        for eps in (1e-3, 1e-5):
            jobs.append((0.0, eps, 0.0, E, dict(rtol=1e-12)))            # rtol control
            jobs.append((0.0, eps, 0.0, 0.0022, {}))                     # tilt sensitivity (+10%)
            jobs.append((0.0, eps, 0.0, 0.0018, {}))                     # (-10%)
            jobs.append((0.0, eps, 0.0, E, dict(atol=1e-10)))            # atol control
    elif grp == 'B':                                                     # nearly circular, eps >= 1e-5
        for delta in (1e-3, 1e-2, 3e-2, 0.1):
            for eps in (1e-3, 1e-4, 1e-5):
                for ph in (0.0, 0.25, 0.5, 0.75):
                    jobs.append((delta, eps, ph, E, {}))
    elif grp == 'C':                                                     # nearly circular, eps = 1e-6
        for delta in (1e-2, 0.1):
            for ph in (0.0, 0.5):
                jobs.append((delta, 1e-6, ph, E, {}))
    out = []
    with Pool(2) as P:
        for r in P.imap_unordered(job, jobs):
            out.append(r)
            print('delta=%-6g eps=%-6g ph=%.2f epst=%.4f %s Jr0/q=%.3e  n-/q=%.4e  Nphi/q=%.4e  W=%.1e q=%.1e  %.1fs'
                  % (r['delta'], r['eps'], r['phase'], r['epst'], r['extra'], r['J_r_over_q_initial'], r['nminus'], r['Nphi'],
                     r['work_rel'], r['qrel'], r['sec']), flush=True)
    C.save_json('circ_transition_%s.json' % grp, dict(note='fx_chain item 3; M1 potential, circular and nearly circular '
                                                           'tau-20-charge orbits; group ' + grp, runs=out))


if __name__ == '__main__':
    main()
