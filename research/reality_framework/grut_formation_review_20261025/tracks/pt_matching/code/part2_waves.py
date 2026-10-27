"""Part 2f: do waves (nonzero modes) switch on their mass adiabatically through the phi condensation?
Linear test modes chi_k'' = -(k^2 + m_w^2 + g^2 phi(s)^2) chi_k driven by the phi(s) of the homogeneous
transition runs (no backreaction; the coherent-X modulation of the wave mass, already present before the
transition, is left out so that only the switch-on is tested).  m_w^2 = 20 is the measured late Hartree mass
[Q, Codex WINDOW_base48_24_34].  Bogoliubov |beta_k|^2 = n_u + n_v - 1/2 from two real solutions started in
the instantaneous adiabatic vacuum.  Reference: the same modes with g^2 phi^2 replaced by the smooth g^2 M(s)^2.
Independent code."""
import sys, os, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *
from part2_floquet import orbit_info
from part2_transition import phase_state

KS = np.array([0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
MW2 = 20.0
CHECK_M_OVER_MPERI = [2.0, 4.0, 8.0]   # M_apo/M_peri = 8.06 for the tau-20 orbit
CHECK_M_OVER_MAPO = [2.0, 5.0]


def job(a):
    eps, ph, wslave = a
    y4 = tau20_state(); info = orbit_info(y4); st = phase_state(y4, info['T_r'], ph)
    Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']
    M_start = 0.5 * Mperi; kappa = eps * info['omega_r']
    k2 = np.repeat(KS ** 2, 2)
    c_set_waves(k2, MW2, wslave)
    par = make_par(mode=1, M0=M_start, kappa=kappa, growth=0)
    w0 = np.sqrt(KS ** 2 + MW2)
    y = np.zeros(8 + 2 * len(k2))
    y[:8] = [st[0], st[1], 1e-6 * Mperi, st[2], st[3], 0.0, 0.0, 0.0]
    for j, w in enumerate(w0):
        y[8 + 4 * j] = 1 / np.sqrt(2 * w)              # u: chi = 1/sqrt(2w), pi = 0
        y[8 + 4 * j + 3] = w / np.sqrt(2 * w)          # v: chi = 0, pi = w/sqrt(2w)
    y = np.ascontiguousarray(y)
    checks = [f * Mperi for f in CHECK_M_OVER_MPERI] + [f * Mapo for f in CHECK_M_OVER_MAPO]
    out = []
    for Mc in checks:
        sc = np.log(Mc / M_start) / kappa
        # average the adiabatic numbers over a short window (~20 wave periods) to remove phi-oscillation ripple
        c_run(par, y, sc, 0.08, 20.0, 0.5, 8)
        acc = np.zeros(len(KS)); nacc = 0
        Tw = 20 * 2 * np.pi / np.sqrt(MW2 + G2 * Mc ** 2)
        for t in np.linspace(sc, sc + Tw, 41)[1:]:
            c_run(par, y, t, 0.08, 20.0, 0.5, 8)
            M = M_start * np.exp(kappa * y[6])
            w = np.sqrt(KS ** 2 + MW2 + G2 * (M ** 2 if wslave else y[2] ** 2))
            chi_u, pi_u, chi_v, pi_v = y[8::4], y[9::4], y[10::4], y[11::4]
            nu = (pi_u ** 2 + w ** 2 * chi_u ** 2) / (2 * w); nv = (pi_v ** 2 + w ** 2 * chi_v ** 2) / (2 * w)
            acc += nu + nv - 0.5; nacc += 1
        out.append(dict(M=Mc, M_over_Mperi=Mc / Mperi, beta2=(acc / nacc).tolist(), phi_over_M=y[2] / (M_start * np.exp(kappa * y[6]))))
    c_set_waves([], 0.0, 0)
    return dict(eps=eps, phase=ph, wslave=wslave, checks=out)


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(eps, ph, ws) for eps in [1e-3, 1e-4] for ws in [0, 1] for ph in np.arange(8) / 8]
    with Pool(2) as p:
        rows = p.map(job, jobs)
    summ = []
    for eps in [1e-3, 1e-4]:
        for ws in [0, 1]:
            sel = [r for r in rows if r['eps'] == eps and r['wslave'] == ws]
            for ic in range(len(sel[0]['checks'])):
                b = np.array([r['checks'][ic]['beta2'] for r in sel])
                ent = dict(eps=eps, background=['actual phi(s)', 'smooth g M(s)'][ws], M_over_Mperi=sel[0]['checks'][ic]['M_over_Mperi'],
                           k=KS.tolist(), beta2_mean=b.mean(0).tolist(), beta2_min=b.min(0).tolist(), beta2_max=b.max(0).tolist())
                summ.append(ent)
                print(eps, ws, f"M/Mperi {ent['M_over_Mperi']:.1f}", ' '.join(f'{x:.2e}' for x in ent['beta2_mean']))
    save_json('waves_switch_on.json', dict(k=KS.tolist(), m_w2=MW2, summary=summ,
                                            note='|beta_k|^2 from the adiabatic vacuum; phase-averaged over 8 initial orbital phases; seed 1e-6'))
    print('elapsed', time.time() - t0)
