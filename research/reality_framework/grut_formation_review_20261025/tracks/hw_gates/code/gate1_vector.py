#!/usr/bin/env python3
"""hw_gates gate 1 (Claude, own code): longitudinal hidden-vector production at the ls_scaling spot point.

Piece A (inflationary Bunch-Davies Goldstone, Graham-Mardon-Rajendran type [R]) with the evolving Higgs norm:
  single-Higgs-direction longitudinal Proca mode S_k = 1/2 int dtau [z^2|B'|^2 - a^2 m^2 |B|^2],
  z^2 = a^2 m^2/(k^2 + a^2 m^2); inflation a = 1/(2 - tau) (H = 1), radiation a = tau (H = 1/a^2), a = 1 at release.
  Pair (B, Pi = z^2 B'): B' = Pi/z^2, Pi' = -a^2 m^2 B; adiabatic vacuum at k/aH = 1000; late n_k from the
  Hamiltonian over the frequency. N_A = int dlnk k^3 n_k/(2 pi^2) in H_i^3 units (comparable with N_X).
  Profiles: P0 constant vacuum mass 6 e h (GMR control); P1 m^2 = m_vac^2 + (m_inf^2 - m_vac^2)/a^2 (the X part of
  the Higgs norm redshifting away; at this point F_inf/F_vac = 1.11 only, because Sigma already dominates F).
  (Same formulation as our sp_vector proca_modes.py; re-implemented with run-time parameters and an e scan.)
Deposit: one vector -> one X Xbar pair at <a>_A = sqrt(pi H_i/(2 Gamma_A)); R_A = N_A (m_A/H)<a>_A/E_X.
e scan: N_A ∝ e^{1/2} (GMR) while the deposit per vector ∝ e^{-1/2}: R_A is e-independent while the vector is light
during inflation (e < 1/F_inf = 3.8e-3); above that the mode is heavy and piece A disappears (computed, not assumed).
Floquet check for a larger e: transverse vector modes A'' + [k^2 + 2e^2|chi_X|^2 + 36 e^2 a^2 h^2]A = 0 over the late
eccentric X orbit (E_X, N_X of the spot point), Sigma mass term frozen at a = 1, 2, 5 (it only detunes).
Usage: python3 -B gate1_vector.py <out.json>
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spotmodel as S

P = S.params()
H = S.H_SPOT
# spot-point release window means (own hwgrel run, relaxed preparation, a = 30-60) [C]
EX, NX, NYH, NYL, NSIG = 38398.22, 2788.161, 0.01725, 0.33446, 3.161e-09
RL = 87.85
XINF = 80.2020754      # aligned inflationary |X| (own)
F2_INF = 2 * (80.2020754**2 + 9 * 3.44245486**2 + 36 * 28.31074753**2)
KGRID = np.geomspace(1e-5, 1.0, 49)


def evolve_mode(k, m_infl, mfun, ends):
    a0 = k / 1000.0
    tau0 = 2.0 - 1.0 / a0

    def a_of(t):
        return 1.0 / (2.0 - t) if t < 1.0 else t

    def mass(t):
        return m_infl if t < 1.0 else mfun(a_of(t))

    def rhs(t, y):
        a = a_of(t); m = mass(t); am2 = (a * m)**2
        z2 = am2 / (k * k + am2)
        return [y[2] / z2, y[3] / z2, -am2 * y[0], -am2 * y[1]]

    a = a0; am2 = (a * m_infl)**2; z2 = am2 / (k * k + am2); om = np.sqrt(k * k + am2)
    B = 1.0 / np.sqrt(2 * z2 * om); Pi = -1j * np.sqrt(z2 * om / 2)
    y = np.array([B.real, B.imag, Pi.real, Pi.imag])
    segs = sorted(set([tau0, 1.0] + list(ends)))
    out = {}
    for t1, t2 in zip(segs[:-1], segs[1:]):
        sol = solve_ivp(rhs, (t1, t2), y, method='DOP853', rtol=1e-10, atol=1e-13 * float(np.max(np.abs(y))))
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        if t2 in ends:
            a = a_of(t2); m = mass(t2); am2 = (a * m)**2; z2 = am2 / (k * k + am2); om = np.sqrt(k * k + am2)
            Bc = y[0] + 1j * y[1]; Pc = y[2] + 1j * y[3]
            out[t2] = (abs(Pc)**2 / (2 * z2) + am2 * abs(Bc)**2 / 2) / om - 0.5
    return out


def mode_job(args):
    k, m_inf, m_vac, prof = args
    if prof == 'P0':
        mfun = lambda a: m_vac
    else:
        mfun = lambda a: np.sqrt(m_vac**2 + (m_inf**2 - m_vac**2) / a**2)
    # n_k freezes once a >> a_* (adiabatic: omega'/omega^2 ~ 1/(a^2 m) << 1 for every k), so cut at 20 and 30 a_*
    a_star = max(1.0 / np.sqrt(m_vac), 1.0)
    ends = [max(20 * a_star, 30.0), max(30 * a_star, 45.0)]
    r = evolve_mode(k, m_inf, mfun, ends)
    n1, n2 = r[ends[0]], r[ends[1]]
    return dict(k=k, n_k=n2, rel_change=abs(n2 - n1) / max(abs(n2), 1e-300))


def spectrum(e, prof, pool):
    m_vac = 6 * e * P['h']; m_inf = e * np.sqrt(F2_INF)
    rows = pool.map(mode_job, [(k, m_inf, m_vac, prof) for k in KGRID])
    lk = np.log(KGRID)
    integ = np.array([kk**3 * r['n_k'] for kk, r in zip(KGRID, rows)]) / (2 * np.pi**2)
    NA = float(np.trapezoid(integ, lk))
    ipk = int(np.argmax(integ))
    return dict(e=e, profile=prof, m_vac_over_H=m_vac, m_inf_over_H=m_inf, a_star=1 / np.sqrt(m_vac), k_star=np.sqrt(m_vac),
                N_A=NA, peak_k=float(KGRID[ipk]), dN_dlnk_peak=float(integ[ipk]),
                max_rel_change=float(max(r['rel_change'] for r in rows)), modes=rows)


def gmr_formula(e):
    """Omega_A/Omega_DM = sqrt(m/6e-6 eV) (H_I/1e14 GeV)^2 [R]; converted to comoving N_A per H_i^3 at instant reheating."""
    m_A = 6 * e * P['h'] * H
    Om = np.sqrt(m_A / 6e-15) * (H / 1e14)**2
    rho_s = 4.4e-10 * Om
    n_s = rho_s / m_A
    s, T = S.entropy_RH(H)
    return dict(Omega_A_over_DM=Om, nA_over_s=n_s, N_A=n_s * s / H**3, T_RH_GeV=T)


def ledger_e(e, NA):
    p = S.params(e=e)
    w = dict(EX=EX, NYH=NYH, NYL=NYL, Nsig=NSIG)
    L = S.ledger(p, w)
    RA = NA * L['deposit_per_A_in_EX']
    return dict(R_Y=L['R_Y'], R_Sigma=L['R_Sigma'], R_cascade=L['R_cascade'], R_A=RA, R_total=L['R_central'] + RA,
                R_L=RL, total_over_RL=(L['R_central'] + RA) / RL, deposit_per_A_in_EX=L['deposit_per_A_in_EX'],
                pairs_per_net_charge=NA / NX, allowance_pairs_per_net=(RL - L['R_central']) / L['deposit_per_A_in_EX'] / NX,
                Gamma_A_GeV=L['Gamma_A_GeV'], T_vis_decay_GeV=(90 / (np.pi**2 * S.GSTAR))**0.25 * np.sqrt(L['Gamma_A_GeV'] * S.MRED))


# ---------------------------------------------------------------- Floquet (transverse modes, late X orbit)
def orbit(E, N, lx=1e-3):
    f = lambda r: E - N**2 / (4 * r * r) - lx * r**4
    rmax = brentq(f, 10.0, 200.0)
    rmin = brentq(f, 0.5, 0.999 * rmax)
    y0 = [rmax, 0.0, 0.0, -N / (2 * rmax)]

    def rhs(t, y):
        c = y[0] + 1j * y[1]
        acc = -2 * lx * abs(c)**2 * c
        return [y[2], y[3], acc.real, acc.imag]

    def ev(t, y):
        return y[0] * y[2] + y[1] * y[3]
    ev.direction = -1
    sol = solve_ivp(rhs, (0, 50), y0, method='DOP853', rtol=1e-12, atol=1e-12, events=ev, dense_output=True)
    T = [t for t in sol.t_events[0] if t > 1e-6][0]
    return sol, T, rmax, rmin


def floquet(sol, T, e, k, a_frozen):
    msig2 = 36 * e * e * a_frozen**2 * P['h']**2

    def rhs(t, y):
        c = sol.sol(t); m2 = 2 * e * e * (c[0]**2 + c[1]**2) + msig2; w2 = k * k + m2
        return [y[1], -w2 * y[0], y[3], -w2 * y[2]]
    s = solve_ivp(rhs, (0, T), [1, 0, 0, 1], method='DOP853', rtol=1e-11, atol=1e-13)
    M = np.array([[s.y[0, -1], s.y[2, -1]], [s.y[1, -1], s.y[3, -1]]])
    return float(np.log(max(abs(np.linalg.eigvals(M)))) / T)


def floquet_job(args):
    e, a_frozen, sol_pack = args
    sol, T = sol_pack
    omega = 2 * np.pi / T
    mA_max = e * np.sqrt(2 * sol.sol(0)[0]**2 + 36 * a_frozen**2 * P['h']**2)
    ks = np.linspace(0.0, 3.0 * max(omega, mA_max), 181)
    mus = [floquet(sol, T, e, k, a_frozen) for k in ks]
    i = int(np.argmax(mus))
    return dict(e=e, a_frozen=a_frozen, mu_max_per_unit_a=float(mus[i]), k_at_max=float(ks[i]),
                a_mA_min_over_omega=float(e * np.sqrt(36 * a_frozen**2 * P['h']**2) / omega),
                growth_exponent_over_10_units_a=float(2 * mus[i] * 10))


def main(out):
    t0 = time.time()
    res = dict(method=__doc__, point=dict(H_GeV=H, F_inf_over_H=float(np.sqrt(F2_INF)), F_vac_over_H=6 * P['h'],
                                           frozen_angle_factor=float(36 * P['h']**2 / F2_INF), EX=EX, NX=NX, R_L=RL))
    with Pool(4) as pool:
        specs = []
        for e in (1e-7,):
            for prof in ('P0', 'P1'):
                specs.append(spectrum(e, prof, pool)); print(prof, e, specs[-1]['N_A'], specs[-1]['peak_k'], time.time() - t0, flush=True)
        for e in (1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3, 1e-2, 3e-2):
            specs.append(spectrum(e, 'P1', pool)); print('P1', e, specs[-1]['N_A'], specs[-1]['peak_k'], time.time() - t0, flush=True)
        res['spectra'] = specs
        rows = []
        for sp in specs:
            g = gmr_formula(sp['e'])
            rows.append(dict(e=sp['e'], profile=sp['profile'], N_A=sp['N_A'], N_A_GMR_formula=g['N_A'], ratio_to_GMR=sp['N_A'] / g['N_A'],
                             m_inf_over_H=sp['m_inf_over_H'], a_star=sp['a_star'], **ledger_e(sp['e'], sp['N_A'])))
        res['ledger_rows'] = rows
        for r in rows:
            print('e=%.0e %s N_A=%.4g (GMR %.4g, ratio %.3f) pairs/net=%.3e allow=%.3e R_A=%.2f R_tot=%.2f (%.3f R_L)'
                  % (r['e'], r['profile'], r['N_A'], r['N_A_GMR_formula'], r['ratio_to_GMR'], r['pairs_per_net_charge'],
                     r['allowance_pairs_per_net'], r['R_A'], r['R_total'], r['total_over_RL']), flush=True)
        # Floquet
        sol, T, rmax, rmin = orbit(EX, NX)
        res['orbit'] = dict(E_X=EX, N_X=NX, r_max=rmax, r_min=rmin, radial_period=T, omega_radial=2 * np.pi / T,
                            X_part_of_mA2_max_over_Sigma_part_at_a1=float(2 * rmax**2 / (36 * P['h']**2)))
        fj = [(e, a, (sol, T)) for e in (1e-7, 1e-3, 3e-3, 5e-3, 1e-2, 3e-2, 1e-1) for a in (1.0, 2.0, 5.0)]
        fl = [floquet_job(a) for a in fj]
        res['floquet'] = fl
        for f in fl:
            print('floquet', f, flush=True)
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print('done', time.time() - t0)


if __name__ == '__main__':
    main(sys.argv[1])
