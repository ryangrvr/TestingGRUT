#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): gravitational production of the longitudinal vector with an evolving mass.

Single-Higgs-direction longitudinal Proca mode (unitary-gauge reduction; Graham-Mardon-Rajendran 2016 action [R]):
  S_k = 1/2 int dtau [ z^2 |B'|^2 - a^2 m(tau)^2 |B|^2 ],   z^2 = a^2 m^2 / (k^2 + a^2 m^2),
in conformal time, comoving k, a = 1 at the end of inflation (= Codex's release), H_i = 1 units.
  Inflation (tau < 1): a = 1/(2 - tau), H = 1.     Radiation (tau >= 1): a = tau, H = 1/a^2.
Canonical pair (B, Pi = z^2 B'):  B' = Pi/z^2,  Pi' = -a^2 m^2 B.  Pi and B are continuous across mass jumps.
Initial state: instantaneous adiabatic vacuum at k/(aH) = 1000 during inflation.
Late number per mode: n_k + 1/2 = [ |Pi|^2/(2 z^2) + a^2 m^2 |B|^2 / 2 ] / omega,  omega^2 = k^2 + a^2 m^2.
Comoving number N_A = int dlnk k^3 n_k / (2 pi^2), in H_i^3 units -> directly comparable with Codex's comoving N_X.

Mass profiles (m in H_i units; m_inf = e F_g(inflation), m_vac = 6 e h):
  P0  GMR control: m = m_vac throughout (constant vacuum mass, as in our earlier gp_leftovers/gp_critic estimate)
  P1  m^2 = m_vac^2 + (m_inf^2 - m_vac^2)/a^2 after release (X-dominated Higgs norm redshifting as 1/a)
  P2  instant drop to m_vac at release (stress test; not physical: the conformal X amplitude does not vanish
      at release, so the X part of the Higgs norm cannot be removed in a step)
  P4  release-measured window-averaged F_g^2(a) from our homogeneous run (vrel output), X part extrapolated ~a^-2
  P5  m = m_inf throughout (GMR with the inflationary mass; for orientation only)
Usage: python3 proca_modes.py ../results/background.json ../results/proca_modes.json [vrel_F2.json]
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def make_profile(name, m_inf, m_vac, F2tab=None, e=1e-7):
    if name == 'P0':
        return (lambda a: m_vac), m_vac, []
    if name == 'P5':
        return (lambda a: m_inf), m_inf, []
    if name == 'P1':
        return (lambda a: np.sqrt(m_vac**2 + (m_inf**2 - m_vac**2) / a**2)), m_inf, []
    if name == 'P2':
        return (lambda a: m_vac), m_inf, []
    if name == 'P4':
        la, lF = np.log(F2tab['a']), np.log(F2tab['F2'])
        F2vac = F2tab['F2vac']; aend = F2tab['a'][-1]; Xend = F2tab['F2X_end']

        def m(a):
            if a <= aend:
                return e * np.sqrt(np.exp(np.interp(np.log(a), la, lF)))
            return e * np.sqrt(F2vac + Xend * (aend / a)**2)
        return m, m_inf, []
    raise ValueError(name)


def evolve_mode(k, mfun, m_inf_infl, tau_end, breaks):
    """Return n_k evaluated at the requested late conformal times."""
    a0 = k / 1000.0
    tau0 = 2.0 - 1.0 / a0

    def a_of(t):
        return 1.0 / (2.0 - t) if t < 1.0 else t

    def mass(t):
        return m_inf_infl if t < 1.0 else mfun(a_of(t))

    def rhs(t, y):
        a = a_of(t); m = mass(t)
        am2 = (a * m)**2
        z2 = am2 / (k * k + am2)
        return [y[2] / z2, y[3] / z2, -am2 * y[0], -am2 * y[1]]

    a = a0; m = m_inf_infl; am2 = (a * m)**2
    z2 = am2 / (k * k + am2); om = np.sqrt(k * k + am2)
    B = 1.0 / np.sqrt(2 * z2 * om); P = -1j * np.sqrt(z2 * om / 2)
    y = np.array([B.real, B.imag, P.real, P.imag])
    scale = np.array([abs(B), abs(B), abs(P), abs(P)])
    segs = [tau0, 1.0] + [b for b in breaks if 1.0 < b < tau_end[-1]] + list(tau_end)
    segs = sorted(set(segs))
    out = {}
    for t1, t2 in zip(segs[:-1], segs[1:]):
        sol = solve_ivp(rhs, (t1, t2), y / 1.0, method='DOP853', rtol=1e-10,
                        atol=1e-13 * np.maximum(np.abs(y), scale * 1e-6))
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        if t2 in tau_end:
            a = a_of(t2); m = mass(t2); am2 = (a * m)**2
            z2 = am2 / (k * k + am2); om = np.sqrt(k * k + am2)
            Bc = y[0] + 1j * y[1]; Pc = y[2] + 1j * y[3]
            Hk = abs(Pc)**2 / (2 * z2) + am2 * abs(Bc)**2 / 2
            wr = (Bc * np.conj(Pc) - np.conj(Bc) * Pc).imag  # Wronskian, should be 1
            # roundoff-limited: report the error relative to 2|B||Pi|
            out[t2] = dict(n=Hk / om - 0.5, wronskian_err=abs(wr - 1) / max(1.0, 2 * abs(Bc) * abs(Pc)))
    return out


def run_profile(name, bg, F2tab, kgrid):
    m_inf = bg['m_A_inf_over_H']; m_vac = bg['m_A_vac_over_H']; a_star = bg['a_star']
    mfun, m_infl, breaks = make_profile(name, m_inf, m_vac, F2tab)
    if name == 'P5':
        a_star_eff = 1 / np.sqrt(m_inf)
    else:
        a_star_eff = a_star
    rows = []
    for k in kgrid:
        if k < 1e-2 or name == 'P5':
            ends = [20 * a_star_eff, 30 * a_star_eff]
        else:
            ends = [700.0 / k, 1000.0 / k]
        ends = [max(e, 1.0 + 1e-9) for e in ends]
        r = evolve_mode(k, mfun, m_infl, ends, breaks)
        n1, n2 = r[ends[0]]['n'], r[ends[1]]['n']
        rows.append(dict(k=k, n_k=n2, n_k_prev=n1, rel_change=abs(n2 - n1) / max(abs(n2), 1e-300),
                         wronskian_err=r[ends[1]]['wronskian_err']))
    lk = np.log(kgrid)
    integrand = np.array([kk**3 * r['n_k'] for kk, r in zip(kgrid, rows)]) / (2 * np.pi**2)
    NA = float(np.trapezoid(integrand, lk))
    ipk = int(np.argmax(integrand))
    cuts = {}
    for acut in (1000.0, 3000.0):  # modes still superhorizon when Sigma dominates the Higgs norm
        sel = kgrid <= 1.0 / acut
        cuts[f'N_A_k_below_1_over_{int(acut)}'] = float(np.trapezoid(integrand[sel], lk[sel]))
    return dict(profile=name, N_A=NA, peak_k=float(kgrid[ipk]), modes=rows, **cuts,
                max_rel_change=float(max(r['rel_change'] for r in rows)),
                max_wronskian_err=float(max(r['wronskian_err'] for r in rows)))


def physical(NA, bg, m_late):
    HI = bg['H_i_GeV']; MRED = 2.435e18; g = 106.75
    TR = (90 * MRED**2 * HI**2 / (np.pi**2 * g))**0.25
    s_i = 2 * np.pi**2 / 45 * g * TR**3
    nA_s = NA * HI**3 / s_i
    rho_s = nA_s * m_late * HI
    RHO_DM_S = 4.37e-10  # GeV per entropy (Omega_c h^2 = 0.12) [R], as in gp_leftovers
    return dict(T_R_GeV=TR, nA_over_s=nA_s, rhoA_over_s_GeV=rho_s, Omega_A_over_DM=rho_s / RHO_DM_S)


def main(bgfile, out, f2file=None):
    bg = json.loads(Path(bgfile).read_text())
    F2tab = json.loads(Path(f2file).read_text()) if f2file else None
    kgrid = np.geomspace(1e-7, 1e-1, 49)
    res = dict(background=dict((k, bg[k]) for k in ('H_i_GeV', 'm_A_inf_over_H', 'm_A_vac_over_H', 'a_star',
                                                    'k_star_over_Hi', 'mass_ratio_sq_vac_over_inf')),
               method='longitudinal Proca mode, adiabatic-vacuum start at k/aH=1000, DOP853 rtol 1e-10',
               profiles={})
    profiles = ['P0', 'P1', 'P2', 'P5'] + (['P4'] if F2tab else [])
    for p in profiles:
        r = run_profile(p, bg, F2tab, kgrid)
        m_late = bg['m_A_inf_over_H'] if p == 'P5' else bg['m_A_vac_over_H']
        r['physical'] = physical(r['N_A'], bg, m_late)
        r['pairs_per_net_charge_actual'] = r['N_A'] / abs(-238.64868)
        r['pairs_per_net_charge_YQnorm'] = r['physical']['nA_over_s'] / 4.8294009076988175e-9
        res['profiles'][p] = r
        print(p, r['N_A'], r['peak_k'], r['max_rel_change'], r['max_wronskian_err'], r['physical'], flush=True)
    P0 = res['profiles']['P0']['N_A']
    for p in profiles:
        res['profiles'][p]['ratio_to_P0'] = res['profiles'][p]['N_A'] / P0
    # GMR analytic formula for comparison [R]: Omega_A/Omega_DM = sqrt(m/6e-6 eV) (H_I/1e14 GeV)^2
    m_vac_GeV = bg['m_A_vac_over_H'] * bg['H_i_GeV']
    res['GMR_formula_Omega_over_DM'] = float(np.sqrt(m_vac_GeV / 6e-15) * (bg['H_i_GeV'] / 1e14)**2)
    res['P0_over_GMR_formula'] = res['profiles']['P0']['physical']['Omega_A_over_DM'] / res['GMR_formula_Omega_over_DM']
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({p: dict(N_A=res['profiles'][p]['N_A'], ratio_to_P0=res['profiles'][p]['ratio_to_P0'],
                              pairs_actual=res['profiles'][p]['pairs_per_net_charge_actual'],
                              pairs_YQ=res['profiles'][p]['pairs_per_net_charge_YQnorm'])
                      for p in profiles}, indent=1))
    print('GMR formula', res['GMR_formula_Omega_over_DM'], 'P0/formula', res['P0_over_GMR_formula'])


if __name__ == '__main__':
    main(*sys.argv[1:])
