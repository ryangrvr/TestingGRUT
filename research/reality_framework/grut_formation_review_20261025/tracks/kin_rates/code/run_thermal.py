#!/usr/bin/env python3
"""kin_rates: in-medium masses, the phase, and the in-medium s width.  Writes results/thermal_masses.json.
usage: python3 -B run_thermal.py <out.json>
Two schemes: one-loop Dolan-Jackiw (primary) and the Parwani daisy resummation (bracket).  For each T: v_T, m_X(T),
m_s(T), lam3(T), lam_s3(T), m_s/(2 m_X), the rest-frame in-medium width Gamma_s(T) = Gamma_vac(m_s(T), m_X(T),
lam3(T)) [1 + 2 f_B(m_s/2)] (mu = 0) and the Bose factor; T_c, the temperature where s -> X Xb closes, and the
deviations of the vacuum masses (1 %, 5 %, 10 % crossings)."""
import json
import sys
import numpy as np
from scipy import optimize
import model as M
import thermal as TH


def row(T, scheme):
    d = TH.thermal_masses(T) if scheme == 'oneloop' else TH.thermal_masses_daisy(T)
    out = dict(T_over_mX=T / M.MX, T_GeV=T, vT_over_v=d['vT'] / M.V, mX_GeV=d['mX'], ms_GeV=d['ms'],
               mX_over_mX0=d['mX'] / M.MX, ms_over_ms0=d['ms'] / M.MS, lam3_GeV=d['lam3'], lam_s3_GeV=d['lam_s3'],
               broken=bool(d['broken']), ms_over_2mX=(d['ms'] / (2 * d['mX']) if d['mX'] > 0 else None),
               PiX_over_mX0sq=d.get('PiX', 0.0) / M.MX**2, ir_s_over_ms0sq=d.get('ir_s', 0.0) / M.MS**2)
    open_ = d['ms'] > 2 * d['mX'] and d['broken']
    if open_:
        g_rest, enh = TH.width_s_medium(d['ms'], T, ms=d['ms'], mx=d['mX'], lam3=d['lam3'])
        out.update(Gamma_s_vacformula_medium_masses_GeV=TH.width_s_vac(d['ms'], d['mX'], d['lam3']),
                   Gamma_s_rest_medium_GeV=g_rest, bose_enh_rest=enh, decay_open=True)
    else:
        out.update(Gamma_s_vacformula_medium_masses_GeV=0.0, Gamma_s_rest_medium_GeV=0.0, bose_enh_rest=None, decay_open=False)
    return out


def crossing(f, lo, hi, n=200):
    xs = np.linspace(lo, hi, n)
    v = [f(x) for x in xs]
    for i in range(n - 1):
        if v[i] * v[i + 1] < 0:
            return float(optimize.brentq(f, xs[i], xs[i + 1], xtol=1e-7))
    return None


def main(out):
    Tgrid = sorted(set(np.round(np.concatenate([np.geomspace(0.05, 30, 49), [0.3, 1, 3, 6, 7.344, 8.118, 10, 12, 13, 14, 15, 16, 17, 18, 18.1, 20, 25]]), 6)))
    res = dict(scope='One-loop thermal masses (Dolan-Jackiw) with the X doublet and the radial s in the loop; daisy (Parwani) '
                     'resummation as a bracket.  Zero-T CW piece absorbed in the vacuum parameters (relative size g^4/16pi^2 ~ 6e-6).  '
                     'mu = 0 for the loop particles (chemical potentials of the gas would add to the tadpoles; stated as a limit).')
    for scheme in ('oneloop', 'daisy'):
        rows = [row(t * M.MX, scheme) for t in Tgrid]
        fm = (lambda T: TH.thermal_masses(T)) if scheme == 'oneloop' else (lambda T: TH.thermal_masses_daisy(T))
        Tc = TH.Tc_oneloop()['T_minimum_disappears'] if scheme == 'oneloop' else TH.Tc_daisy()
        close = crossing(lambda T: (lambda d: d['ms'] - 2 * d['mX'])(fm(T)), 5 * M.MX, 0.999 * Tc)
        dev = {}
        for sp, key in (('mX', M.MX), ('ms', M.MS)):
            for frac in (0.01, 0.05, 0.10):
                dev[f'{sp}_dev_{int(frac*100)}pct_T_over_mX'] = (lambda c: c / M.MX if c else None)(
                    crossing(lambda T: abs(fm(T)[sp] / key - 1) - frac, 0.2 * M.MX, 0.999 * Tc))
        res[scheme] = dict(rows=rows, T_c_GeV=Tc, T_c_over_mX=Tc / M.MX, T_close_s_decay_GeV=close,
                           T_close_s_decay_over_mX=(close / M.MX if close else None), vacuum_mass_deviation_crossings=dev)
        print(scheme, 'Tc/m =', Tc / M.MX, ' s->XXb closes at T/m =', close / M.MX if close else None, dev)
    res['oneloop'].update(TH.Tc_oneloop())
    res['hiT_estimate'] = dict(T_c_over_mX=np.sqrt(M.LP * M.V**2 / (M.LP / 4 + M.GX**2 / 6)) / M.MX,
                               formula='v_T^2 = v^2 - T^2 (lP/4 + g^2/6)/lP;  m_X^2 = g^2 v_T^2 + T^2 (lX/3 + g^2/12)')
    # Bose factor of the in-medium width as a function of the s energy at a few T (mu = 0), one-loop masses
    tab = []
    for t in (0.3, 1.0, 3.0, 7.344, 12.0):
        T = t * M.MX
        d = TH.thermal_masses(T)
        if d['ms'] <= 2 * d['mX']:
            continue
        for Eom in (1.0, 1.2, 1.5, 2, 3, 5, 10, 30, 100):
            E = Eom * d['ms']
            g, enh = TH.width_s_medium(E, T, ms=d['ms'], mx=d['mX'], lam3=d['lam3'])
            tab.append(dict(T_over_mX=t, E_over_ms=Eom, Gamma_lab_GeV=g, bose_enh=enh, Gamma_lab_over_vac=g / M.GAMMA_S_VAC))
    res['width_vs_energy_oneloop'] = tab
    json.dump(res, open(out, 'w'), indent=1, default=float)
    for r in res['oneloop']['rows']:
        if r['T_over_mX'] in (0.3, 1, 3, 6, 7.344, 8.118, 10, 12, 13, 14, 15, 17, 18, 20):
            print(f"T/m={r['T_over_mX']:6.3f} vT/v={r['vT_over_v']:.4f} mX/mX0={r['mX_over_mX0']:.4f} ms/ms0={r['ms_over_ms0']:.4f} ms/2mX={r['ms_over_2mX']} "
                  f"Gamma_rest/vac={r['Gamma_s_rest_medium_GeV']/M.GAMMA_S_VAC:.4f} enh={r['bose_enh_rest']}")


if __name__ == '__main__':
    main(sys.argv[1])
