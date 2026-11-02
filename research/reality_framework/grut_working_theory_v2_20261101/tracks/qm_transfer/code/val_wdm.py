"""qm_transfer validation 2 (own): a pure thermal warm relic (Fermi-Dirac, g = 2 fermion, relativistic at early times)
through the same Boltzmann code (q-binned hierarchy, relativistic eps, background energy e(a) included), against the
Viel et al. 2005 fit T(k) = [1 + (alpha k)^(2 nu)]^(-5/nu), nu = 1.12,
alpha = 0.049 (m/keV)^-1.11 (Omega_x/0.25)^0.11 (h/0.7)^1.22 Mpc/h  (astro-ph/0501562) [R],
with (T_x/T_nu)^3 = Omega_x h^2 (94 eV / m_x) (Bode, Ostriker & Turok 2001; Viel et al. 2005) [R].
Also a Maxwellian relic with the same <q^2> (non-thermal-shape control) and a 1 keV mixed case (f = 0.2).  [C]
usage: val_wdm.py OUT.json
"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from qt_bg import Cosmo
from qt_hidden import Frozen
from qt_boltz import Mode, CDM


def viel_T(m_keV, k, c):
    al = 0.049 * m_keV**(-1.11) * ((c.Odm) / 0.25)**0.11 * (c.h / 0.7)**1.22 / c.h
    return (1 + (al * k)**(2 * 1.12))**(-5 / 1.12)


def Tq_of_m(m_keV, c):
    Tnu = (4 / 11)**(1 / 3) * c.T0 * 8.617333262e-5          # eV
    Tx = Tnu * (c.omc * 94.0 / (m_keV * 1e3))**(1 / 3)
    return Tx / (m_keV * 1e3)


def run_case(m_keV, frac, dist, ks, c0):
    Tq = Tq_of_m(m_keV, c0)
    if dist == 'FD':
        fr = Frozen(c0, fraction=frac, dist='FD', Tq=Tq, Nq=28, Lv=14, a_lo=1e-11)
    else:
        # Maxwellian with the same <q^2> as the FD relic (numerical moment of the FD grid)
        frFD = Frozen(c0, fraction=frac, dist='FD', Tq=Tq, Nq=28, Lv=14, a_lo=1e-11)
        q2 = float(frFD.wq @ (frFD.T['N'][0] * frFD.q**2)) / float(frFD.wq @ frFD.T['N'][0])
        fr = Frozen(c0, fraction=frac, dist='MB', sigma=math.sqrt(q2 / 3), Nq=28, Lv=14, a_lo=1e-11)
    c = Cosmo(h=c0.h, omb=c0.omb, omc=c0.omc, Neff=c0.Neff, dm_efac=fr.efac)
    fr.c = c
    cref = Cosmo(h=c0.h, omb=c0.omb, omc=c0.omc, Neff=c0.Neff)
    rows = []
    for k in ks:
        t0 = time.time()
        d1 = Mode(c, fr, k, K_rsa=200., K_end=400.).run()[1.0]
        d0 = Mode(cref, CDM(), k, K_rsa=200., K_end=400.).run()[1.0]
        P = (d1 / d0)**2
        Pv = viel_T(m_keV, k, c0)**2
        rows.append(dict(k=k, P=P, P_viel=Pv))
        print('m=%g keV f=%g %s k=%g: P=%.4f  Viel=%.4f  (%.0fs)' % (m_keV, frac, dist, k, P, Pv, time.time() - t0), flush=True)
    return dict(m_keV=m_keV, fraction=frac, dist=dist, Tq=Tq, sigma_p=fr.sigma_p, rows=rows)


def main():
    out = sys.argv[1]
    which = sys.argv[2] if len(sys.argv) > 2 else 'all'
    c0 = Cosmo()
    res = json.load(open(out)) if os.path.exists(out) else {}
    cases = {'m3': (3.0, 0.999, 'FD', [2., 5., 10., 15., 20., 30., 50.]),
             'm5.3': (5.3, 0.999, 'FD', [5., 10., 20., 30., 40., 60., 80.]),
             'm1': (1.0, 0.999, 'FD', [1., 2., 4., 6., 10., 15., 25.]),
             'm3MB': (3.0, 0.999, 'MB', [5., 10., 20., 30.]),
             'm1f0.2': (1.0, 0.2, 'FD', [2., 5., 10., 20., 50., 100.])}
    for key in (cases if which == 'all' else which.split(',')):
        res[key] = run_case(*cases[key][:3], cases[key][3], c0)
        json.dump(res, open(out, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
