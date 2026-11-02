"""Controls for vapour_transfer.py (own): (1) a cold gas through the Gilbert branch reproduces CDM; (2) step halving;
(3) marching vs pure ODE; (4) a Maxwellian gas with the 1D velocity dispersion of a 3 keV thermal relic (f = 1) against
the Viel et al. 2005 fit (astro-ph/0501562) [R]; (5) the thermal-WDM curves of the Lyman-alpha and satellite bounds
(5.3 keV Irsic et al. 2017 arXiv:1702.01764; 6.5 keV Nadler et al. 2021 arXiv:2008.00022) at k = 10, 20, 50 /Mpc."""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vapour_transfer as vt
RES = os.path.join(HERE, '..', 'results')


def viel_T2(m_keV, k):
    al = 0.049 * m_keV**(-1.11) * (0.27 / 0.25)**0.11 * (0.674 / 0.7)**1.22 / 0.674
    nu = 1.12
    return ((1 + (al * k)**(2 * nu))**(-5 / nu))**2


out = {}
fake = dict(a_m=1e-6, f_v_mass=0.5, sigma_p_today=0.0, _ceff2=(np.array([1e-12, 1e-6]), np.array([0.0, 0.0])))
c = {}
for k in (3.0, 100.0):
    r1 = vt.growth(k, fake) / vt.growth(k, fake, reference=True)
    r2 = vt.growth(k, fake, dN=0.002) / vt.growth(k, fake, dN=0.002, reference=True)
    f2 = dict(fake); f2['a_m'] = 0.999999
    rr = vt.growth(k, fake, reference=True) / vt.growth(k, f2, reference=True)
    c['k%g' % k] = dict(cold_gas_ratio=r1, cold_gas_ratio_halfstep=r2, marching_vs_ODE=rr)
out['cold_gas_and_step'] = c
sp_ = 7.8e-8 * (1 / 3.)**(4 / 3.)
w = dict(a_m=1e-7, f_v_mass=0.999, sigma_p_today=sp_, _ceff2=(np.array([1e-12, 1e-7]), np.array([0.0, 0.0])))
wd = {}
for k in (5., 10., 20., 50.):
    wd['k%g' % k] = dict(gilbert=(vt.growth(k, w) / vt.growth(k, w, reference=True))**2,
                         fluid=(vt.growth(k, w, 'fluid') / vt.growth(k, w, 'fluid', reference=True))**2, viel_3keV=viel_T2(3.0, k))
out['maxwellian_3keV_vs_viel'] = wd
out['WDM_bound_curves_P'] = {('m%gkeV' % m): {('k%g' % k): viel_T2(m, k) for k in (10, 20, 50)} for m in (5.3, 6.5, 3.0, 1.0)}
json.dump(out, open(os.path.join(RES, 'controls_vt.json'), 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
