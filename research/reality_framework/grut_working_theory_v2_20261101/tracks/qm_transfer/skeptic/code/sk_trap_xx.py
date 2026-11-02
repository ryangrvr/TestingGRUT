"""qm_transfer SKEPTIC: three side estimates (own code).  [C on E/Q inputs]
(a) Interior trapping (tracks/qm_micro section 1f [Q]): an absorbed X becomes one bulk phonon of energy ~Delta + E_k
    (k ~ 0.41); unless it decays (Beliaev, mean free path l ~ 2e7 m_X^-1 at T ~ 0.005-0.01, qm_micro eos_phonon.json)
    before it re-crosses the surface, it is re-emitted as an X of the same energy (a delayed ELASTIC scattering).
    Thermalising capture is s_eff = s_bar(T) P_trap, P_trap = 1 - exp(-8R/3l) (qm_micro's escape-cone estimate), so
    zeta_eff = s_bar P_trap Q^(-1/3) -> 8 s_bar (3/4 pi n_B)^(1/3) / (3 l) for R << l, independent of Q.
(b) Thermostat criterion for the track's 'thermo' re-emission: an object can hold one latent heat below its
    evaporation threshold only if S(Delta) > Delta/T_o (S from phonons + ripplons, as in sk_object_mc.py).
(c) X-X elastic scattering (sigma_XX/m = 6.0e-8 cm^2/g, wt_critic [Q/D]) against vapour-object collisions
    (pi R^2 per object, any sticking): ratio = N_v sigma_XX sqrt(2) / (K Q^(-1/3) C_o^(2/3)).
usage: sk_trap_xx.py OUT.json"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sk_object_mc as OM
TR = os.path.join(HERE, '..', '..')
QM = os.path.join(TR, '..', 'qm_micro', 'results')
NB = 449.31336788
K = math.pi * (3 / (4 * math.pi * NB))**(2 / 3)
eos = json.load(open(os.path.join(QM, 'eos_phonon.json')))
st = json.load(open(os.path.join(QM, 'stick_thermal.json')))
mfp = lambda T: float(np.exp(np.interp(math.log(T), np.log([r['T'] for r in eos['rows']]), np.log([r['mfp'] for r in eos['rows']]))))
sbar = lambda T: float(np.interp(T, [r['T_h'] for r in st['thermal']], [r['s_X_boltz'] for r in st['thermal']]))
out = dict(scope=__doc__.strip().splitlines()[0])
# (a)
rows = []
for T in (0.004, 0.005, 0.009):
    l, s = mfp(T), sbar(T)
    for lQ in range(6, 25, 1):
        Q = 10.0**lQ; R = (3 * Q / (4 * math.pi * NB))**(1 / 3)
        Pt = 1 - math.exp(-8 * R / (3 * l))
        rows.append(dict(T=T, Q=Q, R=R, mfp=l, s_bar=s, P_trap=Pt, zeta_geom=Q**(-1 / 3), zeta_perfect=s * Q**(-1 / 3), zeta_eff=s * Pt * Q**(-1 / 3)))
    zmax = max(r['zeta_eff'] for r in rows if r['T'] == T)
    print('T=%.3f  l=%.2e  s_bar=%.3f  max_Q zeta_eff = %.2e   (R<<l limit %.2e)' % (T, l, s, zmax, 8 * s * (3 / (4 * math.pi * NB))**(1 / 3) / (3 * l)))
out['trapping'] = rows
for r in rows:
    if r['T'] == 0.004 and r['Q'] in (1e8, 1e9, 1e10, 1e12, 1e20):
        print('   Q=%.0e R=%.0f P_trap=%.2e zeta(perfect)=%.2e zeta_eff=%.2e' % (r['Q'], r['R'], r['P_trap'], r['zeta_perfect'], r['zeta_eff']))
# required damping length for zeta_eff = zeta_target at s_bar(0.004): need P_trap = zeta Q^(1/3)/s_bar <= 1
req = []
for zt in (7e-5, 1e-4, 3e-4, 1e-3, 3e-3):
    s = sbar(0.004)
    Qmax = (s / zt)**3                       # largest object that can reach zt even with perfect trapping
    R = (3 * Qmax / (4 * math.pi * NB))**(1 / 3)
    req.append(dict(zeta=zt, Q_max_perfect=Qmax, R=R, l_needed_lt=R, ratio_to_beliaev=mfp(0.004) / R))
    print('zeta=%.0e needs Q <= %.1e (s_bar=%.2f, perfect trapping) and a thermalisation length <~ R = %.0f m_X^-1 = Beliaev/%.1e' % (zt, Qmax, s, R, mfp(0.004) / R))
out['trapping_requirement'] = req
# (b)
crit = []
for D in (0.207, 0.219, 0.224):
    for To in (0.003, 0.0035, 0.004, 0.005, 0.009):
        f = lambda lQ: OM.S_of_T(brentq(lambda lT: math.log(OM.U_of_T(math.exp(lT), 10**lQ)) - math.log(D), -20, 0) and
                                  math.exp(brentq(lambda lT: math.log(OM.U_of_T(math.exp(lT), 10**lQ)) - math.log(D), -20, 0)), 10**lQ) - D / To
        lQc = brentq(f, 4, 20)
        crit.append(dict(Delta=D, T_o=To, Q_c=10**lQc, zeta_c_s1=10**(-lQc / 3)))
print('thermostat threshold Q_c (S(Delta) = Delta/T_o):', [(c['Delta'], c['T_o'], '%.1e' % c['Q_c']) for c in crit])
out['thermostat_threshold'] = crit
# (c)
sig_XX = 6.0e-8 * (5.0677e13**2 / 5.6096e23) * 0.0924089637 * 0.0924089637**2       # cm^2/g -> m_X^-2
G = json.load(open(os.path.join(TR, 'results', 'transfer_grid.json')))['models']
xx = []
for m in G:
    if not m.get('Q'):
        continue
    Nv = m['f_v_charge']
    ratio = Nv * sig_XX * math.sqrt(2) / (K * m['Q']**(-1 / 3) * (1 - Nv)**(2 / 3))
    xx.append(dict(name=m['name'], Q=m['Q'], f_v_charge=Nv, ratio_XX_to_object=ratio))
out['sigma_XX_mX2'] = sig_XX
out['xx'] = xx
for r in xx:
    if r['name'].startswith(('M_', 'P_SQ1.235')):
        print('%-22s Q=%.1e  X-X / object collision rate = %.2e' % (r['name'], r['Q'], r['ratio_XX_to_object']))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
