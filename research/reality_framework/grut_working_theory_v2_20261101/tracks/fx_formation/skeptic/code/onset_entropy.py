"""fx_formation SKEPTIC: how far can the massless 3-dof entropy at onset be off?  Fit an effective mass to each thermal
species from the fx_chain Hartree numbers at onset (n_phi, n_Xbar, net thermal X charge; heat_onset.json, read only),
then compare the massive Bose entropy with the massless 3-dof value used for Y_h.  [E]"""
import json, math, os, sys
from scipy.optimize import brentq, fsolve
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bessel_vapour as bv
Z3 = 1.2020569031595942
HO = json.load(open(os.path.join(HERE, '..', '..', '..', 'fx_chain', 'results', 'heat_onset.json')))
out = {}
for c in HO['chain']:
    if not c['phi_in_equilibrium']:
        continue
    o = c['onset_untilted']
    T = o['T_h_over_mX']; nQ = o['nQ_over_mX3']
    nth = nQ * (1 - c['condensate_charge_fraction'])
    nXb = c['Xbar_per_net_massless'] * nQ; nphi = c['phi_per_net'] * nQ
    n0 = Z3 / math.pi ** 2 * T ** 3
    # phi: mass from n_phi (mu = 0); if n_phi >= n0 take massless
    if nphi < n0:
        mphi = brentq(lambda m: bv.bose(m, 0.0, T)[0] - nphi, 1e-6, 20 * T)
    else:
        mphi = 1e-6
    # X, Xbar: Bose-condensed thermal gas, mu = m (BEC) -> fit m from n_Xbar with mu = m, check net charge
    def f(m):
        return bv.bose(m, -m * 0.999999, T)[0] - nXb
    try:
        mX = brentq(f, 1e-4, 10 * T)
    except ValueError:
        mX = None
    s0 = 3 * 2 * math.pi ** 2 / 45 * T ** 3
    if mX:
        sx = bv.bose(mX, mX * 0.999999, T)[3]; sb = bv.bose(mX, -mX * 0.999999, T)[3]; nxt = bv.bose(mX, mX * 0.999999, T)[0]
    else:
        sx = sb = nxt = float('nan')
    sf = bv.bose(mphi, 0.0, T)[3]
    s = sx + sb + sf
    out[c['label']] = dict(T=T, nQ=nQ, m_phi_eff=mphi, m_X_eff=mX, net_thermal_charge_model=nxt - nXb, net_thermal_charge_fxchain=nth,
                           s_massive_over_s_massless=s / s0, Y_massless=nQ / s0, Y_massive=nQ / s, S_per_Q_massive=s / nQ)
    print('%-40s T=%.3f m_phi=%.3f m_X=%s s/s0=%.4f S/Q: massless %.3f -> massive %.3f (net thermal charge model %.1f vs %.1f)' % (
        c['label'], T, mphi, '%.3f' % mX if mX else '-', s / s0, s0 / nQ, s / nQ, nxt - nXb, nth))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
