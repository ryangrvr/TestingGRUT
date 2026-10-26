#!/usr/bin/env python3
"""kin_thermo (Claude, own code): formation alternatives as callable rules, with what each one supplies.

(a) rule_gas_only      : no objects form.  Charge beyond the Bose capacity sits in a homogeneous zero-momentum
                         condensate at mu = m_X, energy m_X per charge (gp_leftovers 'case H' comparator).
(b) rule_instant       : the old proxies' rule.  Whenever the gas would need mu_+ > omega(Q_ref) the excess charge
                         is condensed into objects of charge Q_ref instantly (and dissolves instantly when
                         undersaturated); energy per condensed charge e_cap = omega (consistent linear comparator)
                         or E/Q (the benchmark's booking).  Variant 'finite_Q' books E(Q_ref) per object, i.e. adds
                         the object-number work (E - Q omega) dn_L on top of omega d(n_L Q).
(c) rule_fragmentation : nucleation/fragmentation bracket from the published Q-ball formation literature, scaled
                         to this model.  Credits (read level in brackets): Kusenko-Shaposhnikov PLB 418 (1998) 46
                         [abstract read; linear fragmentation analysis recalled]; Enqvist-McDonald PLB 425 (1998)
                         309 [abstract]; Kasuya-Kawasaki PRD 61 (2000) 041301 and PRD 62 (2000) 023512 [abstracts:
                         'the typical size is set by the most rapidly growing linearised mode; nearly all charge
                         ends in Q-balls; charge depends linearly on the initial charge density']; Hiramatsu-
                         Kawasaki-Takahashi JCAP 06 (2010) 008 [abstract: 3D peak charge 60% above the earlier
                         estimate]; Laine-Shaposhnikov NPB 532 (1998) 376 [abstract; solitosynthesis context];
                         Griest-Kolb 1989 and Postma 2002 solitosynthesis [recalled].
    Model-specific derivation [D]: with phi integrated out adiabatically (phi^2 = v^2 - 2 g^2 |X|^2/lambda_phi >= 0),
        U_eff(|X|) = m_X^2 |X|^2 - lambda_eff |X|^4,   lambda_eff = g^4/lambda_phi - lambda_X = 2.3333e-3 > 0,
    for |X| <= x_c, and U = lambda_X |X|^4 + V0 (phi = 0) beyond.  U_eff is flatter than quadratic, which is exactly
    why the two-field Q-balls exist; the homogeneous condensate X = X0 e^{i omega t} is then linearly unstable
    (own re-derivation, NOTES.md): for modes e^{alpha t - i k x},
        (alpha^2 + k^2)(alpha^2 + k^2 + U'' - omega^2) + 4 omega^2 alpha^2 = 0,
        unstable band 0 < k^2 < k_max^2 = omega^2 - U''(phi0) = 4 lambda_eff |X0|^2,   omega^2 = m^2 - 2 lambda_eff |X0|^2,
    (phi0 = sqrt2 |X0|; U'' is the second radial derivative).  alpha_max is found numerically.  Lump charge
        Q_lump = shape * n_cond * (2 pi/k_max)^3,  shape in [0.1, 1] (literature O(1) systematics),
    n_cond = 2 omega |X0|^2 the condensate charge density.  Fragmentation time t_frag = N_e/alpha_max with
    N_e in [20, 40] e-folds of linear growth from small seeds [E].  Constraints supplied: stationary objects exist
    only for Q >= Q_min (objects_EQ.json: 5.2e4) and only in the broken phase T < T_c (one-loop high-T estimate
    T_c = sqrt(lambda_phi v^2/(lambda_phi/4 + g^2/6)) = 18.1 m_X [E]; the exact one-loop value is kin_rates' job);
    anti-objects: for an elliptic orbit (epsilon < 1) the literature reports Q-balls of both signs [R, Kasuya-
    Kawasaki 2000]; the fraction f_anti of the symmetric number locked in anti-objects is a bracket [0, 1].
    Nucleation from a supersaturated gas: barrier F_crit(mu) = max_Q [E(Q) - mu Q] = (E - Q omega) on the
    thick-wall branch at omega = mu; exponent F_crit/T from the table.
Nothing here is a measured formation rate.  Every rule returns the parameters it assumed.
"""
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

import kt_common as C
import thermo as TH

RES = Path(__file__).resolve().parents[1] / 'results'
MX, G, LX, LP, V = C.MX, C.GX, C.LX, C.LP, C.V
LAM_EFF = G**4 / LP - LX
TC_HIGH_T = np.sqrt(LP * V * V / (LP / 4 + G * G / 6))


def _table():
    d = json.loads((RES / 'objects_EQ.json').read_text())
    rows = d['rows']
    Om = np.array([r['Omega'] for r in rows]); Q = np.array([r['Q'] for r in rows]); E = np.array([r['E_GeV'] for r in rows])
    return d, Om, Q, E


def U_eff(absX):
    x = np.asarray(absX, float)
    inside = x <= C.XC
    return np.where(inside, MX**2 * x * x - LAM_EFF * x**4, LX * x**4 + C.V0)


def dispersion_alpha_max(absX0):
    """Fastest linear growth rate (GeV) and its k (GeV) for the circular condensate of amplitude |X0| <= x_c."""
    x0 = float(absX0)
    om2 = MX**2 - 2 * LAM_EFF * x0 * x0
    Upp = MX**2 - 6 * LAM_EFF * x0 * x0           # U''(phi) with phi = sqrt2 |X|
    kmax2 = om2 - Upp
    if kmax2 <= 0:
        return dict(k_max=0.0, alpha_max=0.0, k_at_max=0.0, omega=np.sqrt(om2))

    def alpha2(k):
        b = 2 * k * k + 3 * om2 + Upp
        c = k * k * (k * k - kmax2)
        return (-b + np.sqrt(b * b - 4 * c)) / 2
    res = minimize_scalar(lambda k: -alpha2(k), bounds=(0, np.sqrt(kmax2)), method='bounded', options=dict(xatol=1e-10))
    return dict(k_max=float(np.sqrt(kmax2)), alpha_max=float(np.sqrt(-res.fun)), k_at_max=float(res.x),
                omega=float(np.sqrt(om2)), n_cond=float(2 * np.sqrt(om2) * x0 * x0))


# ----------------------------------------------------------------------------------------------- rules
def rule_gas_only(**kw):
    return dict(rule='gas_only', objects_form=False, u_cap=1.0, e_cap=1.0,
                supplies='homogeneous condensate at mu = m_X when the gas saturates; no finite-Q objects',
                anti_objects=False)


def rule_instant(t_hint, nq, pinv, rho, variant='omega', Q_ref=C.Q_REF, with_s=True, antilumps=True):
    """Old-proxy rule as an algebraic constraint on the two-potential state (dimensionless inputs, m_X units).
    Returns the state with condensed charge n_c (objects) and n_ac (anti-objects) and the object energy booking."""
    d, Om, Q, E = _table()
    ref = d['reference_Q200000']
    e_cap = {'omega': C.OMEGA_REF, 'E_over_Q': ref['E_over_QmX'], 'finite_Q': C.OMEGA_REF}[variant]
    st = TH.state_two_potential(nq, pinv, rho, C.OMEGA_REF, e_cap=e_cap, anti_cap=C.OMEGA_REF if antilumps else 0.999999,
                                with_s=with_s, tguess=t_hint)
    n_L = st['n_cond'] / Q_ref
    work = (ref['E_GeV'] - Q_ref * C.OMEGA_REF * MX) * n_L * MX**3 if variant == 'finite_Q' else 0.0   # GeV^4
    return dict(rule='instant_condensation', variant=variant, Q_ref=Q_ref, u_cap=C.OMEGA_REF, e_cap=e_cap,
                n_cond=st['n_cond'], n_anticond=st['n_anticond'], n_objects_m3=n_L, n_antiobjects_m3=st['n_anticond'] / Q_ref,
                object_number_work_GeV4=work, state=st,
                supplies='instant equilibrium condensation at mu_+ = omega(Q_ref); no formation time, no Q distribution')


def rule_fragmentation(n_cond_GeV3=None, absX0_GeV=None, shape=0.3, Ne=30.0, f_cap=0.9, f_anti=None, T_GeV=None):
    """Fragmentation of a homogeneous condensate (either the release at |X0| or a gas condensate of density n_cond).
    Returns Q_lump, time scale, object density and the flags."""
    if absX0_GeV is None:
        # circular orbit: n = 2 omega |X|^2 with omega^2 = m^2 - 2 lam |X|^2 -> solve for |X|
        from scipy.optimize import brentq
        f = lambda x: 2 * np.sqrt(max(MX**2 - 2 * LAM_EFF * x * x, 1e-300)) * x * x - n_cond_GeV3
        absX0_GeV = brentq(f, 1e-12, C.XC, xtol=1e-16)
    disp = dispersion_alpha_max(absX0_GeV)
    n_cond = disp['n_cond']
    lam = 2 * np.pi / disp['k_max'] if disp['k_max'] > 0 else np.inf
    Q_lump = shape * n_cond * lam**3
    d, Om, Q, E = _table()
    Qmin = d['Q_min']['Q']
    t_frag = Ne / disp['alpha_max'] if disp['alpha_max'] > 0 else np.inf
    flags = []
    if Q_lump < Qmin:
        flags.append('Q_lump below Q_min = %.3g: fragments are not stationary objects; they must merge or dissolve' % Qmin)
    if T_GeV is not None and T_GeV > TC_HIGH_T:
        flags.append('T above the one-loop high-T T_c = %.3f GeV: no broken phase, no objects' % TC_HIGH_T)
    sE = CubicSpline(Q[np.argsort(Q)], E[np.argsort(Q)]) if False else None
    # energy per lump from the stable branch (interpolate in Q on the stable branch only)
    stable = np.array([r['dQ_dOmega'] < 0 for r in d['rows']])
    Qs, Es, Os = Q[stable], E[stable], Om[stable]
    o = np.argsort(Qs)
    E_lump = float(np.interp(Q_lump, Qs[o], Es[o])) if Qs.min() <= Q_lump <= Qs.max() else None
    om_lump = float(np.interp(Q_lump, Qs[o], Os[o])) if E_lump is not None else None
    return dict(rule='fragmentation', absX0_GeV=float(absX0_GeV), n_cond_GeV3=n_cond, lambda_eff=LAM_EFF,
                k_max_GeV=disp['k_max'], k_max_over_mX=disp['k_max'] / MX, alpha_max_GeV=disp['alpha_max'],
                alpha_max_over_mX=disp['alpha_max'] / MX, k_at_alpha_max_over_mX=disp['k_at_max'] / MX,
                omega_condensate_over_mX=disp['omega'] / MX, wavelength_GeV_inv=lam, shape=shape, Ne=Ne,
                Q_lump=Q_lump, Q_min_table=Qmin, t_frag_GeV_inv=t_frag, f_cap=f_cap,
                n_objects_GeV3=f_cap * n_cond / Q_lump if Q_lump > 0 else 0.0,
                E_lump_GeV=E_lump, omega_lump_over_mX=om_lump,
                object_number_work_per_lump_GeV=(E_lump - Q_lump * om_lump * MX) if E_lump else None,
                f_anti_bracket=[0.0, 1.0] if f_anti is None else f_anti, flags=flags,
                supplies='k_max, alpha_max, Q_lump(shape), t_frag(Ne), f_cap, f_anti bracket; no lattice run')


def nucleation_barrier():
    """F_crit(mu) = E - Q omega on the thick-wall branch (omega = mu), tabulated from objects_EQ.json."""
    d, Om, Q, E = _table()
    rows = [r for r in d['rows'] if r['dQ_dOmega'] > 0]
    return [dict(mu_over_mX=r['Omega'], Q_crit=r['Q'], F_crit_GeV=r['E_minus_Q_omega_GeV'],
                 F_crit_over_T_at_T=dict((str(t), r['E_minus_Q_omega_GeV'] / (t * MX)) for t in (0.3, 1.0, 3.0, 7.34, 12.0)))
            for r in rows]


if __name__ == '__main__':
    print('lambda_eff', LAM_EFF, 'T_c high-T estimate GeV', TC_HIGH_T, TC_HIGH_T / MX, 'm_X')
    print('U_eff continuity at x_c:', U_eff(C.XC * (1 - 1e-9)), U_eff(C.XC * (1 + 1e-9)))
    r = rule_fragmentation(absX0_GeV=C.XC)
    print({k: v for k, v in r.items() if k != 'supplies'})
    r = rule_fragmentation(n_cond_GeV3=0.7 * C.N_Q0)
    print({k: v for k, v in r.items() if k != 'supplies'})
    nb = nucleation_barrier()
    for x in nb[:3] + nb[-3:]:
        print(x)
    st = rule_instant(7.3, C.N_Q0 / MX**3, 0.0685 / MX**3, C.rho_onset(40.5626) / MX**4)
    print({k: v for k, v in st.items() if k != 'state'}, st['state']['t'], st['state']['up'], st['state']['um'])
