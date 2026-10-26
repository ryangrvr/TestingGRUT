#!/usr/bin/env python3
"""kin_thermo (Claude, own code): capture / evaporation / conversion kernels of the object sector at the hidden
temperature, as callables for the solver, with their detailed-balance reverses and brackets.

1. Conversion (bb211d1's channel, charge-conjugate pair included):
       Xbar(E_-) + core(Q) -> X(E_+) + core(Q-2),   E_+ = E_- + E(Q) - E(Q-2)  (= E_- + 2 omega m_X + O(1/Q))
   Forward sigma v(p_-) on the cold Q = 200000 ground profile: our conversion track's own three-channel solver,
   dense table of 741 momenta p_-/m_X = 0.001 ... 0.708 (results/rate_scan.json, [Q own]; CONFIRMED there against
   Codex's cold 1.13039e-16 and daughter 8.1225e-15 cm^3/s to <= 4e-6).  Above p_- = 0.708 the table does not exist;
   two labelled extrapolations are carried: 'hold' (last value) and 'geometric' (pi r90^2 v, an upper bracket).
   Reverse (reciprocity, bb211d1 MATHEMATICS.md sec. 2 [Q]):  (sigma v)_rev(p_+) = (sigma v)_fwd(p_-) p_- E_-/(p_+ E_+).
   Rates per volume with Bose factors:
       R_fwd = n_L int d^3p_-/(2pi)^3 f_-(p_-) sv_f(p_-) [1 + f_+(p_+)],
       R_rev = n_L int d^3p_+/(2pi)^3 f_+(p_+) sv_r(p_+) [1 + f_-(p_-)],
   which balance exactly when mu_+ - mu_- = E(Q) - E(Q-2)  (i.e. mu_+ = omega at mu_- = -mu_+): checked.
   The hard X carries E_+; nothing is 'erased'.  The core loses E(Q) - E(Q-2) = 2 omega m_X of energy and 2 charges.
2. Surface absorption / evaporation of X (and of Xbar on anti-objects), Langmuir kinetic model with detailed
   balance (Laine-Shaposhnikov 1998 thermodynamics [R]: growth for mu_gas > omega, evaporation for mu_gas < omega):
       dQ/dt|_abs = eta A/4 int d^3p/(2pi)^3 v f_+(p; T, mu_+),      dQ/dt|_em = eta A/4 int ... f_eq(p; T, omega m_X),
   A = 4 pi r90^2 from our profile (r90 = 5.813/m_X); eta in [0, 1] is a labelled bracket: at the cold-core linear
   level the two-body X absorption is closed (no hidden state lighter than m_X - omega), so eta -> 0 for a cold core,
   while a thermally excited surface absorbs with eta = O(1).  Energy: absorbed int E v f, emitted int E v f_eq; the
   solver books the mismatch to the gas ('cold object' option) or to an object temperature (not modelled here).
3. Object/anti-object encounters: geometric bracket sigma = pi (2 r90)^2, rate n_L n_Lbar sigma v_rel.  Not computed.
All thermal integrals use Bose occupation with the chemical potentials passed in; the covered fraction of the
f_- weight inside the table's momentum range is reported with every forward average.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

import kt_common as C

RES = Path(__file__).resolve().parents[1] / 'results'
CONV = Path(__file__).resolve().parents[2] / 'conversion' / 'results'
HBARC_CM = 1.973269804e-14
C_CM_S = 2.99792458e10
GEV2_TO_CM3S = HBARC_CM**2 * C_CM_S          # 1 GeV^-2 (times c) in cm^3/s = 1.16733e-17
MX = C.MX
OM = C.OMEGA_REF


def _load_table():
    d = json.loads((CONV / 'rate_scan.json').read_text())['dense_scan']['table']
    p = np.array([x[0] for x in d]); sv = np.array([x[1] for x in d]) / GEV2_TO_CM3S    # GeV^-2
    # prepend the cold limit (s-wave constant) from the conversion track's spot value
    spot = json.loads((CONV / 'scatter_spot.json').read_text())['rows'][0]
    p = np.r_[spot['p_over_m'], p]; sv = np.r_[spot['sigma_v_cm3_s'] / GEV2_TO_CM3S, sv]
    p, idx = np.unique(np.round(p, 9), return_index=True)      # the dense scan repeats a few nodes
    return p, sv[idx]


P_TAB, SV_TAB = _load_table()
P_MAX = float(P_TAB[-1])
_sv_log = PchipInterpolator(np.log(P_TAB), np.log(SV_TAB))


def profile_geometry():
    d = json.loads((RES / 'objects_EQ.json').read_text())['reference_Q200000']
    return dict(r50=d['r50_over_mX'] / MX, r90=d['r90_over_mX'] / MX, r99=d['r99_over_mX'] / MX,
                E_GeV=d['E_GeV'], dE_two_charges_GeV=2 * OM * MX)


GEO = profile_geometry()


def sv_forward(p_over_m, extrap='hold'):
    """sigma v (GeV^-2) for Xbar + core(Q=2e5) -> X + core(Q-2) at Xbar momentum p/m_X; p is the Xbar momentum."""
    p = np.asarray(p_over_m, float)
    out = np.exp(_sv_log(np.log(np.clip(p, P_TAB[0], P_MAX))))
    if extrap == 'geometric':
        v = p / np.sqrt(1 + p * p)
        geo = np.pi * GEO['r90']**2 * v
        out = np.where(p > P_MAX, geo, out)
    return out


def sv_reverse(p_plus_over_m, extrap='hold'):
    """sigma v (GeV^-2) for X(p_+) + core(Q-2) -> Xbar + core(Q), by reciprocity from the forward table.
    Threshold: E_+ >= m + 2 omega m (p_+ >= 2.532 m)."""
    pp = np.asarray(p_plus_over_m, float)
    Ep = np.sqrt(1 + pp * pp)
    Em = Ep - 2 * OM
    open_ = Em > 1.0
    pm = np.sqrt(np.clip(Em * Em - 1, 0, None))
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(open_, sv_forward(np.where(open_, pm, 1.0), extrap) * pm * Em / (pp * Ep), 0.0)
    return r


def bose(E_over_m, t, u):
    x = (E_over_m - u) / t
    return 1.0 / np.expm1(x) if x < 600 else 0.0


def conversion_rates(t, up, um, extrap='hold', bose_factors=True):
    """Forward and reverse conversion rates per object per unit Xbar density ... returned as thermal averages:
       <sv>_f (GeV^-2) over f_-(T, mu_-), <sv>_r over f_+(T, mu_+), both with/without 1+f factors, plus the
       densities so that R_fwd = n_L n_- <sv>_f etc.  Energies in m_X units: t = T/m, up = mu_+/m, um = mu_-/m."""
    DE = 2 * OM            # E(Q) - E(Q-2) in m_X units (table value 2 omega; O(1/Q) ignored)
    pts = [0.0, 0.05, 0.2, 0.5, 0.708, 1.5, 3.0, 8.0, 30.0, 150.0]

    def integ(f, a, b):
        return sum(quad(f, x, y, epsabs=0, epsrel=1e-9, limit=200)[0] for x, y in zip(pts, pts[1:]) if y > a and x < b
                   for x, y in [(max(x, a), min(y, b))])

    def wf(p):                                   # forward weight p^2 f_-(E_-)[1 + f_+(E_+)]
        Em = np.sqrt(1 + p * p); Ep = Em + DE
        return p * p * bose(Em, t, um) * ((1 + bose(Ep, t, up)) if bose_factors else 1.0)

    def wr(pp):                                  # reverse weight p_+^2 f_+(E_+)[1 + f_-(E_-)]
        Ep = np.sqrt(1 + pp * pp); Em = Ep - DE
        if Em <= 1:
            return 0.0
        return pp * pp * bose(Ep, t, up) * ((1 + bose(Em, t, um)) if bose_factors else 1.0)

    nm = integ(lambda p: p * p * bose(np.sqrt(1 + p * p), t, um), 0, 150)            # ∝ n_- (×T^3/2π^2 → m^3)
    npl = integ(lambda p: p * p * bose(np.sqrt(1 + p * p), t, up), 0, 150)
    num_f = integ(lambda p: wf(p) * sv_forward(p, extrap), 0, 150)
    cov = integ(wf, 0, P_MAX) / max(integ(wf, 0, 150), 1e-300)
    pp_thr = np.sqrt((1 + DE)**2 - 1)
    num_r = integ(lambda pp: wr(pp) * sv_reverse(pp, extrap), pp_thr, 150)
    # energy moments (m_X units): energy removed from the gas per forward event = E_- , returned = E_+
    e_in = integ(lambda p: wf(p) * sv_forward(p, extrap) * np.sqrt(1 + p * p), 0, 150)
    scale = 1.0 / (2 * np.pi**2)                 # densities: n = ∫ p^2 f dp /(2π^2) in m^3 units
    n_minus, n_plus = nm * scale, npl * scale
    Rf = num_f * scale                           # per object per volume: R_fwd = n_L * Rf * m^3 (GeV^-2 * GeV^3)
    Rr = num_r * scale
    return dict(t=t, up=up, um=um, extrap=extrap, bose_factors=bose_factors,
                n_minus_m3=n_minus, n_plus_m3=n_plus,
                sv_forward_avg_GeV2=Rf / n_minus if n_minus > 0 else 0.0,
                sv_reverse_avg_GeV2=Rr / n_plus if n_plus > 0 else 0.0,
                rate_fwd_per_object_GeV=Rf * MX**3, rate_rev_per_object_GeV=Rr * MX**3,
                mean_E_in_over_m=e_in / num_f if num_f > 0 else 0.0,
                mean_E_out_over_m=(e_in / num_f + DE) if num_f > 0 else 0.0,
                covered_weight_fraction_in_table=cov,
                net_charge_change_per_object_GeV=-2 * (Rf - Rr) * MX**3,
                core_energy_change_per_object_GeV=-(Rf - Rr) * MX**3 * DE * MX)


def detailed_balance_residual(t, um, extrap='hold'):
    """Net forward-reverse rate at the equilibrium relation mu_+ = mu_- + 2 omega (requires mu_- < 1 - 2 omega)."""
    up = um + 2 * OM
    if up >= 1:
        raise ValueError('mu_+ would reach the mass')
    r = conversion_rates(t, up, um, extrap)
    f, b = r['rate_fwd_per_object_GeV'], r['rate_rev_per_object_GeV']
    return dict(t=t, um=um, up=up, forward=f, reverse=b, relative_residual=(f - b) / (f + b))


def surface_flux(t, u, radius_GeV=None, eta=1.0, mass_over_m=1.0):
    """Langmuir surface number and energy flux of one bosonic species (T, mu) onto a sphere of radius r:
    number/time = eta * 4 pi r^2 * (1/4) ∫ d^3p/(2pi)^3 v f ;  energy/time with an extra E.  GeV units."""
    r = GEO['r90'] if radius_GeV is None else radius_GeV
    A = 4 * np.pi * r * r
    pts = [0.0, 0.1, 0.5, 2.0, 8.0, 30.0, 150.0]

    def I(k):
        f = lambda p: p * p * (p / np.sqrt(p * p + mass_over_m**2)) * bose(np.sqrt(p * p + mass_over_m**2), t, u) * \
            (np.sqrt(p * p + mass_over_m**2) if k == 1 else 1.0)
        return sum(quad(f, a, b, epsabs=0, epsrel=1e-9)[0] for a, b in zip(pts, pts[1:]))
    nflux = eta * A / 4 * I(0) / (2 * np.pi**2) * MX**3          # GeV
    eflux = eta * A / 4 * I(1) / (2 * np.pi**2) * MX**4          # GeV^2
    return dict(number_per_time_GeV=nflux, energy_per_time_GeV2=eflux, radius_GeV=r, eta=eta)


def surface_exchange(t, up, eta=1.0, omega=OM):
    """Net X number and energy exchange of one object with the gas (absorption from (t, mu_+) minus emission at
    (t, omega)).  Positive = object gains.  Evaporation = negative values at mu_+ < omega."""
    a = surface_flux(t, up, eta=eta)
    e = surface_flux(t, omega, eta=eta)
    return dict(t=t, up=up, omega=omega, eta=eta,
                dQ_dt_GeV=a['number_per_time_GeV'] - e['number_per_time_GeV'],
                absorbed_number_rate_GeV=a['number_per_time_GeV'], emitted_number_rate_GeV=e['number_per_time_GeV'],
                energy_in_GeV2=a['energy_per_time_GeV2'], energy_out_GeV2=e['energy_per_time_GeV2'],
                energy_to_core_at_omega_per_charge_GeV2=(a['number_per_time_GeV'] - e['number_per_time_GeV']) * omega * MX,
                energy_mismatch_returned_to_gas_GeV2=(a['energy_per_time_GeV2'] - e['energy_per_time_GeV2'])
                - (a['number_per_time_GeV'] - e['number_per_time_GeV']) * omega * MX)


def object_collision_bracket(n_Lbar_GeV3, v_rel):
    """Geometric object/anti-object encounter rate per object (GeV): sigma = pi (2 r90)^2."""
    return n_Lbar_GeV3 * np.pi * (2 * GEO['r90'])**2 * v_rel


if __name__ == '__main__':
    print('table', len(P_TAB), 'p range', P_TAB[0], P_MAX, 'cold sv GeV^-2', SV_TAB[0], 'r90 GeV^-1', GEO['r90'])
    for t in (0.03, 0.1, 0.3, 1.0, 3.0, 7.34, 12.0):
        r = conversion_rates(t, 0.0, 0.0)
        print('t=%6.2f <sv>_f=%.4e GeV^-2 (%.3e cm3/s) covered=%.4f <E_in>=%.3f <E_out>=%.3f  <sv>_r=%.3e' % (
            t, r['sv_forward_avg_GeV2'], r['sv_forward_avg_GeV2'] * GEV2_TO_CM3S, r['covered_weight_fraction_in_table'],
            r['mean_E_in_over_m'], r['mean_E_out_over_m'], r['sv_reverse_avg_GeV2']))
    for t, um in ((0.3, -1.0), (1.0, -0.8), (3.0, -0.75)):
        print('detailed balance', detailed_balance_residual(t, um))
    for t in (0.1, 1.0, 7.34):
        print('surface', t, surface_exchange(t, 0.5))
