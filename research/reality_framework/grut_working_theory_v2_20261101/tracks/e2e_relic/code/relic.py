#!/usr/bin/env python3
"""e2e_relic items 2-3 (Claude): the dark-matter relic of M1r/M2 as a callable function of the history's outputs, and its
observables.

    relic(SQ, H_i=None, r_hv=None, e=1e-6, Y_Q=None, trapping='both')  -> dict
      SQ   : hidden entropy per charge after the release (the chain's S/Q(e));
      H_i  : GeV (used only for the hidden-vector mass check m_A = 6 e h = 240 e H_i and echoed);
      r_hv : hidden-to-visible entropy ratio s_h/s_vis after the release; Y_Q = n_Q/s = r_hv/(S/Q);
      Y_Q  : alternatively the charge yield directly; if neither r_hv nor Y_Q is given the abundance-calibrated yield
             Y_Q = D/(M m_X) (D = 4.4e-10 GeV, supplied) is solved self-consistently with M;
      e    : hidden gauge coupling (used only for the screening check; at fixed S/Q the relic does not depend on e [D]).
    observables(rel) -> dict : object mass, radius, local number density, flux, impulse, event rates, self-interaction,
                               P(k) proxy (qm_transfer fit), or the free-X properties.

Physics content (every element a model choice inherited from the named track; nothing re-derived here):
  * Branch: Q-matter iff S/Q < 3.4716 (bulk threshold, fx_chain skeptic/dq_entropy) [Q].  Ripened objects (Q ~ 1e20) are
    bulk-like, so the bulk threshold applies (the Q = 2e5 threshold 2.928 is irrelevant once ripening is included) [E].
  * Free X: p(S/Q, Y_Q) from code/pmap.py (case-L network, relaxed start, gap filled); M = 1 + 2p [C on conditional recipe].
  * Q-matter: <Q> (charge-weighted at capture freeze-out) from qm_micro's coarsening with the qm_micro skeptic's Bose and
    latent-heat corrections (coarsen_heat bose_heat_sphi1 / sphi0); f_v from the skeptic's energy-conserving freeze-out
    (sk_freeze), perfect trapping and Beliaev trapping; M = <E/Q>_obj + Delta f_v with <E/Q> = 0.775866, Delta = 1 - Omega_b
    = 0.224135; p ~ 0 (Xbar destroyed on objects; n_Xbar/n_X ~ e^-170 at T_f) [C on E models, qm_micro + skeptic].
    Y_Q dependence (this track, relic_runs.py, S/Q = 3.0, Y_Q x 0.1 and x 10): <Q> ∝ Y_Q^1.00, f_v ∝ Y_Q^-0.035 (perfect)
    / Y_Q^-0.09 (Beliaev); applied at every S/Q [E: measured at one S/Q].  Reference yield Y_Q,ref = 6.0653e-9.
  * S/Q in (3.2, 3.4716): linear extrapolation of the S/Q = 3.0-3.2 tables [E, flagged].  S/Q < 1.5: clamped [E, flagged];
    on M1r's own heat S/Q < 1.84 is unreachable (e -> infinity floor).
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pmap  # noqa: E402

MX = 0.0924089637375362          # GeV
D_TARGET = 4.4e-10               # GeV per unit entropy [supplied; Planck 2018 4.37e-10]
YQ_M1 = 5.07941423e-9
YQ_REF_QM = YQ_M1 * 0.9374 / 0.785   # qm_micro's visible clock (6.0653e-9)
SQ_BULK = 3.4716191256166837
EQ_OBJ = 0.775866
DELTA = 1 - 0.7758647074595397
NB = 449.31336788454803          # bulk charge density, m_X^3 (qm_micro)
HBARC_GEV_FM = 0.1973269804
GEV_KG = 1.78266192e-27
MPL_KG = 2.176434e-8             # Planck mass
MPL_RED_KG = MPL_KG / math.sqrt(8 * math.pi)
GN = 6.67430e-11

# ---- Q-matter tables at Y_Q,ref (qm_micro skeptic, Bose + heat; this track reproduces the S/Q = 3.0 row)
QT_SQ = [1.5, 2.5, 3.0, 3.2]
QT_Q_SPHI1 = [3.82e19, 7.96e19, 1.267e20, 1.45e20]
QT_Q_SPHI0 = [3.61e19, 7.67e19, 1.17e20, 1.50e20]
QT_FV_PERF = [0.0555, 0.0971, 0.1197, 0.1292]
QT_FV_BEL = [0.0827, 0.1414, 0.1754, 0.1890]
QT_NOBJ = [4.45e-20, 1.75e-20, 9.39e-21, 7.4e-21]      # objects per unit charge (Bose + heat)
QT_QCAP_OVER_Q = 0.794                                   # capture-effective / charge-weighted <Q> (S/Q = 3.0 run)
S_EFF_TF = 0.47                                          # flux-averaged sticking at T_f (qm_micro)
T_F = 0.0089                                             # capture freeze-out T_h/m_X (qm_micro, skeptic control)
# Y_Q exponents measured at S/Q = 3.0 (relic_runs.py): below / above the reference yield
A_Q = (1.0093, 1.0017)
A_FV_PERF = (-0.0372, -0.0324)
A_FV_BEL = (-0.0950, -0.0859)


def _tab(SQ, ys):
    if SQ <= QT_SQ[0]:
        return ys[0]
    if SQ >= QT_SQ[-1]:
        return ys[-1] + (ys[-1] - ys[-2]) / (QT_SQ[-1] - QT_SQ[-2]) * (SQ - QT_SQ[-1])
    return float(np.interp(SQ, QT_SQ, ys))


def _pw(f, a):
    return f ** (a[0] if f < 1 else a[1])


def qmatter(SQ, Y_Q):
    f = Y_Q / YQ_REF_QM
    flags = []
    if SQ > QT_SQ[-1]:
        flags.append('S/Q in (3.2, 3.47): linear extrapolation of the qm_micro tables [E]')
    if SQ < QT_SQ[0]:
        flags.append('S/Q < 1.5: tables clamped [E]; on M1r heat S/Q < 1.84 is unreachable')
    if not (0.1 <= f <= 10):
        flags.append('Y_Q outside the computed x0.1..x10 range: power-law extrapolation [E]')
    Q1, Q0 = _tab(SQ, QT_Q_SPHI1) * _pw(f, A_Q), _tab(SQ, QT_Q_SPHI0) * _pw(f, A_Q)
    fvp = _tab(SQ, QT_FV_PERF) * _pw(f, A_FV_PERF)
    fvb = _tab(SQ, QT_FV_BEL) * _pw(f, A_FV_BEL)
    Mp, Mb = EQ_OBJ + DELTA * fvp, EQ_OBJ + DELTA * fvb
    nobj = _tab(SQ, QT_NOBJ) / _pw(f, A_Q)
    return dict(branch='Q-matter', Q_mean=Q1, Q_mean_range=sorted([Q0, Q1]), Q_capture_eff=Q1 * QT_QCAP_OVER_Q,
                objects_per_charge=nobj, f_v=dict(perfect=fvp, beliaev=fvb), M=dict(perfect=Mp, beliaev=Mb),
                M_central=0.5 * (Mp + Mb), M_range=[Mp, Mb], p=0.0, Xbar_over_X='~exp(-170)', YQ_over_ref=f, flags=flags)


def freeX(SQ, Y_Q):
    f = Y_Q / YQ_M1
    p = pmap.p_of_SQ(SQ, f)
    pa = pmap.p_aref(SQ, f) if SQ >= 4.936 else None
    flags = []
    if SQ < 5 or SQ > 958:
        flags.append('Y_Q exponent beta held constant outside S/Q = 5-958 [E]')
    if SQ > 1307.7:
        flags.append('p extrapolated above the last network row (S/Q = 1308) [E]')
    return dict(branch='free X', p=p, p_aref=pa, M_central=1 + 2 * p, M_range=sorted([1 + 2 * p, 1 + 2 * (pa if pa else p)]),
                Xbar_over_X=p / (1 + p), asym_share=1 / (1 + 2 * p), YQ_over_M1=f, flags=flags)


def relic(SQ, H_i=None, r_hv=None, e=1e-6, Y_Q=None, trapping='both'):
    out = dict(inputs=dict(S_per_Q=SQ, H_i_GeV=H_i, r_hv=r_hv, e=e, Y_Q=Y_Q))
    branch = freeX if SQ > SQ_BULK else qmatter
    if r_hv is not None:
        Y_Q = r_hv / SQ
    calibrated = Y_Q is None
    if calibrated:                                   # abundance-calibrated yield, self-consistent with M
        Y_Q = D_TARGET / (0.80 * MX) if SQ <= SQ_BULK else D_TARGET / MX
        for _ in range(50):
            r = branch(SQ, Y_Q)
            Mx = r['M_central'] if trapping == 'both' else r['M'][trapping]
            Yn = D_TARGET / (Mx * MX)
            if abs(Yn / Y_Q - 1) < 1e-12:
                break
            Y_Q = Yn
    r = branch(SQ, Y_Q)
    out.update(r)
    Muse = r['M_central'] if (trapping == 'both' or r['branch'] == 'free X') else r['M'][trapping]
    out['Y_Q'] = Y_Q
    out['r_hv'] = Y_Q * SQ
    out['Y_Q_calibrated'] = calibrated
    out['Omega_over_target'] = Muse * MX * Y_Q / D_TARGET
    # hidden vector: Higgsed by Sigma, m_A = 6 e h, h = 40 H_i [supplied M1r]; at fixed S/Q e enters nowhere else
    if H_i is not None:
        mA = 240 * e * H_i
        out['m_A_GeV'] = mA
        out['vector_screened'] = mA > 10 * MX
        if not out['vector_screened']:
            out.setdefault('flags', []).append('m_A < 10 m_X: a light hidden vector would change the relic; tables invalid')
    return out


# ------------------------------------------------------------------ observables
def shm_mean_speed(v0=238.0, vE=250.5, vesc=544.0, n=4000):
    """<|v|> in the Earth frame for the truncated Standard Halo Model (km/s); Baxter et al. 2021 conventions [R]."""
    rng = np.random.default_rng(1)
    v = rng.normal(0, v0 / math.sqrt(2), size=(n * 50, 3))
    v = v[np.linalg.norm(v, axis=1) < vesc]
    v[:, 2] -= vE
    return float(np.mean(np.linalg.norm(v, axis=1)))


def sigma_freeX():
    """Tree-level low-energy X X and X Xbar elastic cross sections (cm^2/g) from lambda_X|X|^4 and phi exchange
    (trilinear 2 g^2 v phi|X|^2, g v = m_X, m_phi^2 = 6 m_X^2).  [D at tree level]"""
    lamX, g = 1e-3, 0.1
    tri2 = 4 * g ** 2                      # (2 g^2 v)^2 / m_X^2
    M_XX = -4 * lamX + tri2 * (2 / 6)      # t + u channel at threshold
    M_XXb = -4 * lamX + tri2 * (1 / 2 + 1 / 6)   # s (1/(m_phi^2 - 4 m^2)) + t channel
    gev2_cm2 = 0.3893794e-27
    s = 4.0                                 # s in m_X^2
    sig_XX = M_XX ** 2 / (32 * math.pi * s) / MX ** 2 * gev2_cm2      # identical particles
    sig_XXb = M_XXb ** 2 / (16 * math.pi * s) / MX ** 2 * gev2_cm2
    mg = MX * GEV_KG * 1e3
    return dict(sigma_XX_cm2=sig_XX, sigma_XXbar_cm2=sig_XXb, sigma_over_m_XX=sig_XX / mg, sigma_over_m_XXbar=sig_XXb / mg)


def observables(rel, rho_loc=0.3, v0=238.0, vE=250.5, vesc=544.0, sensor_mass_kg=1e-6, b_m=1e-3):
    vbar = shm_mean_speed(v0, vE, vesc)          # km/s
    v_ms = vbar * 1e3
    yr = 3.15576e7
    o = dict(rho_loc_GeV_cm3=rho_loc, mean_speed_km_s=vbar)
    if rel['branch'] == 'Q-matter':
        Q = rel['Q_mean']
        m_obj_GeV = EQ_OBJ * Q * MX
        m_kg = m_obj_GeV * GEV_KG
        R_mX = (3 * Q / (4 * math.pi * NB)) ** (1 / 3)
        R_nm = R_mX * HBARC_GEV_FM / MX * 1e-6
        frac_obj = 1 - (rel['f_v']['perfect'] * 1 / rel['M']['perfect'])     # mass fraction in objects (vapour mass m_X per charge)
        n_obj_cm3 = rho_loc * frac_obj / m_obj_GeV
        flux_m2_yr = n_obj_cm3 * 1e6 * v_ms * yr
        I = 2 * GN * m_kg * sensor_mass_kg / (b_m * v_ms)                       # impulse on a sensor at impact parameter b
        tau = b_m / v_ms
        sig_geo_cm2 = math.pi * (2 * R_nm * 1e-7) ** 2
        zeta = S_EFF_TF * rel['Q_capture_eff'] ** (-1 / 3)
        fv_mass = rel['f_v']['perfect'] / rel['M']['perfect']
        sp = [0.038 * (zeta * T_F) ** 0.99, 0.081 * (zeta * T_F) ** 0.99]       # qm_transfer track / skeptic converged
        meff = [19.3 * (s / 1e-8) ** -0.49 * (fv_mass / 0.1) ** -0.45 for s in sp]
        Qr = rel['Q_mean_range']
        o.update(object_mass_GeV=m_obj_GeV, object_mass_kg=m_kg, object_mass_ug=m_kg * 1e9,
                 object_mass_range_kg=[EQ_OBJ * q * MX * GEV_KG for q in Qr],
                 object_mass_over_M_Planck=m_kg / MPL_KG, object_mass_over_reduced_M_Planck=m_kg / MPL_RED_KG,
                 radius_nm=R_nm, mass_fraction_in_objects=frac_obj,
                 n_local_per_cm3=n_obj_cm3, n_local_per_m3=n_obj_cm3 * 1e6,
                 flux_per_m2_per_yr=flux_m2_yr, mean_separation_m=(n_obj_cm3 * 1e6) ** (-1 / 3),
                 event_rate_per_yr={'A=1 m^2': flux_m2_yr, 'A=(10 cm)^2': flux_m2_yr * 1e-2,
                                    'A=(3 m)^2': flux_m2_yr * 9},
                 impulse_on_sensor=dict(sensor_mass_kg=sensor_mass_kg, b_m=b_m, impulse_kg_m_s=I, crossing_time_s=tau,
                                        note='I = 2 G M_obj m_s/(b v) (straight-line flyby); detectability is the proposals\' '
                                             'claim, not assessed here'),
                 self_interaction=dict(sigma_geo_cm2=sig_geo_cm2, sigma_over_m_cm2_g=sig_geo_cm2 / (m_kg * 1e3),
                                       note='contact (merger) cross-section pi(2R)^2 per object mass; upper estimate'),
                 vapour=dict(f_v_mass_perfect=fv_mass, zeta=zeta, sigma_p=sp, m_eff_Lya_keV=meff,
                             note='qm_transfer fit m_eff = 19.3 keV (sigma_p/1e-8)^-0.49 (f_v/0.1)^-0.45, sigma_p = C (zeta T_f)^0.99, '
                                  'C = 0.038 (track) / 0.081 (skeptic, converged); fit not portable at <~15% [qm_transfer skeptic R5]; '
                                  'qm_transfer: |dP/P| <= 1e-3 (numerical floor) at k <= 136/Mpc at these zeta'))
    else:
        sx = sigma_freeX()
        p = rel['p']
        xf = p / (1 + 2 * p)                     # Xbar fraction of particles
        so = (1 - xf) ** 2 * sx['sigma_over_m_XX'] + xf ** 2 * sx['sigma_over_m_XX'] + 2 * xf * (1 - xf) * sx['sigma_over_m_XXbar']
        n_cm3 = rho_loc / MX
        o.update(particle_mass_GeV=MX, Xbar_fraction_of_particles=xf, n_local_per_cm3=n_cm3,
                 flux_per_cm2_per_s=n_cm3 * vbar * 1e5, self_interaction=dict(**sx, sigma_over_m_mix_cm2_g=so),
                 annihilation='X Xbar -> phi phi closed (m_phi = 2.45 m_X); no light hidden state; no annihilation signal [D]',
                 P_k='CDM-like (free streaming 0.3-0.5 pc comoving) [Q WORKING_THEORY section 6]',
                 direct_detection='none at zero portals [Q WORKING_THEORY section 6]')
    return o


def regression():
    """Regression against qm_micro's M1r Q-matter benchmark (e = 3e-6, S/Q = 2.99) and the free-X M1r point (S/Q = 4.955)."""
    q = relic(2.99, Y_Q=YQ_REF_QM)
    q3 = relic(3.0, Y_Q=YQ_REF_QM)
    fx = relic(4.955)
    q288 = relic(2.88)
    return dict(
        qmatter_SQ2p99=dict(Q_mean=q['Q_mean'], f_v=q['f_v'], M=q['M'],
                            reference='qm_micro+skeptic: <Q> 1.17-1.27e20, f_v 0.120 (perfect) / 0.175 (Beliaev), M 0.8023-0.8152'),
        qmatter_SQ3p0_at_ref=dict(Q_mean=q3['Q_mean'], f_v=q3['f_v'], M=q3['M'],
                                  reference='this track rerun: <Q> 1.267e20, f_v 0.1197 / 0.1751; skeptic M 0.8027 / 0.8152'),
        qmatter_SQ2p88_calibrated=dict(M_central=q288['M_central'], reference='WORKING_THEORY 4b: M = 0.807 m_X (c0 = 1.40 HFa)'),
        freeX_SQ4p955=dict(M=fx['M_central'], p=fx['p'], reference='WORKING_THEORY 4b: M = 1.002 (M1r, c0 = 1)'))


if __name__ == '__main__':
    import json
    print(json.dumps(regression(), indent=1, default=float))
