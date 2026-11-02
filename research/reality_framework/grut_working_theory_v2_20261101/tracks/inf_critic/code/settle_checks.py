#!/usr/bin/env python3
"""inf_critic (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): quick checks that settle the
contradictions between the build, chain and skeptic tracks, plus the per-benchmark derived quantities the final verdict
needs (T_RH, leptogenesis by-product, portal window, Higgs extrapolations, retrospective n_s pulls).

Inputs, all READ-ONLY: tracks/inf_chain/results/summary.json and chain_*.json; tracks/inf_chain/skeptic/results/
sk_higgs_TB.json; tracks/inf_higgs/code (ih_inflation, ih_common, smpot via higgs_reheat) for one inflation check.
Every input is a model choice; nothing here is a validation or an exclusion.

Checks (labels):
  K1 [C]  S/Q from the dq_entropy closed form at inf_chain's W for every row (ledger re-derivation of the final step).
  K2 [C]  secant fixed points of the one-step calibration iterations (HT6t, T1) in ln H.
  K3 [C/E] decomposition of the Gamma_cal disagreement inf_reheat (1e-6) / its skeptic (1e-7) / inf_chain (3.4e-9).
  K4 [C]  T_RH = (90/pi^2 g*)^{1/4} (Gamma M_P)^{1/2}; T_max from a perturbative-decay fluid ODE.
  K5 [E]  non-thermal leptogenesis from phi_I -> N N for TB1/TB2; Davidson-Ibarra fraction; inverse-decay washout.
  K6 [E]  R-even portal window lam_min = (20/3)(m_I/M_P)^2 at TB; S1 induced loop terms; S2 decay constant f.
  K7 [C]  Higgs during inflation on the TB1 history with and without lam_IH = lam_min (inf_higgs criterion, tree).
  K8 [E]  post-inflation gravitational-only M_t,crit extrapolated in ln H_i from the inf_higgs skeptic's N6 table.
  K9 [R/C] n_s pulls: Planck 2018 0.9649 +- 0.0042; ACT DR6 P-ACT-LB 0.9743 +- 0.0034 (arXiv:2503.14452, as quoted by
          the inf_potential skeptic; retrospective).
Output ../results/settle_checks.json
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]
RES = HERE.parent / 'results'
MP = 2.435e18
GS = 106.75
H_M1 = 1.342842735e11
OUT = {}

S = json.load(open(TRK / 'inf_chain' / 'results' / 'summary.json'))
ROWS = S['rows']


def sq_closed(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


# ---------------------------------------------------------------- K1
k1 = {}
for k, r in ROWS.items():
    W = r['heat_e1e6']['W']
    k1[k] = dict(W=W, SQ_closed=sq_closed(W), SQ_track=r['S_per_Q_e1e6'], rel=sq_closed(W) / r['S_per_Q_e1e6'] - 1,
                 W_check=r['release']['E_over_Ec'] * (1 + r['heat_e1e6']['R_total']))
OUT['K1_SQ_closed_form'] = k1

# ---------------------------------------------------------------- K2
k2 = {}
for k, (Hpot0, Hreq0, Hreq1) in S['iteration'].items():
    Hpot1 = ROWS[k + '_cal1']['H_i_potential_GeV']
    x0, f0 = math.log(Hpot0), math.log(Hreq0 / Hpot0)
    x1, f1 = math.log(Hpot1), math.log(Hreq1 / Hpot1)
    xs = x1 - f1 * (x1 - x0) / (f1 - f0)
    k2[k] = dict(mismatch_ln=[f0, f1], contraction=f1 / f0, secant_fixed_point_GeV=math.exp(xs),
                 slope_dlnHreq_dlnHpot=(math.log(Hreq1) - math.log(Hreq0)) / (x1 - x0),
                 r_at_fixed_point=ROWS[k]['r'] * (math.exp(xs) / Hpot0) ** 2)
OUT['K2_iteration_secant'] = k2

# ---------------------------------------------------------------- K3
tb1 = json.load(open(TRK / 'inf_chain' / 'results' / 'chain_TB1_g1.000e-06.json'))
Y6 = tb1['Y_over_YQ_at_HM1']                       # chain yield at Gamma = 1e-6 H_i(S1)
cY = Y6 / 1e-3                                     # Y = cY Gamma^{1/2}
H_TB1 = ROWS['TB1']['H_i_potential_GeV']; mI = ROWS['TB1']['m_I_over_Hi']; M_TB1 = ROWS['TB1']['M']
Yreq1 = 0.9374 * (H_M1 / H_TB1) ** 1.5             # M = 1 requirement at H_i(potential)
G_M1 = (Yreq1 / cY) ** 2
# inf_reheat static formula Y = 2.206 (Gamma/H_end)^{1/2} at its own H_end = H_i convention; B1 H_end = 5.10e12 GeV
Hend_B1 = 5.10e12; conv = H_TB1 / Hend_B1
Y_reheat_static_at_same_phys = 2.206 * (1e-6 * conv) ** 0.5 / conv ** 1.5   # re-expressed in H(S1) units at Gamma(S1)=1e-6
OUT['K3_Gamma_decomposition'] = dict(
    chain_Y_over_YQ_at_Gamma1e6=Y6, static_formula_same_physical_Gamma=Y_reheat_static_at_same_phys,
    rolling_release_yield_ratio=Y6 / Y_reheat_static_at_same_phys,
    Gamma_cal_if_M1_over_Hi=G_M1, Gamma_cal_if_M1_over_mI=G_M1 / mI,
    Gamma_cal_chain_over_Hi=ROWS['TB1']['reheating']['Gamma_cal_over_Hi'], M=M_TB1,
    M_squared=M_TB1 ** 2, predicted_Gamma_cal=G_M1 / M_TB1 ** 2,
    inf_reheat_B1_Gamma_over_m=1.12e-6, inf_reheat_skeptic_B1_Gamma_over_m=3.8e-7,
    note='Gamma_cal = (Y_req/(cY M))^2: the symmetric relic (M^2 ~ 4.7e3) dominates; the computed rolling-release yield '
         'per unit Gamma^{1/2} is below the static oscillator formula, which pushes the other way.')

# ---------------------------------------------------------------- K4
def T_RH(Gam_GeV):
    return (90 / (math.pi ** 2 * GS)) ** 0.25 * math.sqrt(Gam_GeV * MP)


def T_max(H0, Gam):
    """perturbative decay from rho_phi = 3 H0^2 M_P^2 (GeV units); max of T(rho_r)."""
    def rhs(t, y):
        rp, rr = np.exp(y)
        H = math.sqrt((rp + rr) / 3) / MP
        return [-3 * H - Gam, -4 * H + Gam * rp / rr]
    rp0 = 3 * H0 ** 2 * MP ** 2
    y0 = [math.log(rp0), math.log(rp0 * 1e-30)]
    ts = np.logspace(math.log10(1e-3 / H0), math.log10(30 / Gam), 4000)
    sol = solve_ivp(rhs, (0, ts[-1]), y0, t_eval=ts, rtol=1e-9, atol=1e-12, method='LSODA')
    rr = np.exp(sol.y[1])
    T = (30 * rr / (math.pi ** 2 * GS)) ** 0.25
    return float(T.max())


k4 = {}
for k, r in ROWS.items():
    Hi = r['H_i_potential_GeV']
    G = r['Gamma_I_over_Hi'] * Hi
    d = dict(Gamma_GeV=G, T_RH_GeV=T_RH(G), m_I_GeV=r['m_I_over_Hi'] * Hi)
    if r['Gamma_I_over_Hi'] < 1e-3:
        d['T_max_GeV'] = T_max(Hi, G)
    k4[k] = d
OUT['K4_reheating_temperatures'] = k4

# ---------------------------------------------------------------- K5
YB = 8.7e-11                                     # Planck 2018 eta_B = 6.1e-10 -> Y_B = n_B/s [R]
YL_need = YB * 79 / 28
v = 174.0; m3 = 0.050e-9                          # GeV; sqrt(Delta m^2_atm) [R]
k5 = {}
for k in ('TB1', 'TB2'):
    r = ROWS[k]; Hi = r['H_i_potential_GeV']; mI = r['m_I_over_Hi'] * Hi
    TR = k4[k]['T_RH_GeV']
    eps_need = YL_need / (1.5 * TR / mI)
    rows = []
    for MN in (1e12, 2e12, 4e12, 0.45 * mI):
        z = MN / TR
        eps_DI = 3 / (8 * math.pi) * MN * m3 / v ** 2
        Kw = 0.05e-9 / 1.08e-12                  # m_tilde / m_* with m_tilde = 0.05 eV, m_* = 1.08e-3 eV [R]
        # inverse-decay washout exponent from z_RH (BDP approximation W_ID = Kw z^3 K1(z)/4), integrated
        from scipy.special import k1 as K1
        zz = np.linspace(z, z + 60, 6000)
        WID = 0.25 * Kw * zz ** 3 * K1(zz)
        wash = float(np.trapezoid(WID, zz))
        GN = Kw * 1.08e-12 * 0 + (m3 * MN / v ** 2) * MN / (8 * math.pi)   # Gamma_N = (y^2) M/8pi, y^2 = m_tilde M/v^2
        rows.append(dict(M_N_GeV=MN, z_RH=z, eps_DI=eps_DI, eps_need_over_DI=eps_need / eps_DI, washout_exponent=wash,
                         eps_need_with_washout_over_DI=eps_need * math.exp(min(wash, 700)) / eps_DI,
                         Gamma_N_GeV=GN, Gamma_N_over_Gamma_I=GN / (r['Gamma_I_over_Hi'] * Hi)))
    k5[k] = dict(T_RH_GeV=TR, T_max_GeV=k4[k]['T_max_GeV'], m_I_GeV=mI, T_RH_over_mI=TR / mI, Y_L_needed=YL_need,
                 eps1_needed_nonthermal=eps_need, by_M_N=rows)
OUT['K5_leptogenesis'] = k5

# ---------------------------------------------------------------- K6
k6 = {}
for k in ('TB1', 'TB2', 'HT6t', 'T1'):
    r = ROWS[k]; Hi = r['H_i_potential_GeV']; mI = r['m_I_over_Hi'] * Hi
    lam_min = 20 / 3 * (mI / MP) ** 2
    d = dict(m_I_GeV=mI, lam_min=lam_min)
    if k.startswith('TB'):
        G = r['Gamma_I_over_Hi'] * Hi
        y = r['reheating']['y_majorana_N1']
        # S1 (scalar Yukawa) induced terms at M_N = 4e12 GeV and phi_* ~ 4 M_P (T-model alpha = 0.3/1) [E]
        MN = 4e12; ynu2 = m3 * MN / v ** 2
        phis = 4.1 * MP if k == 'TB1' else 5.0 * MP
        Hs = r['H_star_GeV']
        A_ind = y * ynu2 * MN / (8 * math.pi ** 2)
        Mphi = MN + y * phis
        dm2_loop = ynu2 / (8 * math.pi ** 2) * Mphi ** 2         # times [ln(M^2/mu^2) - 1]
        d.update(y_S1=y, A_ind_GeV=A_ind, A_phi_star_over_Hstar2=A_ind * phis / Hs ** 2,
                 S1_loop_mass2_over_Hstar2_per_log=dm2_loop / Hs ** 2,
                 portal_mass2_over_Hstar2_at_lam_min=lam_min * phis ** 2 / 2 / Hs ** 2,
                 min_m2_at_zero_crossing_over_Hi2=-(A_ind ** 2) / (2 * lam_min) / Hi ** 2,
                 # S2: y_p = M_N/f, Majorana pseudoscalar width Gamma = y_p^2 m beta/(16 pi)
                 S2_f_GeV_at_MN=dict((('%.0e' % M), M * math.sqrt(math.sqrt(1 - 4 * M ** 2 / mI ** 2) / (16 * math.pi * G / mI)))
                                     for M in (1e12, 4e12)))
    k6[k] = d
OUT['K6_portal_and_channels'] = k6

# ---------------------------------------------------------------- K7: Higgs during inflation on TB1, lam_IH on/off
sys.path.insert(0, str(TRK / 'inf_higgs' / 'code'))
sys.path.insert(0, str(TRK / 'inf_chain' / 'code'))
try:
    import ih_inflation as IHI                      # read-only
    import ch_common as cc                          # read-only
    ic = cc.ic
    b = cc.tmodel_bench(0.3, Nstar=ROWS['TB1']['N_star'], name='TB1')
    Hi = b.Hi; A = b.V0 / b.eta2
    Nb = np.unique(np.concatenate([np.linspace(0, 0.3, 16), np.linspace(0.3, 1.2, 13), np.linspace(1.2, 6.0, 12)]))
    # Nb counted from eps_H = 1 (inf_higgs convention); S1 is Nb = 1
    rows = []
    for x in Nb:
        f, fN = ic.bg_at(b.bg, x)
        eps = 0.5 * fN * fN
        H = math.sqrt(A * b.m.v(f) / (3 - eps))
        rows.append((x, H * Hi, 6 * H * H * (2 - eps) * Hi ** 2, eps, f))
    f, fN = ic.bg_at(b.bg, 50.0)
    eps = 0.5 * fN * fN
    Hp = math.sqrt(A * b.m.v(f) / (3 - eps)) * Hi
    hist = dict(rows=np.array(rows), H_plateau=Hp, phiI_plateau=f, H_GeV=Hi)
    lam_min = k6['TB1']['lam_min']
    k7 = dict(H_i_GeV=Hi, H_plateau_GeV=Hp, lam_min=lam_min, N0=[58.8, 60.0], by_Mt={})
    for Mt in (172.57, 172.86, 173.15, 173.75, 174.6):
        sm = IHI.SMfull(Mt)
        rr = {}
        for N0 in (58.8, 60.0):
            l0, i0 = IHI.log10_bubbles(hist, sm, 0.1, N0, False)
            l1, i1 = IHI.log10_bubbles(hist, sm, 0.1, N0, False, lamIH=lam_min)
            l2, i2 = IHI.log10_bubbles(hist, sm, 0.1, N0, False, lamIH=0.1 * lam_min)
            rr['N0=%g' % N0] = dict(log10N_no_portal=l0, B_end_no_portal=i0['B_end'], log10N_lam_min=l1,
                                    B_end_lam_min=i1['B_end'], log10N_0p1_lam_min=l2, B_end_0p1=i2['B_end'])
        k7['by_Mt']['%.2f' % Mt] = rr
        print('K7', Mt, json.dumps(rr), flush=True)
    sk = json.load(open(TRK / 'inf_chain' / 'skeptic' / 'results' / 'sk_higgs_TB.json'))['rows']['TB1']['by_Mt']
    k7['regression_vs_inf_chain_skeptic_log10N_xi01'] = {m: sk[m].get('log10N_at_xi01') for m in sk}
    OUT['K7_Higgs_inflation_TB1'] = k7
except Exception as exc:                            # noqa: BLE001
    OUT['K7_Higgs_inflation_TB1'] = dict(error=repr(exc))
    print('K7 failed', repr(exc))

# ---------------------------------------------------------------- K8
Ht = np.array([1.29e11, 1.5e11, 2e11, 3e11, 5e11, 1e12]); Mc = np.array([172.917, 172.873, 172.790, 172.684, 172.561, 172.447])
lnH = np.log(Ht)
lo = np.polyfit(lnH[:3], Mc[:3], 1); hi = np.polyfit(lnH[-2:], Mc[-2:], 1)
k8 = dict(table_H=Ht.tolist(), table_Mtcrit=Mc.tolist(), slope_low=lo[0], slope_high=hi[0])
for k in ('HT6t', 'T1', 'T1_cal1', 'TB1', 'TB2'):
    r = ROWS[k]
    Hreq = r['H_i_req_GeV'] if not k.startswith('TB') else r['H_i_potential_GeV']
    if Hreq < Ht[0]:
        est = float(np.polyval(lo, math.log(Hreq)))
    else:
        est = float(np.polyval(hi, math.log(Hreq)))
    k8[k] = dict(H_i_consistent_GeV=Hreq, Mtcrit_grav_only_extrapolated=est)
OUT['K8_Higgs_post_inflation_extrapolation'] = k8

# ---------------------------------------------------------------- K9
k9 = {}
for k, r in ROWS.items():
    ns = r['n_s']
    k9[k] = dict(n_s=ns, pull_Planck2018=(ns - 0.9649) / 0.0042, pull_ACT_DR6=(ns - 0.9743) / 0.0034, r=r['r'])
OUT['K9_ns_pulls'] = k9

OUT['method'] = __doc__
json.dump(OUT, open(RES / 'settle_checks.json', 'w'), indent=1, default=float)
for key in OUT:
    if key != 'method':
        print('##', key, json.dumps(OUT[key], default=float)[:2000])
