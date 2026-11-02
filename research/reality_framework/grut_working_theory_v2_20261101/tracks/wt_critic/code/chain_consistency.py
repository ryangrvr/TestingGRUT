"""wt_critic (Claude, own code): consistency checks across fx_holder, fx_chain, fx_formation, fx_observables.

1. Release-state invariant.  For the pre-onset X sector (massless quartic X + phi, conformal), the relaxed state at fixed
   comoving charge N and conformal energy E depends only on E / E_c(N), with E_c the circular-orbit energy of charge N.
   The engines record eps = N_X / N_c(E_X) with N_c = c E^{3/4}; we check that c is the same constant in the kin_source
   spot point and in every fx_holder run, so E/E_c = |eps|^{-4/3} compares releases without any normalisation of a.
2. Heat budget per charge in units of E_c:  W = |eps|^{-4/3} (1 + R_med) + nu * W_vec,spot(e),
   R_med from fx_holder's own decay ledger for that release (Y_L + cascade + sigma, at the same e),
   W_vec,spot = |eps_spot|^{-4/3} R_vec,spot(e) (kin_source; vectors are inflationary GMR longitudinal modes, so they do
   not depend on the X release) and the transfer factor nu = (H/H_spot)^2 (Y_spot/Y)^{4/3} [D/E, derived in NOTES].
3. S/Q from the copied Hartree solver (fx_chain skeptic hartree_indep.py) at R_eff = W / |eps_spot|^{-4/3} - 1, against
   the survival thresholds sigma_min of fx_chain skeptic survival.py (bulk 3.4716, Q=2e5 2.9284, binding 2.5531).
4. p(S/Q): fx_chain's gas-only R scan mapped through the same Hartree S/Q(R); log-log interpolation.
5. H_i ledger for the repaired chain.
Labels: [C] computed here, [E] estimate (nu transfer), [D] arithmetic."""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hartree_copy as HI   # unchanged copy of fx_chain/skeptic/code/hartree_indep.py

BASE = '/home/user/TestingGRUT/research/reality_framework'
WT = BASE + '/grut_working_theory_v2_20261101/tracks'
TR = BASE + '/grut_formation_review_20261025/tracks'
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'results', 'chain_consistency.json')

src = json.load(open(TR + '/kin_source/results/source_interface.json'))
FP = src['meta']['fixed_point']
VS = src['vector_sources']; MD = src['mediator_debris']['totals']
spot = dict(NX=abs(FP['N_X']), EX=FP['E_X'], eps=abs(FP['eps']), H=FP['H_GeV'], Y=FP['YQ'])
C_spot = spot['NX'] / (spot['eps'] * spot['EX'] ** 0.75)
W0_spot = spot['eps'] ** (-4.0 / 3.0)

# ---- scalings in the interface [C]
scal = dict(vector_R_ratio_1e6_over_1e7=VS['e1e-6_band']['per_sign_totals']['R'] / VS['e1e-7_band']['per_sign_totals']['R'],
            cascade_R_ratio_1e6_over_1e7=src['mediator_debris']['cascade_vectors']['e1e-6']['R'] / src['mediator_debris']['cascade_vectors']['e1e-7']['R'])
scal['vector_exponent'] = math.log(scal['vector_R_ratio_1e6_over_1e7']) / math.log(10.0)
scal['cascade_exponent'] = math.log(scal['cascade_R_ratio_1e6_over_1e7']) / math.log(10.0)

Rvec = dict(band=VS['e1e-6_band']['per_sign_totals']['R'], broad=VS['e1e-6_broad']['per_sign_totals']['R'])  # at e = 1e-6
Wvec_spot = {k: W0_spot * v for k, v in Rvec.items()}   # E_c units per charge at e = 1e-6; scales as 1/e

# ---- releases from fx_holder (its own engine; our code)
TRI = json.load(open(WT + '/fx_holder/results/tri_runs.json'))['runs']
rel = {}
for hist, runs in TRI.items():
    for name, r in runs.items():
        if not isinstance(r, dict) or 'EX' not in r:
            continue
        NX, EX, eps = abs(r['NX_400_1000']), r['EX'], abs(r['eps'])
        L6 = r.get('ledger_e1e-6', {}); L7 = r.get('ledger_e1e-7', {})
        rel[hist + '/' + name] = dict(NX=NX, EX=EX, eps=eps, c=NX / (eps * EX ** 0.75), W0=eps ** (-4.0 / 3.0),
                                       R_med_e1e6=L6.get('R_central'), R_Y=L6.get('R_Y'), R_casc_e1e6=L6.get('R_cascade'),
                                       R_med_e1e7=L7.get('R_central'), NYL=r.get('NYL'), yield_over_YQ=r.get('yield_over_YQ'))
rel['spot_instant_release (kin_source; used by fx_chain/fx_formation)'] = dict(
    NX=spot['NX'], EX=spot['EX'], eps=spot['eps'], c=C_spot, W0=W0_spot,
    R_med_e1e6=MD['e1e-6']['R_mediator_total'], R_Y=MD['e1e-6']['R_Y'], R_casc_e1e6=MD['e1e-6']['R_cascade'],
    R_med_e1e7=MD['e1e-7']['R_mediator_total'], NYL=FP['N_YL'])

# ---- Hartree S/Q at an effective R (spot normalisation)
_g = [None]
def SQ_of_W(W):
    Reff = W / W0_spot - 1.0
    st, g = HI.solve(Reff, True, 0.0, _g[0])
    if (not st['success']) or max(abs(x) for x in st['resid']) > 1e-9:
        st, g = HI.solve(Reff, True, 0.0, None)
    if (not st['success']) or max(abs(x) for x in st['resid']) > 1e-9:
        raise RuntimeError('Hartree solve failed at R_eff = %g' % Reff)
    _g[0] = g
    a = HI.onset_a(st); mc = HI.MX * a / HI.H
    return dict(R_eff=Reff, S_per_Q=st['S'] / HI.NX, T_on_over_mX=st['T'] / mc, nQ_on=HI.NX / mc ** 3,
                Xbar_per_Q_on=st['n_minus'] / HI.NX, cond_frac=st['cond_frac'], resid=max(abs(x) for x in st['resid']))

SURV = json.load(open(WT + '/fx_chain/skeptic/results/survival_threshold.json'))['thresholds']
sig = {k: SURV[k]['sigma_min'] for k in ('bulk_Q_inf', 'Q_3e7', 'Q_2e5', 'Q_bind_71263')}

# W at the thresholds (invert S/Q(W))
def W_of_SQ(s):
    return brentq(lambda W: SQ_of_W(W)['S_per_Q'] - s, W0_spot * 0.9, W0_spot * 8.0, xtol=1e-10)
W_thr = {k: W_of_SQ(v) for k, v in sig.items()}

# ---- p(S/Q) from fx_chain's gas-only scan [C mapping; E that p depends on S/Q only, fx_chain skeptic]
KIN = json.load(open(WT + '/fx_chain/results/kin_m1.json'))
scan = [(R, p) for R, p in KIN['R_scan'] if abs(R - 4.96654) > 1e-6]
sq_scan = []
for R, p in scan:
    st, _ = HI.solve(R, True, 0.0, None)
    sq_scan.append((R, st['S'] / HI.NX, p))
lsq = np.log([x[1] for x in sq_scan]); lp = np.log([x[2] for x in sq_scan])
def p_of_SQ(s):
    return float(np.exp(np.interp(math.log(s), lsq, lp))) if lsq[0] <= math.log(s) <= lsq[-1] else None
dlnp_dlnsq = [(lp[i + 1] - lp[i]) / (lsq[i + 1] - lsq[i]) for i in range(len(lp) - 1)]

# ---- H_i ledger [D]
D_GEV, MX = 4.4e-10, 0.0924089637375362
YQ_ref = 5.07941423e-9            # D / (0.9374 m_X)
H_M1 = 1.342842735e11
yield_tri05 = 0.990733            # fx_holder FX11 (skeptic 0.990733), Y/Y_Q at H_M1, portals zero at M_Y
yield_tri05_planck = 1.30606      # fx_holder skeptic N1, portals zero at M_P
def H_for(yield_ratio, M_per_charge):
    """H_i giving Y = D/(M m_X) with Y = yield_ratio * YQ_ref * (H/H_M1)^{3/2}."""
    Yt = D_GEV / (M_per_charge * MX)
    return H_M1 * (yield_ratio * YQ_ref / Yt) ** (-2.0 / 3.0)

def nu_transfer(H, Y):
    return (H / spot['H']) ** 2 * (spot['Y'] / Y) ** (4.0 / 3.0)

# ---- assemble per release at e = 1e-6 (and e scan)
cases = {}
for name, r in rel.items():
    if r['R_med_e1e6'] is None:
        continue
    Wmed = r['W0'] * (1 + r['R_med_e1e6'])
    row = dict(eps=r['eps'], c=r['c'], E_over_Ec=r['W0'], R_med_e1e6=r['R_med_e1e6'], W_no_vectors=Wmed,
               SQ_no_vectors=SQ_of_W(Wmed)['S_per_Q'])
    for lab in ('band', 'broad'):
        W = Wmed + Wvec_spot[lab]           # nu = 1
        s = SQ_of_W(W)
        row['nu1_' + lab] = dict(W=W, **s, p_gas=p_of_SQ(s['S_per_Q']))
        # nu needed for bulk / Q=2e5 survival at e = 1e-6, and e_min at nu = 1 (vector part ~ 1/e; cascade inside R_med held)
        row['nu_crit_bulk_' + lab] = (W_thr['bulk_Q_inf'] - Wmed) / Wvec_spot[lab]
        row['nu_crit_Q2e5_' + lab] = (W_thr['Q_2e5'] - Wmed) / Wvec_spot[lab]
        row['e_min_bulk_nu1_' + lab] = 1e-6 * Wvec_spot[lab] / (W_thr['bulk_Q_inf'] - Wmed) if W_thr['bulk_Q_inf'] > Wmed else None
    cases[name] = row

# ---- the repaired point (trilinear, lambda = 0.5, (3,0.5)) with the transferred nu
m1r_key = 'lam0.5_b3_g0.5/trilinear'
r = rel[m1r_key]
Wmed = r['W0'] * (1 + r['R_med_e1e6'])
# iterate: p -> M -> H -> nu -> S/Q -> p
p = 9.1e-4
for it in range(6):
    M = 1 + 2 * p
    H = H_for(yield_tri05, M); Y = D_GEV / (M * MX); nu = nu_transfer(H, Y)
    s_band = SQ_of_W(Wmed + nu * Wvec_spot['band']); p_new = p_of_SQ(s_band['S_per_Q'])
    p = p_new
s_broad = SQ_of_W(Wmed + nu * Wvec_spot['broad'])
m1r = dict(release=r, W_med=Wmed, nu=nu, H_i_GeV=H, Y_target=Y, M_per_charge=1 + 2 * p,
           band=dict(**s_band, p_gas=p_of_SQ(s_band['S_per_Q'])), broad=dict(**s_broad, p_gas=p_of_SQ(s_broad['S_per_Q'])),
           margin_bulk_band=s_band['S_per_Q'] / sig['bulk_Q_inf'], margin_bulk_broad=s_broad['S_per_Q'] / sig['bulk_Q_inf'],
           nu_crit_bulk_band=(W_thr['bulk_Q_inf'] - Wmed) / Wvec_spot['band'],
           e_crit_bulk_band=1e-6 * nu * Wvec_spot['band'] / (W_thr['bulk_Q_inf'] - Wmed),
           e_crit_Q2e5_band=1e-6 * nu * Wvec_spot['band'] / (W_thr['Q_2e5'] - Wmed),
           e_crit_bulk_broad=1e-6 * nu * Wvec_spot['broad'] / (W_thr['bulk_Q_inf'] - Wmed))
# Planck boundary condition for the induced portals: release state not recorded by the skeptic run; vector side only [E]
H_pl = H_for(yield_tri05_planck, 1 + 2 * p); nu_pl = nu_transfer(H_pl, D_GEV / ((1 + 2 * p) * MX))
m1r['planck_bc'] = dict(H_i_GeV=H_pl, nu=nu_pl, SQ_band_same_release=SQ_of_W(Wmed + nu_pl * Wvec_spot['band'])['S_per_Q'],
                        note='release state (E_X, eps) of the Planck-BC run is not recorded; same release state assumed [E]')

# ---- H_i ledger
ledger = dict(
    M1_reference=dict(H=H_M1, note='Codex M1 (Option 1A, lambda=1), Q=2e5 object relic'),
    trilinear_lam05_Q2e5_relic=dict(H=H_for(yield_tri05, 0.9374), note='fx_holder FX11 (assumes the Q=2e5 relic)'),
    trilinear_lam05_freeX_relic=dict(H=H_for(yield_tri05, m1r['M_per_charge']), M=m1r['M_per_charge'],
                                     factor=(0.9374 / m1r['M_per_charge']) ** (2 / 3)),
    trilinear_lam05_freeX_planckBC=dict(H=H_pl),
    trilinear_lam05_Qmatter_relic=dict(H_range=[H_for(yield_tri05, 0.802), H_for(yield_tri05, 0.782)], M_range=[0.782, 0.802]),
    g_I_freeX=8.306742069498e-7 * H_for(yield_tri05, m1r['M_per_charge']) / H_M1)

out = dict(scalings=scal, spot=dict(**spot, c=C_spot, E_over_Ec=W0_spot), Wvec_spot_e1e6=Wvec_spot, sigma_min=sig, W_thresholds=W_thr,
           releases=rel, cases_e1e6_nu1=cases, p_scan_SQ=[dict(R=a, S_per_Q=b, p=c) for a, b, c in sq_scan],
           dlnp_dlnSQ=dlnp_dlnsq, m1r=m1r, H_ledger=ledger)
json.dump(out, open(OUT, 'w'), indent=1)

print('c constant: spot %.5f; fx_holder runs %.5f..%.5f' % (C_spot, min(v['c'] for v in rel.values()), max(v['c'] for v in rel.values())))
print('scalings', {k: round(v, 5) for k, v in scal.items()})
print('W thresholds (E_c units):', {k: round(v, 3) for k, v in W_thr.items()}, ' spot W0 = %.4f' % W0_spot)
for k, v in cases.items():
    print('%-70s eps %.4f E/Ec %6.2f Rmed %.3f  SQ0 %.3f  SQ(nu=1 band) %.3f  nu_crit(bulk,band) %.3f  e_min %.3g' % (
        k, v['eps'], v['E_over_Ec'], v['R_med_e1e6'], v['SQ_no_vectors'], v['nu1_band']['S_per_Q'], v['nu_crit_bulk_band'],
        v['e_min_bulk_nu1_band'] or float('nan')))
print('p scan:', [(round(a, 3), round(b, 4), '%.3g' % c) for a, b, c in sq_scan], ' slopes', [round(x, 2) for x in dlnp_dlnsq])
print('M1r:', json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if not isinstance(vv, (dict, list))}) for k, v in m1r.items() if k != 'release'}, indent=0))
print('H ledger:', json.dumps(ledger, indent=0))
