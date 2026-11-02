#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude): hidden-heat re-ledger of the inf_potential benchmarks, S/Q(e).

Inputs (read-only): dq_heat/results/heat_table.json (M1r ledger rows per e), track results exit_factor.json,
yy_harmonics.json, engine_test.json; skeptic results sk_engine.json (our re-run of rpinf.c MODE 4 from S1).
S/Q closed form of dq_entropy: S/Q = 0.190255 (W - 0.97494)^0.75013 - 0.00428, W = E_total/E_c, E_c = (N_X/c)^{4/3},
c = 6.977568 [Q dq_entropy / wt_critic].  Ledger of dq_heat [Q]: each Y_L parent deposits m_L <a_L>, <a_L> =
(pi K/2 Gamma_L)^{1/2}; each Y_H (and 2 per sigma) also yields one A_L that decays with the primordial vectors and deposits
m_A <a_A> (∝ e^{-1/2} K^{1/2}); primordial vectors N_A ∝ K^{3/2} f_eff, deposit ∝ K^{1/2} e^{-1/2}.
Scenarios:
 A  track R6/R7 as stated (M1r release, exit factor only, M1r K for the amplitude)          -> reproduction
 B  A with the K^2 amplitude of R_vec (and K^1/2 of the Y deposits) applied consistently with the k-shift, K = K_end
 C  Send (R7) + the track's own R_YY(Send, Gamma = 0.5) on the M1r release
 D  S1 (R6) + R_YY(S1) on the M1r release (the track's own 'S/Q 7-8' statement)
 F  the track's rough re-ledger on its smoke release (reproduction; omits N_YH, N_sigma)
 E  the same smoke release WITH the homogeneous N_YH and N_sigma the engine reports (end normalisation)
 E0 E without R_YY (as if the harmonic/fundamental Y-pair channel were absent)
Output ../results/sk_reledger.json.  Labels: [C] arithmetic on [C]/[E] inputs; the ledger itself is [E] (dq_heat).
"""
import json, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]; WT = TRK.parent
RES = TRK / 'results'; OUT = HERE.parent / 'results'
C_HART = 6.977568
SQTHR = 3.4716191256166837


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


hb = json.load(open(WT / 'dq_heat' / 'results' / 'heat_table.json'))['budgets']['portals_vanish_at_MY_benchmark']
EXM, NXM, KM, W0 = hb['EX'], hb['NX'], hb['K'], hb['W0']
ML, GL = hb['mediators']['m_L'], hb['mediators']['Gamma_L']
rows = {r['e']: r for r in hb['rows']}
ES = sorted(rows)
ef = json.load(open(RES / 'exit_factor.json'))['results']
yy = json.load(open(RES / 'yy_harmonics.json'))['runs']
eng = json.load(open(RES / 'engine_test.json'))['smoke']
ske = json.load(open(OUT / 'sk_engine.json'))['runs']
KB = {n: json.load(open(RES / ('bench_%s.json' % n)))['by_Hi'][0]['K_end_gamma05'] for n in ('HT6t', 'T1', 'E1')}


def depYL(K):
    return ML * math.sqrt(math.pi * K / (2 * GL))


def depA(e, K):
    return rows[e]['deposit_per_vector_EX'] * EXM * math.sqrt(K / KM)


def feff(name, e, conv):
    return [x for x in ef['%s_%s' % (name, conv)]['by_e'] if x['e'] == e][0]['f_eff']


def yyN(name, lab, gam=0.5):
    r = [x for x in yy if x['benchmark'] == name and x['gamma_over_Hi'] == gam][0]['res'][lab]
    t = r['tail_NY_estimate'] or 0.0
    return r['N_YL'] + r['N_YH'] + t, r['N_YH'] + 0.5 * t      # (Y_L parents, A_L from Y_H)


def e_crit(fun):
    vals = [(e, fun(e)) for e in ES]
    for (e1, s1), (e2, s2) in zip(vals, vals[1:]):
        if (s1 - SQTHR) * (s2 - SQTHR) < 0:
            x = math.log(e1) + (math.log(e2) - math.log(e1)) * (s1 - SQTHR) / (s1 - s2)
            return math.exp(x)
    return None


out = dict(method=__doc__, benchmarks={})
for n in ('HT6t', 'T1', 'E1'):
    K = KB[n]; res = {}

    def mlr(e, conv, kamp=False, add_yy=None):
        r = rows[e]; f = feff(n, e, conv)
        kv = (K / KM) ** 2 if kamp else 1.0; ky = math.sqrt(K / KM) if kamp else 1.0
        Rv = (r['R_vec'] - r['R_vec_k_le_1'] + f * r['R_vec_k_le_1']) * kv
        Rt = Rv + (r['R_Y_hom'] + r['R_Y_k'] + r['R_cascade']) * ky
        if add_yy:
            L, A = yyN(n, add_yy)
            Rt += (L * depYL(KM) + A * depA(e, KM)) / EXM
        return SQ(W0 * (1 + Rt))
    res['A_S1'] = {str(e): mlr(e, 'S1') for e in ES}; res['A_S1_ecrit'] = e_crit(lambda e: mlr(e, 'S1'))
    res['A_Send'] = {str(e): mlr(e, 'Send') for e in ES}; res['A_Send_ecrit'] = e_crit(lambda e: mlr(e, 'Send'))
    res['B_S1_Kamp'] = {str(e): mlr(e, 'S1', kamp=True) for e in ES}; res['B_S1_Kamp_ecrit'] = e_crit(lambda e: mlr(e, 'S1', kamp=True))
    res['C_Send_plusYY'] = {str(e): mlr(e, 'Send', add_yy='Send') for e in ES}; res['C_ecrit'] = e_crit(lambda e: mlr(e, 'Send', add_yy='Send'))
    res['D_S1_plusYY'] = {str(e): mlr(e, 'S1', add_yy='S1') for e in ES}; res['D_ecrit'] = e_crit(lambda e: mlr(e, 'S1', add_yy='S1'))

    # smoke release (track's numbers) -- F reproduces the track's rough re-ledger at e = 1e-6
    sm = eng[n]; NX, EX = abs(sm['NX_endnorm']), sm['EX_endnorm']
    Ec = (NX / C_HART) ** (4 / 3)
    sk = ske['%s_t40' % n]['late_mean'] if '%s_t40' % n in ske else ske['%s_t60' % n]['late_mean']
    NYL, NYH, NS = sm['NYL'] / math.e ** 3, sm['NYH'] / math.e ** 3, sk['Nsig']

    def smoke(e, with_YH=True, with_YY=True, yl=1.0):
        r = rows[e]; f = feff(n, e, 'S1')
        NA = (r['N_A_k_le_1'] * f + r['N_A_k_ge_1']) * (K / KM) ** 1.5
        D = NA * depA(e, K)
        D += NYL * depYL(K) * yl
        if with_YH:
            D += NYH * depYL(K) * yl + (NYH + 2 * NS) * depA(e, K)
        if with_YY:
            L, A = yyN(n, 'S1')
            D += L * depYL(K) * yl + A * depA(e, K)
        return SQ((EX + D) / Ec), D / EX
    Fv = smoke(1e-6, with_YH=False)
    res['F_track_rough_e1e-6'] = dict(SQ=Fv[0], R_own=Fv[1], note='track summary.json gives %.3f; our F reproduces the structure' %
                                      json.load(open(RES / 'summary.json'))['benchmarks'][n]['rough_reledger_E']['SQ_all'])
    res['E_full'] = {str(e): dict(zip(('SQ', 'R_own_EX'), smoke(e))) for e in ES}
    res['E0_noYY'] = {str(e): dict(zip(('SQ', 'R_own_EX'), smoke(e, with_YY=False))) for e in ES}
    res['E1_noYY_noYLdeposit'] = {str(e): dict(zip(('SQ', 'R_own_EX'), smoke(e, with_YY=False, yl=0.0))) for e in ES}
    # e -> infinity: vectors and A_L deposits vanish; only Y_L parents
    Dinf = (NYL + NYH) * depYL(K) + yyN(n, 'S1')[0] * depYL(K)
    res['E_full_e_inf'] = SQ((EX + Dinf) / Ec)
    res['E0_noYY_e_inf'] = SQ((EX + (NYL + NYH) * depYL(K)) / Ec)
    res['inputs'] = dict(K_end=K, NX_end=NX, EX_end=EX, Ec=Ec, W_X=EX / Ec, NYL_end=NYL, NYH_end=NYH, Nsig_end=NS,
                         NYH_over_M1r=NYH / hb['NYH_hom'], NYL_over_M1r=NYL / hb['NYL_hom'], Nsig_over_M1r=NS / hb['Nsig_hom'])
    out['benchmarks'][n] = res
    print('%s: A(S1,1e-6)=%.3f ecrit=%s | B(K^2,S1)=%.3f ecrit=%s | A(Send)=%.3f | C(Send+YY)=%.1f | D(S1+YY)=%.2f | F=%.2f (R_own=%.1f) | E(1e-6)=%.1f E(1e-5)=%.1f E(inf)=%.1f | E0(1e-6)=%.1f E0(inf)=%.1f | E1(noYY,noYLdep,1e-6)=%.1f (1e-5)=%.1f | NYH/M1r=%.0f' % (
        n, res['A_S1']['1e-06'], res['A_S1_ecrit'], res['B_S1_Kamp']['1e-06'], res['B_S1_Kamp_ecrit'], res['A_Send']['1e-06'],
        res['C_Send_plusYY']['1e-06'], res['D_S1_plusYY']['1e-06'], Fv[0], Fv[1], res['E_full']['1e-06']['SQ'],
        res['E_full']['1e-05']['SQ'], res['E_full_e_inf'], res['E0_noYY']['1e-06']['SQ'], res['E0_noYY_e_inf'],
        res['E1_noYY_noYLdeposit']['1e-06']['SQ'], res['E1_noYY_noYLdeposit']['1e-05']['SQ'], res['inputs']['NYH_over_M1r']))
json.dump(out, open(OUT / 'sk_reledger.json', 'w'), indent=1)
