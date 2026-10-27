"""Skeptic: recompute the frozen line, tension and instability scale with an independent pipeline (sk_rge.py),
and quantify the sensitivity of P3a's line to a Planck-suppressed |H|^6 operator (GRUT's L_higher, ACTION.md).

Run: python3 -B sk_line.py   (~1 min).  Writes ../results/sk_line.json.
"""
import json
import os

import numpy as np
from scipy.stats import norm

import sk_rge as R

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'sk_line.json')

MH, SMH = 125.13, 0.11            # PDG 2026 S126M (pdglive, read 2026-10-27 by this skeptic)
AS, SAS = 0.1180, 0.0009          # PDG 2026 QCD review (1 June 2026), world average (Fig. 9.4 / sec. 9.4.8)
MT_D, SMT_D = 172.60, 0.27        # PDG 2026 Q007TP (OUR AVERAGE, S = 1.5)
MT_X, SMT_X = 172.1, 0.6          # PDG 2026 Q007TP4 (cross-section pole mass)
MT_AC, SMT_AC = 172.52, 0.33      # ATLAS+CMS 2402.08713
AS_FLAG, SAS_FLAG = 0.1183, 0.0007
AS_EW, SAS_EW = 0.1221, 0.0027    # PDG 2026 QCD review Table 9.1, 'electroweak' subfield (unweighted)

res = {'inputs_checked_at_primary_source': {
    'Mh': 'PDG 2026 S126M 125.13 +- 0.11 (S=1.5) [Q pdglive]',
    'Mt_direct': 'PDG 2026 Q007TP 172.60 +- 0.27 (S=1.5) [Q pdglive]',
    'Mt_xsec': 'PDG 2026 Q007TP4 172.1 +- 0.6 [Q pdglive]',
    'alpha_s': 'PDG 2026 QCD review (dated 1 June 2026, rev. Aug 2025) 0.1180 +- 0.0009; FLAG 2024 0.1183 +- 0.0007 '
               '(eq. 9.23); electroweak subfield 0.1221 +- 0.0027 (Table 9.1) [Q pdftotext]'}}

# ---------------- validation against the primary sources
sol = R.run(173.34, 125.15, 0.1184, 'B13')
y = sol.sol(np.log(R.MPL ** 2))
res['validation'] = {
    'Buttazzo central, lambda(M_Pl)': [float(y[4]), 'paper: -0.0143 (TeX l.1296)'],
    'Buttazzo central, y_t(M_Pl)': [float(np.sqrt(y[3])), 'paper: 0.3825 (TeX l.1294)'],
    'B15 at Bednyakov inputs (125.7, 0.1185)': [R.mt_crit(125.7, 0.1185, 'B15')[0], 'paper: 171.44'],
    'B13 at Buttazzo inputs (125.15, 0.1184)': [R.mt_crit(125.15, 0.1184, 'B13')[0],
                                                'paper eq. inverted: 171.115; paper quoted: 171.53'],
}

# ---------------- the line at PDG 2026
line = {}
for sch in ('B13', 'B15'):
    for L in (2, 3):
        mt, mus = R.mt_crit(MH, AS, sch, L)
        line[f'{sch}_loop{L}'] = {'Mt_crit': mt, 'mu_star': mus}
    m0 = line[f'{sch}_loop3']['Mt_crit']
    line[f'{sch}_dMt_dMh'] = (R.mt_crit(MH + 0.5, AS, sch)[0] - R.mt_crit(MH - 0.5, AS, sch)[0]) / 1.0
    line[f'{sch}_dMt_per_0.0009_as'] = (R.mt_crit(MH, AS + 0.0009, sch)[0] - R.mt_crit(MH, AS - 0.0009, sch)[0]) / 2
    line[f'{sch}_Mpl_reduced_cut'] = R.mt_crit(MH, AS, sch, 3, mu_hi=R.MPLR)[0]
# matching-error propagation (Buttazzo's own +-0.00030 lambda, +-0.00050 y_t)
mtl = R.mt_crit(MH, AS, 'B13', dlam=3e-4)[0] - line['B13_loop3']['Mt_crit']
mty = R.mt_crit(MH, AS, 'B13', dyt=5e-4)[0] - line['B13_loop3']['Mt_crit']
line['B13_dMt_from_dlam_3e-4'] = mtl
line['B13_dMt_from_dyt_5e-4'] = mty
line['B13_matching_error_quadrature'] = float(np.hypot(mtl, mty))
line['track_SMDR_value'] = 170.997
line['spread_B13_B15_SMDR'] = [line['B13_loop3']['Mt_crit'], line['B15_loop3']['Mt_crit'], 170.997]
res['line_PDG2026'] = line

# ---------------- tension (same recipe as the track, with each matching)
def ztab(mtc, sl_h, sl_a, sth):
    out = {}
    cases = [('direct', MT_D, SMT_D, AS, SAS, 0.0), ('direct+0.3MC', MT_D, SMT_D, AS, SAS, 0.3),
             ('direct+0.6MC', MT_D, SMT_D, AS, SAS, 0.6), ('ATLAS+CMS', MT_AC, SMT_AC, AS, SAS, 0.0),
             ('xsec', MT_X, SMT_X, AS, SAS, 0.0), ('direct, FLAG as', MT_D, SMT_D, AS_FLAG, SAS_FLAG, 0.0),
             ('direct, EW-fit as', MT_D, SMT_D, AS_EW, SAS_EW, 0.0)]
    for name, mt, smt, a, sa, mc in cases:
        crit = mtc + sl_a * (a - AS) / 0.0009
        s = np.sqrt(smt ** 2 + (sl_h * SMH) ** 2 + (sl_a * sa / 0.0009) ** 2 + sth ** 2 + 0.032 ** 2 + mc ** 2)
        z = (mt - crit) / s
        out[name] = {'Delta': mt - crit, 'sigma': s, 'z': z, 'p_one_sided': float(norm.sf(z))}
    return out


res['tension'] = {}
for sch in ('B13', 'B15'):
    res['tension'][sch] = ztab(line[f'{sch}_loop3']['Mt_crit'], line[f'{sch}_dMt_dMh'],
                               line[f'{sch}_dMt_per_0.0009_as'], 0.113)
# Bednyakov's own asymmetric scale-variation error on M_t^cri (+0.17/-0.36, prl.tex l.375/589): upward 0.17
res['tension']['B15_with_Bednyakov_upward_th_0.17'] = ztab(line['B15_loop3']['Mt_crit'], line['B15_dMt_dMh'],
                                                         line['B15_dMt_per_0.0009_as'], 0.17)

# ---------------- instability scale Lambda_lambda (MSbar lambda = 0)
sc = {}
for sch in ('B13', 'B15'):
    sc[sch] = {f'{mt}': R.lambda_zero(mt, MH, AS, sch) for mt in (172.33, 172.60, 172.87, 172.1)}
res['Lambda_lambda'] = sc

# ---------------- Planck-suppressed operator: V = h^4/4 [lambda(h) + c6 h^2/(2 Mbar^2)], cut at h <= cut
c6s = [-1.0, -0.1, -0.01, -0.001, 0.0, 0.01, 0.1, 0.3, 1.0, 3.0, 4 * np.pi]
ops = {}
for cut_name, cut in (('Mbar_Pl', R.MPLR), ('M_Pl', R.MPL)):
    ops[cut_name] = {}
    for c6 in c6s:
        r = R.mt_crit(MH, AS, 'B15', 3, mu_hi=cut, c6=c6, lo=150.0, hi=180.0)
        ops[cut_name][f'{c6:.4g}'] = None if r is None else {'Mt_crit': r[0], 'h_star': r[1]}
res['c6_operator_B15'] = ops
res['c6_note'] = ('V_eff ~ h^4/4 [lambda(mu=h) + c6 h^2/(2 Mbar^2)], c6|H|^6/Mbar^2 with |H|^2 = h^2/2; RG-improved '
                  'tree level, operator anomalous dimension neglected [E]. P3a criterion: bracket >= 0 for h <= cut.')
# for the central metastable point: how large must c6 be to remove the deeper region?
sol = R.run(MT_D, MH, AS, 'B15', 3)
mus = np.logspace(10, np.log10(R.MPLR), 400)
lam = R.lam_of_mu(sol, mus)
need = np.max(np.where(lam < 0, -lam * 2 * R.MPLR ** 2 / mus ** 2, 0))
res['c6_needed_to_stabilise_PDG_central'] = float(need)

json.dump(res, open(OUT, 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
