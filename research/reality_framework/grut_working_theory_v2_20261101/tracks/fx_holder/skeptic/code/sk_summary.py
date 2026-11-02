#!/usr/bin/env python3
"""fx_holder_skeptic (Claude; own code): collect the skeptic runs (../results/runs/*.json), compare with the track's
results (read-only from ../../results/*.json), and derive the new sensitivities.  Output: ../results/sk_summary.json
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent / 'results' / 'runs'
TR = HERE.parent.parent / 'results'
H_M1 = 134284273491.30815
LN_PL = math.log(2.435e18 / (30 * H_M1))

R = {p.stem: json.load(open(p)) for p in sorted(RUNS.glob('*.json')) if not p.stem.startswith('test')}
tr = json.load(open(TR / 'tri_runs.json'))['runs']
ls = {round(r['lam'], 3): r for r in json.load(open(TR / 'lam_scan.json'))['rows']}


def NX(k):
    return R[k]['NX_400_1000']


out = dict(note='skeptic = scipy DOP853 on the prescribed background, holder written as u = mu + (10i - mu) I(t)/I(0); track = rptri.c',
           reproduction={}, derived={})
pairs = {
    'orig_l1_b3': ('lam1_b3_g0.5', 'original'), '1A_l1_b3': ('lam1_b3_g0.5', 'option1A'), 'tri_l1_b3': ('lam1_b3_g0.5', 'trilinear'),
    'opp_l1_b3': ('lam1_b3_g0.5', 'trilinear_opposite_branch'), 'tri_l1_b3_s1.01': ('lam1_b3_g0.5', 'trilinear_slice_1.01'),
    'tri_l1_b3_s0.99': ('lam1_b3_g0.5', 'trilinear_slice_0.99'),
    'orig_l1_b5': ('lam1_b5_g1', 'original'), '1A_l1_b5': ('lam1_b5_g1', 'option1A'), 'tri_l1_b5': ('lam1_b5_g1', 'trilinear'),
    'orig_l05_b3': ('lam0.5_b3_g0.5', 'original'), '1A_l05_b3': ('lam0.5_b3_g0.5', 'option1A'), 'tri_l05_b3': ('lam0.5_b3_g0.5', 'trilinear'),
    'orig_l05_b5': ('lam0.5_b5_g1', 'original'), '1A_l05_b5': ('lam0.5_b5_g1', 'option1A'), 'tri_l05_b5': ('lam0.5_b5_g1', 'trilinear'),
}
for k, (a, b) in pairs.items():
    if k in R:
        t = tr[a][b]
        out['reproduction'][k] = dict(skeptic_NX=NX(k), track_NX=t['NX_400_1000'], rel_diff=NX(k) / t['NX_400_1000'] - 1,
                                      skeptic_Y_over_YQ=R[k]['Y_over_YQ_400_1000'], track_Y_over_YQ=t['yield_over_YQ_400_1000'],
                                      skeptic_H=R[k]['H_new_GeV'], track_H=t['H_new_GeV_400_1000'],
                                      skeptic_NYL=R[k]['NYL_400_1000'], track_NYL=t['NYL'])
for lam, k in ((0.25, 'tri_l025_b3'), (0.75, 'tri_l075_b3')):
    if k in R:
        out['reproduction'][k] = dict(skeptic_Y_over_YQ=R[k]['Y_over_YQ_400_1000'], track_Y_over_YQ=ls[lam]['trilinear']['yield_over_YQ'],
                                      rel_diff=R[k]['Y_over_YQ_400_1000'] / ls[lam]['trilinear']['yield_over_YQ'] - 1)
if 'orig_l025_b3' in R:
    out['reproduction']['orig_l025_b3'] = dict(skeptic_Y_over_YQ=R['orig_l025_b3']['Y_over_YQ_400_1000'], track_Y_over_YQ=ls[0.25]['original']['yield_over_YQ'])
out['reproduction']['max_abs_rel_diff_NX'] = max(abs(v['rel_diff']) for v in out['reproduction'].values() if isinstance(v, dict) and 'rel_diff' in v)
out['reproduction']['M1_convention_check_orig_l1_640_800_Y_over_YQ'] = R['orig_l1_b3']['Y_over_YQ_640_800']

d = out['derived']
d['ratios'] = {f'{h}_{l}': dict(tri_over_orig=NX(f'tri_{l}_{h}') / NX(f'orig_{l}_{h}'), A1_over_orig=NX(f'1A_{l}_{h}') / NX(f'orig_{l}_{h}'))
               for h in ('b3', 'b5') for l in ('l1', 'l05')}
d['opposite_over_reference'] = NX('opp_l1_b3') / NX('tri_l1_b3')
d['NYL_1A_over_tri'] = {k: R[f'1A_{k}']['NYL_400_1000'] / R[f'tri_{k}']['NYL_400_1000'] for k in ('l1_b3', 'l05_b3', 'l1_b5', 'l05_b5')}
d['NY_std_1A_over_tri'] = {k: R[f'1A_{k}']['NY_std_late'] / R[f'tri_{k}']['NY_std_late'] for k in ('l1_b3', 'l05_b3', 'l1_b5', 'l05_b5')}


def clock(kp, km, dl=math.log(1.01 / 0.99)):
    return math.log(abs(NX(kp)) / abs(NX(km))) / dl


d['clock'] = dict(tri_lam1=clock('tri_l1_b3_s1.01', 'tri_l1_b3_s0.99'), tri_lam05=clock('tri_l05_b3_s1.01', 'tri_l05_b3_s0.99'),
                  A1_lam05=clock('1A_l05_b3_s1.01', '1A_l05_b3_s0.99'), A1_lam1_track=1.3372986684790071,
                  H_exponent_tri_lam05=-2 / 3 * clock('tri_l05_b3_s1.01', 'tri_l05_b3_s0.99'),
                  note='d ln|N_X| / d ln(holder coupling) at fixed history, re-prepared; equivalently d ln|N_X|/d ln phi_I(slice) for the trilinear')
# portal counterterms cport (|u|^2 - mu^2)(|S|^2 + 2|Y|^2), one unit = 1/(4 pi^2) = one unit of ln Q of the one-loop running
cp = {}
for l in ('l1', 'l05'):
    n0, npl, nmi = NX(f'tri_{l}_b3'), NX(f'tri_{l}_b3_cpP'), NX(f'tri_{l}_b3_cpM')
    cp[f'tri_{l}'] = dict(dlnNX_per_unit_lnQ=(math.log(abs(npl)) - math.log(abs(nmi))) / 2, plus=npl / n0 - 1, minus=nmi / n0 - 1,
                          dlnH_per_unit_lnQ=-2 / 3 * (math.log(abs(npl)) - math.log(abs(nmi))) / 2)
    cp[f'1A_{l}'] = dict(dlnNX_per_unit_lnQ_one_sided=math.log(abs(NX(f'1A_{l}_b3_cpP')) / abs(NX(f'1A_{l}_b3'))))
cp['planck_bc'] = dict(units=-LN_PL, ln_MPl_over_MY=LN_PL,
                       tri_lam05=dict(NX_ratio=NX('tri_l05_b3_cpPlanck') / NX('tri_l05_b3'), Y_over_YQ=R['tri_l05_b3_cpPlanck']['Y_over_YQ_400_1000'],
                                      H_GeV=R['tri_l05_b3_cpPlanck']['H_new_GeV'], NYL=R['tri_l05_b3_cpPlanck']['NYL_400_1000']),
                       tri_lam05_reversed_sign=dict(Y_over_YQ=R['tri_l05_b3_cpPlanckRev']['Y_over_YQ_400_1000'], H_GeV=R['tri_l05_b3_cpPlanckRev']['H_new_GeV']),
                       tri_lam1=dict(NX_ratio=NX('tri_l1_b3_cpPlanck') / NX('tri_l1_b3'), Y_over_YQ=R['tri_l1_b3_cpPlanck']['Y_over_YQ_400_1000'],
                                     H_GeV=R['tri_l1_b3_cpPlanck']['H_new_GeV']),
                       note='c(M_Y) = -ln(M_P/M_Y)/(4 pi^2): the induced portals vanish at M_P and run to M_Y with the track\'s own beta functions')
d['portal_counterterms'] = cp
# history dependence at fixed Lagrangian vs at fixed prepared state
hd = {}
for l in ('l1', 'l05'):
    hd[l] = dict(H_b3=R[f'tri_{l}_b3']['H_new_GeV'], H_b5_retuned_gI=R[f'tri_{l}_b5']['H_new_GeV'], H_b5_fixed_gI=R[f'tri_{l}_b5_fixg']['H_new_GeV'],
                 u0_b5_fixed_gI=R[f'tri_{l}_b5_fixg']['u0'], prep_delta_b5_fixed=R[f'tri_{l}_b5_fixg']['prep_delta'],
                 H_1A_b3=R[f'1A_{l}_b3']['H_new_GeV'], H_1A_b5=R[f'1A_{l}_b5']['H_new_GeV'],
                 dlnH_tri_retuned=math.log(R[f'tri_{l}_b5']['H_new_GeV'] / R[f'tri_{l}_b3']['H_new_GeV']),
                 dlnH_tri_fixed=math.log(R[f'tri_{l}_b5_fixg']['H_new_GeV'] / R[f'tri_{l}_b3']['H_new_GeV']),
                 dlnH_1A=math.log(R[f'1A_{l}_b5']['H_new_GeV'] / R[f'1A_{l}_b3']['H_new_GeV']),
                 dlnH_unheld=math.log(R[f'orig_{l}_b5']['H_new_GeV'] / R[f'orig_{l}_b3']['H_new_GeV']))
d['history'] = hd
d['checks'] = dict(max_gauss=max(r['max_gauss'] for r in R.values()), min_corr_NX_6NS=min(r['corr_NX_6NS'] for r in R.values()))
json.dump(out, open(HERE.parent / 'results' / 'sk_summary.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
