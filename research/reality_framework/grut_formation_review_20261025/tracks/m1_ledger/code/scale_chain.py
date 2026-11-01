#!/usr/bin/env python3
"""m1_ledger (Claude): the mass-normalisation scale chain.  Item 2 of the brief.

Part A [D]: reconstruct Codex's chain from the numbers quoted in d4747cb MATCHED_RESULTS.json /
DEFORMED_MATCHED_RESULTS.json and 1fcc6b7 ABUNDANCE_LEDGER.json; test the H_i^{3/2} exponent implied by its own rows.
Part B [C, own engine]: run our rp_reheating C engine (tracks/rp_reheating/code/rp.py, Claude's own code) on the
(3, 0.5) and (5, 1) histories at three coordinated scales lambda = H_i/H_old with all heavy scales co-scaled
(dimensionless M_Y/H_i, h/H_i, mu/H_i, m_I/H_i, Gamma_I/H_i fixed; only eta^2 = (H_i/M_red)^2 and the unit
conversion change) and measure d ln Y / d ln H_i directly.
Usage: OMP_NUM_THREADS=1 python3 -B scale_chain.py            (single process, ~1 min)
Output: ../results/scale_chain.json
"""
import json, os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACKS = HERE.parents[1]
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/m1_ledger')
os.environ.setdefault('RP_SCRATCH', str(SP / 'rp'))
sys.dont_write_bytecode = True
sys.path.insert(0, str(TRACKS / 'rp_reheating' / 'code'))
import rp                                         # Claude's own engine (rp_reheating track)

OUT = HERE.parent / 'results'

# ---------------- Part A: Codex's own numbers [Q] -----------------------------------------------------------------
YQ_old = 4.8294009076988175e-9
D = 4.4e-10
eps_old = 0.0853226872719292            # GeV, old object E/Q      [Q MATCHED_RESULTS]
eps_new = 0.08662416173210513           # GeV, M1 object E/Q       [Q DEFORMED_MATCHED_RESULTS]
H_old = 133349619398.33566              # GeV                      [Q]
codex = {
    'b3_g0.5_baseline': dict(scale=1.0, Y=5.026475516056739e-09, yield_over_YQ=1.040807257903098, T=7702499481901.287, H=133349619398.33566),
    'b3_g0.5_matched_oldobj': dict(scale=1.0172235041614452, Y=5.156893366462146e-09, yield_over_YQ=1.0678122328260746, T=7768548315222.136, H=135646367122.97003),
    'b5_g1_baseline': dict(scale=1.0, Y=5.033450222058231e-09, yield_over_YQ=1.0422514755472314, T=7700753422824.31, H=133349619398.33566),
    'M1_b3_g0.5_matched': dict(scale=1.0070090495735164, Y=5.079414232732486e-09, yield_over_YQ=1.051769014379218, T=7729445947446.704, H=134284273491.30815),
    'M1_b5_g1_matched': dict(scale=1.006078579112091, Y=5.079414232730545e-09, yield_over_YQ=1.051769014378816, T=7724122783045.604, H=134160195609.41566),
}
A = {}
A['targets'] = dict(YQ_old=YQ_old, target_old_object=D / eps_old, target_M1_object=D / eps_new,
                    old_object_fraction_of_D=YQ_old * eps_old / D)
# yield_over_YQ is relative to the OLD supplied yield 4.8294e-9 in every row (check)
A['yield_over_YQ_is_relative_to_YQ_old'] = {k: v['Y'] / YQ_old / v['yield_over_YQ'] for k, v in codex.items()}
# exponent implied by each matched row relative to its own baseline
pairs = [('b3_g0.5_matched_oldobj', 'b3_g0.5_baseline'), ('M1_b3_g0.5_matched', 'b3_g0.5_baseline'), ('M1_b5_g1_matched', 'b5_g1_baseline')]
A['implied_exponents'] = {}
for m, b in pairs:
    ry = codex[m]['Y'] / codex[b]['Y']; rs = codex[m]['scale']; rT = codex[m]['T'] / codex[b]['T']; rH = codex[m]['H'] / codex[b]['H']
    A['implied_exponents'][m] = dict(Y_ratio=ry, scale=rs, H_ratio=rH, exponent_Y=np.log(ry) / np.log(rs), exponent_T=np.log(rT) / np.log(rs),
                                     scale_predicted_from_3_2=ry ** (2 / 3))
# the naive reading in the brief: 5.07941423/4.8294 = 1.0518 would need +3.42% if the OLD POINT delivered 4.8294e-9
A['naive_reading'] = dict(ratio_to_YQ_old=(D / eps_new) / YQ_old, scale_if_old_point_gave_YQ_old=((D / eps_new) / YQ_old) ** (2 / 3),
                          resolution='the old heavy point on the (3,0.5) history already yields 1.040807 x 4.8294e-9 = 5.0265e-9, '
                                     'so the needed factor is 5.07941/5.02648 = 1.010532 = 1.00700905^{3/2}; the 1.722% / 1.628% repairs '
                                     'are the OLD-object targets 5.15689e-9 = 1.0678 x 4.8294e-9 (3,0.5)/(5,1), superseded by M1')
A['chain'] = [
    '1. supplied Y_Q,old = 4.8294e-9 with v-calibrated mass per net charge D/Y_Q,old = 0.98592 m_X (old target, not an object E/Q)',
    '2. old object E/Q = 0.923316 m_X delivers 0.936494 D at Y_Q,old (1fcc6b7); needs x1.0678 -> 5.15689e-9',
    '3. (3,0.5) history at H_old = 1.3335e11 actually yields 1.040807 Y_Q,old = 5.02648e-9 (ours at tau=80: 1.04542; see Part B); (5,1): 1.042251',
    '4. old-object repair: scale = (5.15689/5.02648)^{2/3} = 1.017224 [1.722%]; (5,1): (5.15689/5.03345)^{2/3} = 1.016277 [1.628%]',
    '5. M1 object E/Q = 0.9374 m_X -> target 5.07941e-9; scale = (5.07941/5.02648)^{2/3} = 1.0070090 [0.7009%]; (5,1): 1.0060786 [0.6079%]',
    '6. Y propto H_i^{3/2}: n_Q = N_X H_i^3/a^3 at fixed dimensionless N_X, a; rho_r = r M_red^2 H_i^2 -> T propto (H_i M_red)^{1/2}, '
    's propto (H_i M_red)^{3/2}; dimensionless history depends on the scale only through eta^2 = (H_i/M_red)^2 (hidden back-reaction ~1e-9)',
]

# ---------------- Part B: own engine ---------------------------------------------------------------------------------
def yield_at(res, lam):
    """Recompute Y from the engine's dimensionless end state at H_i = lam * H_old (rp.run's own conversion assumes lam = 1)."""
    H = lam * rp.H0
    a, Rr = res['final_a'], res['Rr_end']
    rr = Rr / a ** 4
    rho = rr * rp.MRED ** 2 * H ** 2
    T = (30 * rho / (np.pi ** 2 * rp.GSTAR)) ** 0.25
    s = 2 * np.pi ** 2 / 45 * rp.GSTAR * T ** 3
    return abs(res['late']['NX']) * H ** 3 / (a ** 3 * s), T

B = dict(engine='tracks/rp_reheating/code/rpreheat.c + rp.py (Claude, own code); rtol 1e-10, atol 1e-13; late = tau-weighted last 20%',
         runs=[])
lams = [1.0, 1.0070090495735164, 1.0172235041614452, 1.05]
t0 = time.time()
for b, g in ((3.0, 0.5), (5.0, 1.0)):
    for end in (80.0, 800.0):
        if end == 80.0 and b != 3.0:
            continue
        for lam in (lams if end == 800.0 else [1.0]):
            eta2 = (lam * rp.H0 / rp.MRED) ** 2
            res, S = rp.run('m1_b%g_g%g_l%.5f_t%g' % (b, g, lam, end), mode=0, b=b, gam=g, end=end, eta2=eta2, nsamp=2000, save=False)
            Y, T = yield_at(res, lam)
            B['runs'].append(dict(b=b, gamma=g, end=end, lam=lam, eta2=eta2, Y=Y, Y_over_YQ_old=Y / YQ_old, Tvis_end_GeV=T,
                                  NX=res['late']['NX'], eps=res['late']['eps'], NYL=res['late']['NYL'], NYH=res['late']['NYH'],
                                  final_a=res['final_a'], Rr_end=res['Rr_end'], dilution=res['dilution'],
                                  half_window_charge_change=res['half_window_charge_change'], checks=res['checks'], seconds=res['seconds']))
            print('b %g g %g end %g lam %.6f  Y %.6e  Y/YQold %.6f  NX %.4f eps %.5f a %.4f  %.1fs' %
                  (b, g, end, lam, Y, Y / YQ_old, res['late']['NX'], res['late']['eps'], res['final_a'], res['seconds']), flush=True)

def find(b, end, lam):
    return next(r for r in B['runs'] if r['b'] == b and r['end'] == end and abs(r['lam'] - lam) < 1e-12)

B['exponents'] = {}
for b in (3.0, 5.0):
    base = find(b, 800.0, 1.0)
    for lam in lams[1:]:
        r = find(b, 800.0, lam)
        B['exponents']['b%g_lam%.5f' % (b, lam)] = dict(Y_ratio=r['Y'] / base['Y'], exponent=np.log(r['Y'] / base['Y']) / np.log(lam),
                                                       NX_ratio=r['NX'] / base['NX'], a_ratio=r['final_a'] / base['final_a'],
                                                       Rr_ratio=r['Rr_end'] / base['Rr_end'])
b3 = find(3.0, 800.0, 1.0); b5 = find(5.0, 800.0, 1.0); b3s = find(3.0, 80.0, 1.0)
B['comparison_with_codex'] = dict(
    b3_g0_5_tau800_ours=b3['Y'], codex=codex['b3_g0.5_baseline']['Y'], rel=b3['Y'] / codex['b3_g0.5_baseline']['Y'] - 1,
    b5_g1_tau800_ours=b5['Y'], codex_b5=codex['b5_g1_baseline']['Y'], rel_b5=b5['Y'] / codex['b5_g1_baseline']['Y'] - 1,
    b3_g0_5_tau80_ours=b3s['Y'] / YQ_old, rp_reheating_extra_json=1.0454205564913566,
    tau80_over_tau800=b3s['Y'] / b3['Y'], NX_tau80_over_tau800=b3s['NX'] / b3['NX'],
    M1_scale_from_our_tau800=(D / eps_new / b3['Y']) ** (2 / 3), codex_scale=1.0070090495735164,
    M1_scale_b5_from_our_tau800=(D / eps_new / b5['Y']) ** (2 / 3), codex_scale_b5=1.006078579112091,
    oldobj_scale_from_our_tau800=(D / eps_old / b3['Y']) ** (2 / 3), codex_oldobj_scale=1.0172235041614452)
# M1 matched row in our engine: run at Codex's scale and compare Y to the target
m1 = find(3.0, 800.0, 1.0070090495735164)
B['comparison_with_codex']['our_Y_at_codex_M1_scale'] = m1['Y']
B['comparison_with_codex']['our_Y_at_codex_M1_scale_over_target'] = m1['Y'] / (D / eps_new)
B['seconds_total'] = time.time() - t0

OUT.mkdir(exist_ok=True)
(OUT / 'scale_chain.json').write_text(json.dumps(dict(method=__doc__, partA_codex_numbers=A, partB_own_engine=B), indent=1, default=float))
print(json.dumps(A['implied_exponents'], indent=1, default=float))
print(json.dumps(B['exponents'], indent=1, default=float))
print(json.dumps(B['comparison_with_codex'], indent=1, default=float))
