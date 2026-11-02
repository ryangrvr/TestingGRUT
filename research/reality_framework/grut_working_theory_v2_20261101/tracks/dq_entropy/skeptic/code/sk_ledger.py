#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): re-derivation of dq_entropy's release invariants, debris ledger, R_crit, the two
readings of the provisional heat, the H_i ledger, nu, and the branch map with dq_heat's direct R(e), using
  * the independent late-time frozen release of sk_late.c (results/sk_late.json),
  * the independent Bessel-series Hartree solver sk_hartree_bessel.py (S/Q at given W),
  * dq_entropy's own increments (+0.004..+0.0068) only where stated.
Every formula is re-typed here from its definition (not imported from dq_entropy's code).
Usage: sk_ledger.py ../results   (reads sk_late.json, sk_hartree_bessel.json, sk_homfunc.json, sk_portal.json; writes sk_ledger.json)"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_hartree_bessel as HB

RES = Path(sys.argv[1])
TR = Path(__file__).resolve().parents[3]
DQ = json.load(open(TR / 'dq_entropy' / 'results' / 'summary.json'))
FR = json.load(open(TR / 'dq_entropy' / 'results' / 'freeze.json'))
HEAT = json.load(open(TR / 'dq_heat' / 'results' / 'summary.json'))
LATE = json.load(open(RES / 'sk_late.json'))
HBJ = json.load(open(RES / 'sk_hartree_bessel.json'))

LX = 1e-3
CC = 2 * math.sqrt(2 * LX) * (3 * LX) ** -0.75
H_M1, YQ, D_GEV, MX = 1.3428427349130815e11, 5.07941423e-9, 4.4e-10, 0.0924089637375362
SIG = dict(bulk=3.4716191256166837, Q_2e5=2.9283778774213727, binding=2.5530624453282833)


def e_over_ec(N, E):
    return E / (N / CC) ** (4 / 3)


def SQ_of_W(W, release=(HB.NX_M1R, HB.EX_M1R)):
    N, E = release
    w0 = e_over_ec(N, E)
    return HB.solve(N, E * W / w0)['S'] / N


def main():
    out = {}
    # ---------------------------------------------------------------- 1. frozen release, independent integrator
    fz = LATE['frozen']
    NXs, EXs = abs(fz['NX']), fz['Econf']
    eps = NXs / (CC * EXs ** 0.75)
    out['release'] = dict(N_X_sk=NXs, E_X_sk=EXs, eps_sk=eps, E_over_Ec_sk=eps ** (-4 / 3), E_per_charge=EXs / NXs,
                          dq=dict(N_X=8047.840302004927, E_X=188737.70931829128, eps=0.1273741466807691, E_over_Ec=15.603605078079994),
                          rel_dN=NXs / 8047.840302004927 - 1, rel_dE=EXs / 188737.70931829128 - 1,
                          c_definition=CC)
    out['release']['E_over_Ec_rel_diff'] = out['release']['E_over_Ec_sk'] / 15.603605078079994 - 1
    # yield at the M1 H_i: Y/Y_Q scales with |N_X| at the same end-of-run entropy (fx_holder convention)
    y_dq = FR['base_1e5']['yield_over_YQ_frozen']
    y_sk = y_dq * NXs / 8047.840302004927
    out['yield'] = dict(Y_over_YQ_sk=y_sk, Y_over_YQ_dq=y_dq)
    # ---------------------------------------------------------------- 2. debris ledger from the definitions
    K = 1.9858052754037523
    M, h, mu, kap = 30.0, 40.0, 961402368973.181 / H_M1, 0.0034586847211042523
    M = 4028528204739.244 / H_M1; h = 5371370939652.326 / H_M1
    mL = math.sqrt(M * M - math.sqrt(2) * mu * h)
    GL = 3 * kap ** 2 * mL / (256 * math.pi ** 3)          # Y_L -> XXX + XbarXbarXbar (derived in SKEPTIC.md)
    aL = math.sqrt(math.pi * K / (2 * GL))                 # <a> at decay, radiation a^2 H = K
    deb = FR['base_1e5']['debris_late']
    RY = (deb['NYL'] + deb['NYH']) * mL * aL / EXs
    e = 1e-6; mA = 6 * e * h; GA = e ** 2 * mA / (48 * math.pi); aA = math.sqrt(math.pi * K / (2 * GA))
    Rc = (deb['NYH'] + 2 * deb['Nsig']) * mA * aA / EXs
    out['ledger'] = dict(m_L=mL, Gamma_L=GL, a_mean_L=aL, R_Y=RY, R_casc_plus_Sigma_e1e6=Rc,
                         dq=dict(R_Y=0.3785839760987034, R_casc_plus_Sigma=0.018343106488297886 + 7.078020914503917e-05),
                         dq_heat=dict(R_Y_hom=0.37883994745548166, R_cascade=0.01855662421470789,
                                      NYH_window=0.010377131381878653, NYH_late_dq=deb['NYH']),
                         note='dq_entropy takes N_YH, N_sigma from a >= 315 windows; Y_H decays at <a> = 9.4 and sigma at a = 2.2, '
                              'so the ledger should use their values at those epochs (not recorded); effect <= 1e-3 in R at e = 1e-6')
    # ---------------------------------------------------------------- 3. R_crit from the independent Hartree W_crit
    w_sk = out['release']['E_over_Ec_sk']
    out['R_crit'] = {k: dict(W_crit_sk=v['W_crit'], R_crit_dq_release=v['R_crit'], R_crit_sk_release=v['W_crit'] / w_sk - 1,
                             dq=DQ['benchmark']['R_crit'].get({'bulk': 'bulk_Q_inf', 'Q_2e5': 'Q_2e5', 'binding': 'Q_bind_71263'}.get(k, k)))
                     for k, v in HBJ['thresholds'].items() if k in ('bulk', 'Q_2e5', 'binding')}
    out['W_crit_dq_heat'] = HEAT['e_crit']['benchmark']['bulk_Q_inf']['W_threshold']
    # ---------------------------------------------------------------- 4. dq_heat's direct R(e) -> S/Q (independent Hartree)
    # dq_heat's R is in units of its 400-1000 E_X = 188706.29 at N_X = 8045.485 (window); convert the deposited energy to
    # the frozen release: W = (E_X,frozen + R E_X,win) / E_c(N_frozen).
    EXw, NXw = HEAT['release']['EX'], HEAT['release']['NX']
    Ec_fz = (NXs / CC) ** (4 / 3)
    rows = []
    for r in HEAT['table']['benchmark_MY']:
        W_fz = (EXs + r['R_total'] * EXw) / Ec_fz
        W_win = HEAT['release']['E_over_Ec'] * (1 + r['R_total'])
        s_fz = SQ_of_W(W_fz)
        rows.append(dict(e=r['e'], R_total_dq_heat=r['R_total'], W_frozen=W_fz, W_window=W_win, SQ_hartree_frozen=s_fz,
                         SQ_hartree_window=SQ_of_W(W_win), SQ_final_lo=s_fz + 0.0040, SQ_final_hi=s_fz + 0.0068,
                         branch='free X' if s_fz + 0.004 > SIG['bulk'] else 'Q-matter',
                         margin=(s_fz + 0.004) / SIG['bulk']))
        print('e=%.2e R=%.4f W_fz=%.3f W_win=%.3f S/Q=%.4f (window %.4f) %s' % (r['e'], r['R_total'], W_fz, W_win, s_fz, rows[-1]['SQ_hartree_window'], rows[-1]['branch']), flush=True)
    out['direct_heat_branch'] = rows
    # e_crit with the direct heat (power-law interpolation of R_vec and R_cascade in e; R_Y const), threshold incl. increments
    tab = HEAT['table']['benchmark_MY']
    le = np.log([r['e'] for r in tab]); lv = np.log([r['R_vec'] for r in tab]); lc = np.log([r['R_cascade'] for r in tab])
    RYh = tab[0]['R_Y_hom'] + tab[0]['R_Y_k']

    def Rtot(e):
        x = math.log(e)
        return RYh + math.exp(np.interp(x, le, lv)) + math.exp(np.interp(x, le, lc))
    ecr = {}
    for k, s in SIG.items():
        for lab, inc in (('hartree', 0.0), ('with_increments_hi', 0.0068), ('with_increments_lo', 0.0040)):
            f = lambda e: SQ_of_W((EXs + Rtot(e) * EXw) / Ec_fz) + inc - s
            try:
                ecr[k + '_' + lab] = brentq(f, 1.0e-6, 9.9e-6, xtol=1e-11)
            except ValueError:
                ecr[k + '_' + lab] = None
    out['e_crit_direct_heat'] = ecr
    out['e_crit_dq_heat'] = {k: v['e_crit'] for k, v in HEAT['e_crit']['benchmark'].items()}
    print('e_crit (direct heat):', ecr, flush=True)
    # ---------------------------------------------------------------- 5. readings A and B (re-typed) at e = 1e-6 and e_crit
    w0 = 15.603605078079994; WSP = 13.041178                       # spot E/E_c (wt_critic)
    H_SPOT, Y_SPOT = 1.3334961939833566e11, 4.8294009076988175e-09
    Mfx = 1.00114
    H_fx = H_M1 * (0.9374 / (Mfx * y_dq)) ** (2 / 3)
    Y_fx = D_GEV / (Mfx * MX)
    nu = (H_fx / H_SPOT) ** 2 * (Y_SPOT / Y_fx) ** (4 / 3)
    RB = lambda e, rv: RY + Rc * (1e-6 / e) ** 0.5 + nu * (WSP / w0) * rv * (1e-6 / e)
    RA = lambda e: 0.8056 + 4.0538 * (1e-6 / e) + 0.10716 * (1e-6 / e) ** 0.5
    out['readings'] = dict(nu=nu, H_freeX=H_fx, Y_freeX=Y_fx,
                           W_vec_transfer_band=nu * WSP * 4.0538, W_vec_transfer_broad=nu * WSP * 4.7905,
                           W_vec_direct_dq_heat=HEAT['release']['E_over_Ec'] * tab[1]['R_vec'],
                           SQ_B_band_1e6=SQ_of_W(w0 * (1 + RB(1e-6, 4.0538))), SQ_A_1e6=SQ_of_W(w0 * (1 + RA(1e-6))),
                           SQ_B_band_2e6=SQ_of_W(w0 * (1 + RB(2e-6, 4.0538))), SQ_B_broad_2e6=SQ_of_W(w0 * (1 + RB(2e-6, 4.7905))))
    out['readings']['direct_over_transfer_band'] = out['readings']['W_vec_direct_dq_heat'] / out['readings']['W_vec_transfer_band']
    out['readings']['direct_over_transfer_broad'] = out['readings']['W_vec_direct_dq_heat'] / out['readings']['W_vec_transfer_broad']
    # ---------------------------------------------------------------- 6. H_i ledger (re-typed)
    out['H_ledger'] = dict(freeX=H_fx, Qmatter=[H_M1 * (0.9374 / (M_ * y_sk)) ** (2 / 3) for M_ in (0.802, 0.782)],
                           Q2e5=H_M1 * (0.9374 / (0.9374 * y_sk)) ** (2 / 3))
    # ---------------------------------------------------------------- 7. analytic ideal-gas law vs Hartree
    out['ideal_gas_law'] = dict(formula='S/Q = (4/3)(3 pi^2/30)^(1/4)(3/4)^(3/4) lX^(1/4) (W-1)^(3/4)',
                                prefactor=(4 / 3) * (3 * math.pi ** 2 / 30) ** 0.25 * 0.75 ** 0.75 * LX ** 0.25,
                                W_crit_bulk_ideal=1 + (SIG['bulk'] / ((4 / 3) * (3 * math.pi ** 2 / 30) ** 0.25 * 0.75 ** 0.75 * LX ** 0.25)) ** (4 / 3))
    # ---------------------------------------------------------------- 8. spinodal increment re-typed
    hf = json.load(open(RES / 'sk_homfunc.json'))
    fc_, Ton, nQ = 0.9356634, 7.4155351, 154.098955
    ncon = fc_ * nQ
    af = (ncon / hf['spinodal']['n_sp']) ** (1 / 3)
    out['spinodal_increment'] = dict(dE=hf['latent_heat_T0']['dE_per_charge'], a_sp_over_a_on=af, T_sp=Ton / af,
                                     dSQ=fc_ * hf['latent_heat_T0']['dE_per_charge'] / (Ton / af))
    json.dump(out, open(RES / 'sk_ledger.json', 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ('direct_heat_branch',)}, indent=1, default=float))


if __name__ == '__main__':
    main()
