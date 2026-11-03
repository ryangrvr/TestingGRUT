#!/usr/bin/env python3
"""e2e_chain item 2 (Claude): self-consistent points from the holder scans.

For a closure (GM: the mf_hybrid skeptic's box-mean Pauli saturation on the lattice, which the kinetic code's
homogeneous exact-mode blocking agrees with at t_sw; P0: local per-site blocking) two final histories exist, at
H_i = 8.4e10 (lattice H84_*) and 1.3e11 GeV (lattice design), each lattice -> kinetic moments hand-over at t_sw = 3.
At fixed holder (c0, delta), f(ln H) = ln(H_req/H_i) is interpolated linearly in ln H between the two histories (the
V0 rescaling of mf_chain iterate.py; the dimensionless history moves with H only through the early stage and K) and
its root is the self-consistent H_i [C at the two nodes, E in between; extrapolation flagged].  S/Q, e_crit, M, K are
interpolated in ln H to the root.  N_* is re-solved at the root in cmb.py.
Usage: python3 -B fixpoint.py  -> results/fixpoints.json
"""
import json, math
import numpy as np
import e2c_common as E

CL = 'cL0_vanish_at_MY'
SETS = {'GM': ('hH84_GM', 'hHFa_GM'), 'P0': ('hH84_P0', 'hHFa_P0')}


def runs(tag):
    return E.load(E.RES / ('scan_%s.json' % tag))['runs']


def node(o, cl=CL):
    b = o['by_e'][cl]
    return dict(H=o['H_i_GeV'], K=o['K'], f=math.log(b['e1e-6']['H_req_over_H_i']), SQ6=b['e1e-6']['S_per_Q'], SQ3=b['e3e-6']['S_per_Q'],
                branch6=b['e1e-6']['branch'], branch3=b['e3e-6']['branch'], M6=b['e1e-6']['M'], ecrit=b['e_crit_bulk'],
                ecritQ=b['e_crit_Q2e5'], N_YL=o['N_YL'], N_YH=o['N_YH'], eps=o['eps'], EEc=o['E_over_Ec'],
                Rv=b['e1e-6']['R_vec'], RY=b['e1e-6']['R_Y'], Rc=b['e1e-6']['R_casc'], Rl=b['e1e-6']['R_leak'], Rg=b['e1e-6']['dR_grav'],
                SQinf=o['S_per_Q_e_inf'])


def solve(lo, hi):
    x0, x1 = math.log(lo['H']), math.log(hi['H'])
    s = (hi['f'] - lo['f']) / (x1 - x0)
    x = x0 - lo['f'] / s
    w = (x - x0) / (x1 - x0)
    lin = lambda a, b: None if (a is None or b is None) else a + w * (b - a)
    loglin = lambda a, b: None if (a is None or b is None) else math.exp(math.log(a) + w * (math.log(b) - math.log(a)))
    out = dict(H_star_GeV=math.exp(x), w=w, extrapolated=not (0 <= w <= 1), slope_dlnratio_dlnH=s,
               K=lin(lo['K'], hi['K']), SQ_1e6=lin(lo['SQ6'], hi['SQ6']), SQ_3e6=lin(lo['SQ3'], hi['SQ3']),
               e_crit_bulk=loglin(lo['ecrit'], hi['ecrit']), e_crit_Q2e5=loglin(lo['ecritQ'], hi['ecritQ']),
               M_1e6=lin(lo['M6'], hi['M6']), N_YL=lin(lo['N_YL'], hi['N_YL']), N_YH=lin(lo['N_YH'], hi['N_YH']), eps=lin(lo['eps'], hi['eps']),
               E_over_Ec=lin(lo['EEc'], hi['EEc']), R_vec=lin(lo['Rv'], hi['Rv']), R_Y=lin(lo['RY'], hi['RY']), R_casc=lin(lo['Rc'], hi['Rc']),
               R_leak=lin(lo['Rl'], hi['Rl']), dR_grav=lin(lo['Rg'], hi['Rg']), SQ_inf=lin(lo['SQinf'], hi['SQinf']),
               nodes=dict(lo=lo, hi=hi))
    sq = out['SQ_1e6']
    out['branch_1e6'] = 'free X' if sq > 3.4716191 else ('Q-matter' if sq < 2.9283779 else 'Q-matter (bulk only)')
    out['branch_3e6'] = 'free X' if out['SQ_3e6'] > 3.4716191 else 'Q-matter'
    out['node_branches_agree'] = lo['branch6'][:6] == hi['branch6'][:6]
    return out


def main():
    res = {}
    for cl in (CL, 'cL_1/16pi^2'):
        for nm, (tlo, thi) in SETS.items():
            rl, rh = runs(tlo), runs(thi)
            for key in sorted(set(rl) & set(rh)):
                if key.startswith('1.000'):
                    continue
                o = solve(node(rl[key], cl), node(rh[key], cl))
                res.setdefault(cl, {}).setdefault(nm, {})[key] = o
                print('%-12s %s %s: H* %.4e (w %.2f%s) K %.2f S/Q 1e-6 %.3f [%s, nodes agree %s] 3e-6 %.3f e_crit %s M %.4f' % (
                    cl, nm, key, o['H_star_GeV'], o['w'], ' EXTRAP' if o['extrapolated'] else '', o['K'], o['SQ_1e6'], o['branch_1e6'],
                    o['node_branches_agree'], o['SQ_3e6'], ('%.3g' % o['e_crit_bulk']) if o['e_crit_bulk'] else '-', o['M_1e6']))
    E.dump(E.RES / 'fixpoints.json', res)


if __name__ == '__main__':
    main()
