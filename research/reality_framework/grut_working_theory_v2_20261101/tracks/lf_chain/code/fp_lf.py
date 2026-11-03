#!/usr/bin/env python3
"""lf_chain step 3 (Claude): branch-aware self-consistency on the LF closure, windows, representative points.

Nodes (final histories, all LF closure, HFa): computed lattice blocks at 8.4e10 (KFmid_H84_nl) and 1.3e11 (KFmid);
two extra nodes built with e2e_chain's rescaling rule (lattice block nearest in H, kinetic network at that H) [A]:
5.5e10 (8.4e10 lattice block) and 2.7e11 (1.3e11 lattice block).
Rule (the e2e_chain skeptic's, sk_static.py D, re-typed unchanged): for branch b in {free X, Q-matter},
f_b(node) = ln[(H_req/H_i)_node (M_b(S/Q_node)/M_node)^(-2/3)], interpolated linearly in ln H between nodes; a root counts
iff the interpolated S/Q(1e-6) is on branch b.  PRIMARY = the computed pair (8.4e10, 1.3e11), extrapolation flagged by w,
exactly as e2e_chain/its skeptic; ALL = piecewise over every node of that key (roots inside a pair, or outside the end
pairs, flagged).  The e2e_chain rule (fixpoint.solve, read-only) is recorded for the computed pair.
Windows [model choice, stated]: M2-Q = a self-consistent Q-matter root at e = 1e-6; M2-X = a self-consistent free-X root
at 1e-6 that has an e_crit <= 1e-3 (M1r's two-branch picture), "M2-X proper" with 1e-6 < e_crit <= 4e-6 (e2e_chain's
"e_crit ~ 2e-6").  Output: results/fixpoints_LF.json.  Labels: [C] at nodes, [E] between, [E/A] with the [A] nodes.
"""
import json, math
import lfc_common as L
import pmap as PM
import relic as RL
import fixpoint as FP                                          # e2e_chain, read-only (node(), solve())

SQB = 3.4716191256166837
CL = 'cL0_vanish_at_MY'
NODES = [('hLF55_mid', 'A'), ('hLF84_mid', 'C'), ('hLF_mid', 'C'), ('hLF270_mid', 'A')]


def M_branch(b, sq):
    if b == 'free X':
        p, _ = PM.p_selfconsistent(max(sq, SQB + 1e-6))
        return 1 + 2 * p
    return RL.relic(min(sq, SQB - 1e-6))['M_central']


def nd(o):
    b = o['by_e'][CL]
    return dict(H=o['H_i_GeV'], K=o['K'], SQ=b['e1e-6']['S_per_Q'], SQ3=b['e3e-6']['S_per_Q'], M=b['e1e-6']['M'],
                R=b['e1e-6']['H_req_over_H_i'], Hreq=b['e1e-6']['H_i_req_GeV'], ecrit=b['e_crit_bulk'], ecritQ=b['e_crit_Q2e5'],
                NYL=o['N_YL'], eps=o['eps'], EEc=o['E_over_Ec'], branch=b['e1e-6']['branch'], SQinf=o['S_per_Q_e_inf'],
                r=b['e1e-6']['r_at_req'])


def fb(n, b):
    return math.log(n['R'] * (M_branch(b, n['SQ']) / n['M']) ** (-2.0 / 3.0))


def root(lo, hi, b, rule='lnH'):
    x = math.log if rule == 'lnH' else (lambda h: h)
    f0, f1 = fb(lo, b), fb(hi, b)
    if f1 == f0:
        return None
    w = -f0 / (f1 - f0)
    lin = lambda k: lo[k] + w * (hi[k] - lo[k])
    sq = lin('SQ')
    Hs = math.exp(x(lo['H']) + w * (x(hi['H']) - x(lo['H']))) if rule == 'lnH' else x(lo['H']) + w * (x(hi['H']) - x(lo['H']))
    ec = None
    if lo['ecrit'] and hi['ecrit']:
        ec = math.exp(math.log(lo['ecrit']) + w * (math.log(hi['ecrit']) - math.log(lo['ecrit'])))
    elif lo['ecrit'] or hi['ecrit']:
        ec = 'one node only: %.3g at %.3g' % ((lo['ecrit'] or hi['ecrit']), (lo['H'] if lo['ecrit'] else hi['H']))
    return dict(branch=b, H_star=Hs, w=w, extrapolated=not (0 <= w <= 1), SQ=sq, SQ3=lin('SQ3'), K=lin('K'), eps=lin('eps'),
                NYL=lin('NYL'), self_consistent=bool((sq > SQB) == (b == 'free X')), e_crit=ec, M=M_branch(b, sq),
                pair=[lo['H'], hi['H']])


def main():
    scans = {t: L.load(L.RES / ('scan_%s.json' % t), {}).get('runs', {}) for t, _ in NODES}
    keys = sorted({k for t in scans for k in scans[t]} - {'1.000:+0'}, key=lambda k: (float(k.split(':')[0]), float(k.split(':')[1])))
    out = {}
    for k in keys:
        ns = [(t, q, nd(scans[t][k])) for t, q in NODES if k in scans[t]]
        rec = dict(nodes={t: dict(n, quality=q) for t, q, n in ns})
        cp = [n for t, q, n in ns if q == 'C']
        if len(cp) == 2:
            lo, hi = cp
            rec['primary'] = {b: root(lo, hi, b) for b in ('free X', 'Q-matter')}
            rec['primary_linearH'] = {b: root(lo, hi, b, 'H') for b in ('free X', 'Q-matter')}
            try:
                e2c = FP.solve(FP.node(scans['hLF84_mid'][k]), FP.node(scans['hLF_mid'][k]))
                rec['e2e_chain_rule'] = dict(H_star=e2c['H_star_GeV'], SQ=e2c['SQ_1e6'], branch=e2c['branch_1e6'], M=e2c['M_1e6'], w=e2c['w'])
            except Exception as ex:                                 # noqa: BLE001
                rec['e2e_chain_rule'] = str(ex)
        # piecewise over all nodes
        allr = []
        nn = [n for _, _, n in ns]
        for b in ('free X', 'Q-matter'):
            for i in range(len(nn) - 1):
                r = root(nn[i], nn[i + 1], b)
                if r is None:
                    continue
                inside = 0 <= r['w'] <= 1
                edge_lo, edge_hi = i == 0 and r['w'] < 0, i == len(nn) - 2 and r['w'] > 1
                if inside or edge_lo or edge_hi:
                    r['uses_A_node'] = any(q == 'A' for t, q, n in ns[i:i + 2])
                    allr.append(r)
        rec['all_nodes'] = allr
        def classify(rs):
            cls = []
            for r in rs:
                if not r or not r['self_consistent']:
                    continue
                if r['branch'] == 'Q-matter':
                    cls.append('M2-Q')
                else:
                    ec = r['e_crit']
                    if isinstance(ec, float) and ec <= 1e-3:
                        cls.append('M2-X proper' if 1e-6 < ec <= 4e-6 else 'M2-X')
                    elif isinstance(ec, str):
                        cls.append('M2-X (e_crit at one node only)')
                    else:
                        cls.append('free X, no e_crit')
            return cls
        rec['class_primary'] = classify(list(rec.get('primary', {}).values()))
        rec['class_all'] = classify(allr)
        rec['n_roots_primary'] = len(rec['class_primary']); rec['n_roots_all'] = len(rec['class_all'])
        out[k] = rec
        pr = rec.get('primary', {})
        s = ' | '.join('%s H*%.3g w%.2f SQ%.3f %s ec %s' % (b[:6], r['H_star'], r['w'], r['SQ'], 'OK' if r['self_consistent'] else '--',
                       ('%.2g' % r['e_crit']) if isinstance(r['e_crit'], float) else ('1n' if r['e_crit'] else '-')) for b, r in pr.items() if r)
        sa = '; '.join('%s %.3g%s%s' % (r['branch'][:6], r['H_star'], '*' if r['uses_A_node'] else '', 'x' if r['extrapolated'] else '') for r in allr if r['self_consistent'])
        print('%-9s %s || all: %s || %s' % (k, s, sa, rec['class_all']))
    L.dump(L.RES / 'fixpoints_LF.json', out)


if __name__ == '__main__':
    main()
