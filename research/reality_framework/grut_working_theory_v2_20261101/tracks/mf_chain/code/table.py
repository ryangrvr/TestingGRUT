#!/usr/bin/env python3
"""mf_chain: print/collect the chain results (results/engine/chain/chain_<tag>.json) as one table.
Usage: python3 -B table.py [TAG ...]   (default: all)"""
import json, sys, glob, os
import mfc_common as M


def row(tag):
    o = json.load(open(M.ENG / 'chain' / ('chain_%s.json' % tag)))
    rel = json.load(open(M.ENG / 'release' / ('release_%s.json' % tag)))['result']
    a, b = o['e1e-6'], o['e3e-6']
    return dict(tag=tag, K=o['K'], eps=o['eps'], E_over_Ec=o['E_over_Ec'], N_X=o['N_X'], E_X=o['E_X'], N_YL=o['N_YL'], N_YH=o['N_YH'],
                N_sigma=o['N_sigma'], N_A_1e6=o['N_A_1e6'], R_vec_1e6=a['R_vec'], R_Y=a['R_Y'], R_casc_1e6=a['R_casc'],
                SQ_1e6=a['S_per_Q'], branch_1e6=a['branch'], SQ_3e6=b['S_per_Q'], branch_3e6=b['branch'], SQ_inf=o['S_per_Q_e_inf'],
                e_crit_bulk=o['e_crit_bulk'], e_crit_Q2e5=o['e_crit_Q2e5'], p_1e6=a.get('p'), M_1e6=a['M'], M_3e6=b['M'],
                H_req_1e6=a['H_i_req_GeV'], H_req_3e6=b['H_i_req_GeV'], r_1e6=a['r'], Y_over_YQ=o['Y_over_YQ_at_HM1'],
                N_X_err=rel['N_X_err'], gauss=rel['gauss'], work=rel['hidden_work'], Hcons=rel['H_consistency'])


def main(tags):
    if not tags:
        tags = sorted(os.path.basename(p)[6:-5] for p in glob.glob(str(M.ENG / 'chain' / 'chain_*.json')))
    rows = []
    for t in tags:
        try:
            rows.append(row(t))
        except FileNotFoundError:
            pass
    hdr = '%-16s %6s %6s %6s %7s %7s %7s %7s %6s %7s %7s %7s %-6s %7s %-6s %6s %8s %6s %9s'
    print(hdr % ('tag', 'K', 'eps', 'E/Ec', 'N_YL', 'N_YH', 'N_sig', 'N_A', 'Rvec', 'R_Y', 'Rcasc', 'SQ1e-6', 'br', 'SQ3e-6', 'br', 'SQinf',
                 'e_crit', 'M', 'H_req'))
    for r in rows:
        print('%-16s %6.3f %6.4f %6.2f %7.3f %7.4f %7.2g %7.3f %6.2f %7.3f %7.3f %7.2f %-6s %7.2f %-6s %6.2f %8s %6.3f %9.3e' % (
            r['tag'], r['K'], r['eps'], r['E_over_Ec'], r['N_YL'], r['N_YH'], r['N_sigma'], r['N_A_1e6'], r['R_vec_1e6'], r['R_Y'],
            r['R_casc_1e6'], r['SQ_1e6'], r['branch_1e6'][:6], r['SQ_3e6'], r['branch_3e6'][:6], r['SQ_inf'],
            ('%.2e' % r['e_crit_bulk']) if r['e_crit_bulk'] else '-', r['M_1e6'], r['H_req_1e6']))
    return rows


if __name__ == '__main__':
    main(sys.argv[1:])
