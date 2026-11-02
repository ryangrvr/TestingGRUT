#!/usr/bin/env python3
"""mf_critic check K10 (Claude): the holder amplitude c0 on M1r's OWN supplied history.

Question.  The mf_chain skeptic found that, on the hybrid histories, the branch at e = 1e-6 is decided by the holder
amplitude c0 (mu_eff(slice) = mu + c0 (10i - mu); c0 = 1 is M1r's 10i, inherited from Codex's Option 1A and kept by
fx_holder's g_I calibration).  Is the same true on M1r's own (m_I, Gamma_I) = (3, 0.5) static-slice history?  If so,
"recovering M1r's picture" on a hybrid at c0 != 1 is the same holder re-tune applied to a different history, and the
decisive pass's e_crit = 2.06e-6 is a c0 = 1 number.

Method.  mf_ideal's R2 bridge table (M1r's own recorded background and holder c = I/I0, cdot = I'/I0; read-only) is
copied into this track's scratch with the c and cdot columns multiplied by c0, so that the table-following holder
(MFI_STAB = 1) gives u(t) = mu + (10i - mu) c0 I(t)/I0.  The preparation is re-solved at u = mu + c0 (10i - mu), R = 12
(release_mfi does this from c0).  Then mf_ideal's release (dense + long, tau <= 1e4), vectors on the recorded background
(f_exit = 1) and chain (dq_heat ledger, closed form, p map, H_i), all read-only, outputs in scratch.
c0 = 1 must reproduce the decisive pass: S/Q(1e-6) = 4.955, e_crit = 2.06e-6.
Usage: python3 -B m1r_c0.py c0 [c0 ...]     -> results/m1r_c0.json (merged)
Labels: [C] given the model choices.
"""
import json, sys, time
import numpy as np
import crc_common as C


def table_for(c0):
    path = C.SP / 'tab' / ('table_m1r_c%.3f.txt' % c0)
    if path.exists():
        return path
    lines = open(C.M1R_TABLE).read().split('\n')
    out = [lines[0]]
    for ln in lines[1:]:
        if not ln.strip():
            continue
        v = [float(x) for x in ln.split()]
        v[4] *= c0; v[5] *= c0
        out.append(' '.join('%.17g' % x for x in v))
    open(path, 'w').write('\n'.join(out) + '\n')
    return path


def run(c0):
    tag = 'M1r_c%.3f' % c0
    t0 = time.time()
    path = table_for(c0)
    C.rm.release(tag, 1.986, 3.0, 0.0, gs=0.5, sstar=1e300, c0=c0, bg='tables:%s' % path, Tl=1.0e4, bcap=3.0)
    C.vm.main(tag)
    o = C.cm.analyse(tag)
    C.dump(C.ENG / 'chain' / ('chain_%s.json' % tag), o)
    s = C.short(o)
    s['c0'] = c0; s['seconds'] = time.time() - t0
    a, b = s['e1e-6'], s['e3e-6']
    print('%s: eps %.4f E/Ec %.3f N_YL %.4g N_YH %.4g | 1e-6: Rv %.3g RY %.3g Rc %.3g S/Q %.4f %s M %.4g H_req %.3e | 3e-6: S/Q %.4f %s | inf %.3f | e_crit %s / %s (%.0f s)' % (
        tag, s['eps'], s['E_over_Ec'], s['N_YL'], s['N_YH'], a['R_vec'], a['R_Y'], a['R_casc'], a['S_per_Q'], a['branch'], a['M'],
        a['H_i_req_GeV'], b['S_per_Q'], b['branch'], s['S_per_Q_e_inf'], s['e_crit_bulk'], s['e_crit_Q2e5'], s['seconds']), flush=True)
    return s


if __name__ == '__main__':
    C.mc.build()
    outp = C.RES / 'm1r_c0.json'
    allr = json.load(open(outp)) if outp.exists() else {}
    for c in sys.argv[1:]:
        r = run(float(c))
        allr['%.3f' % float(c)] = r
        C.dump(outp, dict(sorted(allr.items())))
