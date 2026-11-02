#!/usr/bin/env python3
"""mf_ideal regression (Claude, own code).  Three checks, all on M1r's (m_I, Gamma_I) = (3, 0.5) release (lambda = 0.5,
trilinear holder, prepared at u = 10i, R = 12; portals vanishing at M_Y):
  R1  MODE 0 of our copy rpmfi.c against dq_entropy's rpdq.c and inf_chain's rpchain.c (same arguments, tau <= 1000):
      byte-for-byte identical stdout, and N_X(400-1000), N_YL, N_YH, E_X against the quoted fx_holder / dq_entropy digits.
  R2  the MODE 5 table interface ("bridge"): M1r's own recorded background (t, ln a, H, R) and its holder c = I/I_0
      (with cdot = I'/I_0) written as a bench-like table, fed to MODE 5 with MFI_STAB = 1, to tau = 1000.  This is the
      path mf_chain will use for bench_*.json files.  Reported: relative differences of the late means.
  R3  MODE 5 fluid + evolved spectator in the M1r-equivalent limit is NOT M1r (different background); only a sanity run.
Usage: python3 -B regress.py   -> ../results/regression.json
"""
import json, math, subprocess, sys
import numpy as np
import mfi_common as mc

QUOTED = dict(NX_400_1000=-8045.485140302082, NYL=0.566640911258934, NYH=0.010377131381878653, EX=188706.29275778163)


def m1r_prep():
    st = mc.prep(mc.U_PREP, 12.0)
    G, I0, ph0 = mc.fc.trilinear_G(st['z'], mc.P05, 3.0)
    return st, G


def late(out, lo=400.0, hi=1000.0):
    B = out['B']
    return {k: mc.bin_mean(B, k, lo, hi) for k in ('NX', 'EX', 'NYL', 'NYH', 'Nsig', 'eps')}


def write_table(path, D, I0, ms2=9.0):
    """bench-like table from MODE 0 dump rows D: t y[0..18] H R."""
    t = D[:, 0]; la = D[:, 1]; I = D[:, 2]; J = D[:, 3]; H = D[:, 20]; R = D[:, 21]
    with open(path, 'w') as fh:
        fh.write('%d 7\n' % len(t))
        for i in range(len(t)):
            fh.write('%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n' % (t[i], la[i], H[i], R[i], I[i] / I0, J[i] / I0, ms2))
    return len(t)


def main():
    mc.build()
    sp = mc.SP / 'build'
    for src, name, flags in ((mc.TRK / 'dq_entropy' / 'code' / 'rpdq.c', 'rpdq_ref', ['-O3', '-std=c11']),
                             (mc.TRK / 'inf_chain' / 'code' / 'rpchain.c', 'rpchain_ref', ['-O2'])):
        subprocess.run(['gcc'] + flags + ['-o', str(sp / name), str(src), '-lm'], check=True)
    st, G = m1r_prep()
    z = st['z']
    out = dict(prep=dict(amps=[float(x) for x in np.abs(z)], eig_min=float(sorted(st['eig'])[1]), grad=st['grad'], G=[G.real, G.imag],
                         F_s=mc.F_of(z)))
    # ---- R1
    runs = {}
    for name in ('rpmfi', 'rpdq_ref', 'rpchain_ref'):
        runs[name] = mc.run('reg_m0_' + name, 0, 3.0, 0.5, z, G, 1000.0, 5000, binw=20.0, save=False, binary=sp / name)
    same_dq = runs['rpmfi']['stdout'] == runs['rpdq_ref']['stdout']
    same_ch = runs['rpmfi']['stdout'] == runs['rpchain_ref']['stdout']
    L = late(runs['rpmfi'])
    r1 = dict(byte_identical_vs_rpdq=same_dq, byte_identical_vs_rpchain=same_ch,
              mine=dict(NX_400_1000=L['NX'], NYL=L['NYL'], NYH=L['NYH'], EX=L['EX']), quoted=QUOTED,
              rel_diff={k: (dict(NX_400_1000=L['NX'], NYL=L['NYL'], NYH=L['NYH'], EX=L['EX'])[k] / v - 1) for k, v in QUOTED.items()},
              END=runs['rpmfi']['END'], seconds=runs['rpmfi']['seconds'])
    print('R1 byte-identical: rpdq %s rpchain %s; NX(400-1000) %.15g (quoted %.15g) NYL %.15g EX %.15g' % (
        same_dq, same_ch, L['NX'], QUOTED['NX_400_1000'], L['NYL'], L['EX']), flush=True)
    out['R1'] = r1
    # ---- R2 bridge: dense MODE 0 dump -> table -> MODE 5 table run
    res = {}
    for dt in (0.002, 0.001):
        nd = int(round(100.0 / dt))
        d0 = mc.run('reg_m0_dump', 0, 3.0, 0.5, z, G, 100.0, nd, save=False, dump=True)
        D = d0['D']; I0e = d0['I0']['I0']; ph0 = d0['I0']['phihat0']
        # prepend the t = 0 row (the dump starts at the first sample); row layout: t y[0..18] H R
        S0 = d0['S'][0]
        first = np.r_[0.0, 0.0, I0e, 0.0, np.zeros(16), S0[mc.COLS.index('H')], S0[mc.COLS.index('R')]]
        D = np.vstack([first, D])
        path = mc.SP / ('table_m1r_dt%g.txt' % dt)
        n = write_table(path, D, I0e)
        Geff = G * ph0
        env = dict(MFI_TABLE=str(path), MFI_STAB=1, MFI_MS2=9.0, MFI_GS=0.5, MFI_BR=0.0)
        m5 = mc.run('reg_m5_bridge_dt%g' % dt, 5, 3.0, 0.0, z, Geff, 1000.0, 5000, env=env, binw=20.0, save=False)
        L5 = late(m5)
        rel = {k: L5[k] / L[k] - 1 for k in ('NX', 'EX', 'NYL', 'NYH', 'Nsig', 'eps')}
        res['dt%g' % dt] = dict(rows=n, rel_diff_vs_mode0=rel, late=L5, END=m5['END'], seconds=m5['seconds'],
                                max_lna_drift=float(np.max(np.abs(mc.col(m5['S'], 'Hcheck')))), u0=m5['I0']['u0'])
        print('R2 dt=%g (%d rows): rel diffs %s | END %s' % (dt, n, {k: '%.2e' % v for k, v in rel.items()}, m5['END']), flush=True)
    out['R2'] = res
    json.dump(out, open(mc.RES / 'regression.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
