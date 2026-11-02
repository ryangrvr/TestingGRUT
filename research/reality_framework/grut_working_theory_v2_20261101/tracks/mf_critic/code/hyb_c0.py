#!/usr/bin/env python3
"""mf_critic check K11 (Claude): is the strong-holder point of the mf_chain skeptic robust to the reheating bracket?

The skeptic's Q-matter point (Pauli-blocked HFa, Gamma_th = 0, perturbative tails free, H_i = 8.40e10 GeV, K = 6.50,
c0 = 1.40: S/Q(1e-6) = 2.875) was computed on ONE corner of the Pauli bracket.  Here the same history construction
(mf_chain's hist_build.py, read-only, run as a subprocess with its scratch redirected to ours) is rebuilt with the other
end of the skeptic's thermalisation bracket, Gamma_th = 10 H_i (K = 3.19 at 8.40e10 GeV), and the full chain is run at
c0 = 1 and at the skeptic's c0 values.  A regression run reproduces the skeptic's c0 = 1.40 point on its own P0 history
(its bench file in its scratch, read-only).

Steps per history: mf_ideal's bench_run (MODE 5 table, spectator evolved with the bench's m_s^2(t), Gamma_s from the
bench), then for each c0 a release with MFI_SM2TAB = 1 on the same table (as the skeptic's sk_c0chain.py), vectors and
chain.  All of mf_ideal's drivers are read-only; outputs go to this track's scratch.
Usage: python3 -B hyb_c0.py BENCH.json TAG c0 [c0 ...]    -> results/hyb_c0.json (merged)
Labels: [C] given the model choices (the Pauli fluid and its Gamma_th are the mf_hybrid skeptic's [E] model).
"""
import json, sys, time
import crc_common as C
import bench_run as br                                   # mf_ideal, read-only


def base(bench, tag):
    if not (C.ENG / 'chain' / ('chain_%s.json' % tag)).exists():
        sys.argv = ['bench_run.py', bench, '--tag', tag, '--R', 'raw', '--tlong', '1e4']
        br.main()
    return json.load(open(C.ENG / 'chain' / ('chain_%s.json' % tag)))


def c0run(tag, c0):
    rel0 = json.load(open(C.ENG / 'release' / ('release_%s.json' % tag)))
    sp = rel0['spec']
    orig = C.rm.spec_env

    def spec_env(*a, **kw):
        env, gam = orig(*a, **kw)
        env['MFI_SM2TAB'] = '1'
        return env, gam
    C.rm.spec_env = spec_env
    t2 = '%s_C%.3f' % (tag, c0)
    try:
        C.rm.release(t2, sp['K_target'], sp['m_s'], 0.0, gs=sp['Gamma_s'], sstar=sp['s_star_GeV'], c0=c0, bg=sp['bg'],
                     Tl=rel0['result']['T_long'], bcap=sp['step_cap_b'])
    finally:
        C.rm.spec_env = orig
    C.vm.main(t2)
    o = C.cm.analyse(t2)
    C.dump(C.ENG / 'chain' / ('chain_%s.json' % t2), o)
    return o


def line(tag, s):
    a, b = s['e1e-6'], s['e3e-6']
    return '%s: K %.3f eps %.4f E/Ec %.3f N_YL %.4g N_YH %.4g | 1e-6: Rv %.3g RY %.3g Rc %.3g S/Q %.4f %s M %.4g H_req %.4e | 3e-6: S/Q %.4f %s | inf %.3f | e_crit %s / %s' % (
        tag, s['K'], s['eps'], s['E_over_Ec'], s['N_YL'], s['N_YH'], a['R_vec'], a['R_Y'], a['R_casc'], a['S_per_Q'], a['branch'], a['M'],
        a['H_i_req_GeV'], b['S_per_Q'], b['branch'], s['S_per_Q_e_inf'], s['e_crit_bulk'], s['e_crit_Q2e5'])


if __name__ == '__main__':
    argv = list(sys.argv)
    bench, tag = argv[1], argv[2]
    C.mc.build()
    outp = C.RES / 'hyb_c0.json'
    allr = json.load(open(outp)) if outp.exists() else {}
    d = json.load(open(bench))
    meta = dict(bench=bench, H_i_GeV=d['params'].get('H_i_GeV'), reheating={k: d.get('reheating', {}).get(k) for k in
                ('K', 'Gamma_th', 'tails', 't_sw', 'frad_t03', 'frad_t06', 'frad_t3', 'frad_t100', 't_SM_1em3', 't_SM_half', 't_rad_99')})
    t0 = time.time()
    s = C.short(base(bench, tag)); s['c0'] = 1.0; s['seconds'] = time.time() - t0
    print(line(tag, s), flush=True)
    allr.setdefault(tag, dict(meta=meta, runs={}))['runs']['1.000'] = s
    C.dump(outp, allr)
    for c in argv[3:]:
        t0 = time.time()
        s = C.short(c0run(tag, float(c))); s['c0'] = float(c); s['seconds'] = time.time() - t0
        print(line('%s_C%s' % (tag, c), s), flush=True)
        allr[tag]['runs']['%.3f' % float(c)] = s
        C.dump(outp, allr)
