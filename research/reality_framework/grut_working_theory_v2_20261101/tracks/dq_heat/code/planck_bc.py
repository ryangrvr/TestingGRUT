#!/usr/bin/env python3
"""dq_heat (Claude, own code): the M1r release with the induced inflaton portals vanishing at M_P instead of at M_Y.

V_c = c (|u|^2 - mu^2)(|Sigma|^2 + 2|Y|^2), c = -ln(M_P/M_Y)/(4 pi^2) = -0.3371983 (fx_holder SKEPTIC N1), in rptri_cp.c.
Preparation: conditional minimum at u = 10 i with the c term, i.e. M_Y^2 -> M_Y^2 + 2c', h^2 -> h^2 - c'/lambda_Sigma,
c' = c(|10 i|^2 - mu^2) (exact completion of the square at fixed u).  Controls: c = 0 in rptri_cp.c reproduces rptri_dq.c;
N_X(400-1000) is compared with the skeptic's independent code (-10606.160).  Writes ../results/background_planck.json and
the scratch npz in the format of background.py (for vector_modes.py with DQ_BG=planck).
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_model as fc
from background import series, Q

RES = Path(__file__).resolve().parent.parent / 'results'
B, GAM, LAM = 3.0, 0.5, 0.5
CPL = -0.33719832310397263


def run_cp(tag, c, p, z0, G, end, nsamp):
    os.environ['DQ_CPORT'] = '%.17g' % c
    old = fc.BIN
    fc.BIN = fc.SP / 'build' / 'rptri_cp'
    fc.SRC = Path(__file__).resolve().parent / 'rptri_cp.c'
    try:
        return fc.run(tag, p, B, GAM, z0=z0, G=G, end=end, nsamp=nsamp, binw=20.0)
    finally:
        fc.BIN = old; fc.SRC = Path(__file__).resolve().parent / 'rptri_dq.c'
        os.environ.pop('DQ_CPORT', None)


def main():
    t0 = time.time()
    p = fc.params(LAM)
    st0 = fc.stationary(fc.U_PREP, p); z00 = st0['z']
    G, I0, ph0 = fc.trilinear_G(z00, p, B)
    # control: c = 0 through rptri_cp
    ctrl = run_cp('cp0_control', 0.0, p, z00, G, 1000.0, 5000)
    ref = json.load(open(RES / 'background.json'))['validation']
    ctrl_nx = fc.window_mean(ctrl, 'NX', 400.0, 1000.0)
    # Planck boundary condition
    cprime = CPL * (abs(fc.U_PREP) ** 2 - p['mu'] ** 2)
    pe = dict(p); pe['M'] = np.sqrt(p['M'] ** 2 + 2 * cprime); pe['h'] = np.sqrt(p['h'] ** 2 - cprime / p['ls'])
    st = fc.stationary(fc.U_PREP, pe); zp = st['z']
    rc = run_cp('planck_coarse', CPL, p, zp, G, 2000.0, 10000)
    rd = run_cp('planck_dense', CPL, p, zp, G, 200.0, 20000)
    val = dict(NX_400_1000=fc.window_mean(rc, 'NX', 400.0, 1000.0), checks=rc['END'])
    for k in ('NYL', 'NYH', 'Nsig', 'EX', 'eps', 'NS'):
        val[k] = fc.window_mean(rc, k, 400.0, 1000.0)
    Y, ratio, Rr = fc.yield_over_target(rc, val['NX_400_1000'])
    val.update(yield_over_YQ=ratio, Rr_end=Rr, K_late=float(np.sqrt(Rr / 3.0)), skeptic_NX_400_1000=-10606.160254392727,
               skeptic_NYL=1.407731893097742, skeptic_NYH=0.02343329962747681, skeptic_Y_over_YQ=1.3060583336393292)
    val['rel_diff_NX_vs_skeptic'] = val['NX_400_1000'] / val['skeptic_NX_400_1000'] - 1
    sd, sc = series(rd['D']), series(rc['D'])
    np.savez(fc.SP / 'runs' / 'planck_background.npz', **{'d_' + k: v for k, v in sd.items()}, **{'c_' + k: v for k, v in sc.items()})
    sel = sc['t'] > 1600
    out = dict(method=__doc__, c=CPL, c_prime_at_slice=cprime, prep_effective=dict(M=pe['M'], h=pe['h']),
               prep=dict(z=[[c.real, c.imag] for c in zp], amps=np.abs(zp).tolist(), hessian_eig=st['eig'].tolist(),
                         F_inf=float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(zp) ** 2))), skeptic_amps=[96.47722349207996, 8.606647396734143, 29.016616720317433]),
               control_c0=dict(NX_400_1000=ctrl_nx, rptri_dq_NX_400_1000=ref['NX_400_1000'], rel=ctrl_nx / ref['NX_400_1000'] - 1),
               validation=val,
               late=dict(a_end=float(sc['a'][-1]), H_end=float(sc['H'][-1]), eta_end=float(sc['eta'][-1]), K_end=float(sc['a'][-1] ** 2 * sc['H'][-1]),
                         chi2_mean_late=float(np.mean((sc['a'][sel] * np.abs(sc['z'][sel, 0])) ** 2)), absS_late=float(np.mean(np.abs(sc['z'][sel, 2])))),
               seconds=time.time() - t0)
    out['prep']['F_vac'] = 6 * p['h']
    json.dump(out, open(RES / 'background_planck.json', 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != 'method'}, indent=1, default=float))


if __name__ == '__main__':
    main()
