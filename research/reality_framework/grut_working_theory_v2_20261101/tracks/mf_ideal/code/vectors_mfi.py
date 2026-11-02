#!/usr/bin/env python3
"""mf_ideal item 2 (Claude, own code): longitudinal hidden vectors on each idealised history's own recorded background.

Method [C] (inf_chain vectors.py, adapted; vmodes.c copied unchanged from tracks/inf_chain/code):
  1. Background table for t >= 0 (a = 1, H = 1 at the slice): the dense MODE 5 run (dt = 0.005, tau <= 100) gives a, H and
     F^2 = 2(|X|^2 + 9|Y|^2 + 36|Sigma|^2); afterwards the fluid is continued exactly (dust x_m e^{-gam t}/a^3 decaying
     into radiation, from the run's last rho_m, rho_r) until gam t = 60, then exact radiation; F^2 -> F_v^2 + 2 chi^2/a^2.
     For table backgrounds the continuation follows the table to its end, then exact radiation.
  2. Modes (vmodes.c): exact massless de Sitter Goldstone start at k/aH = 50 with the prepared F_s (constant-H inflation,
     so the start is exact), integrated through the recorded release, read out as the adiabatic number once Omega/H > 100
     (checked at 200).  Single-Higgs-direction Proca action (dq_heat and its skeptic).
  3. Horizon-exit factor: H is exactly constant and s frozen before the slice, so the frozen-angle factor of the dq_heat
     skeptic, f = [h_r F_s/F(exit)]^2, is exactly 1 for every mode (the static slice).
  4. N_A = int dlnk k^3 n_k / 2 pi^2; deposit m_A <a_A>, <a_A> = (pi K/2 Gamma_A)^{1/2}, Gamma_A = e^2 m_A/48 pi, m_A = e F_v;
     R_vec = N_A m_A <a_A>/E_X (E_X of the history's own release).
Usage: python3 -B vectors_mfi.py TAG [TAG ...]    Output: ../results/vectors/vectors_<TAG>.json
"""
import json, math, os, subprocess, sys
import numpy as np
from scipy.integrate import solve_ivp
import mfi_common as mc

VSRC = mc.HERE / 'vmodes.c'
VBIN = mc.SP / 'build' / 'vmodes'
ES = tuple(float(x) for x in os.environ.get('MFI_ES', '1e-7,1e-6,3e-6,1e-5,1e-4').split(','))
KLO = np.geomspace(1e-4, 1.0, 41)
KHI = np.geomspace(1.0, 30.0, 61)[1:]
KG = np.concatenate([KLO, KHI])
OUT = mc.RES / 'vectors'


def vbuild():
    if not VBIN.exists() or VBIN.stat().st_mtime < VSRC.stat().st_mtime:
        subprocess.run(['gcc', '-O2', '-o', str(VBIN), str(VSRC), '-lm'], check=True)


def fluid_late(t0, a0, xm0, rr0, gam, tmax=1e10):
    """dust (rho_m a^3 = xm, decaying at gam) + radiation (rho_r a^4 = rr) from t0; log-t grid."""
    def rhs(lt, y):
        t = math.exp(lt); la, xm, rr = y
        a = math.exp(la)
        H = math.sqrt((xm / a ** 3 + rr / a ** 4) / 3)
        return [t * H, -t * gam * xm, t * gam * xm * a]
    t1 = min(tmax, max(t0 * 10, 60.0 / gam if gam > 0 else t0 * 10))
    lts = np.linspace(math.log(t0), math.log(t1), 3000)
    s = solve_ivp(rhs, (lts[0], lts[-1]), [math.log(a0), xm0, rr0], t_eval=lts, method='DOP853', rtol=1e-11, atol=1e-14)
    t = np.exp(s.t); la = s.y[0]; xm = s.y[1]; rr = s.y[2]
    a = np.exp(la)
    H = np.sqrt((xm / a ** 3 + rr / a ** 4) / 3)
    return t, la, H


def build_table(tag):
    rel = json.load(open(mc.RES / 'release' / ('release_%s.json' % tag)))
    S = np.load(mc.RUNS / (tag + '_dense.npz'))['S']
    C = mc.COLS
    t, a, H = S[:, 0], S[:, 1], S[:, 2]
    X, Y, Sv = S[:, C.index('absX')], S[:, C.index('absY')], S[:, C.index('S_over_vev')]
    Sg = Sv * mc.P05['h'] / math.sqrt(2)
    F2 = 2 * (X ** 2 + 9 * Y ** 2 + 36 * Sg ** 2)
    Td = float(t[-1])
    late = t >= 0.75 * Td
    chi2 = float(np.mean(X[late] ** 2 * a[late] ** 2)); Fv2 = float(np.mean(72 * Sg[late] ** 2))
    rows = [(float(tt), math.log(aa), float(hh), float(ff)) for tt, aa, hh, ff in zip(t, a, H, F2)]
    ad, Hd = float(a[-1]), float(H[-1])
    gam = rel['spec']['gam_bg']
    if rel['spec']['bg'] == 'fluid' and gam > 0 and S[-1, C.index('rho_I')] > 0:
        xm0 = float(S[-1, C.index('rho_I')]) * ad ** 3; rr0 = float(S[-1, C.index('rho_r')]) * ad ** 4
        tl, la, Hl = fluid_late(Td, ad, xm0, rr0, gam)
        for tt, l, hh in zip(tl[1:], la[1:], Hl[1:]):
            rows.append((float(tt), float(l), float(hh), Fv2 + 2 * chi2 / math.exp(2 * l)))
        te, ae, He = float(tl[-1]), math.exp(la[-1]), float(Hl[-1])
        cont = 'fluid to gam t = 60, then radiation'
    elif rel['spec']['bg'].startswith('table'):
        # follow the MODE 5 table beyond the dense run, then exact radiation from the table's last row
        path = rel['spec']['bg'].split(':', 1)[1]
        T = np.loadtxt(path, skiprows=1, usecols=(0, 1, 2))
        T = T[T[:, 0] > Td * 1.0001]
        for tt, l, hh in T[::max(1, len(T) // 6000)]:
            rows.append((float(tt), float(l), float(hh), Fv2 + 2 * chi2 / math.exp(2 * l)))
        if len(T):
            te, ae, He = float(T[-1, 0]), math.exp(float(T[-1, 1])), float(T[-1, 2])
            if abs(rows[-1][0] - te) > 1e-12:
                rows.append((te, math.log(ae), He, Fv2 + 2 * chi2 / ae ** 2))
        else:
            te, ae, He = Td, ad, Hd
        cont = 'table to its end, then radiation'
    else:
        te, ae, He = Td, ad, Hd
        cont = 'radiation'
    K = ae * ae * He
    tl = np.geomspace(te * 1.0001, 1e10, 4000)
    al = np.sqrt(ae ** 2 + 2 * K * (tl - te))
    for tt, aa in zip(tl, al):
        rows.append((float(tt), math.log(aa), K / aa ** 2, Fv2 + 2 * chi2 / aa ** 2))
    path = mc.SP / ('vtable_%s.txt' % tag)
    with open(path, 'w') as fh:
        fh.write('%d\n' % len(rows))
        for r in rows:
            fh.write('%.17g %.17g %.17g %.17g\n' % r)
    return dict(K=K, K_release=rel['result']['K'], continuation=cont, FS=rel['prep']['F'], Fv2=Fv2, Fv=math.sqrt(Fv2), chi2=chi2,
                path=str(path), Td=Td, E_X=rel['result']['E_X'])


def run_modes(info, ks, e, AD=100.0, nproc=2):
    TR = info['Td']
    chunks = [ks[i::nproc] for i in range(nproc)]
    procs = []
    for c in chunks:
        args = [str(VBIN), info['path'], '%.17g' % info['FS'], '%.17g' % TR, '%.17g' % AD, 'D', '0']
        for k in c:
            args += ['%.17g' % k, '%.6g' % e]
        procs.append(subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
    out = {}
    for pr in procs:
        o, er = pr.communicate()
        if pr.returncode != 0:
            raise RuntimeError(er[:300])
        for line in o.split('\n'):
            if line.startswith('M '):
                v = line.split()
                out[float(v[1])] = dict(n1=float(v[3]), n2=float(v[4]), t1=float(v[5]), t2=float(v[6]), steps=int(v[7]))
    return [out[min(out, key=lambda x: abs(x - k))] for k in ks]


def assemble(ks, modes, e, info):
    n = np.array([m['n2'] for m in modes]); n1 = np.array([m['n1'] for m in modes])
    dN = ks ** 3 * n / (2 * math.pi ** 2)
    lk = np.log(ks)
    lo = ks <= 1.0; hi = ks >= 1.0
    NA_lo = float(np.trapezoid(dN[lo], lk[lo])); NA_hi = float(np.trapezoid(dN[hi], lk[hi]))
    mA = e * info['Fv']; GA = e * e * mA / (48 * math.pi)
    aA = math.sqrt(math.pi * info['K'] / (2 * GA))
    dep = mA * aA
    NA = NA_lo + NA_hi
    return dict(e=e, N_A=NA, N_A_k_le_1=NA_lo, N_A_k_ge_1=NA_hi, f_exit=1.0, m_A=mA, Gamma_A=GA, a_mean=aA, deposit_per_vector=dep,
                R_vec=NA * dep / info['E_X'],
                max_readout_change=float(np.max(np.abs(n - n1) / np.maximum(np.abs(n), 1e-300))),
                spectrum=[dict(k=float(k), n_k=float(x)) for k, x in zip(ks[::4], n[::4])])


def main(tag):
    vbuild()
    info = build_table(tag)
    out = dict(method=__doc__, tag=tag, background={k: v for k, v in info.items() if k != 'path'}, by_e=[])
    for e in ES:
        modes = run_modes(info, list(KG), e)
        r = assemble(KG, modes, e, info)
        r['steps_total'] = int(sum(m['steps'] for m in modes))
        out['by_e'].append(r)
        print('%s e=%.0e: N_A=%.5g (k<=1 %.5g; k>=1 %.4g) dep/E_X=%.4g R_vec=%.5g K=%.5f readout chg %.1e' % (
            tag, e, r['N_A'], r['N_A_k_le_1'], r['N_A_k_ge_1'], r['deposit_per_vector'] / info['E_X'], r['R_vec'], info['K'],
            r['max_readout_change']), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(OUT / ('vectors_%s.json' % tag), 'w'), indent=1, default=float)
    try:
        os.remove(info['path'])
    except OSError:
        pass


if __name__ == '__main__':
    for t in sys.argv[1:]:
        if (OUT / ('vectors_%s.json' % t)).exists() and not os.environ.get('MFI_RERUN'):
            print('skip', t); continue
        main(t)
