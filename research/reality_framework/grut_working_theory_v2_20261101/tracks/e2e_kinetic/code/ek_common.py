#!/usr/bin/env python3
"""e2e_kinetic (Claude, end-to-end pass of the GRUT working theory, GitHub issue #13 context): common module.

Every coupling, rate normalisation, renormalisation scale, start state and switch time below is a MODEL CHOICE; none
follows from closure.  Labels in comments: [Q] quoted (track or literature, arXiv numbers in NOTES.md), [C] computed by
this track, [E] estimate, [A] argued.

Physics constants and rates (units H_i = 1, the H at sigma = sigma_c; a = 1 there):
  * seesaw (mf_hybrid's completion, 4 singlets at m~ = 0.05 eV) [Q mf_hybrid]: h_j^2 = sum_i |h_ji|^2 = m~_j M_j / v^2,
    Gamma_rest,j = h_j^2 M_j / (8 pi), v = 174.1 GeV.  e2e_sm's completion enters ONLY through m~_j (the rates are linear
    in h_j^2), so every result is also reported against m~ (scan in ek_runs.py).
  * in-medium N damping [Q Besak-Bodeker arXiv:1202.1288, leading order]: total production rate
    Gamma_prod/(h^2 T^4) = G12 + G22, G22 = {y_t^2 c_Q + (3 g^2 + g'^2) [ln 1/(3 g^2 + g'^2) + c_V]}/(1536 pi),
    c_Q = 2.52, c_V = 3.17; G12 = 4.3e-4 (their Table 2 falls from 5.4e-4 at 1e8 GeV to 5.0e-4 at 6.7e9 GeV; extrapolated
    to 1e14-1e15 GeV, model choice +-10%).  Damping-rate form Gamma(k) = c_BB h^2 Tcal^2 / omega_k [model choice for the
    k dependence: t-channel / collinear ~ m_inf^2/k], normalised so that 2 int d^3k/(2pi)^3 Gamma f_FD = Gamma_prod.
  * Tcal^2 = T_th^2 + (12/73) n_h^2/rho_h : the medium scale (thermal part plus hard SM gas, = T^2 in equilibrium,
    73 = g_b + g_f/2 of the SM).  Asymptotic masses m_H^2 = kH Tcal^2, m_L^2 = kL Tcal^2 [Q standard HTL coefficients]
    close N -> L H when M < m_L + m_H.
  * SM couplings at mu = 1e15 GeV (~ 2 pi T) from higgs_reheat's 3-loop SM RGE (Buttazzo et al. arXiv:1307.3536 matching),
    run read-only at M_t = 172.57, alpha_s = 0.1180 [C, ek_common.couplings()].
  * SM thermalisation [Q Kurkela-Lu arXiv:1405.6318, AMY effective kinetic theory with LPM]: under-occupied
    t_eq = (34 + 21 ln(Q/T)) / (1 + 0.037 ln 1/lam) (Q/T)^{1/2} / (lam^2 T), lam = N_c g^2 (central: SU(2), the charge of
    the primaries L, H; variant SU(3)).  Used as the relaxation rate of the hard SM energy into the thermal bath.
"""
import json, math, os, subprocess, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRK = TRACK.parent
RES = TRACK / 'results'
SP = Path(os.environ.get('EK_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_kinetic'))
SP.mkdir(parents=True, exist_ok=True)
BIN = SP / 'ek_modes'
BENCH_DIR = TRK / 'mf_hybrid' / 'results'

MP_GEV = 2.435e18
H_REF = 1.30e11
VEW = 174.1
MNU_REF = 0.05e-9          # GeV, m~ of mf_hybrid's seesaw
GSTAR = 106.75
GSUM_F = 73.0

# ------------------------------------------------------------------ SM couplings (3-loop RGE of higgs_reheat, read-only)
_COUP = None


def couplings(mu=1.0e15):
    global _COUP
    if _COUP is not None and _COUP[0] == mu:
        return _COUP[1]
    sys.path.insert(0, str(TRK / 'higgs_reheat' / 'code'))
    import sm_rge as S                                    # read-only import (python -B: no bytecode written)
    c = S.matching(172.57, 125.20, 0.1180, MW=80.369)
    r = S.run(c, 172.57, mu_max=1e19)
    k = S.couplings_at(r, mu)
    out = dict(mu=mu, g1=float(k['g1']), gY=float(math.sqrt(3 / 5) * k['g1']), g2=float(k['g2']), g3=float(k['g3']),
               yt=float(k['yt']), lam=float(k['lam']))
    _COUP = (mu, out)
    return out


def sm_coeffs(mu=1.0e15, G12=4.3e-4, lam_group='SU2'):
    c = couplings(mu)
    g2, gY, yt, lam = c['g2'], c['gY'], c['yt'], c['lam']
    s = 3 * g2 ** 2 + gY ** 2
    G22 = (yt ** 2 * 2.52 + s * (math.log(1 / s) + 3.17)) / (1536 * math.pi)
    neq, inv_k = 2 * 0.75 * 1.2020569 / math.pi ** 2, (math.pi ** 2 / 12) / (1.5 * 1.2020569)
    cBB = (G12 + G22) / (neq * inv_k)
    kH = 3 * g2 ** 2 / 16 + gY ** 2 / 16 + yt ** 2 / 4 + lam / 2
    kL = s / 16
    lam_eff = {'SU2': 2 * g2 ** 2, 'SU3': 3 * c['g3'] ** 2}[lam_group]
    return dict(couplings=c, G12=G12, G22=G22, cBB=cBB, kH=kH, kL=kL, lam_eff=lam_eff, lam_group=lam_group,
                mean_damping_over_h2T=(G12 + G22) / neq)


# ------------------------------------------------------------------ benchmark model (mf_hybrid bench, read-only)
def load_bench(name):
    return json.load(open(BENCH_DIR / ('bench_%s.json' % name)))


def bench_state(name, t0=0.19):
    """state at t0 from the exported homogeneous history (fractions sum to 1)."""
    d = load_bench(name)
    h = d['history']
    T = np.array(h['t'], float)
    g = lambda k: float(np.interp(t0, T, np.array([np.nan if x is None else x for x in h[k]], float)))
    return dict(t=t0, a=g('a'), H=g('H'), f_sig=g('f_sig'), f_chi=g('f_chi'), f_rad=g('f_rad'), f_SM=g('f_SM'))


def species(name, H_GeV=H_REF, mnu=MNU_REF):
    """4 axial channels [Q mf_hybrid params]: sigma -> N1 and chi -> N2, N3, N4. mnu: m~ (GeV), scalar or 4-list."""
    pr = load_bench(name)['params']
    mn = mnu if isinstance(mnu, (list, tuple)) else [mnu] * 4
    out = []
    chans = [(0, c) for c in pr['axial_sigma']] + [(1, c) for c in pr['axial_chi']]
    mpar = [pr['m_sig_over_Hi'], pr['m_chi_over_Hi']]
    for j, (par, c) in enumerate(chans):
        M, f, m = c['M_over_Hi'], c['f_over_Hi'], mpar[par]
        h2 = mn[j] * M * H_GeV / VEW ** 2
        Grest = h2 * M / (8 * math.pi)
        beta2 = 1 - 4 * M * M / (m * m)
        kres = 0.5 * math.sqrt(beta2) if beta2 > 0 else 0.0
        Gp = (M / f) ** 2 * m * math.sqrt(beta2) / (16 * math.pi) if beta2 > 0 else 0.0
        out.append(dict(par=par, M=M, f=f, m=m, Mh=M / m, f_over_m=f / m, h2=h2, Grest=Grest, kres_m=kres, Gpert=Gp,
                        mnu=mn[j]))
    return out, mpar


# ------------------------------------------------------------------ p tables (ek_ptable.py writes SP/ptable_<tag>.npz)
def load_ptable(Mh):
    files = sorted(SP.glob('ptable_*.npz'))
    if not files:
        raise FileNotFoundError('no ptable in %s: run ek_ptable.py first' % SP)
    fn = min(files, key=lambda f: abs(float(f.stem.split('_')[1]) - Mh))
    if abs(float(fn.stem.split('_')[1]) - Mh) > 1e-3 * Mh:
        raise FileNotFoundError('no ptable for M/m = %.5f (nearest %s): run ek_ptable.py %.5f' % (Mh, fn.name, Mh))
    z = np.load(fn)
    return z['lA'], z['pp']


# ------------------------------------------------------------------ C engine input / output
def compile_engine():
    src = HERE / 'ek_modes.c'
    if not BIN.exists() or BIN.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-march=native', '-o', str(BIN), str(src), '-lm'], check=True)


def grid_linear(qmax, nq):
    q = (np.arange(nq) + 0.5) * (qmax / nq)
    w = q ** 2 * (qmax / nq) / (2 * math.pi ** 2)
    return q, w


def grid_log(qlo, qhi, nq):
    e = np.geomspace(qlo, qhi, nq + 1)
    q = np.sqrt(e[1:] * e[:-1])
    w = (e[1:] ** 3 - e[:-1] ** 3) / 3 / (2 * math.pi ** 2)
    return q, w


def write_input(path, st, spec, mpar, sm, flags, grids, occ, bloch=None, frozen=None, tables=None, Apert=1.0, cdt=0.04,
                ckcut=1.3, KLscale=1.0, winv=1.0, dshell=0.0):
    """st: dict t, t_end, a, mP, rho_sig, rho_chi, rho_h, n_h, rho_th (physical, H_i units)."""
    L = []
    L.append('%d %d %d %d' % (flags.get('coherent', 1), flags.get('noexp', 0), flags.get('fixamp', 0), flags.get('nocoll', 0)))
    L.append('%.17g %.17g %.17g %.17g' % (st['mP'], st['t'], st['t_end'], st['a']))
    ev = st.get('events', (st['t'], st['t'], -1.0, -1.0, -1.0, -1.0))          # tph0 tph1 tnext0 tnext1 tlastc0 tlastc1
    L.append(' '.join('%.17g' % x for x in ev))
    L.append('%.17g %.17g %.17g %.17g' % (mpar[0], mpar[1], st['rho_sig'], st['rho_chi']))
    L.append('%.17g %.17g %.17g' % (st['rho_h'], st['n_h'], st['rho_th']))
    L.append('%.17g %.17g %.17g %.17g %.17g %.17g' % (sm['lam_eff'], GSTAR, sm['kH'], sm['kL'], sm['cBB'], GSUM_F))
    L.append('%.17g %.17g %.17g %.17g %.17g' % (34.0, 21.0, 0.037, KLscale, winv))
    L.append('%.17g %.17g %.17g %.17g' % (Apert, cdt, ckcut, dshell))
    for j, s in enumerate(spec):
        q, w = grids[j]
        L.append('%d %.17g %.17g %.17g %.17g %.17g %.17g' % (s['par'], s['M'], s['f'], s['Grest'], s['h2'], s['kres_m'], s['Gpert']))
        L.append('%d' % len(q))
        L.extend('%.17g %.17g' % (a_, b_) for a_, b_ in zip(q, w))
        L.extend('%.12g' % x for x in occ[j])
        if bloch is not None and bloch[j] is not None:
            L.append('1'); L.extend('%.12g' % x for x in bloch[j])
        else:
            L.append('0')
        fz = frozen[j] if frozen is not None and frozen[j] is not None else np.zeros(2 * len(q), int)
        L.extend('%d' % x for x in fz)
        if tables is not None and tables[j] is not None:
            lA, pp = tables[j]
            L.append('%d %d' % (len(lA), pp.shape[1]))
            L.extend('%.12g' % x for x in lA)
            L.extend('%.8g' % x for x in pp.ravel())
        else:
            L.append('0 0')
    open(path, 'w').write('\n'.join(L) + '\n')


REC_COLS = ['t', 'a', 'H', 'f_sig', 'f_chi', 'f_N', 'f_hard', 'f_th', 'T_th', 'Tcal2', 'T_SM', 'tau_KL']
for _j in range(4):
    REC_COLS += ['fN%d' % (_j + 1), 'Gband%d' % (_j + 1), 'Gall%d' % (_j + 1), 'nband%d' % (_j + 1), 'A%d' % (_j + 1), 'nN%d' % (_j + 1)]
REC_COLS += ['n_h', 'rho', 'rho_a4']


def read_rec(prefix):
    a = np.loadtxt(str(prefix) + '_rec.txt', ndmin=2)
    return {c: a[:, i] for i, c in enumerate(REC_COLS)}


def read_state(prefix, nqs):
    L = open(str(prefix) + '_state.txt').read().split('\n')
    head = list(map(float, L[0].split()))
    st = dict(zip(['t', 'a', 'H', 'rho_sig', 'rho_chi', 'rho_h', 'n_h', 'rho_th', 'Ecoll_com', 'Eprod_com', 'rho_tot'], head))
    i = 1
    occ, fz = [], []
    for j in range(4):
        nq = int(L[i]); i += 1
        arr = np.array([L[i + k].split() for k in range(2 * nq)], float); i += 2 * nq
        occ.append(arr[:, 0]); fz.append(arr[:, 1].astype(int))
    st['occ'], st['frozen'] = occ, fz
    return st


def read_bloch(prefix, nqs):
    L = open(str(prefix) + '_bloch.txt').read().split()
    ev = tuple(float(x) for x in L[:6])
    vals = np.array(L[6:], float)
    out, i = [], 0
    for nq in nqs:
        out.append(vals[i:i + 6 * nq]); i += 6 * nq
    return ev, out


def run_engine(inp, prefix, timeout=590):
    compile_engine()
    r = subprocess.run([str(BIN), 'run', str(inp), str(prefix)], capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return read_rec(prefix)


def rebin(q_src, w_src, n_src, q_dst_edges):
    """occupation re-binning (both helicities stacked [h*nq + i]), number-conserving within each destination bin."""
    nq = len(q_src)
    out = []
    w_dst = (q_dst_edges[1:] ** 3 - q_dst_edges[:-1] ** 3) / 3 / (2 * math.pi ** 2)
    dq = np.diff(np.r_[0.0, 0.5 * (q_src[1:] + q_src[:-1]), q_src[-1] + 0.5 * (q_src[-1] - q_src[-2])])
    edges_src = np.r_[0.0, np.cumsum(dq)]
    for h in range(2):
        nn = n_src[h * nq:(h + 1) * nq]
        cum = np.r_[0.0, np.cumsum(w_src * nn)]                 # number below edges_src
        Nd = np.interp(q_dst_edges, edges_src, cum, right=cum[-1])
        out.append(np.diff(Nd) / w_dst)
    return np.concatenate(out)


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(path, 'w'), indent=1, default=float)
