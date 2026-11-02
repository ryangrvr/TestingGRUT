#!/usr/bin/env python3
"""mf_ideal track (Claude, own code; issue #13 multi-field inflaton pass): common module and driver of rpmfi.c.

rpmfi.c is our copy of tracks/inf_chain/code/rpchain.c (<- inf_potential rpinf.c <- dq_entropy rpdq.c <- fx_holder rptri.c)
with MODE 5 (external background + holder on a spectator s; see the header of rpmfi.c).

Reuse, stated (READ-ONLY imports; nothing in those tracks is modified or executed for writing):
  * tracks/dq_heat/skeptic/code/dq_model_copy.py : M1r hidden parameters (lambda_Y = lambda_Sigma = 0.5, M_Y = 30, h = 40,
    mu = 7.159456, kappa, xi_X = -5/6) and the exact stationary-point solver of the held preparation (fc.stationary),
    fc.trilinear_G (the MODE 0 inflaton trilinear).
Every history, ramp shape, decay rate and slice convention below is a MODEL CHOICE; none follows from closure.

Units: H_i = H at the slice = 1, a = 1 at the slice, tau = H_i t, hidden fields in H_i.  In MODE 5 the holder is
u = mu + G_eff c with c = s/s_*, G_eff = u_slice - mu (= 10i - mu for the M1r calibration).
"""
import json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]                                   # .../tracks
sys.path.insert(0, str(TRK / 'dq_heat' / 'skeptic' / 'code'))
import dq_model_copy as fc                              # noqa: E402  read-only

SP = Path(os.environ.get('MFI_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_ideal'))
BIN = SP / 'build' / 'rpmfi'
SRC = HERE / 'rpmfi.c'
RES = HERE.parent / 'results'
RUNS = SP / 'runs'
P05 = fc.params(0.5)
MU = P05['mu']
U_PREP = 10j
COLS = fc.COLS
ETA2 = fc.ETA2                                          # (H_i/M_P)^2 at the M1 H_i; enters only at O(1e-15)
MP = 2.435e18
H_M1 = 1.342842735e11                                   # M1 reference H_i [Q]
YQ_TARGET = 5.07941423e-9                               # [Q] M1 calibration target at E/Q = 0.9374 m_X
GSTAR = 106.75
C_HART = 6.977568                                       # N_c/E^{3/4} [Q dq_entropy]
Q = np.array([1.0, 3.0, 6.0])


def F_of(z):
    """hidden Higgs norm F = (2 sum q_i^2 |z_i|^2)^{1/2}."""
    z = np.asarray(z, complex)
    return float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(z) ** 2)))


def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O2', '-o', str(BIN), str(SRC), '-lm'], check=True)


def prep(u=U_PREP, R=12.0, guess=None):
    """held preparation: stationary point of V + xi R |X|^2 at fixed u (constant-H inflation: R = 12 H_i^2)."""
    st = fc.stationary(u, P05, R=R, guess=guess)
    if st['grad'] > 1e-8:
        st = fc.stationary(u, P05, R=R, guess=None)
    return st


# ------------------------------------------------------------------------------------------------ engine driver
def args_for(mode, b, gam, G, end, nsamp, z, binw=0.0, cbin=0.0, rtol=1e-10, atol=1e-13, cport=0.0, p=P05, eta2=ETA2):
    g = lambda v: '%.17g' % float(v)
    a = [str(BIN), str(mode), g(b), g(gam), g(p['xi']), g(eta2), g(p.get('cH', 0) if mode == 2 else 0), g(p['lx']), g(p['ly']),
         g(p['ls']), g(p['kap']), g(p['M']), g(p['h']), g(p['mu']), '0', '0', g(G.real), g(G.imag), g(cport), g(end), g(rtol),
         g(atol), str(int(nsamp)), g(binw), g(cbin), 'F']
    for c in z:
        a += [g(c.real), g(c.imag)]
    a += ['0'] * 6
    return a


def parse(txt):
    S, B, C, D, W, Y, END, I0, ETA = [], [], [], [], {}, None, None, None, None
    for line in txt.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('B '):
            B.append([float(x) for x in line.split()[2:]])
        elif line.startswith('C '):
            C.append([float(x) for x in line.split()[2:]])
        elif line.startswith('D '):
            D.append([float(x) for x in line.split()[1:]])
        elif line.startswith('W') or line.startswith('A3040'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('Y '):
            Y = [float(v) for v in line.split()[1:]]
        elif line.startswith('ETA '):
            ETA = float(line.split()[1])
        elif line.startswith('I0 '):
            v = [float(x) for x in line.split()[1:]]; I0 = dict(I0=v[0], phihat0=v[1], u0=[v[2], v[3]], R0=v[4])
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]),
                                        max_hidden_work=float(v[5]), min_mL2=float(v[6]))
    return dict(S=np.array(S), B=np.array(B), C=np.array(C), D=np.array(D), W=W, Y=Y, END=END, I0=I0, eta_end=ETA)


def run(tag, mode, b, gam, z, G, end, nsamp, env=None, binw=0.0, cbin=0.0, capT=None, rtol=1e-10, atol=1e-13, dump=False,
        save=True, binary=None, cport=0.0):
    build()
    e = dict(os.environ)
    for k in [k for k in e if k.startswith('MFI_') and k != 'MFI_SCRATCH']:
        del e[k]
    e.update({k: str(v) for k, v in (env or {}).items()})
    if capT:
        e['RPCH_CAPT'] = '%.17g' % capT
    if dump:
        e['MFI_DUMP'] = '1'
    a = args_for(mode, b, gam, G, end, nsamp, z, binw=binw, cbin=cbin, rtol=rtol, atol=atol, cport=cport)
    if binary:
        a[0] = str(binary)
    t0 = time.time()
    r = subprocess.run(a, capture_output=True, text=True, env=e)
    if r.returncode != 0:
        raise RuntimeError('%s: rc %d %s' % (tag, r.returncode, r.stderr[-800:]))
    out = parse(r.stdout)
    out.update(tag=tag, seconds=time.time() - t0, stdout=r.stdout if not save else None, stderr=r.stderr, end=end, gam=gam,
               env={k: v for k, v in (env or {}).items()})
    if save:
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RUNS / (tag + '.npz'), S=out['S'], B=out['B'], C=out['C'])
        json.dump(dict(tag=tag, END=out['END'], I0=out['I0'], W=out['W'], Y=out['Y'], seconds=out['seconds'], end=end, gam=gam,
                       env=out['env'], eta_end=out['eta_end']), open(RUNS / (tag + '_meta.json'), 'w'), indent=1, default=float)
    return out


def col(S, name):
    return S[:, COLS.index(name)]


def bin_mean(B, name, lo, hi, off=0):
    """mean over physical-time bins [t0, t1] inside [lo, hi] (B rows: t0 t1 <NQ means>)."""
    t0, t1 = B[:, 0], B[:, 1]; v = B[:, 2 + COLS.index(name)]; dt = t1 - t0
    s = (t0 >= lo - 1e-9) & (t1 <= hi + 1e-9)
    return float(np.sum(v[s] * dt[s]) / np.sum(dt[s]))


# ------------------------------------------------------------------------------------------------ idealised backgrounds
def fluid_K(gam, x_m=1.0, tmax_fac=60.0):
    """K = lim a^2 H for dust (fraction x_m of rho = 3 at a = 1, H = 1) decaying into radiation at rate gam
    (3H^2 = rho_m + rho_r, M_P = H_i = 1).  gam = inf or x_m = 0 -> pure radiation, K = 1."""
    if x_m == 0 or not np.isfinite(gam):
        return 1.0

    # independent variable N = ln a; y = (t, ln xm, rr) with xm = rho_m a^3, rr = rho_r a^4
    def rhs(N, y):
        t, lxm, rr = y
        a = math.exp(N); xm = math.exp(lxm)
        H = math.sqrt((xm / a ** 3 + rr / a ** 4) / 3)
        return [1 / H, -gam / H, gam * xm * a / H]

    def done(N, y):
        return gam * y[0] - tmax_fac
    done.terminal = True
    s = solve_ivp(rhs, (0, 200.0), [0.0, math.log(3 * x_m), 3 * (1 - x_m)], method='DOP853', rtol=1e-12, atol=1e-14, events=done)
    N = s.t[-1]; t, lxm, rr = s.y[:, -1]
    a = math.exp(N); xm = math.exp(lxm)
    H = math.sqrt((xm / a ** 3 + rr / a ** 4) / 3)
    return float(a * a * H)


def gam_for_K(K, x_m=1.0):
    """decay rate of the dust that gives K (K = 1: pure radiation, returned as inf)."""
    if K <= 1.0 + 1e-12:
        return float('inf')
    return float(math.exp(brentq(lambda lg: fluid_K(math.exp(lg), x_m) - K, math.log(1e-3), math.log(1e2), xtol=1e-12)))


def history_env(K, ms, ton, gs=0.0, br=0.0, c0=1.0, mi2=0.0):
    """MODE 5 environment and the engine 'gamma' (dust decay) for an idealised fast history."""
    g = gam_for_K(K)
    env = dict(MFI_MS2='%.17g' % (ms * ms), MFI_MI2='%.17g' % mi2, MFI_TON='%.17g' % ton, MFI_GS='%.17g' % gs,
               MFI_BR='%.17g' % br, MFI_C0='%.17g' % c0)
    if not np.isfinite(g):
        env['MFI_FM'] = '0'                              # pure radiation from the slice
        g = 0.0
    return env, g


def br_of(s_star_GeV, Hi=H_M1):
    """back-reaction coefficient (H_i/s_*)^2 of the holder on s."""
    return (Hi / s_star_GeV) ** 2
