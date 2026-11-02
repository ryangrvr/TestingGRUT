#!/usr/bin/env python3
"""inf_reheat track (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): shared constants and the driver
of rprh.c (our engine: copy of tracks/dq_entropy/code/rpdq.c with change RH1, an optional T-model plateau).

Reuse, stated: model constants, preparation (fx_common.stationary), the trilinear coupling convention
(fx_common.trilinear_G), the yield convention (fx_common.yield_over_target) and the decay ledger (fx_common.ledger) are
imported READ-ONLY from tracks/fx_holder/code/fx_common.py (through tracks/dq_entropy/code/dq_common.py for the portal
helper).  run() is a copy of dq_common.run with the extra engine argument wT.  Labels: [C] for engine output.
"""
import math, os, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACKS = HERE.parent.parent
sys.path.insert(0, str(TRACKS / 'fx_holder' / 'code'))
import fx_common as fc  # noqa: E402  (read-only reuse)

SP = Path(os.environ.get('IR_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_reheat'))
BIN = SP / 'build' / 'rprh'
SRC = HERE / 'rprh.c'
RES = HERE.parent / 'results'
COLS = fc.COLS

# ---------------------------------------------------------------------------------------------- constants (GeV)
MP = 2.435e18                  # reduced Planck mass (engine convention)
H_M1 = fc.H_GEV                # 1.3428e11, the engine's H_i (dimensionless histories are H_i-independent)
H_CAL_FREEX = 1.291e11         # free-X calibration on (3, 0.5) [Q dq_critic hi_calibration: 1.290-1.292e11]
H_CAL_QM = 1.49e11             # Q-matter calibration [Q dq_critic: 1.483-1.499e11]
YQ_REF_3_05 = 0.991023         # Y/Y_Q of the (3, 0.5) release at H_M1, frozen value [Q dq_entropy]
NX_REF_3_05 = 8047.84          # |N_X| frozen, (3, 0.5) [Q dq_entropy]
EX_REF_3_05 = 188737.7         # conformal E_X, (3, 0.5) [Q dq_entropy]
K_REF = 1.9858052754013398     # a^2 H in radiation, (3, 0.5) [Q dq_heat, skeptic sk_K]
AS = 2.10e-9                   # Planck 2018 [R]
R_BOUND = 0.036                # BICEP/Keck 2021 (2110.00483) [R]
GSTAR = 106.75
VEW = 174.1                    # <H> = v/sqrt(2) [R]
MNU = 0.05e-9                  # atmospheric-scale neutrino mass, GeV [R]
KAPPA_T = 0.11                 # Higgs thermal mass^2 / T^2 [Q higgs_reheat]
NB_OBS = 8.7e-11               # n_B/s [R, Planck 2018]


def r_of_H(H):
    """tensor ratio for H at CMB exit: r = 2 H^2 / (pi^2 A_s M_P^2)."""
    return 2 * H ** 2 / (math.pi ** 2 * AS * MP ** 2)


def H_of_r(r):
    return MP * math.sqrt(r * math.pi ** 2 * AS / 2)


def T_of_rho(rho):
    return (30 * rho / (math.pi ** 2 * GSTAR)) ** 0.25


def T_RH(Gamma_GeV):
    """T at H = Gamma (3 M_P^2 Gamma^2 = pi^2 g T^4 / 30): 0.546 sqrt(Gamma M_P) for g = 106.75."""
    return (90 / (math.pi ** 2 * GSTAR)) ** 0.25 * math.sqrt(Gamma_GeV * MP)


# ---------------------------------------------------------------------------------------------- engine driver
def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-o', str(BIN), str(SRC), '-lm'], check=True)


def inflaton_I0(z, u, p, b, wT=0.0):
    """I0/M_r at rest with H = H_i: V(I0) = need = 3 f - eta2 V_h (as fx_common.inflaton_I0, plus the RH1 plateau)."""
    V, _ = fc.potential(z, u, p)
    f = 1 - 2 * p['xi'] * fc.ETA2 * abs(z[0]) ** 2
    need = 3 * f - fc.ETA2 * V
    if wT <= 0:
        return float(np.sqrt(2 * need) / b)
    s = np.sqrt(2 * need) / (b * wT)
    if s >= 1:
        raise ValueError('plateau below the release energy')
    return float(wT * np.arctanh(s))


def trilinear_G(z, p, b, wT=0.0, u_target=fc.U_PREP):
    """G = g_I e^{i alpha_I} with mu + G phihat(0) = 10 i at the slice (g_I re-fixed for each inflaton history)."""
    I0 = inflaton_I0(z, u_target, p, b, wT)
    ph0 = I0 / np.sqrt(fc.ETA2)
    return (u_target - p['mu']) / ph0, I0, ph0


def prepare(lam=0.5, b=3.0, wT=0.0):
    """M1r held preparation at u = 10 i (lambda_Y = lambda_Sigma = 0.5, portals vanish at M_Y), G re-fixed."""
    p = fc.params(lam)
    st = fc.stationary(fc.U_PREP, p)
    G, I0, ph0 = trilinear_G(st['z'], p, b, wT)
    return p, st, G, dict(I0=I0, phihat0=ph0, g_I=abs(G), alpha_I=float(np.angle(G)))


def run(tag, p, b, gam, z0, G=0j, cport=0.0, end=1000.0, rtol=1e-10, atol=1e-13, nsamp=5000, binw=20.0, cbin=0.0, wT=0.0, v0=None):
    build()
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), '0', g(b), g(gam), g(p['xi']), g(fc.ETA2), g(p['cH']), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']),
            g(p['h']), g(p['mu']), g(0.0), g(0.0), g(G.real), g(G.imag), g(cport), g(end), g(rtol), g(atol), str(int(nsamp)), g(binw), g(cbin), g(wT)]
    v0 = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)
    args += ['F'] + [g(v) for c in z0 for v in (c.real, c.imag)] + [g(v) for c in v0 for v in (c.real, c.imag)]
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(tag + ': ' + out.stderr)
    S, Bn, W, Y, END, I0 = [], [], {}, None, None, None
    for line in out.stdout.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('B '):
            Bn.append([float(x) for x in line.split()[2:]])
        elif line.startswith('W') or line.startswith('A3040'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('Y '):
            Y = [float(v) for v in line.split()[1:]]
        elif line.startswith('I0 '):
            v = [float(x) for x in line.split()[1:]]; I0 = dict(I0=v[0], phihat0=v[1], u0=[v[2], v[3]], R0=v[4])
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]),
                                        max_hidden_work=float(v[5]), min_mL2=float(v[6]))
    S = np.array(S); Bn = np.array(Bn)
    (SP / 'runs').mkdir(parents=True, exist_ok=True)
    np.savez_compressed(SP / 'runs' / (tag + '.npz'), S=S, bins=Bn)
    return dict(tag=tag, S=S, bins=Bn, W=W, Y=Y, END=END, I0=I0, seconds=time.time() - t0, b=b, gam=gam, end=end, wT=wT,
                G=[G.real, G.imag], rtol=rtol)


def bin_col(rr, name):
    return rr['bins'][:, 2 + COLS.index(name)]


def window_mean(rr, name, lo, hi=None):
    B = rr['bins']; hi = rr['end'] if hi is None else hi
    t0, t1 = B[:, 0], B[:, 1]; vals = bin_col(rr, name); dt = t1 - t0
    sel = (t0 >= lo - 1e-9) & (t1 <= hi + 1e-9)
    return float(np.sum(vals[sel] * dt[sel]) / np.sum(dt[sel]))


def last(rr):
    return dict(zip(COLS, rr['S'][-1]))


def dilution(rr):
    """D = (rho_r a^4 / 3)^{3/4} at the end of the run (comoving entropy relative to an instant release) and K = a^2 H."""
    L = last(rr)
    return float((L['rho_r'] * L['a'] ** 4 / 3) ** 0.75), float(L['a'] ** 2 * L['H'])
