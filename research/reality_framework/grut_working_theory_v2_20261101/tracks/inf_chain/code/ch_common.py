#!/usr/bin/env python3
"""inf_chain track (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): common module.

Benchmarks, slice state, re-solved M1r preparation at the slice, holder re-fixing, and the driver of the engine copy
rpchain.c (MODE 4, explicit inflaton potential).

Reuse, stated (all READ-ONLY imports; nothing in those tracks is executed for writing or modified):
  * tracks/inf_potential/code/inf_common.py, calibrate.py : potentials, background ODE, Mukhanov-Sasaki spectra, N_*.
  * tracks/dq_heat/skeptic/code/dq_model_copy.py          : M1r hidden parameters (lambda = 0.5) and the exact
                                                            stationary-point solver of the held preparation.
  * tracks/inf_potential/results/bench_*.json              : the HT6t and T1 benchmark records (shape parameters, V0).
Every benchmark, slice convention, holder calibration point and width is a MODEL CHOICE; none follows from closure.

Units: H_i = H at the slice = 1, a = 1 at the slice, tau = H_i t; hidden fields in H_i; inflaton I in reduced-Planck
units M_P = 2.435e18 GeV; eta2 = (H_i/M_P)^2.
"""
import json, math, os, subprocess, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]                               # .../tracks
sys.path.insert(0, str(TRK / 'inf_potential' / 'code'))
sys.path.insert(0, str(TRK / 'dq_heat' / 'skeptic' / 'code'))
import inf_common as ic                              # noqa: E402  read-only
import calibrate as cal                              # noqa: E402  read-only
import dq_model_copy as fc                           # noqa: E402  read-only

SP = Path(os.environ.get('CH_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_chain'))
BIN = SP / 'build' / 'rpchain'
SRC = HERE / 'rpchain.c'
RES = HERE.parent / 'results'
RUNS = SP / 'runs'
MP = ic.MPL_GEV
Q = np.array([1.0, 3.0, 6.0])
P05 = fc.params(0.5)                                 # M1r hidden sector, lambda_Y = lambda_Sigma = 0.5
MU = P05['mu']
U_PREP = 10j
COLS = fc.COLS
H_M1 = 1.342842735e11                                # M1 reference H_i [Q]
YQ_TARGET = 5.07941423e-9                            # [Q] M1 calibration target at E/Q = 0.9374 m_X
GSTAR = 106.75
C_HART = 6.977568                                    # N_c / E^{3/4} [Q dq_entropy]


def F_of(z):
    """hidden Higgs norm F = (2 sum q_i^2 |z_i|^2)^{1/2} (q = 1, 3, 6)."""
    z = np.asarray(z, complex)
    return float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(z) ** 2)))


def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O2', '-o', str(BIN), str(SRC), '-lm'], check=True)


# ------------------------------------------------------------------------------------------------ benchmarks
class Bench:
    """A benchmark inflaton (model choice) with its slice.  kind: 'HT', 'T', 'E' (inf_common.Model)."""

    def __init__(self, name, model, V0, Nslice=1.0, Ncal=None, gam=0.5, note=''):
        self.name, self.m, self.V0, self.Nslice = name, model, V0, Nslice
        self.Ncal = Nslice if Ncal is None else Ncal       # holder calibration point (mu_eff = 10i there)
        self.gam = gam
        self.note = note
        self.bg = ic.background(model, 75.0)
        f, fN = ic.bg_at(self.bg, Nslice)
        h, eps = ic.h_of(model, f, fN)
        self.phi_s, self.phiN_s, self.h_s, self.eps_s = f, fN, h, eps
        self.Hi = math.sqrt(V0) * h * MP                    # H at the slice (GeV)
        self.eta2 = (self.Hi / MP) ** 2
        self.pmin = model.phi_min()
        self.phiI_s = f - self.pmin                         # holder field at the slice (M_P)
        fc_, _ = ic.bg_at(self.bg, self.Ncal)
        self.phiI_cal = fc_ - self.pmin
        # holder: mu + G phihat_I(cal) = 10 i  ->  at the slice u = mu + (10i - mu) phi_I(s)/phi_I(cal)
        self.G = (U_PREP - MU) / (self.phiI_cal / math.sqrt(self.eta2))
        self.u_s = MU + (U_PREP - MU) * self.phiI_s / self.phiI_cal
        self.R_s = 6.0 * (2.0 - eps)                        # Ricci scalar at the slice (H_i units)
        self.mI = math.sqrt(model.m2_min() * V0) * MP / self.Hi   # m_I / H_i at the minimum
        self.gI = abs(self.G)

    def kind_code(self):
        return {'T': 1, 'E': 2, 'HT': 3}[self.m.kind]

    def env(self):
        e = dict(os.environ, RPINF_KIND=str(self.kind_code()), RPINF_AMP='%.17g' % (self.V0 / self.eta2),
                 RPINF_PHIMIN='%.17g' % self.pmin, RPINF_I0='%.17g' % self.phiI_s,
                 RPINF_SIGN='%d' % (1 if self.phiN_s > 0 else -1))
        if self.m.kind == 'HT':
            e.update(RPINF_P1='%.17g' % self.m.vev, RPINF_P2='%d' % self.m.pw, RPINF_P3='%.17g' % self.m.beta,
                     RPINF_FTOP='%.17g' % self.m.ftop)
        else:
            e.update(RPINF_P1='%.17g' % self.m.M)
        return e

    def prep(self, guess=None):
        """held preparation re-solved at the slice: stationary point at u = u_s and R = 6 (2 - eps_H) H_i^2."""
        st = fc.stationary(self.u_s, P05, R=self.R_s, guess=guess)
        if st['grad'] > 1e-8:
            st = fc.stationary(self.u_s, P05, R=self.R_s, guess=None)
        return st

    def history(self, Ngrid):
        """pre-slice path: s = phi_I/phi_I(cal), h_r = H/H_i, eps_H at N e-folds before the slice."""
        out = []
        for N in Ngrid:
            f, fN = ic.bg_at(self.bg, N + self.Nslice)
            h, eps = ic.h_of(self.m, f, fN)
            out.append(((f - self.pmin) / self.phiI_cal, h / self.h_s, eps))
        return out

    def summary(self):
        return dict(name=self.name, kind=self.m.kind, params=self.m.p, V0_MP4=self.V0, H_i_GeV=self.Hi, eta2=self.eta2,
                    N_slice_before_end=self.Nslice, N_cal_before_end=self.Ncal, eps_H_slice=self.eps_s, R_slice=self.R_s,
                    phi_slice=self.phi_s, phiI_slice=self.phiI_s, phiI_cal=self.phiI_cal, phiN_slice=self.phiN_s,
                    G=[self.G.real, self.G.imag], g_I=self.gI, alpha_I=math.atan2(self.G.imag, self.G.real),
                    u_slice=[self.u_s.real, self.u_s.imag], m_I_over_Hi=self.mI, Gamma_I_over_Hi=self.gam, note=self.note)


def bench_from_inf_potential(name, Hi_index=0, gam=0.5, Ncal=None):
    rec = json.load(open(TRK / 'inf_potential' / 'results' / ('bench_%s.json' % name)))['by_Hi'][Hi_index]
    kind = {'T1': 'T', 'E1': 'E', 'HT6t': 'HT'}[name]
    m = ic.Model(kind, **rec['params'])
    tag = name + ('' if Hi_index == 0 else '_H%.2f' % (rec['H_i_GeV'] / 1e11))
    b = Bench(tag, m, rec['V0_MP4'], Nslice=1.0, Ncal=Ncal, gam=gam, note='inf_potential bench_%s.json by_Hi[%d]' % (name, Hi_index))
    b.rec = rec
    return b


def tmodel_bench(alpha, Nstar=None, Nslice=1.0, gam=0.0, name=None):
    """alpha-attractor T-model (inf_reheat's stand-in family) with V0 fixed by A_s at the pivot; H_i := H(slice)."""
    m = ic.Model('T', alpha=alpha)
    bg = ic.background(m, 75.0)
    Ns = 53.4 if Nstar is None else Nstar
    sp = ic.spectrum_at_pivot(m, bg, Ns)
    V0 = ic.AS_TARGET / sp['Pz']
    b = Bench(name or ('Ta%g' % alpha), m, V0, Nslice=Nslice, gam=gam, note='T-model alpha=%g, A_s at N_* = %.3f' % (alpha, Ns))
    b.spec = sp; b.Nstar = Ns
    return b


# ------------------------------------------------------------------------------------------------ engine driver
def engine_args(b, gam, end, nsamp, cbin=0.0, rtol=1e-10, atol=1e-13, binw=0.0, z=None, G=None, mode=4, bcap=None):
    g = lambda v: '%.17g' % float(v)
    G = b.G if G is None else G
    args = [str(BIN), str(mode), g(b.mI if bcap is None else bcap), g(gam), g(P05['xi']), g(b.eta2), '0', g(P05['lx']), g(P05['ly']),
            g(P05['ls']), g(P05['kap']), g(P05['M']), g(P05['h']), g(MU), '0', '0', g(G.real), g(G.imag), '0', g(end), g(rtol),
            g(atol), str(int(nsamp)), g(binw), g(cbin), 'F']
    for c in z:
        args += [g(c.real), g(c.imag)]
    args += ['0'] * 6
    return args


def parse(txt):
    S, C, END, I0, W = [], [], None, None, {}
    for line in txt.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('C '):
            C.append([float(x) for x in line.split()[2:]])
        elif line.startswith('W') or line.startswith('A3040'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('I0 '):
            v = [float(x) for x in line.split()[1:]]; I0 = dict(I0=v[0], phihat0=v[1], u0=[v[2], v[3]], R0=v[4])
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]),
                                        max_hidden_work=float(v[5]), min_mL2=float(v[6]))
    return np.array(S), np.array(C), END, I0, W


def run_engine(b, tag, gam, end, nsamp, cbin=0.0, capT=None, z=None, G=None, rtol=1e-10, save=True):
    build()
    if z is None:
        z = b.prep()['z']
    env = b.env()
    if capT:
        env['RPCH_CAPT'] = '%.17g' % capT
    args = engine_args(b, gam, end, nsamp, cbin=cbin, z=z, G=G, rtol=rtol)
    import time
    t0 = time.time()
    r = subprocess.run(args, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError('%s: %s' % (tag, r.stderr[:500]))
    S, C, END, I0, W = parse(r.stdout)
    out = dict(tag=tag, S=S, C=C, END=END, I0=I0, W=W, seconds=time.time() - t0, gam=gam, end=end)
    if save:
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RUNS / (tag + '.npz'), S=S, C=C)
        json.dump(dict(tag=tag, END=END, I0=I0, W=W, seconds=out['seconds'], gam=gam, end=end, bench=b.summary()),
                  open(RUNS / (tag + '_meta.json'), 'w'), indent=1, default=float)
    return out


def col(S, name):
    return S[:, COLS.index(name)]
