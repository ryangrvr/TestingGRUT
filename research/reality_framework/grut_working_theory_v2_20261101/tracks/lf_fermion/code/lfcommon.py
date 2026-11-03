"""lfcommon: COPY of tracks/e2e_lattice/code/latcommon.py for the lf_fermion track (Claude; issue #13 context).
Changed only: scratch SP, binary and source (lf3d.c), extra run.par lines (extra=...), and the extended diag of lf3d.
Original docstring follows.
e2e_lattice common module (Claude, end-to-end pass; GitHub issue #13 context).

Builds lattice runs of the HF waterfall (lat3d.c): model parameters, the homogeneous attractor up to just before
sigma_c, the local transfer-rate tables, and a Gaussian vacuum realisation of the sub-cutoff modes.

REUSE (read-only, nothing written in other tracks; Python must run with -B):
  * tracks/mf_chain/code/hist_build.py: model(name, H_GeV) (mf_hybrid HFa/HFb potential at a design scale, with
    mf_chain's rescaling rule), V2/grad2/hess2/valley/M2c, channels(P) (the mf_hybrid skeptic's Chan objects: the
    half-period production rates eps(A) of mf_hybrid's axial_tables.json, the skeptic's Pauli saturation F_sat(A) from
    exact multi-period mode runs, the seesaw rest widths rescaled with H), stage1() (homogeneous regression target) and
    stage2() (the skeptic's per-species Pauli fluid, used to continue the lattice end state to K convergence).
    Imported unchanged; MFC_SCRATCH points mf_chain's scratch dirs at this track's scratch.
  * Nothing is copied except the slow-roll start recipe of hist_build.stage1 (re-typed in `homogeneous_start`, with
    the friction off before the slice exactly as there).
Every potential term, closure, cut-off and seed is a model choice; nothing follows from closure.
"""
import json, math, os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRK = TRACK.parent
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/lf_fermion')
os.environ['MFC_SCRATCH'] = str(SP / 'mfc')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
sys.path.insert(0, str(TRK / 'mf_chain' / 'code'))
import hist_build as HB                                   # noqa: E402  (read-only)
from scipy.integrate import solve_ivp                     # noqa: E402
from scipy.optimize import brentq                         # noqa: E402

BIN = SP / 'build' / 'lf3d'
SRC = HERE / 'lf3d.c'
RUNS = SP / 'runs'
RES = TRACK / 'results'
MP_GEV = 2.435e18
H_REF = 1.30e11
CLOSURES = {'hubble': 0, 'track': 1, 'pauli': 2}


def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if (not BIN.exists()) or BIN.stat().st_mtime < SRC.stat().st_mtime:
        import subprocess
        subprocess.run(['gcc', '-O3', '-march=native', '-fopenmp', '-o', str(BIN), str(SRC), '-lm'], check=True)  # into scratch
        BIN.touch()
    return BIN


def model(name, H_GeV=H_REF):
    d, P = HB.model(name, H_GeV)
    P['s_star'] = P['s_star_MP'] * P['mP']
    return d, P


def homogeneous_start(P, t0):
    """Slow-roll attractor 1 e-fold before sigma_c (hist_build.stage1 recipe), friction off, integrated to the slice;
    returns the state at t0 (< 0) relative to the slice, with a(t0) normalised so that a = 1 at the slice."""
    mP = P['mP']

    def Nsr(z):
        g = np.linspace(P['sc'] * (1 + 1e-6), z, 2001)
        f = [HB.valley(P, q)[0] / (mP ** 2 * HB.valley(P, q)[1]) for q in g]
        return float(np.trapezoid(f, g))
    hi = P['sc'] + 0.002 * mP
    s_start = brentq(lambda z: Nsr(z) - 1.0, P['sc'] * (1 + 1e-4), hi, xtol=1e-9 * mP)
    Vst, V1st = HB.valley(P, s_start)
    sd0 = -V1st / math.sqrt(3 * Vst) * mP
    x0 = P['cb'] / HB.M2c(P, s_start)

    def rhs(t, y):
        s, sd, x, xd, la = y
        H = math.sqrt(max(0.5 * (sd * sd + xd * xd) + HB.V2(P, s, x), 1e-300) / (3 * mP ** 2))
        Vs, Vx = HB.grad2(P, s, x)
        return [sd, -3 * H * sd - Vs, xd, -3 * H * xd - Vx, H]
    ev = lambda t, y: y[0] - P['sc']
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0, 50), [s_start, sd0, x0, 0.0, 0.0], method='DOP853', rtol=1e-11, atol=1e-13, events=ev,
                    max_step=0.01, dense_output=True)
    tc = sol.t_events[0][0]
    yc = sol.y_events[0][0]
    y0 = sol.sol(tc + t0) if t0 < 0 else yc
    s, sd, x, xd, la = y0
    a0 = math.exp(la - yc[4])
    H0 = math.sqrt((0.5 * (sd * sd + xd * xd) + HB.V2(P, s, x)) / (3 * mP ** 2))
    return dict(t0=t0, s=s, sd=sd, x=x, xd=xd, a=a0, H=H0, slice=dict(s=yc[0], sd=yc[1], x=yc[2], xd=yc[3]))


def rate_tables(P, path, ntab=700):
    """Per channel: production rate Gamma_j(rho), saturation F_sat,j(rho), band unblocking rate Gamma_band,j(rho)
    and injection energy max(3/4 mu5, M) on a grid in u = ln(rho / v^2) (rho in design H^4 units)."""
    allc = HB.channels(P)
    v = P['v']
    V0 = P['V0']
    u0, u1 = math.log(1e-16 * V0 / v ** 2), math.log(30.0 * V0 / v ** 2)
    us = np.linspace(u0, u1, ntab)
    du = us[1] - us[0]
    lines = ['%d %.17g %.17g' % (ntab, u0, du)]
    meta = []
    for c, w in allc:
        lines.append('%d %.17g %.17g %.17g' % (w, c.Gp, c.Grest, c.M))
        for u in us:
            rho = v * v * math.exp(u)
            mu5 = c.mu5(rho)
            lines.append('%.10g %.10g %.10g %.10g' % (c.rate(rho), c.fsat(rho), c.Grest * c.M / max(c.M, 0.75 * mu5),
                                                       max(0.75 * mu5, c.M)))
        meta.append(dict(w=w, m=c.m, M=c.M, f=c.f, Gp=c.Gp, Grest=c.Grest))
    Path(path).write_text('\n'.join(lines) + '\n')
    return meta


def vacuum_field(N, L, omega2_of_k, a, kcut, rng):
    """Gaussian Wigner realisation of the vacuum for |k| < kcut (comoving): <|phi_k|^2> = 1/(2 w a^3),
    <|phidot_k|^2> = w/(2 a^3), w^2 = k^2/a^2 + M^2.  Returns phi(x), phidot(x) (physical field units)."""
    kk = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing='ij')
    K = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    w2 = np.maximum(K ** 2 / a ** 2 + omega2_of_k, 1e-12)
    w = np.sqrt(w2)
    mask = (K < kcut) & (K > 0)
    out = []
    for amp in (1.0 / (2 * w * a ** 3), w / (2 * a ** 3)):
        W = rng.standard_normal((N, N, N))
        F = np.fft.fftn(W) * np.sqrt(N ** 3 * amp / L ** 3) * mask
        out.append(np.fft.ifftn(F).real)
    return out[0], out[1]


def make_run(tag, name='HFa', H_GeV=H_REF, closure='pauli', Gth=0.0, N=64, L=0.42, dt=None, kcut_mchi=1.2, seed=1,
             t0=-0.02, t_end=3.0, wall=540.0, ndiag=None, snaps=None, gw_snaps=(), pert_free=1, Bglobal=0, loop=1,
             fluct=True, s_fluct=True, tau_av=0.05, gw=1, extra=()):
    build()
    d, P = model(name, H_GeV)
    rd = RUNS / tag
    rd.mkdir(parents=True, exist_ok=True)
    for f in rd.glob('snap_*.bin'):
        f.unlink()
    for f in ('restart.bin', 'done.flag', 'diag.txt'):
        if (rd / f).exists():
            (rd / f).unlink()
    st = homogeneous_start(P, t0)
    v, mP = P['v'], P['mP']
    dx = L / N
    if dt is None:
        dt = dx / 20.0
    if ndiag is None:
        ndiag = max(1, int(round(0.002 / dt)))
    tables = rd / 'tables.txt'
    meta = rate_tables(P, tables)
    rng = np.random.default_rng(seed)
    a = st['a']
    shape = (N, N, N)
    sg = np.full(shape, st['s'] / v); psg = np.full(shape, a ** 3 * st['sd'] / v)
    ch = np.full(shape, st['x'] / v); pch = np.full(shape, a ** 3 * st['xd'] / v)
    sa = np.ones(shape); psa = np.zeros(shape); sb = np.ones(shape); psb = np.zeros(shape)
    kcut = kcut_mchi * P['m_chi']
    if fluct and kcut > 0:
        Vss, Vxx, Vsx = HB.hess2(P, st['s'], st['x'])
        f1, d1 = vacuum_field(N, L, Vss, a, kcut, rng)
        f2, d2 = vacuum_field(N, L, Vxx, a, kcut, rng)
        sg += f1 / v; psg += a ** 3 * d1 / v
        ch += f2 / v; pch += a ** 3 * d2 / v
        if s_fluct:
            for arr, parr in ((sa, psa), (sb, psb)):
                f3, d3 = vacuum_field(N, L, 0.0, a, kcut, rng)
                arr += f3 / P['s_star']; parr += a ** 3 * d3 / P['s_star']
    init = rd / 'init.bin'
    with open(init, 'wb') as fh:
        np.array([a, st['H']], dtype=np.float64).tofile(fh)
        for arr in (sg, psg, ch, pch, sa, psa, sb, psb):
            np.ascontiguousarray(arr, dtype=np.float64).tofile(fh)
    if snaps is None:
        snaps = sorted(set([round(x, 4) for x in list(np.arange(0.05, 1.0001, 0.025)) + list(np.arange(1.1, t_end + 1e-9, 0.1))]))
    par = dict(N=N, L=L, dt=dt, t0=t0, t_end=t_end, wall_limit=wall, ndiag=ndiag,
               Lv2=P['lam'] * v * v, Gv2=P['g'] ** 2 * v * v, sct=P['sc'] / v, cbv=P['cb'] / v, c1v=P['c1'] / v,
               m2c=P['m2'], mu2=P['mu2'], Vc=P['Vc'] / v ** 2, q=(v / mP) ** 2, srel=(P['s_star'] / v) ** 2,
               s0t=P['s0'] / v, vH=v, m0b2=9.0, xib=-0.75, gsa2=(P['g_s'] * v) ** 2, loop=loop,
               fric=CLOSURES[closure], pert_free=pert_free, Bglobal=Bglobal, Gth=Gth,
               own_s=1 if P['smech'] == 'a' else 2, gw=gw, tau_av=tau_av, t_fric_on=0.0, restart=1,
               rundir=str(rd), init=str(init), tables=str(tables))
    lines = ['%s %s' % (k, repr(v_) if isinstance(v_, float) else v_) for k, v_ in par.items()]
    lines += list(extra)
    for ts in snaps:
        lines.append('snap %.6f %d' % (ts, 1 if any(abs(ts - g) < 1e-6 for g in gw_snaps) else 0))
    (rd / 'run.par').write_text('\n'.join(lines) + '\n')
    info = dict(tag=tag, name=name, H_GeV=H_GeV, closure=closure, Gth=Gth, N=N, L=L, dx=dx, dt=dt, kcut_mchi=kcut_mchi,
                seed=seed, t0=t0, t_end=t_end, gw=gw, pert_free=pert_free, Bglobal=Bglobal, fluct=fluct, ndiag=ndiag,
                start=st, channels=meta, gw_snaps=list(gw_snaps), tau_av=tau_av,
                P={k: (float(v_) if isinstance(v_, (int, float, np.floating)) else v_) for k, v_ in P.items()
                   if isinstance(v_, (int, float, str, np.floating))})
    (rd / 'info.json').write_text(json.dumps(info, indent=1, default=float))
    return rd


def run(tag, budget=560.0, threads=2):
    """Run (or resume) until done or until about `budget` seconds of wall time have been used in this call."""
    import subprocess, time
    rd = RUNS / tag
    t_start = time.time()
    while not (rd / 'done.flag').exists():
        left = budget - (time.time() - t_start)
        if left < 30:
            break
        par = (rd / 'run.par').read_text().splitlines()
        par = [l for l in par if not l.startswith('wall_limit')] + ['wall_limit %.1f' % (left - 15)]
        (rd / 'run.par').write_text('\n'.join(par) + '\n')
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        r = subprocess.run([str(BIN), str(rd / 'run.par')], capture_output=True, text=True, env=env)
        print(r.stdout.strip(), r.stderr.strip()[-500:], flush=True)
        if r.returncode not in (0, 3):
            raise RuntimeError('lat3d failed %d' % r.returncode)
    return (rd / 'done.flag').exists()


COLS = ('t a H Hdot eps w rho fsig fchi fgrad fN1 fN2 fN3 fN4 fSM fs R_mean R_rms chi_mean chi_std chiav_mean sig_mean '
        'sig_std fminus fminus_av area_av area_inst sa_mean sa_std sb_mean sb_std ms2a_mean ms2a_std ms2b_mean ms2b_std '
        'B1 B2 B3 B4 cap1 cap2 cap3 cap4 Gam_sig Gam_chi budget hdot_uu fcoh xdot_std chiabs_mean Gam_chi_w rN1 rN2 rN3 '
        'rN4 rSM nN1 nN2 nN3 nN4 Ks Kx Gs Gx Vm').split()


def load_diag(tag):
    rd = RUNS / tag
    D = np.loadtxt(rd / 'diag.txt')
    if D.ndim == 1:
        D = D[None, :]
    # drop duplicate times after restarts
    _, idx = np.unique(D[:, 0], return_index=True)
    D = D[np.sort(idx)]
    return {c: D[:, i] for i, c in enumerate(COLS)}


def load_snap(path):
    with open(path, 'rb') as f:
        N, gw = np.fromfile(f, dtype=np.int32, count=2)
        t, a, L = np.fromfile(f, dtype=np.float64, count=3)
        NT = int(N) ** 3
        names = ['chi', 'chiav', 'sig', 'sa', 'sb', 'R', 'chidot']
        out = dict(N=int(N), t=float(t), a=float(a), L=float(L))
        for nm in names:
            out[nm] = np.fromfile(f, dtype=np.float32, count=NT).reshape(N, N, N)
        if gw:
            out['hdot'] = [np.fromfile(f, dtype=np.float32, count=NT).reshape(N, N, N) for _ in range(6)]
    return out


def dump(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=1, default=float))
