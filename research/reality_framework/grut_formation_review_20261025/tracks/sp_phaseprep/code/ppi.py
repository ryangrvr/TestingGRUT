#!/usr/bin/env python3
"""sp_phaseprep: a minimal supplied phase-preparation interaction (PPI) of Dine-Randall-Thomas type, tested (own code).

NOT part of the stated action. Supplied, tagged, test-only inflationary terms (H = H_i = 1 units):
    V_PPI = c_S H^2 |Sigma|^2 + c_Y H^2 |Y|^2 + [ muP Sigma^* Y^2 + c.c. ]     (muP complex, |muP| ~ a_H H)
i.e. Hubble-induced masses and a Hubble-induced A-term on the same gauge-invariant operator as mu, with a CP phase
relative to mu (cf. Dine, Randall & Thomas 1995/96 [R]).  Z2_phi even, gauge invariant, dimension <= 4.

'Tuned' design for a given c_S: muP and c_Y are fixed so that the inflationary minimum is EXACTLY Codex's supplied
state (aligned radii 23.826064 / 1.893762 / 0.430868 H and source phase delta = arg(mu + muP) = pi/2):
    |mu + muP| = [2 lam_S s (s^2 - h^2/2) + c_S s] / y^2,   c_Y = 2 (|mu + muP| - mu) s.
'Untuned' design: c_Y = 0, |mu + muP| fixed by |Sigma| = 0.430868 H only (Y and X radii then shift).
c_S = 0 gives muP = mu (i - 1), c_Y = 0: the supplied state is exactly the minimum of a pure A-term rotation.

Tests:
 1. minimum (radial elimination + 6D multistart check), 6x6 spectrum, light phase mass and f_eff;
 2. release with instant removal (Codex's assumption) from that minimum;
 3. finite switch-off: g = exp(-(a-1)/Dl) for the PPI alone, and for the PPI and the X Hubble mass together;
    and DRT-like persistence g_A = H/H_i = a^-2, g_m = a^-4 in radiation;
 4. Gaussian ensembles of all six components drawn with the exact Bunch-Davies coarse-grained variance of each
    eigenmode of the PPI minimum (modes with k/aH < x_c at the end of inflation; spectrum.sigma2_coarse), zero
    velocities, instant removal: what net charge fraction survives the PPI's own residual fluctuations.

Usage: python3 ppi.py ../results/ppi.json
"""
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, brentq

HERE = Path(__file__).resolve().parent
TR = HERE.parents[1]
sys.path.insert(0, str(TR / 'sp_point' / 'code'))
sys.path.insert(0, str(HERE))
import spmodel as spm  # noqa: E402  (our sp_point model: potential, aligned minimum, output parser)
from spectrum import sigma2_coarse  # noqa: E402

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/sp_phaseprep')
BIN = SCR / 'pprel'
H = spm.H_STAR
M_, h_, mu_ = spm.pars(H)
CPAR = 0.025
XA, YA, SA = spm.minimum(H, 1.0)


def build():
    src = HERE / 'pprel.c'
    if not BIN.exists() or BIN.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-o', str(BIN), str(src), '-lm'], check=True)


def V_full(z, cS, muP, cY=0.0):
    x, y, s = z
    return float(spm.V(z, H, 1.0) + cS * abs(s)**2 + cY * abs(y)**2 + 2 * np.real(muP * np.conj(s) * y * y))


def radial_min(cS, mutot_abs, cY=0.0):
    """Aligned-phase radial minimum of V_inf + V_PPI with |mu + muP| = mutot_abs (phases at their minimum)."""
    def Vr(r):
        x, y, s = r
        return (spm.LX * x**4 - x * x + spm.LY * y**4 + (M_**2 + cY) * y * y - 2 * spm.KAP * y * x**3
                + spm.LS * (s * s - h_ * h_ / 2)**2 + cS * s * s - 2 * mutot_abs * s * y * y)
    def g(r):
        x, y, s = r
        return np.array([4 * spm.LX * x**3 - 2 * x - 6 * spm.KAP * y * x * x,
                         4 * spm.LY * y**3 + 2 * (M_**2 + cY) * y - 2 * spm.KAP * x**3 - 4 * mutot_abs * s * y,
                         4 * spm.LS * s * (s * s - h_ * h_ / 2) + 2 * cS * s - 2 * mutot_abs * y * y])
    best = None
    for s0 in (SA, 0.5 * SA, 2 * SA, 0.05):
        r = minimize(Vr, [XA, YA, s0], jac=g, method='BFGS', options=dict(gtol=1e-14, maxiter=10000))
        if best is None or r.fun < best.fun:
            best = r
    return best.x, best.fun


def design(cS, tuned=True, delta=np.pi / 2, s_target=SA):
    """muP (and c_Y if tuned) giving source phase delta and the required radii at the PPI minimum."""
    if tuned:
        m_abs = (2 * spm.LS * SA * (SA * SA - h_ * h_ / 2) + cS * SA) / YA**2
        cY = 2 * (m_abs - mu_) * SA
    else:
        cY = 0.0
        f = lambda m: radial_min(cS, m)[0][2] - s_target
        m_abs = brentq(f, 1e-4, 5.0, xtol=1e-15, rtol=1e-14)
    r, v = radial_min(cS, m_abs, cY)
    mutot = m_abs * np.exp(1j * delta)
    muP = mutot - mu_
    x, y, s = r
    z0 = np.array([x, -y, -s * np.exp(1j * delta)], complex)
    return dict(cS=cS, cY=cY, tuned=tuned, mutot_abs=m_abs, muP=muP, z0=z0, radii=r)


def spectrum6(z0, cS, muP, cY=0.0):
    e = 1e-4
    def Vu(u):
        z = z0 + (u[0::2] + 1j * u[1::2]) / np.sqrt(2)
        return V_full(z, cS, muP, cY)
    Hs = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            di = np.zeros(6); dj = np.zeros(6); di[i] = e; dj[j] = e
            Hs[i, j] = Hs[j, i] = (Vu(di + dj) - Vu(di - dj) - Vu(-di + dj) + Vu(-di - dj)) / (4 * e * e)
    # project out the exact gauge zero mode dz_j = i q_j z_j (charges 1, 3, 6)
    gz = 1j * np.array([1, 3, 6]) * z0
    g = np.sqrt(2) * np.ravel([[q.real, q.imag] for q in gz]); g /= np.linalg.norm(g)
    P = np.eye(6) - np.outer(g, g)
    w, U = np.linalg.eigh(P @ Hs @ P)
    gauge = int(np.argmax(np.abs(U.T @ g)))
    w[gauge] = 0.0
    # source phase psi2 = 2 thY - thS with z = R e^{-i th}: d th_j = -Im(dz_j / z_j)
    rows = []
    for k in range(6):
        v = U[:, k]
        dz = (v[0::2] + 1j * v[1::2]) / np.sqrt(2)
        dth = -np.imag(dz / z0)
        rows.append(dict(m2=float(w[k]), dpsi2=float(2 * dth[1] - dth[2]), dpsi1=float(3 * dth[0] - dth[1]),
                         dlnS=float(np.real(dz[2] / z0[2]))))
    phase_like = [r for r in rows if abs(r['dpsi2']) > 0.5 and r['m2'] > 1e-6]
    light = min(phase_like, key=lambda r: r['m2'])
    return dict(mass2=w.tolist(), modes=rows, light_phase=dict(mass_over_H=float(np.sqrt(light['m2'])),
                f_eff=float(1 / abs(light['dpsi2'])), dpsi2=light['dpsi2'], dlnS=light['dlnS'])), U, w


def run(z0, cS=0.0, muP=0j, gmode=0, Dl=1.0, cX=0.0, cY=0.0, a_end=1000.0):
    if not BIN.exists():
        build()
    env = dict(os.environ, PPI_CS=repr(cS), PPI_CY=repr(cY), PPI_CX=repr(cX), PPI_MUR=repr(float(np.real(muP))),
               PPI_MUI=repr(float(np.imag(muP))), PPI_GMODE=str(gmode), PPI_DL=repr(Dl))
    windows = (1, 300, 600, 1000)
    args = [str(BIN), '8', repr(CPAR), repr(M_), repr(h_), repr(mu_), repr(a_end), '4', '0', str(len(windows))]
    args += [repr(float(w)) for w in windows]
    for zz in z0:
        args += [repr(float(zz.real)), repr(float(zz.imag))]
    p = subprocess.run(args, capture_output=True, text=True, env=env)
    if p.returncode:
        raise RuntimeError(p.stderr)
    out = {}
    for line in p.stdout.splitlines():
        f = line.split()
        if f[0] == 'W':
            v = list(map(float, f[1:]))
            d = dict(zip(spm.QN, v[3:3 + len(spm.QN)]))
            out['%g-%g' % (v[0], v[1])] = d
        if f[0] == 'END':
            out['gauss'] = float(f[3])
    w = out['600-1000']
    return dict(eps=w['epsX'], NX=w['NX'], Ncirc=w['NcircX'], JY=w['JY'], JYeig=w['JYeig'],
                YX_over_YQ=w['NX'] * spm.yield_factor(H) / spm.YQ, pairs=1.5 * w['JY'] / abs(w['NX']),
                eps300=out['300-600']['epsX'], gauss=out['gauss'])


def _run_kw(kw):
    return run(**kw)


def ensemble_states(d, U, w, xc, n, seed):
    rng = np.random.default_rng(seed)
    sig = np.array([np.sqrt(sigma2_coarse(m2, xc)) if m2 > 0 else 0.0 for m2 in w])
    zs = []
    for _ in range(n):
        u = U @ (rng.normal(size=6) * sig)
        zs.append(d['z0'] + (u[0::2] + 1j * u[1::2]) / np.sqrt(2))
    return zs, sig


def net(rows):
    nx = np.array([r['NX'] for r in rows]); nc = np.array([r['Ncirc'] for r in rows])
    e = nx.mean() / nc.mean()
    sem = np.std(nx - e * nc, ddof=1) / (np.sqrt(len(nx)) * nc.mean())
    eps = nx / nc
    jy = np.array([r['JY'] for r in rows])
    return dict(N=len(rows), eps_net=float(e), sem=float(sem), frac_abs_ge_019=float(np.mean(np.abs(eps) >= 0.19)),
                frac_same_sign=float(np.mean(np.sign(nx) == np.sign(nx.mean()))),
                pairs_per_net_charge=float(1.5 * jy.mean() / abs(nx.mean())))


def tolerance_scan(pool):
    """Release from Codex's supplied state with one radius rescaled: how precisely must a preparation place radii."""
    base = np.array([XA, -YA, -1j * SA])
    jobs, keys = [], []
    for k, eps_list in ((0, (-1e-3, -3e-4, -1e-4, 1e-4, 3e-4, 1e-3)), (1, (-5e-3, -2e-3, -1e-3, -3e-4, -1e-4, 1e-4, 3e-4, 1e-3, 2e-3, 5e-3)),
                        (2, (-0.05, -0.03, -0.01, -3e-3, 3e-3, 0.01, 0.03, 0.05))):
        for e in eps_list:
            z = base.copy(); z[k] = z[k] * (1 + e)
            jobs.append(dict(z0=z)); keys.append(('XYS'[k], e))
    R = list(pool.map(_run_kw, jobs))
    out = {}
    for (f, e), r in zip(keys, R):
        out.setdefault(f, []).append(dict(rel_shift=e, eps=r['eps'], pairs=r['pairs']))
    return out


def main(out):
    build()
    t0 = time.time()
    res = dict(H_GeV=H, aligned_radii=[XA, YA, SA], mu_over_H=mu_)
    designs = {}
    for tuned in (True, False):
        for cS in ((0.0, 0.5, 1.0, 2.0) if tuned else (0.5, 1.0, 2.0)):
            d = design(cS, tuned)
            sp, U, w = spectrum6(d['z0'], cS, d['muP'], d['cY'])
            d['spectrum'] = sp
            d['U'], d['w'] = U, w
            designs[(tuned, cS)] = d
            print(f"tuned={tuned} cS={cS} cY={d['cY']:.5f} |mutot|={d['mutot_abs']:.6f} muP={d['muP']:.6f} "
                  f"dR/R={(d['radii']/np.array([XA,YA,SA])-1)} light={sp['light_phase']} "
                  f"m={np.sqrt(np.clip(w,0,None)).round(4)}", flush=True)
    d1 = designs[(True, 1.0)]
    rng = np.random.default_rng(3)
    best = np.inf
    for _ in range(16):
        u0 = rng.normal(size=6) * np.repeat(np.abs(d1['z0']) * 1.5 + 0.1, 2)
        r = minimize(lambda u: V_full(u[0::2] + 1j * u[1::2], 1.0, d1['muP'], d1['cY']), u0, method='BFGS',
                     options=dict(gtol=1e-12, maxiter=20000))
        best = min(best, r.fun)
    res['multistart_check_tuned_cS1'] = dict(V_design=V_full(d1['z0'], 1.0, d1['muP'], d1['cY']), V_best=float(best))
    print('multistart', res['multistart_check_tuned_cS1'], flush=True)

    jobs, keys = [], []
    for (tuned, cS), d in designs.items():
        tag = ('tuned' if tuned else 'untuned') + f'_cS{cS}'
        jobs.append(dict(z0=d['z0'])); keys.append((tag, 'instant'))
        if not tuned:
            continue
        kw = dict(z0=d['z0'], cS=cS, muP=d['muP'], cY=d['cY'])
        for Dl in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
            jobs.append(dict(kw, gmode=1, Dl=Dl)); keys.append((tag, f'exp_PPI_only_Dl{Dl}'))
            jobs.append(dict(kw, gmode=1, Dl=Dl, cX=1.0)); keys.append((tag, f'exp_all_Hubble_Dl{Dl}'))
        jobs.append(dict(kw, gmode=2)); keys.append((tag, 'H_tracking_PPI_only'))
        jobs.append(dict(kw, gmode=2, cX=1.0)); keys.append((tag, 'H_tracking_all_Hubble'))
    jobs.append(dict(z0=np.array([XA, -YA, -1j * SA]))); keys.append(('codex_supplied', 'instant'))
    for Dl in (0.003, 0.01, 0.03, 0.1):
        jobs.append(dict(z0=np.array([XA, -YA, -1j * SA]), gmode=1, Dl=Dl, cX=1.0)); keys.append(('codex_supplied', f'exp_XHubble_only_Dl{Dl}'))
    with ProcessPoolExecutor(max_workers=4) as pool:
        R = list(pool.map(_run_kw, jobs))
        rel = {}
        for (tag, k), r in zip(keys, R):
            rel.setdefault(tag, {})[k] = r
        res['releases'] = rel
        for tag in rel:
            print(tag, {k: round(v['eps'], 4) for k, v in rel[tag].items()}, flush=True)
        print('t', time.time() - t0, flush=True)
        res['radius_tolerance'] = tolerance_scan(pool)
        for f, rows in res['radius_tolerance'].items():
            print('tol', f, [(r['rel_shift'], round(r['eps'], 4)) for r in rows], flush=True)
        ens = {}
        for cS, xc, seed in ((0.0, 1 / 40, 11), (0.5, 1 / 40, 12), (1.0, 1 / 40, 13), (2.0, 1 / 40, 14),
                             (0.0, 1 / 500, 15), (1.0, 1 / 500, 16)):
            d = designs[(True, cS)]
            zs, sig = ensemble_states(d, d['U'], d['w'], xc, 128, seed)
            rows = list(pool.map(_run_kw, [dict(z0=z) for z in zs]))
            dpsi = [float(np.angle(np.conj(z[2]) * z[1]**2) - np.angle(np.conj(d['z0'][2]) * d['z0'][1]**2)) for z in zs]
            k = f'tuned_cS{cS}_xc1over{1/xc:.0f}'
            ens[k] = dict(mode_mass_over_H=np.sqrt(np.clip(d['w'], 0, None)).tolist(), mode_sigma_canonical=sig.tolist(),
                          psi2_rms_drawn=float(np.std(dpsi)), **net(rows),
                          members=[dict(eps=r['eps'], NX=r['NX'], dpsi2=p) for r, p in zip(rows, dpsi)])
            print(k, {kk: v for kk, v in ens[k].items() if kk != 'members'}, time.time() - t0, flush=True)
        res['ensembles'] = ens
    for (tuned, cS), d in designs.items():
        res.setdefault('designs', {})[('tuned' if tuned else 'untuned') + f'_cS{cS}'] = dict(
            cS=cS, cY=d['cY'], mutot_abs_over_H=d['mutot_abs'], mutot_over_mu=d['mutot_abs'] / mu_,
            muP_re=float(d['muP'].real), muP_im=float(d['muP'].imag), abs_muP_over_H=float(abs(d['muP'])),
            arg_muP=float(np.angle(d['muP'])), abs_muP_GeV=float(abs(d['muP']) * H), radii=d['radii'].tolist(),
            rel_shift_radii_vs_aligned=(d['radii'] / np.array([XA, YA, SA]) - 1).tolist(),
            spectrum=d['spectrum'], mass_over_H=np.sqrt(np.clip(d['w'], 0, None)).tolist())
    Path(out).write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
