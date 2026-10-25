#!/usr/bin/env python3
"""sp_phaseprep: stochastic (Fokker-Planck / Langevin) distribution of the source phase during inflation (own code).

H_i = 1.1211167304108865e13 GeV, c_H = 1, stated couplings; H = 1 units.  Three treatments:

 A. 1D compact source phase psi2 at fixed (aligned) radii, V = B (1 - cos psi2), canonical phi = f_eff psi2
    (f_eff and B from spectrum.json).  Overdamped Starobinsky FP
        dP/dN = (1/f^2) [ d/dpsi (V' P / 3) + (1/8 pi^2) d^2P/dpsi^2 ],
    finite-volume on a periodic grid, spectrum of the FP operator (relaxation rates) and evolution from
    (i) uniform, (ii) a narrow peak at the aligned minimum, (iii) a narrow peak at the required psi2 = pi/2.
 B. 2D complex Sigma plane (the light phase and the Sigma radial mode, 0.61 H and 1.05 H, are not hierarchical):
    equilibrium P ~ exp(-8 pi^2 V / 3) d^2 sigma_c with X, Y held at the aligned minimum (B1) or re-minimised at each
    Sigma (B2).  Probabilities that one Hubble patch lies in the tolerated band around +-pi/2, and in the
    'supplied box' (phase band AND |Sigma| within 5% of its aligned value, as needed by the release: sp_point).
 C. Second-order Langevin in e-fold time for Sigma (X, Y fixed), noise H/2pi per sqrt(e-fold) on each canonical
    component [E: standard stochastic-with-momentum prescription], started at Codex's supplied state: how fast does
    a preparation that is switched off N e-folds before the end of inflation decay?

All of this is the overdamped / light-field stochastic formalism applied at m ~ 0.6-1 H; corrections O(m^2/3H^2)
are expected [E].  Usage: python3 fokker_planck.py ../results/fokker_planck.json
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm, eig

HERE = Path(__file__).resolve().parent
TR = HERE.parents[1]
sys.path.insert(0, str(TR / 'sp_point' / 'code'))
import spmodel as spm  # noqa: E402  (our sp_point radial minimum and potential)

H = spm.H_STAR
SPEC = json.load(open(HERE.parent / 'results' / 'spectrum.json'))['H_star']
F_EFF = SPEC['light']['f_eff_over_H']
B = SPEC['B_over_H4']
K_SY = 8 * np.pi**2 / 3
WIDTHS = [0.03, 0.09, 0.18]   # tolerated half-widths of the band around +-pi/2 [rad]; 0.09 ~ the coherent band
                              # 1.46-1.64 rad of sp_point; 0.03 Codex's noise level; 0.18 generous


def band_mask(t, w):
    return (np.abs(t - np.pi / 2) < w) | (np.abs(t + np.pi / 2) < w)


def part_A():
    n = 1024
    t = -np.pi + (np.arange(n) + 0.5) * 2 * np.pi / n
    dt = 2 * np.pi / n
    # finite-volume (Scharfetter-Gummel-like symmetric) discretisation of d/dpsi[ D (dP/dpsi + K V' P) ]
    D = 1 / (8 * np.pi**2 * F_EFF**2)
    V = B * (1 - np.cos(t))
    L = np.zeros((n, n))
    for i in range(n):
        for j in ((i + 1) % n, (i - 1) % n):
            # flux i<-j with equilibrium-preserving weights exp(-K (V_i - V_j)/2)
            w = D / dt**2 * np.exp(-K_SY * (V[i] - V[j]) / 2)
            L[i, j] += w
            L[i, i] -= D / dt**2 * np.exp(-K_SY * (V[j] - V[i]) / 2)
    ev = np.sort(np.real(eig(L, right=False)))[::-1]
    peq = np.exp(-K_SY * V); peq /= peq.sum()
    res = dict(f_eff=F_EFF, B=B, kappa=K_SY * B, relaxation_rates_per_efold=(-ev[1:5]).tolist(),
               relaxation_efolds=float(-1 / ev[1]))
    res['equilibrium'] = {f'w{w}': float(peq[band_mask(t, w)].sum()) for w in WIDTHS}
    res['equilibrium']['mean_cos'] = float((peq * np.cos(t)).sum())
    res['equilibrium']['density_ratio_pi2_to_0'] = float(np.exp(-K_SY * B))
    # time evolution
    Ns = [0, 0.5, 1, 2, 4, 8, 16, 30, 60]
    starts = dict(uniform=np.ones(n) / n,
                  aligned=np.exp(-t**2 / (2 * 0.03**2)),
                  required=np.exp(-(t - np.pi / 2)**2 / (2 * 0.03**2)))
    evo = {}
    for name, p0 in starts.items():
        p0 = p0 / p0.sum()
        rows = []
        for N in Ns:
            p = expm(L * N) @ p0
            rows.append(dict(N=N, **{f'P_band_w{w}': float(p[band_mask(t, w)].sum()) for w in WIDTHS},
                             P_within_0p09_of_plus_pi2=float(p[np.abs(t - np.pi / 2) < 0.09].sum()),
                             mean_cos=float((p * np.cos(t)).sum()), mean_sin=float((p * np.sin(t)).sum())))
        evo[name] = rows
    res['evolution'] = evo
    return res


def _Vs(S, x, y):
    """Inflationary potential as a function of complex Sigma with X = x, Y = -y fixed (H units)."""
    M, h, mu = spm.pars(H)
    return (spm.LX * x**4 + spm.LY * y**4 + spm.LS * (np.abs(S)**2 - h * h / 2)**2 + M * M * y * y
            - 2 * spm.KAP * y * x**3 + 2 * mu * y * y * np.real(S) - x * x)


def relaxed_xy(s, cos_mu, x0, y0, iters=40):
    """Vectorised Newton for the X, Y radial balance at given Sigma modulus s and cos_mu (cos_mu=-1 aligned)."""
    M, h, mu = spm.pars(H)
    x = np.full_like(s, x0); y = np.full_like(s, y0)
    for _ in range(iters):
        gx = 4 * spm.LX * x**3 - 6 * spm.KAP * y * x * x - 2 * x
        gy = 4 * spm.LY * y**3 + 2 * M * M * y - 2 * spm.KAP * x**3 + 4 * mu * s * y * cos_mu
        jxx = 12 * spm.LX * x * x - 12 * spm.KAP * y * x - 2
        jxy = -6 * spm.KAP * x * x
        jyy = 12 * spm.LY * y * y + 2 * M * M + 4 * mu * s * cos_mu
        det = jxx * jyy - jxy * jxy
        dx = (-gx * jyy + gy * jxy) / det
        dy = (-gy * jxx + gx * jxy) / det
        x += dx; y += dy
    return x, y


def part_B():
    xa, ya, sa = spm.minimum(H, 1.0)
    M, h, mu = spm.pars(H)
    ns, nt = 1200, 1440
    s = (np.arange(ns) + 0.5) * (1.6 / ns)            # |Sigma| in H units, 0..1.6 (aligned value 0.4309)
    th = -np.pi + (np.arange(nt) + 0.5) * 2 * np.pi / nt  # angle relative to the aligned direction
    S_, T_ = np.meshgrid(s, th, indexing='ij')
    Sig = -S_ * np.exp(1j * T_)                        # aligned direction is Sigma = -|Sigma|
    out = {}
    for label in ('B1_XY_fixed', 'B2_XY_relaxed'):
        if label == 'B1_XY_fixed':
            Vg = _Vs(Sig, xa, ya)
        else:
            x, y = relaxed_xy(S_, -np.cos(T_), xa, ya)
            Vg = (spm.LX * x**4 + spm.LY * y**4 + spm.LS * (S_**2 - h * h / 2)**2 + M * M * y * y
                  - 2 * spm.KAP * y * x**3 - 2 * mu * S_ * y * y * np.cos(T_) - x * x)
        V0 = Vg.min()
        w = S_ * np.exp(-K_SY * (Vg - V0))            # measure d^2 sigma_c ~ s ds dtheta
        Z = w.sum()
        Pth = w.sum(axis=0) / Z
        box_r = np.abs(S_ / sa - 1) < 0.05
        r = dict(V_min_minus_aligned=float(V0 - _Vs(-sa + 0j, xa, ya)) if label == 'B1_XY_fixed' else None,
                 mean_cos=float((Pth * np.cos(th)).sum()), rms_theta=float(np.sqrt((Pth * th**2).sum())),
                 mean_abs_sigma=float((w * S_).sum() / Z), rms_abs_sigma=float(np.sqrt((w * S_**2).sum() / Z)))
        for wd in WIDTHS:
            mb = band_mask(T_, wd)
            pb = float(w[mb].sum() / Z)
            pbox = float(w[mb & box_r].sum() / Z)
            r[f'w{wd}'] = dict(P_phase_band=pb, P_supplied_box=pbox, P_radius_ok_given_band=pbox / pb if pb else None)
        # conditional radius distribution at theta = pi/2
        k = np.argmin(np.abs(th - np.pi / 2))
        col = w[:, k] / w[:, k].sum()
        r['at_theta_pi2'] = dict(mean_abs_sigma=float((col * s).sum()), mode_abs_sigma=float(s[np.argmax(col)]),
                                 P_abs_sigma_within_5pct_of_aligned=float(col[np.abs(s / sa - 1) < 0.05].sum()))
        r['density_ratio_supplied_to_min'] = float(np.exp(-K_SY * (_Vs(-1j * sa, xa, ya) - _Vs(-sa + 0j, xa, ya))))
        out[label] = r
    out['aligned_radii'] = [xa, ya, sa]
    out['dV_supplied'] = float(_Vs(-1j * sa, xa, ya) - _Vs(-sa + 0j, xa, ya))
    return out


def part_C(nreal=20000, dN=1e-3, Nmax=3.0, seed=20261025):
    xa, ya, sa = spm.minimum(H, 1.0)
    M, h, mu = spm.pars(H)
    rng = np.random.default_rng(seed)
    # canonical components sigma_c = sqrt2 (Re Sigma, Im Sigma); supplied state Sigma = -i sa
    u = np.zeros((nreal, 2)); u[:, 1] = -np.sqrt(2) * sa
    v = np.zeros((nreal, 2))
    def grad(u):
        S = (u[:, 0] + 1j * u[:, 1]) / np.sqrt(2)
        dVdS = 2 * spm.LS * (np.abs(S)**2 - h * h / 2) * S + mu * ya * ya   # dV/dS* ; dV/dsig = sqrt2 Re/Im(dV/dS*)
        return np.sqrt(2) * np.stack([dVdS.real, dVdS.imag], axis=1)
    def stats(u, v, N):
        S = (u[:, 0] + 1j * u[:, 1]) / np.sqrt(2)
        thv = np.angle(-S)  # relative to aligned direction -|S|
        rad_ok = np.abs(np.abs(S) / sa - 1) < 0.05
        d = dict(N=N, mean_abs_sigma=float(np.abs(S).mean()), mean_cos=float(np.cos(thv).mean()),
                 rms_speed=float(np.sqrt((v**2).sum(1).mean())))
        for wd in WIDTHS:
            mb = np.abs(np.abs(thv) - np.pi / 2) < wd
            d[f'P_band_w{wd}'] = float(mb.mean()); d[f'P_box_w{wd}'] = float((mb & rad_ok).mean())
        return d
    rows = [stats(u, v, 0.0)]
    marks = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0]
    N = 0.0; k = 0
    nsteps = int(round(Nmax / dN))
    for it in range(1, nsteps + 1):
        # velocity Verlet for the deterministic part plus noise kick on the field (Euler-Maruyama)
        a = -3 * v - grad(u)
        v_half = v + 0.5 * dN * a
        u = u + dN * v_half + rng.normal(size=u.shape) * np.sqrt(dN) / (2 * np.pi)
        a2 = -3 * v_half - grad(u)
        v = v_half + 0.5 * dN * a2
        N = it * dN
        if k < len(marks) and N >= marks[k] - 1e-12:
            rows.append(stats(u, v, round(N, 6))); k += 1
    # deterministic (noise-free) control
    ud = np.array([[0.0, -np.sqrt(2) * sa]]); vd = np.zeros((1, 2)); det = []
    k = 0
    for it in range(1, nsteps + 1):
        a = -3 * vd - grad(ud); vh = vd + 0.5 * dN * a; ud = ud + dN * vh; vd = vh + 0.5 * dN * (-3 * vh - grad(ud))
        N = it * dN
        if k < len(marks) and N >= marks[k] - 1e-12:
            S = (ud[0, 0] + 1j * ud[0, 1]) / np.sqrt(2)
            det.append(dict(N=round(N, 6), abs_sigma=float(abs(S)), theta=float(np.angle(-S)))); k += 1
    return dict(nreal=nreal, dN=dN, noisy=rows, deterministic=det,
                note='Sigma only, X and Y fixed at the aligned minimum; second-order Langevin; estimate.')


def main(out):
    res = dict(H_GeV=H, A_fixed_radius_1D=part_A())
    print('A', {k: v for k, v in res['A_fixed_radius_1D'].items() if k != 'evolution'}, flush=True)
    for name, rows in res['A_fixed_radius_1D']['evolution'].items():
        print('  ', name, [(r['N'], round(r['P_band_w0.09'], 4), round(r['mean_cos'], 3)) for r in rows])
    res['B_sigma_plane_equilibrium'] = part_B()
    print('B', json.dumps(res['B_sigma_plane_equilibrium'], indent=1), flush=True)
    res['C_langevin_from_supplied'] = part_C()
    for r in res['C_langevin_from_supplied']['noisy']:
        print('  C', r)
    print('  det', res['C_langevin_from_supplied']['deterministic'])
    Path(out).write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
