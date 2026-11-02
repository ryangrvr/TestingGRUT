"""Waterfall, defects and fragmentation with coupled (delta sigma, delta chi) modes and three friction models [C on E].

Own integrator (cosmic time, H_ref units, a = 1 at sigma = sigma_c), own potential (sk_common), the track's axial
production tables for the friction rates (validated in sk_axial.py). Differences from the track's waterfall.py:
- linear modes are the coupled 2x2 system (delta sigma, delta chi) with the full tree Hessian (the track evolves delta chi
  alone with M_chi^2 only; the off-diagonal 2 g^2 sigma chi is ~ m_chi^2 during the waterfall);
- two independent vacuum solutions (a quantum in delta sigma or in delta chi at the slice), summed in quadrature;
- the wall test uses the whole time series: ratio(t) = |chi_cl|/delta chi_rms at every step, minimised over
  t in [0.02, 0.6], so the inner turning points of the chi oscillation are included;
- three friction models for the background: 'track' (the track's unblocked half-period rates), 'pauli' (the same rates
  switched off channel by channel once the band saturates, F_sat from sk_axial), 'hubble' (no decay friction).
delta chi_rms counts only the amplified part: sum over k of (|u_k|^2 - |u_k^vac|^2), with |u^vac|^2 = 1/(2 omega_k) at
the slice -- the unamplified vacuum is not a classical displacement (the track adds the full vacuum, which is
conservative for its ratio; both are reported).
Usage: python3 -B sk_waterfall.py HFa -> ../results/sk_waterfall_HFa.json
"""
import math, sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import erfc
import sk_common as C
import sk_cmb as CMB
import sk_pauli as SP

mP = C.mP


def make(name):
    d, P = C.load_bench(name)
    s0, x0 = C.vacuum(P)
    Vc = -C.V2(P, s0, x0, 0.0, loop=False)
    return d, P, s0, x0, Vc


def run(name, friction='track', t_end=0.6, ks=(0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3), rtol=1e-8):
    d, P, s0v, x0v, Vc = make(name)
    cs = [SP.Chan('N1', P['m_sig'], 250., 1200.)]
    cc = [SP.Chan('N2', P['m_chi'], 400., 2000.), SP.Chan('N3', P['m_chi'], 200., 2000.), SP.Chan('N4', P['m_chi'], 120., 2000.)]
    allc = [(c, 0) for c in cs] + [(c, 1) for c in cc]
    V, V1, V2 = CMB.setup(P)
    # start 1 e-fold before sigma_c on the slow-roll attractor (single-field valley, chi at its biased minimum)
    lo, hi = P['sc'] * (1 + 1e-4), P['sc'] + 0.002 * mP
    from scipy.optimize import brentq
    s_start = brentq(lambda z: CMB.N_sr(P, V, V1, z) - 1.0, lo, hi, xtol=1e-6)
    sd = -V1(s_start) / math.sqrt(3 * V(s_start)) * mP
    x_start = P['cb'] / C.M2(P, s_start)
    nk = len(ks)
    kk = np.array(ks) * P['m_chi']
    nC = len(allc)

    def energies(s, sdt, x, xdt):
        Es = 0.5 * sdt * sdt + 0.5 * P['g'] ** 2 * x * x * (s - s0v) ** 2
        Et = 0.5 * (sdt * sdt + xdt * xdt) + C.V2(P, s, x, Vc)
        return Es, max(Et - Es, 0.0)

    def bg(t, y, on):
        s, sdt, x, xdt, la, ur = y[:6]
        Eex = y[6:6 + nC]
        a = math.exp(la)
        rr = max(ur, 0) / a ** 4
        Vv = C.V2(P, s, x, Vc)
        Kin = 0.5 * (sdt * sdt + xdt * xdt)
        H = math.sqrt((Kin + Vv + rr) / (3 * mP ** 2))
        Vs, Vx = C.grad2(P, s, x)
        Gs = Gc = 0.0
        dE = [0.0] * nC
        if on and friction != 'hubble':
            Es, Ec = energies(s, sdt, x, xdt)
            for j, (c, w) in enumerate(allc):
                E = Es if w == 0 else Ec
                G = c.rate(E)
                if friction == 'pauli':
                    cap = c.fsat(E) * E
                    G *= max(0.0, 1.0 - max(Eex[j], 0) / a ** 4 / cap) if cap > 0 else 0.0
                if w == 0:
                    Gs += G
                    dE[j] = G * sdt * sdt * a ** 4
                else:
                    Gc += G
                    dE[j] = G * xdt * xdt * a ** 4
        return [sdt, -(3 * H + Gs) * sdt - Vs, xdt, -(3 * H + Gc) * xdt - Vx, H,
                (Gs * sdt * sdt + Gc * xdt * xdt) * a ** 4] + dE, H

    def rhs(t, y, on, modes):
        out, H = bg(t, y, on)
        if not modes:
            return out + [0.0] * (len(y) - len(out))
        s, x, la = y[0], y[2], y[4]
        a = math.exp(la)
        Vss, Vxx, Vsx = C.hess2(P, s, x)
        m = y[6 + nC:]
        dm = []
        for i in range(nk):
            k2 = (kk[i] / a) ** 2
            for sol in range(2):          # two vacuum solutions
                b = (i * 2 + sol) * 8
                sr, si, xr, xi, dsr, dsi, dxr, dxi = m[b:b + 8]
                dm += [dsr, dsi, dxr, dxi,
                       -3 * H * dsr - (k2 + Vss) * sr - Vsx * xr, -3 * H * dsi - (k2 + Vss) * si - Vsx * xi,
                       -3 * H * dxr - (k2 + Vxx) * xr - Vsx * sr, -3 * H * dxi - (k2 + Vxx) * xi - Vsx * si]
        return out + dm

    y0 = [s_start, sd, x_start, 0.0, 0.0, 0.0] + [0.0] * nC + [0.0] * (16 * nk)
    ev = lambda t, y, *a: y[0] - P['sc']
    ev.terminal, ev.direction = True, -1
    s1 = solve_ivp(rhs, (0, 50), y0, args=(False, False), method='DOP853', rtol=1e-10, atol=1e-12, events=ev, max_step=0.01)
    tc = s1.t_events[0][0]
    yc = list(s1.y_events[0][0])
    lac = yc[4]
    yc[4] = 0.0                            # a = 1 at the slice
    # vacuum initial data at the slice: omega_k^2 = k^2 + diag Hessian (floor 1e-6 m_chi^2 for the k = 0 mode)
    Vss, Vxx, Vsx = C.hess2(P, yc[0], yc[2])
    mm = []
    for i in range(nk):
        for sol in range(2):
            m2 = Vss if sol == 0 else Vxx
            w = math.sqrt(max(kk[i] ** 2 + m2, (0.05 * P['m_chi']) ** 2))
            amp = 1 / math.sqrt(2 * w)
            if sol == 0:
                mm += [amp, 0.0, 0.0, 0.0, 0.0, -w * amp, 0.0, 0.0]
            else:
                mm += [0.0, 0.0, amp, 0.0, 0.0, 0.0, 0.0, -w * amp]
    yc = yc[:6 + nC] + mm
    tt = np.linspace(0, t_end, 6001)
    s2 = solve_ivp(rhs, (0, t_end), yc, args=(True, True), method='DOP853', rtol=rtol, atol=1e-9, t_eval=tt,
                   max_step=0.02 / P['m_chi'])
    Y = s2.y
    chi = Y[2]
    a = np.exp(Y[4])
    # delta chi_rms: amplified part, and the track's 'full vacuum' version
    w0 = []
    for i in range(nk):
        w0.append([math.sqrt(max(kk[i] ** 2 + Vss, (0.05 * P['m_chi']) ** 2)), math.sqrt(max(kk[i] ** 2 + Vxx, (0.05 * P['m_chi']) ** 2))])
    P_k = np.zeros((nk, len(tt)))
    for i in range(nk):
        tot = np.zeros(len(tt))
        for sol in range(2):
            b = 6 + nC + (i * 2 + sol) * 8
            tot += Y[b + 2] ** 2 + Y[b + 3] ** 2
        P_k[i] = tot
    vac_chi = np.array([1 / (2 * w0[i][1]) for i in range(nk)])[:, None]
    amp_k = np.maximum(P_k - vac_chi, 0.0)
    kphys = kk[:, None] / a[None, :]
    # integrate d^3k/(2pi)^3 over the comoving grid (k = 0 point stands for the IR bin [0, k_1/2])
    kc = kk.copy()
    kc[0] = 0.5 * kk[1]
    edges = np.concatenate([[0], 0.5 * (kc[1:] + kc[:-1]), [kc[-1] + 0.5 * (kc[-1] - kc[-2])]])
    shell = (edges[1:] ** 3 - edges[:-1] ** 3) / (6 * math.pi ** 2)
    dchi2_amp = np.sum(shell[:, None] * amp_k, axis=0) / a ** 3
    dchi2_full = np.sum(shell[:, None] * P_k, axis=0) / a ** 3
    ratio_amp = np.abs(chi) / np.sqrt(np.maximum(dchi2_amp, 1e-30))
    ratio_full = np.abs(chi) / np.sqrt(dchi2_full)
    m = (tt >= 0.02)
    j_amp = int(np.argmin(np.where(m, ratio_amp, np.inf)))
    j_full = int(np.argmin(np.where(m, ratio_full, np.inf)))
    # inner turning points of chi after the first fall (t > 0.2)
    tp = [(float(tt[j]), float(chi[j] / P['v'])) for j in range(1, len(tt) - 1)
          if tt[j] > 0.2 and chi[j] < chi[j - 1] and chi[j] <= chi[j + 1]]
    H = np.array([bg(t, Y[:, j], True)[1] for j, t in enumerate(tt)])
    eps = -np.gradient(H, tt) / H ** 2
    j1 = int(np.argmax(eps >= 1.0))
    rad = np.maximum(Y[5], 0) / a ** 4
    tot = rad + np.array([0.5 * (Y[1, j] ** 2 + Y[3, j] ** 2) + C.V2(P, Y[0, j], Y[2, j], Vc) for j in range(len(tt))])
    frad = rad / tot
    res = dict(name=name, friction=friction, t_c_from_start=float(tc), t_eps1=float(tt[j1]), H_eps1=float(H[j1]),
               R_min=float(np.min(6 * H ** 2 * (2 - eps))), frad_t03=float(np.interp(0.3, tt, frad)),
               frad_t06=float(frad[-1]),
               ratio_amp_min=float(ratio_amp[j_amp]), t_ratio_amp_min=float(tt[j_amp]),
               ratio_full_min=float(ratio_full[j_full]), t_ratio_full_min=float(tt[j_full]),
               dchi_rms_over_v_t06=float(math.sqrt(dchi2_amp[-1]) / P['v']),
               max_dchi_rms_over_v=float(np.sqrt(dchi2_amp.max()) / P['v']),
               chi_inner_turning_points_over_v=tp[:12],
               P_wrong_side_at_min=float(0.5 * erfc(ratio_amp[j_amp] / math.sqrt(2))),
               growth_k=[float(np.sqrt(P_k[i].max() / P_k[i][0])) for i in range(nk)], ks_over_mchi=list(ks))
    # first times at which the amplified spread reaches 0.1 v and 0.3 v, and at which ratio(t) < 5 and < 2
    dr = np.sqrt(dchi2_amp) / P['v']
    def first(mask):
        return float(tt[int(np.argmax(mask))]) if mask.any() else None
    res.update(t_dchi_0p1v=first(dr >= 0.1), t_dchi_0p3v=first(dr >= 0.3), t_ratio_lt5=first(m & (ratio_amp < 5)),
               t_ratio_lt2=first(m & (ratio_amp < 2)),
               series=dict(t=[float(x) for x in tt[::200]], chi_over_v=[float(x) for x in (chi / P['v'])[::200]],
                           dchi_rms_over_v=[float(x) for x in dr[::200]], frad=[float(x) for x in frad[::200]]))
    return res


if __name__ == '__main__':
    out = {}
    for nm in (sys.argv[1:] or ['HFa']):
        for fr in ('track', 'pauli', 'hubble'):
            r = run(nm, fr)
            out['%s_%s' % (nm, fr)] = r
            print(json.dumps({k: v for k, v in r.items() if k not in ('chi_inner_turning_points_over_v', 'series')}), flush=True)
            print('   inner turning points chi/v:', r['chi_inner_turning_points_over_v'][:8], flush=True)
    C.dump('sk_waterfall_%s.json' % '_'.join(sys.argv[1:] or ['HFa']), dict(label='C on E (skeptic)', note=__doc__, **out))
