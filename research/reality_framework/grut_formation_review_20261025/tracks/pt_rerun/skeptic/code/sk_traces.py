"""Skeptic of pt_rerun: all eight published histories from FIELD_TRACES.npz (own code, Claude).

(1) Coherent (k=0) Hann negative-frequency action regenerated for all 8 histories x 2 windows from the float64 spatial
    means q0,q1 of the published traces (pt_rerun regenerated only IR12_s1_w1).
(2) Averaged-export coherent outer product <v0 v0^T> recomputed from the traces for all 8 histories and both rules;
    coherent share of charge/energy; what squaring the window-mean would keep.
(3) Matching ambiguity of the LATE coherent component (pt_rerun's qualification guessed it from the 1.5% gap of the
    tau=20 orbit). Hann^2-weighted window means of several counts of the zero mode, per original net charge:
      lab   : the Hann negative-frequency action (the packet's statistic);
      Jq    : instantaneous radial action of the zero mode alone in the bare quartic potential (E0 = |p0|^2/2+lam|psi0|^4/4);
      JH    : instantaneous radial action with the fluctuation Hartree mass M^2 = Hartree_m2 - 2 lam |psi0|^2 added
              (U(1)-symmetric Gaussian fluctuations; the anomalous <dpsi^2> psi0* term is NOT central and is omitted);
      n(m)  : equal-time harmonic count ((|p0|^2/m + m|psi0|^2)/2 - Q0)/2 at m^2 = 7, 12, 20, 26 and the instantaneous
              Hartree m^2 (this is also the sudden-switch value at that mass).
    Time variation of Jq/JH inside a window measures how far the zero mode is from an adiabatic central orbit.
Usage: python3 -B sk_traces.py <FIELD_TRACES.npz> <published packet dir> <out.json>
"""
import sys, json, glob, os
import numpy as np
from scipy.integrate import quad

tr_path, pk, out = sys.argv[1:4]
Z = np.load(tr_path); names = [str(x) for x in Z['names']]; col = {n: i for i, n in enumerate(names)}
LX, G2 = 1e-3, 1e-2
a_m, C = 8.869218738806, 1.9857962750634526
tags = ['base48', 'base48_half', 'IR12_s1_w1', 'IR12_s2_w1', 'IR12_s1_w3', 'IR12_fine_s1_w1', 'IR24_s1_w1', 'IR24_s2_w1']


def Jr_u(E, Q, M2, lam):
    """radial action of H = p^2/2 + M2 r^2/2 + lam r^4/4 at energy E, angular momentum Q>0 (u = r^2 form)."""
    c = [-2 * lam, -4 * M2, 8 * E, -4 * Q * Q]
    rt = np.roots(c); rr = np.sort(rt.real[abs(rt.imag) < 1e-6 * (1 + abs(rt.real))])
    pos = rr[rr > 0]
    if len(pos) < 2: return np.nan
    u2, u3 = pos[-2], pos[-1]
    P = lambda u: max(-2 * lam * u ** 3 - 4 * M2 * u * u + 8 * E * u - 4 * Q * Q, 0.)
    # u = u2 + (u3-u2) sin^2 t  removes both square-root endpoints
    g = lambda t: np.sqrt(P(u2 + (u3 - u2) * np.sin(t) ** 2)) / (4 * (u2 + (u3 - u2) * np.sin(t) ** 2)) * 2 * (u3 - u2) * np.sin(t) * np.cos(t)
    return quad(g, 0, np.pi / 2, epsabs=0, epsrel=1e-11, limit=200)[0] / np.pi


res = {}
for tag in tags:
    s = Z[tag + '_s']; V = Z[tag + '_values']; ds = s[1] - s[0]
    q0, q1, p0, p1 = (V[:, col[n]] for n in ('q0', 'q1', 'p0', 'p1'))
    psi0 = q0 + 1j * q1; Q0 = q0 * p1 - q1 * p0; nQ = V[:, col['nQ']]; Qnet = float(nQ[0])
    r2 = q0 ** 2 + q1 ** 2; pp = p0 ** 2 + p1 ** 2; hm2 = V[:, col['Hartree_m2']]
    M2 = hm2 - 2 * LX * r2
    E0 = .5 * pp + LX / 4 * r2 * r2
    EH = E0 + .5 * M2 * r2
    out_t = {}
    for lo, hi in [(14., 24.), (24., 34.)]:
        ix = np.flatnonzero((s >= lo - 1e-9) & (s <= hi + 1e-9)); M = len(ix); h = np.hanning(M); w = 2 * np.pi * np.fft.fftfreq(M, ds)
        X = np.fft.fft(h * psi0[ix]); P = abs(X) ** 2 / (M * np.sum(h * h))
        lab = float(np.sum((-w * P)[w < 0]) / Qnet); labp = float(np.sum((w * P)[w > 0]) / Qnet)
        wq = h * h / np.sum(h * h)
        win = json.load(open(os.path.join(pk, f'WINDOW_{tag}_{lo:g}_{hi:g}.json')))
        # dominant coherent frequencies
        neg = w < 0; posm = w > 0
        wn = float(w[neg][np.argmax(P[neg])]); wp = float(w[posm][np.argmax(P[posm])])
        wmean_neg = float(np.sum((w * P)[neg]) / np.sum(P[neg])); wmean_pos = float(np.sum((w * P)[posm]) / np.sum(P[posm]))
        Jq = np.array([Jr_u(E0[i], Q0[i], 0., LX) for i in ix]); JH = np.array([Jr_u(EH[i], Q0[i], M2[i], LX) for i in ix])
        harm = {}
        for m2 in [7., 12., 20., 26.]:
            m = np.sqrt(m2); n = ((pp[ix] / m + m * r2[ix]) / 2 - Q0[ix]) / 2
            harm[f'm2_{m2:g}'] = float(wq @ n / Qnet)
        mH = np.sqrt(hm2[ix]); nH = ((pp[ix] / mH + mH * r2[ix]) / 2 - Q0[ix]) / 2
        harm['m2_Hartree_instantaneous'] = float(wq @ nH / Qnet)
        for lab_m, mm in (('m_eq_mean_abs_neg_freq', abs(wmean_neg)), ('m_eq_mean_pos_freq', wmean_pos)):
            n = ((pp[ix] / mm + mm * r2[ix]) / 2 - Q0[ix]) / 2; harm[lab_m] = float(wq @ n / Qnet)
        out_t[f'{lo:g}_{hi:g}'] = dict(
            lab_coherent_negative_trace=lab, lab_published=win['coherent_negative_frequency_action'], lab_minus_published=lab - win['coherent_negative_frequency_action'],
            lab_coherent_positive_minus_negative=labp - lab, Q0_over_Qnet_hann=float(wq @ Q0[ix] / Qnet),
            peak_neg_freq=wn, peak_pos_freq=wp, mean_neg_freq=wmean_neg, mean_pos_freq=wmean_pos,
            Jr_quartic_hann=float(wq @ Jq / Qnet), Jr_quartic_range=[float(np.nanmin(Jq) / Qnet), float(np.nanmax(Jq) / Qnet)],
            Jr_Hartree_hann=float(wq @ JH / Qnet), Jr_Hartree_range=[float(np.nanmin(JH) / Qnet), float(np.nanmax(JH) / Qnet)],
            Hartree_m2_hann=float(wq @ hm2[ix]), fluct_M2_hann=float(wq @ M2[ix]), lam_psi0sq_hann=float(wq @ (LX * r2[ix])),
            harmonic_counts=harm)
        # (2) export check with trapezoid weights (as in pt_production)
        trap = np.ones(M); trap[0] = trap[-1] = .5
        v = np.stack([q0, q1, p0, p1, V[:, col['qphi']], V[:, col['pphi']]], axis=1)[ix]
        A = np.load(os.path.join(pk, f'AVERAGE_{tag}.npz'))
        for rule in ('conformal', 'cosmic'):
            ww = trap * (1. if rule == 'conformal' else a_m + C * s[ix]); ww = ww / ww.sum()
            outer = np.einsum('t,ti,tj->ij', ww, v, v); mean = ww @ v; sq = np.outer(mean, mean)
            key = f'{lo:g}_{hi:g}_{rule}'; Mx = A[key + '_mean_outer']; Qa = float(A[key + '_Q']); rho = float(A[key + '_rho'])
            out_t[f'{lo:g}_{hi:g}'][f'export_{rule}'] = dict(
                mean_outer_vs_trace_max_rel=float(np.max(abs(Mx - outer)) / np.max(abs(Mx))),
                coherent_charge_share=float((Mx[0, 3] - Mx[1, 2]) / Qa), coherent_kinetic_energy_share=float(.5 * (Mx[2, 2] + Mx[3, 3] + Mx[5, 5]) / rho),
                coherent_quartic_energy_share=float(LX / 4 * (wq @ (r2[ix] ** 2)) / rho),
                squared_mean_charge_share=float((sq[0, 3] - sq[1, 2]) / Qa), squared_mean_kinetic_share=float(.5 * (sq[2, 2] + sq[3, 3] + sq[5, 5]) / rho))
    res[tag] = out_t
    print(tag, json.dumps({k: {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items() if not kk.startswith('export')} for k, v in out_t.items()}), flush=True)
json.dump(res, open(out, 'w'), indent=1)
