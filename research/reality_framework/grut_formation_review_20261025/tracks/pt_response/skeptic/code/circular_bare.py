"""SKEPTIC (pt_response_skeptic), INDEPENDENT code.  Three checks of pt_response's Bogoliubov section.

(a) Circular condensate.  Parameters rebuilt from the recorded background (hessian_s20_26.8.npz, our own reduction).
    Branches from the EIGENVALUES of the 4x4 rotating-frame first-order matrix (not the closed-form quartic), and the
    conjugate/normal RMS ratio from a LAB-FRAME integration with H_lab(s) = R(w0 s) diag(Hrr,Htt) R(w0 s)^T at several
    kick phases (pt_response used expm in the rotating frame).  Radial-branch antiparticle content, amplitude-squared
    and |omega|-weighted.
(b) Why the bare (condensate-only) closure is unstable.  The measured coherent orbit is held by more than the bare
    force: we measure the radial force deficit Sigma(s) = -(psi_c'' . psi_c)/|psi_c|^2 - lam |psi_c|^2 from the
    recorded orbit.  With Sigma > 0, the bare linearisation violates the Goldstone identity: the phase-rotation mode
    i psi_c obeys L_bare(i psi_c) = -i Sigma psi_c, a tachyonic tangential mass^2 = -Sigma.  Circular estimate of the
    unstable band and rate.
(c) Growth rates by DIRECT long-time integration (no Floquet construction) of d'' + k^2 d + H(s) d = 0 with
      (c1) bare Hessian on the measured orbit (s 20-26.8),
      (c2) bare Hessian on a SELF-CONSISTENT bare orbit psi'' = -lam|psi|^2 psi with the measured r_max, r_min,
      (c3) a consistent 'Hartree' closure: orbit and linearisation both in V = Sigma0 r^2/2 + lam r^4/4, with Sigma0 set so
           that the radial period equals the measured 0.90 at the measured r_max/r_min.
Usage: python3 -I circular_bare.py <hessian npz> <packet dir> <out json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.signal import argrelmax, argrelmin
from scipy.optimize import brentq

hz = np.load(sys.argv[1]); P = sys.argv[2].rstrip('/') + '/'; OUT = sys.argv[3]
LX = 1e-3
S = hz['s']; m = hz['mean']; psi = m[:, 0] + 1j * m[:, 1]
r = np.abs(psi); th = np.unwrap(np.angle(psi))
mx = argrelmax(r)[0]; mn = argrelmin(r)[0]
w_mean = (th[mx[-1]] - th[mx[0]]) / (S[mx[-1]] - S[mx[0]])
Trad = (S[mx[-1]] - S[mx[0]]) / (len(mx) - 1)
g_coh = LX * np.mean(r ** 2)
# rotating-frame orbit-averaged mean Hessian
Hb = hz['Hbar'][:, :2, :2]
c, s_ = np.cos(th), np.sin(th)
Rm = np.array([[c, -s_], [s_, c]]).transpose(2, 0, 1)
Hrot = np.einsum('tji,tjk,tkl->til', Rm, Hb, Rm)  # R^T H R
sel = (S >= S[mx[0]]) & (S <= S[mx[-1]])
Hrot_mean = Hrot[sel].mean(0)
hf = hz['hf'][sel].mean()
res = dict(method=__doc__, orbit=dict(radial_period=Trad, mean_angular_velocity=w_mean, lam_mean_r2_coh=g_coh, hf_mean=hf,
                                      Hrot_mean=Hrot_mean.tolist(), r_max=float(r[mx].mean()), r_min=float(r[mn].mean())))
print('orbit: T_rad %.4f  w_mean %.4f  lam<|psi_c|^2> %.4f  Hrot %s  hf %.3f  r_max %.2f r_min %.2f' % (
    Trad, w_mean, g_coh, np.round(Hrot_mean, 3).tolist(), hf, r[mx].mean(), r[mn].mean()))

models = {'C1': dict(Hrr=Hrot_mean[0, 0], Htt=Hrot_mean[1, 1], w0=np.sqrt(Hrot_mean[1, 1])),
          'C2': dict(Hrr=w_mean ** 2 + 2 * g_coh, Htt=w_mean ** 2, w0=w_mean),
          'C3': dict(Hrr=hf + g_coh, Htt=hf - g_coh, w0=np.sqrt(hf - g_coh))}

def rot_matrix(k2, Hrr, Htt, w0):
    return np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-(k2 + Hrr - w0 ** 2), 0, 0, 2 * w0], [0, -(k2 + Htt - w0 ** 2), -2 * w0, 0]])

def lab_ratio(k2, Hrr, Htt, w0, phase0):
    D = np.diag([Hrr, Htt])
    def H(t):
        a = w0 * t + phase0; R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]); return R @ D @ R.T
    def f(t, y):
        Y = y.reshape(4, 2); q = Y[:2]; p = Y[2:]
        return np.concatenate([p, -(k2 * q + H(t) @ q)]).ravel()
    t = np.arange(151) * .04
    y0 = np.zeros((4, 2)); y0[2:] = np.eye(2)
    sol = solve_ivp(f, (0, 6), y0.ravel(), t_eval=t, method='DOP853', rtol=1e-11, atol=1e-13)
    Rr = sol.y.reshape(4, 2, -1)[:2]  # [a, b, t]
    gN = .5 * (Rr[0, 0] + Rr[1, 1] + 1j * (Rr[1, 0] - Rr[0, 1])); gA = .5 * (Rr[0, 0] - Rr[1, 1] + 1j * (Rr[1, 0] + Rr[0, 1]))
    return float(np.linalg.norm(gA) / np.linalg.norm(gN))

shells = [x for x in json.load(open(P + 'WINDOW_base48_14_24.json'))['shells'] if x['shell'] > 0]
circ = {}
for name, p in models.items():
    rows = []
    for sh in shells:
        k2 = sh['kfd2_mean']
        ev = np.linalg.eigvals(rot_matrix(k2, **p)); Om = np.sort(np.abs(ev.imag))[::-1][[0, 2]]  # radial, phonon
        maxre = float(np.max(ev.real))
        br = []
        for O in Om:
            beta = 2 * p['w0'] * O / (k2 + p['Htt'] - p['w0'] ** 2 - O ** 2)
            pu, pv = (1 - beta) ** 2, (1 + beta) ** 2
            lu, lv = p['w0'] - O, p['w0'] + O
            # amplitude^2 fraction in negative lab lines, and |omega|-weighted fraction
            neg_amp = (pu * (lu < 0) + pv * (lv < 0)) / (pu + pv)
            neg_act = (abs(lu) * pu * (lu < 0) + abs(lv) * pv * (lv < 0)) / (abs(lu) * pu + abs(lv) * pv)
            br.append(dict(Om=float(O), beta=float(beta), v_over_u=float(np.sqrt(pv / pu)), lab=[float(lu), float(lv)],
                           neg_amp_fraction=float(neg_amp), neg_action_fraction=float(neg_act)))
        rr = [lab_ratio(k2, p['Hrr'], p['Htt'], p['w0'], ph) for ph in (0., 1.1, 2.3)] if sh['shell'] in (1, 2, 3, 5, 10, 15) else None
        rows.append(dict(shell=sh['shell'], k2=k2, max_real_eig=maxre, radial=br[0], phonon=br[1], lab_ratio_phases=rr))
    circ[name] = dict(params={k: float(v) for k, v in p.items()}, rows=rows)
    print(name, {k: round(float(v), 3) for k, v in p.items()})
    for rw in rows:
        if rw['lab_ratio_phases'] is not None:
            print('  sh%2d k2 %.2f radial Om %.3f lab %.2f v/u %.4f negfrac(amp) %.4f negfrac(|w|) %.4f | ratio(lab, 3 phases) %s' % (
                rw['shell'], rw['k2'], rw['radial']['Om'], rw['radial']['lab'][0], rw['radial']['v_over_u'], rw['radial']['neg_amp_fraction'],
                rw['radial']['neg_action_fraction'], np.round(rw['lab_ratio_phases'], 4)))
res['circular'] = circ
vu = [rw['radial']['v_over_u'] for nm in circ for rw in circ[nm]['rows']]
na = [rw['radial']['neg_amp_fraction'] for nm in circ for rw in circ[nm]['rows']]
nw = [rw['radial']['neg_action_fraction'] for nm in circ for rw in circ[nm]['rows']]
res['circular_summary'] = dict(v_over_u_range=[min(vu), max(vu)], neg_amp_fraction_range=[min(na), max(na)], neg_action_fraction_range=[min(nw), max(nw)])
print('radial v/u %.4f-%.4f ; antiparticle fraction amp^2 %.4f-%.4f ; |w|-weighted %.4f-%.4f' % (min(vu), max(vu), min(na), max(na), min(nw), max(nw)))

# ---------------- (b) Goldstone violation of the bare closure on the measured orbit
spq = CubicSpline(S, m[:, :2], axis=0)
acc = spq(S, 2)  # psi_c''
Sigma = -(acc[:, 0] * m[:, 0] + acc[:, 1] * m[:, 1]) / r ** 2 - LX * r ** 2
torque = (m[:, 0] * acc[:, 1] - m[:, 1] * acc[:, 0])  # d/ds of coherent charge density
inner = (S > S[0] + .1) & (S < S[-1] - .1)
SigOrb = float(np.mean(Sigma[sel & inner]))
res['bare_goldstone'] = dict(Sigma_orbit_mean=SigOrb, Sigma_min_max=[float(Sigma[inner].min()), float(Sigma[inner].max())],
                             Sigma_at_apo=float(np.mean(Sigma[mx])), Sigma_at_peri=float(np.mean(Sigma[mn])),
                             hf_minus_2lam_r2c_mean=float(np.mean((hz['hf'] - 2 * LX * r ** 2)[sel])))
# circular estimate: bare closure, Hrr = 3 g_bare, Htt = g_bare, orbit rotating at w0 with w0^2 = g_bare + Sigma
gb = g_coh; w0b = np.sqrt(gb + SigOrb)
band = []
for k2 in np.linspace(0, 12, 121):
    ev = np.linalg.eigvals(rot_matrix(k2, 3 * gb, gb, w0b)); band.append((float(k2), float(max(ev.real))))
unst = [b for b in band if b[1] > 1e-8]
res['bare_goldstone']['circular_estimate'] = dict(g_bare=gb, w0=w0b, rate_k0=band[0][1], unstable_k2_max=max(b[0] for b in unst) if unst else 0.,
                                                  rates=band[::10])
print('bare closure: Sigma(orbit mean) %.3f (apo %.3f, peri %.3f); circular estimate: rate(k=0) %.3f, unstable for k2 < %.2f' % (
    SigOrb, res['bare_goldstone']['Sigma_at_apo'], res['bare_goldstone']['Sigma_at_peri'], band[0][1], res['bare_goldstone']['circular_estimate']['unstable_k2_max']))

# ---------------- (c) growth rates by direct integration
def growth(Hfun, s0, s1, k2, nren=40):
    """max Lyapunov rate of d'' + k2 d + H(s) d = 0 over [s0, s1] (2x2 field), with periodic renormalisation."""
    def f(t, y):
        Y = y.reshape(4, 4); q = Y[:2]; p = Y[2:]
        return np.concatenate([p, -(k2 * q + Hfun(t) @ q)]).ravel()
    Y = np.eye(4); edges = np.linspace(s0, s1, nren + 1); logs = 0.; lg = []
    for a, b in zip(edges[:-1], edges[1:]):
        sol = solve_ivp(f, (a, b), Y.ravel(), method='DOP853', rtol=1e-10, atol=1e-12)
        Y = sol.y[:, -1].reshape(4, 4)
        Q_, R_ = np.linalg.qr(Y); logs += np.log(abs(R_[0, 0])); Y = Q_  # leading direction
        lg.append(logs)
    # fit slope over the second half to remove transients
    t = edges[1:]; h = len(t) // 2
    return float(np.polyfit(t[h:], np.array(lg)[h:], 1)[0])

Hc_sp = CubicSpline(S, hz['Hcoh'][:, :2, :2].reshape(len(S), 4), axis=0)
Hc = lambda t: Hc_sp(t).reshape(2, 2)
k2list = [sh['kfd2_mean'] for sh in shells[:9]]
c1 = [growth(Hc, 20.0, 26.7, k2) for k2 in k2list]
print('c1 bare on measured orbit (direct, s 20-26.7):', np.round(c1, 3))

def orbit_rates(Sig, rmax, rmin, nper=40, k2s=k2list, kfine=None):
    # turning points of V_eff = q^2/(2r^2) + Sig r^2/2 + lam r^4/4
    q2 = (Sig / 2 * (rmax ** 2 - rmin ** 2) + LX / 4 * (rmax ** 4 - rmin ** 4)) / (.5 * (1 / rmin ** 2 - 1 / rmax ** 2))
    q = np.sqrt(q2)
    def fo(t, y):
        x1, x2, p1, p2 = y; rr = x1 * x1 + x2 * x2; F = Sig + LX * rr
        return [p1, p2, -F * x1, -F * x2]
    y0 = [rmax, 0., 0., q / rmax]
    # radial period
    ev = lambda t, y: y[0] * y[2] + y[1] * y[3]
    ev.direction = 1.
    sol = solve_ivp(fo, (0, 50), y0, method='DOP853', rtol=1e-12, atol=1e-10, dense_output=True, events=ev, max_step=.01)
    tmin = sol.t_events[0]
    # pr crosses 0 upward at r minima; period = spacing of minima
    Tr = float(np.mean(np.diff(tmin))) if len(tmin) > 1 else np.nan
    Tend = nper * Tr
    sol = solve_ivp(fo, (0, Tend + 1), y0, method='DOP853', rtol=1e-12, atol=1e-10, dense_output=True, max_step=.01)
    def H(t):
        x1, x2 = sol.sol(t)[:2]; v = np.array([x1, x2]); rr = x1 * x1 + x2 * x2
        return (Sig + LX * rr) * np.eye(2) + 2 * LX * np.outer(v, v)
    ks = k2s if kfine is None else kfine
    rates = [growth(H, 0., Tend, k2, nren=4 * nper) for k2 in ks]
    angle = np.unwrap(np.arctan2(sol.sol(np.linspace(0, Tend, 20000))[1], sol.sol(np.linspace(0, Tend, 20000))[0]))
    return dict(Sigma=Sig, q=float(q), radial_period=Tr, mean_angular_velocity=float((angle[-1] - angle[0]) / Tend),
                k2=list(map(float, ks)), rates=rates)

rmax, rmin = float(r[mx].mean()), float(r[mn].mean())
c2 = orbit_rates(0., rmax, rmin)
print('c2 self-consistent bare orbit: T_rad %.3f w %.3f ; rates %s' % (c2['radial_period'], c2['mean_angular_velocity'], np.round(c2['rates'], 4)))
# c3: Sigma0 such that the radial period equals the measured one
f3 = lambda Sg: orbit_rates(Sg, rmax, rmin, nper=2, k2s=[1.])['radial_period'] - Trad
Sg3 = brentq(f3, 0., 30., xtol=1e-4)
c3 = orbit_rates(Sg3, rmax, rmin)
print('c3 consistent Hartree-like orbit Sigma0 %.3f: T_rad %.3f w %.3f ; rates %s' % (Sg3, c3['radial_period'], c3['mean_angular_velocity'], np.round(c3['rates'], 4)))
# fine k2 scan for the self-consistent bare orbit (where are its genuine parametric bands?)
kf = list(np.round(np.linspace(0.05, 12., 49), 4))
c2f = orbit_rates(0., rmax, rmin, nper=30, kfine=kf)
c3f = orbit_rates(Sg3, rmax, rmin, nper=30, kfine=kf)
print('c2 fine: max rate %.4f at k2 %.2f' % (max(c2f['rates']), kf[int(np.argmax(c2f['rates']))]))
print('c3 fine: max rate %.4f at k2 %.2f' % (max(c3f['rates']), kf[int(np.argmax(c3f['rates']))]))
res['growth'] = dict(k2=k2list, c1_bare_measured_orbit=c1, c2_selfconsistent_bare_orbit=c2, c3_consistent_hartree_orbit=c3,
                     c2_fine=c2f, c3_fine=c3f, trackfloquet_bare_rates_20p09=[0.816, 0.703, 0.411, 0.194, 0.251, 0, 0, 0, 0])
json.dump(res, open(OUT, 'w'), indent=1, default=float)
