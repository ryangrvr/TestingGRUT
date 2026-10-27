"""INDEPENDENT code (Claude, pt_response): relativistic Bogoliubov theory of a CIRCULAR condensate and its
consequences for (i) the conjugate/normal response ratio and (ii) particle counting.

Lab field psi = chi1 + i chi2, condensate psi_c = R e^{i w0 s} (e^{+i w s} = positive charge, Codex convention).
Rotating frame eta = e^{-i w0 s} delta psi = a + i b (a radial, b tangential).  With a homogeneous mean-field
(Hartree-Fock-Bogoliubov) Hessian whose rotating-frame components are H_rr, H_tt:
    a'' - 2 w0 b' + (k^2 + H_rr - w0^2) a = 0,      b'' + 2 w0 a' + (k^2 + H_tt - w0^2) b = 0.
For a consistent circular HFB condensate H_tt = w0^2 (Goldstone) and H_rr - H_tt = 2 g, g = lam R^2 (+ pair term).
Modes (a, b) = A (1, i beta) e^{-i Om s}:
    (k^2 + H_rr - w0^2 - Om^2)(k^2 + H_tt - w0^2 - Om^2) = 4 w0^2 Om^2,
    beta = 2 w0 Om / (k^2 + H_tt - w0^2 - Om^2),
lab lines: delta psi_k ~ A(1-beta) e^{i(w0-Om)s} [u],  delta psi_{-k} ~ A*(1+beta) e^{i(w0+Om)s} [v];
rotating-frame energy E = |A|^2 [Om^2(1+beta^2) + (k^2+H_rr-w0^2) + (k^2+H_tt-w0^2) beta^2], action J = E/Om.
Free limit (g -> 0, H_rr = H_tt = w0^2): Om_+ = E_k + w0 (antiparticle, beta=-1), Om_- = E_k - w0 (particle, beta=+1).
Usage: python3 -I bogoliubov_circular.py <hessian npz> <packet dir> <out json>
"""
import sys, json
import numpy as np
from scipy.linalg import expm

hz = np.load(sys.argv[1]); P = sys.argv[2].rstrip('/') + '/'; OUT = sys.argv[3]
LX = .001
s = hz['s']; m = hz['mean']
th = np.unwrap(np.arctan2(m[:, 1], m[:, 0])); r2 = m[:, 0] ** 2 + m[:, 1] ** 2
w_mean = np.polyfit(s, th, 1)[0]
c, sn = np.cos(th), np.sin(th)
Hrot = np.array([np.array([[ci, si], [-si, ci]]) @ H[:2, :2] @ np.array([[ci, -si], [si, ci]]) for ci, si, H in zip(c, sn, hz['Hbar'])]).mean(0)
hf = hz['hf'].mean(); g_coh = LX * r2.mean()

models = {
    # C1: HFB-consistent circular condensate from the orbit-averaged rotating-frame mean Hessian (Goldstone enforced)
    'C1_HFB_rotframe': dict(Hrr=Hrot[0, 0], Htt=Hrot[1, 1], w0=np.sqrt(Hrot[1, 1])),
    # C2: bare condensate amplitude lam<|psi_c|^2>, measured mean angular velocity; Sigma = w0^2 - g (consistency)
    'C2_bare_amp_measured_w0': dict(Hrr=w_mean ** 2 + 2 * g_coh, Htt=w_mean ** 2, w0=w_mean),
    # C3: bare amplitude, isotropic Hartree normal mass = <hf>: Sigma = hf - 2g, w0^2 = g + Sigma
    'C3_bare_amp_isotropic_hartree': dict(Hrr=hf + g_coh, Htt=hf - g_coh, w0=np.sqrt(hf - g_coh)),
}

def branches(k2, Hrr, Htt, w0):
    A_ = k2 + Hrr - w0 ** 2; B_ = k2 + Htt - w0 ** 2
    # Om^4 - Om^2 (A+B+4w0^2) + A B = 0
    bb = A_ + B_ + 4 * w0 ** 2; disc = bb * bb - 4 * A_ * B_
    out = []
    for sign in (+1, -1):
        Om2 = (bb + sign * np.sqrt(disc)) / 2; Om = np.sqrt(max(Om2, 0.))
        beta = 2 * w0 * Om / (B_ - Om2) if abs(B_ - Om2) > 1e-14 else np.inf
        E = Om2 * (1 + beta ** 2) + A_ + B_ * beta ** 2
        J = E / Om if Om > 0 else np.nan
        pu = (1 - beta) ** 2; pv = (1 + beta) ** 2
        lines = [(w0 - Om, pu), (w0 + Om, pv)]
        Aminus = sum(-w * p for w, p in lines if w < 0); Aplus = sum(w * p for w, p in lines if w > 0)
        out.append(dict(Om=Om, beta=beta, v_over_u=float(np.sqrt(pv / pu)) if pu > 0 else np.inf,
                        lab_u=w0 - Om, lab_v=w0 + Om, r_minus=Aminus / J, r_plus=Aplus / J,
                        charge_per_action=(Aplus - Aminus) / J))
    return out  # [radial (+), phonon (-)]

def ratio_rms(k2, Hrr, Htt, w0, T=6., n=151):
    Amat = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-(k2 + Hrr - w0 ** 2), 0, 0, 2 * w0], [0, -(k2 + Htt - w0 ** 2), -2 * w0, 0]])
    t = np.linspace(0, T, n); E1 = expm(Amat * (t[1] - t[0])); Y = np.eye(4); gN = []; gA = []
    for _ in t:
        r = Y[:2, 2:]  # response of (a,b) to unit (a',b') kicks
        gN.append(.5 * (r[0, 0] + r[1, 1] + 1j * (r[1, 0] - r[0, 1]))); gA.append(.5 * (r[0, 0] - r[1, 1] + 1j * (r[1, 0] + r[0, 1])))
        Y = E1 @ Y
    gN = np.array(gN); gA = np.array(gA)
    return float(np.linalg.norm(gA) / np.linalg.norm(gN)), gN, gA

shells = json.load(open(P + 'WINDOW_base48_14_24.json'))['shells']
res = dict(method=__doc__, orbit=dict(mean_angular_velocity=w_mean, mean_r2=float(r2.mean()), g_coh=g_coh, mean_hf=float(hf),
                                       rotating_frame_mean_Hessian=Hrot.tolist()), models={})
for name, p in models.items():
    rows = []
    for sh in shells:
        if sh['shell'] == 0: continue
        k2 = sh['kfd2_mean']
        br = branches(k2, **p); rr, _, _ = ratio_rms(k2, **p)
        rows.append(dict(shell=sh['shell'], k=sh['k_center'], kfd2=k2, ratio_rms_0_6=rr, radial=br[0], phonon=br[1]))
    # Counting inversion with the measured lab actions (window s 14-24, per net charge)
    inv = []; Jp = Jm = Am = 0.
    for sh, row in zip([x for x in shells if x['shell']], rows):
        M = np.array([[row['radial']['r_minus'], row['phonon']['r_minus']], [row['radial']['r_plus'], row['phonon']['r_plus']]])
        y = np.array([sh['negative_frequency_action_per_net'], sh['positive_frequency_action_per_net']])
        J = np.linalg.solve(M, y)
        inv.append(dict(shell=sh['shell'], A_minus_lab=y[0], A_plus_lab=y[1], J_radial=J[0], J_phonon=J[1]))
        Jp += J[0]; Jm += J[1]; Am += y[0]
    res['models'][name] = dict(params={k: float(v) for k, v in p.items()}, rows=rows, counting=dict(per_shell=inv,
                               waves_A_minus_lab=Am, waves_J_radial=Jp, waves_J_phonon=Jm))
    print(name, {k: round(float(v), 3) for k, v in p.items()})
    for r in rows[:15]:
        print(f"  sh{r['shell']:2d} k2={r['kfd2']:.2f} ratio={r['ratio_rms_0_6']:.3f} | rad Om={r['radial']['Om']:.2f} lab {r['radial']['lab_u']:.2f}/{r['radial']['lab_v']:.2f} v/u={r['radial']['v_over_u']:.3f} r-={r['radial']['r_minus']:.3f}"
              f" | ph Om={r['phonon']['Om']:.2f} lab {r['phonon']['lab_u']:.2f}/{r['phonon']['lab_v']:.2f} v/u={r['phonon']['v_over_u']:.3f} r-={r['phonon']['r_minus']:.3f}")
    print('  waves: lab A- =', round(Am, 4), ' J_radial =', round(Jp, 4), ' J_phonon =', round(Jm, 4))
json.dump(res, open(OUT, 'w'), indent=1, default=float)
