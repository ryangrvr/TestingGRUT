"""P7: what closure forces in the GRUT dark-matter history (own code).

(a) Preparation theorem and its assumptions.
  a1 [C]: three-mode neutral toy of Codex e160356 (charges +1, +3, -6; cutoffs (24,12,9)), own implementation in the
      neutral basis. Reproduce the torque -1.3576 at delta = +pi/2. Then: the K-symmetric cat (psi_+ + psi_-), and the
      record-correlated whole (psi_+ |0> + psi_- |1>)/sqrt2, which is invariant under K (x) X_R: global torque 0,
      branch-conditional torques +-1.3576. K-even data (energy, occupations, entropies) of the branches agree.
  a2 [D, C]: spontaneous CP violation from real couplings. V(A,B) = -a cos A - b cos B + c cos(A+B) (the two
      gauge-invariant angles of kappa Y*X^3, mu Sigma*Y^2 and a real Sigma*Y X^3 term of 'frustrated' sign). The
      CP-conserving point A = B = 0 is a minimum iff c < ab/(a+b); above it there are two degenerate minima +-(A*,B*).
      At every minimum the gradient (the torque) vanishes.
  a3 [C]: rephasing count of physical CP phases among {kappa, mu0, mu1} (rank of the charge matrix).
  a4 [E]: number of release-Hubble patches in today's Hubble volume (random-sign residual ~ N^{-1/2}).
(c) [C]: exact Fourier decomposition of the homogeneous quartic orbit; reproduce the tau = 20 antiparticle content
    1.7387 (en_particles) with own code; the content as a function of E/E_min(Q); what closure-type rules select.
(d) [C]: internal-clock representatives of the inflationary -10 H^2 X mass, on an oscillating-inflaton matter era.
Usage: python3 -B p7_dm.py ../results/p7_dm.json
"""
import sys
import itertools
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import gammaln
from common import dump


# ---------------- a1: neutral three-mode toy ----------------
def build_toy(cut=(24, 12, 9)):
    Na, Nb, Nc = cut
    basis = [(na, nb, nc) for na in range(Na + 1) for nb in range(Nb + 1) for nc in range(Nc + 1)
             if na + 3 * nb - 6 * nc == 0]
    idx = {s: i for i, s in enumerate(basis)}
    d = len(basis)
    H = np.zeros((d, d))
    for i, (na, nb, nc) in enumerate(basis):
        H[i, i] = na + 1.7 * nb + 2.3 * nc
        # b^dag a^3 : (na, nb, nc) -> (na-3, nb+1, nc)
        t = (na - 3, nb + 1, nc)
        if na >= 3 and t in idx:
            amp = np.sqrt(na * (na - 1) * (na - 2) * (nb + 1))
            H[idx[t], i] += 0.025 * amp
            H[i, idx[t]] += 0.025 * amp
        # c b^2 : (na, nb, nc) -> (na, nb-2, nc-1)
        t = (na, nb - 2, nc - 1)
        if nb >= 2 and nc >= 1 and t in idx:
            amp = np.sqrt(nb * (nb - 1) * nc)
            H[idx[t], i] += 0.04 * amp
            H[i, idx[t]] += 0.04 * amp
    Ld = np.array([na + 3 * nb for (na, nb, nc) in basis], float)
    return basis, H, Ld


def proj_coherent(basis, delta):
    al = (np.sqrt(6.0), np.sqrt(2.0), np.sqrt(2.0) * np.exp(1j * delta))
    v = np.array([np.prod([al[k]**n[k] * np.exp(-0.5 * gammaln(n[k] + 1)) for k in range(3)]) for n in basis], complex)
    return v / np.linalg.norm(v)


def toy_observables(psi, basis, H, Ld):
    T_op = 1j * (H * Ld[None, :] - Ld[:, None] * H)       # i[H, L], L diagonal
    E = np.real(np.vdot(psi, H @ psi))
    T = np.real(np.vdot(psi, T_op @ psi))
    p = np.abs(psi)**2
    occ = [float(np.sum(p * np.array([n[k] for n in basis]))) for k in range(3)]
    ent = []
    for k in range(3):
        vals = {}
        for pi, n in zip(p, basis):
            vals[n[k]] = vals.get(n[k], 0) + pi
        q = np.array(list(vals.values()))
        q = q[q > 1e-300]
        ent.append(float(-np.sum(q * np.log(q))))
    return dict(E=E, torque=T, occ=occ, mode_shannon_entropies=ent)


def part_a1():
    basis, H, Ld = build_toy()
    pp = proj_coherent(basis, np.pi / 2)
    pm = proj_coherent(basis, -np.pi / 2)
    K_check = float(np.abs(np.vdot(pm, pp.conj())))      # K psi_+ = psi_- (up to phase)
    ov = np.vdot(pp, pm)
    cat = (pp + pm)
    cat /= np.linalg.norm(cat)
    T_op = 1j * (H * Ld[None, :] - Ld[:, None] * H)
    # record-correlated whole: Psi = (psi_+ |0> + psi_- |1>)/sqrt2 ; hidden torque expectation and conditionals
    Tpp = np.real(np.vdot(pp, T_op @ pp))
    Tmm = np.real(np.vdot(pm, T_op @ pm))
    rec_global = 0.5 * (Tpp + Tmm)
    # reduced hidden state of the record-correlated whole
    rho_h = 0.5 * (np.outer(pp, pp.conj()) + np.outer(pm, pm.conj()))
    w = np.linalg.eigvalsh(rho_h)
    w = w[w > 1e-15]
    S_h = float(-np.sum(w * np.log(w)))
    # invariance of Psi under K (x) X_R: K maps psi_+ -> psi_-, X_R swaps the record
    return dict(neutral_dimension=len(basis), K_maps_plus_to_minus_overlap=K_check,
                overlap_plus_minus=[ov.real, ov.imag],
                psi_plus=toy_observables(pp, basis, H, Ld), psi_minus=toy_observables(pm, basis, H, Ld),
                cat=toy_observables(cat, basis, H, Ld),
                record_whole=dict(global_hidden_torque=rec_global, torque_given_record0=Tpp, torque_given_record1=Tmm,
                                  hidden_entropy_nats=S_h),
                exact_formula_torque=-12 * 0.04 * 2 * np.sqrt(2))


def part_a1b():
    """[C] Scope of the finite K-theorem: K is antiunitary (time-reversal type), so for a K-invariant state at t = 0
    the torque is odd in t and the light charge L (K-even) is even in t; neither vanishes at t != 0."""
    from scipy.linalg import eigh
    basis, H, Ld = build_toy()
    T_op = 1j * (H * Ld[None, :] - Ld[:, None] * H)
    w, Vv = eigh(H)
    out = {}
    for name, dl in (('delta=0', 0.0), ('delta=pi (energy minimum, torque-free at t=0)', np.pi)):
        psi = proj_coherent(basis, dl)
        rows = []
        for t in (-1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0):
            pt = Vv @ (np.exp(-1j * w * t) * (Vv.T @ psi))
            rows.append(dict(t=t, L=float(np.real(np.vdot(pt, Ld * pt))), torque=float(np.real(np.vdot(pt, T_op @ pt)))))
        out[name] = dict(max_imag_of_state=float(np.abs(psi.imag).max()), history=rows)
    return out


# ---------------- a2: frustrated phase potential ----------------
def _V2(x, a, b, c):
    A, B = x
    return -a * np.cos(A) - b * np.cos(B) + c * np.cos(A + B)


def _g2(x, a, b, c):
    A, B = x
    return np.array([a * np.sin(A) - c * np.sin(A + B), b * np.sin(B) - c * np.sin(A + B)])


def global_min2(a, b, c, ng=181):
    g = np.linspace(-np.pi, np.pi, ng)
    AA, BB = np.meshgrid(g, g, indexing='ij')
    VV = _V2((AA, BB), a, b, c)
    i, j = np.unravel_index(np.argmin(VV), VV.shape)
    best = minimize(lambda x: _V2(x, a, b, c), [AA[i, j], BB[i, j]], jac=lambda x: _g2(x, a, b, c),
                    method='BFGS', options=dict(gtol=1e-13))
    x = (best.x + np.pi) % (2 * np.pi) - np.pi
    return x, float(best.fun)


def canted(a, b, c):
    """[D] reciprocal triangle inequalities: the frustrated three-cosine potential has a CP-violating (canted)
    global minimum iff 1/a, 1/b, 1/c satisfy all three triangle inequalities."""
    ia, ib, ic = 1 / a, 1 / b, 1 / c
    return (ia < ib + ic) and (ib < ia + ic) and (ic < ia + ib)


def is_cp_violating(x, tol=1e-5):
    return any(min(abs(v), abs(abs(v) - np.pi)) > tol for v in x)


def part_a2(seed=3):
    rows = []
    cases = [(1, 1, 0.4), (1, 1, 0.5), (1, 1, 0.6), (1, 1, 1.0), (0.1, 1, 13), (1, 1, 13), (10, 1, 13), (100, 1, 13)]
    for a, b, c in cases:
        x, v = global_min2(a, b, c)
        xm = minimize(lambda y: _V2(y, a, b, c), -x, jac=lambda y: _g2(y, a, b, c), method='BFGS',
                      options=dict(gtol=1e-13)).x
        rows.append(dict(a=a, b=b, c=c, criterion_canted=canted(a, b, c), minimum=x.tolist(), V_min=v,
                         V_mirror=float(_V2(xm, a, b, c)), grad_norm_at_min=float(np.linalg.norm(_g2(x, a, b, c))),
                         CP_violating_minimum=is_cp_violating(x)))
    # validation of the criterion on random log-uniform triples (excluding a 2% band around the boundary)
    rng = np.random.default_rng(seed)
    agree = n = 0
    for _ in range(600):
        a, b, c = 10 ** rng.uniform(-1.5, 1.5, 3)
        ia, ib, ic = 1 / a, 1 / b, 1 / c
        margin = min(ib + ic - ia, ia + ic - ib, ia + ib - ic)
        if abs(margin) < 0.02 * max(ia, ib, ic):
            continue
        x, _ = global_min2(a, b, c)
        agree += (is_cp_violating(x) == canted(a, b, c))
        n += 1
    # GRUT-like case: c/b = 13 (gp_phase natural ratio of Sigma*Y X^3 to the mu term); window in a/b
    ratios = np.linspace(0.8, 1.2, 4001)
    win = ratios[[canted(r, 1.0, 13.0) for r in ratios]]
    return dict(cases=rows, validation=dict(n=n, agree=agree),
                window_c_over_b_13=dict(a_over_b_min=float(win.min()), a_over_b_max=float(win.max()),
                                        analytic=[13 / 14, 13 / 12]))


# ---------------- a3: physical phases ----------------
def part_a3():
    # rephasing vectors (theta_X, theta_Y, theta_Sigma) of the coefficient phases
    ops = {'kappa Y*X^3': (3, -1, 0), 'mu0 Sigma*Y^2': (0, 2, -1), 'mu1 phi Sigma*Y^2': (0, 2, -1)}
    out = {}
    for name, keys in (('Z2_phi exact (kappa, mu0)', ['kappa Y*X^3', 'mu0 Sigma*Y^2']),
                       ('Z2_phi broken (kappa, mu0, mu1)', ['kappa Y*X^3', 'mu0 Sigma*Y^2', 'mu1 phi Sigma*Y^2'])):
        Mx = np.array([ops[k] for k in keys], float)
        r = np.linalg.matrix_rank(Mx)
        out[name] = dict(complex_couplings=len(keys), rank=int(r), physical_phases=len(keys) - int(r))
    return out


# ---------------- a4: patch count ----------------
def part_a4():
    T0 = 2.7255 * 8.617333e-14          # GeV
    H0 = 67.4 * 1e3 / 3.0857e22 * 6.582119569e-25   # GeV (h = 0.674)
    gs0, gs = 3.938, 110.0
    out = {}
    for TRH in (2e14, 5e14):
        Hr = 1.3335e11                   # release at ~ H_i (heavy spot point), GeV
        a_ratio = (gs0 / gs)**(1 / 3) * T0 / TRH
        k_ratio = a_ratio * Hr / H0
        N = k_ratio**3
        out[f'T_RH={TRH:.0e}'] = dict(aH_release_over_aH_today=k_ratio, patches=N, residual_fraction=N**-0.5)
    return out


# ---------------- c: orbit antiparticle content ----------------
LX = 1e-3


def orbit_content(y0, lam=LX, nmax=400):
    """Exact line decomposition of psi(t) for H = |pi|^2/2 + (lam/4)|chi|^4 (2D central force)."""
    def rhs(t, y):
        r2 = y[0]**2 + y[1]**2
        return [y[2], y[3], -lam * r2 * y[0], -lam * r2 * y[1]]

    def rdot(t, y):
        return y[0] * y[2] + y[1] * y[3]
    rdot.direction = 1
    # find two successive r-minima (rdot crosses 0 upward)
    sol = solve_ivp(rhs, (0, 50.0), y0, events=rdot, rtol=1e-13, atol=1e-11, method='DOP853', dense_output=True)
    te = sol.t_events[0]
    t1, t2 = te[0], te[1]
    Tr = t2 - t1
    th1 = np.arctan2(sol.sol(t1)[1], sol.sol(t1)[0])
    # unwrap angle along the period
    ts = np.linspace(t1, t2, 2**14 + 1)
    Y = sol.sol(ts)
    psi = Y[0] + 1j * Y[1]
    theta = np.unwrap(np.angle(psi))
    dtheta = theta[-1] - theta[0]
    Om = dtheta / Tr
    wr = 2 * np.pi / Tr
    f = psi * np.exp(-1j * Om * (ts - t1))
    # Fourier coefficients of the periodic f (drop the duplicated endpoint)
    fc = np.fft.fft(f[:-1]) / (len(ts) - 1)
    n = np.fft.fftfreq(len(ts) - 1, d=1.0 / (len(ts) - 1))
    om_n = Om + n * wr
    P = np.abs(fc)**2
    Q = float(np.sum(om_n * P))
    npl = float(np.sum((om_n * P)[om_n > 0]))
    nmi = float(-np.sum((om_n * P)[om_n < 0]))
    Qd = y0[0] * y0[3] - y0[1] * y0[2]
    E = 0.5 * (y0[2]**2 + y0[3]**2) + lam / 4 * (y0[0]**2 + y0[1]**2)**2
    Emin = 0.75 * lam**(1 / 3) * abs(Qd)**(4 / 3)
    r = np.abs(psi)
    return dict(Q_lines=Q, Q_direct=Qd, n_plus=npl, n_minus=nmi, nbar_per_Q=nmi / Qd, E=E, E_over_Emin=E / Emin,
                r_max_over_r_min=float(r.max() / r.min()), radial_period=Tr, Omega=Om, omega_r=wr)


def part_c():
    y_tau20 = [-57.31774205507421, 63.02971554304654, -351.9822243388189, 245.42938010389813]
    tau20 = orbit_content(np.array(y_tau20))
    tau20['quoted_en_particles_n_minus_per_Q0'] = 1.7377883530390552
    # scan in E/Emin at fixed Q: start at apocentre-like point r0 on x axis with tangential speed set by Q = r0 v
    Q = 1000.0
    rows = []
    rc = (Q / np.sqrt(LX))**(1 / 3)
    for s in np.concatenate([np.linspace(1.0, 3.0, 9), np.linspace(3.5, 12, 18)]):
        r0 = rc * s                       # apocentre radius (s >= 1 means beyond the circular radius)
        y0 = np.array([r0, 0.0, 0.0, Q / r0])
        if abs(s - 1.0) < 1e-12:
            Emin = 0.75 * LX**(1 / 3) * Q**(4 / 3)
            rows.append(dict(s=s, E_over_Emin=1.0, nbar_per_Q=0.0, r_max_over_r_min=1.0))
            continue
        oc = orbit_content(y0)
        rows.append(dict(s=float(s), E_over_Emin=oc['E_over_Emin'], nbar_per_Q=oc['nbar_per_Q'],
                         r_max_over_r_min=oc['r_max_over_r_min'], Q_lines_over_Q=oc['Q_lines'] / oc['Q_direct']))
    # harmonic check: psi = A e^{-iwt} + B e^{iwt}: nbar/Q = B^2/(A^2 - B^2); eccentricity 8.06
    e = 8.06
    BA = (e - 1) / (e + 1)
    harm = BA**2 / (1 - BA**2)
    return dict(tau20=tau20, scan=rows, harmonic_ecc_8p06=harm,
                rules=dict(energy_minimum_at_fixed_Q=0.0, stationary_density=0.0,
                           max_entropy_equilibrated_quoted='0.10-0.12 (REVIEW_TIMING / REVIEW_EARLY)',
                           purity='any value: a coherent state |alpha>_X |beta>_Xbar is pure for every alpha, beta'))


# ---------------- d: internal-clock representatives of the -10 H^2 mass ----------------
def part_d(eps_sr=0.01):
    # inflaton chi with V = m^2 chi^2 / 2, units m = 1, M_P(reduced) = 1; start at end of inflation chi ~ sqrt(2)
    def rhs(t, y):
        chi, dchi = y
        rho = 0.5 * dchi**2 + 0.5 * chi**2
        H = np.sqrt(rho / 3)
        return [dchi, -3 * H * dchi - chi]
    sol = solve_ivp(rhs, (0, 200), [1.0, 0.0], rtol=1e-11, atol=1e-13, dense_output=True, max_step=0.01)
    t = np.linspace(100, 200, 200001)
    chi, dchi = sol.sol(t)
    K = 0.5 * dchi**2
    V = 0.5 * chi**2
    rho = K + V
    H2 = rho / 3
    R = 4 * V - 2 * K
    reps = {
        'xi R, xi = -5/6': -5 / 6 * R / H2,
        'c rho/M^2, c = -10/3': -10 / 3 * rho / H2,
        'c V/M^2, c = -10/3': -10 / 3 * V / H2,
        f'c K/M^2, c = -10/eps (eps = {eps_sr})': -10 / eps_sr * K / H2,
    }
    rad = {'xi R, xi = -5/6': 0.0, 'c rho/M^2, c = -10/3': -10.0, 'c V/M^2, c = -10/3': 0.0,
           f'c K/M^2, c = -10/eps (eps = {eps_sr})': 0.0}
    out = {}
    for k, v in reps.items():
        out[k] = dict(de_Sitter=-10.0, matter_era_mean=float(np.mean(v)), matter_era_min=float(v.min()),
                      matter_era_max=float(v.max()), radiation_era=rad[k])
    return out


def part_d2():
    """[E] THM 2.14 athermality budget versus the relative entropy of a charge asymmetry.
    Second order in the chemical potential: D(rho_mu || rho_0) ~ (beta mu)^2 Var(Q)/2 = n_Q^2/(2 chi T), with
    chi T ~ (g/3) T^3 for relativistic bosons, i.e. D/S ~ (eta^2/2) (s/T^3)^2 / (g/3) per entropy unit."""
    rho_c_h2 = 1.05375e-5    # GeV cm^-3
    s0 = 2891.2              # cm^-3
    rho_dm_over_s = 0.120 * rho_c_h2 / s0
    m_X = 0.092              # GeV, Candidate 1 target
    eta = rho_dm_over_s / m_X
    gstar = 110.0
    s_over_T3 = 2 * np.pi**2 / 45 * gstar
    g_hidden = 2.0
    D_over_S = 0.5 * eta**2 * s_over_T3**2 / (g_hidden / 3) / s_over_T3
    return dict(rho_dm_over_s_GeV=rho_dm_over_s, eta_X_if_all_DM_is_mX=eta, D_per_entropy_unit=D_over_S,
                initial_athermality_per_entropy_unit_order='O(1) (pure condensate vs thermal state of equal energy)',
                orders_of_margin=float(-np.log10(D_over_S)))


def main(out):
    res = dict(a1=part_a1(), a1b=part_a1b(), a2=part_a2(), a3=part_a3(), a4=part_a4(), c=part_c(), d=part_d(), d2=part_d2())
    dump(res, out)
    import json
    print(json.dumps(res, indent=1, default=float)[:9000])


if __name__ == '__main__':
    main(sys.argv[1])
