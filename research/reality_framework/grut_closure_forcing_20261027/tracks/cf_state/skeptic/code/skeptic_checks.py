"""Skeptic checks of track cf_state (key cf_state_skeptic, 27 Oct 2026). Own code, independent of ../../code.

Usage: python3 -B skeptic_checks.py ../results/skeptic_checks.json

Parts (labels as in SKEPTIC.md):
  A [D, C]  collision-model ("closed unitary whole") check: the track's 8.9e-16 compares the collision channel with a
            map that is algebraically the same channel; here the closed dilation is compared with the EXACT Lindblad
            solution as the step shrinks (matched and unmatched theta, Lie and Strang splitting, one full dilation).
  B [C]     CSL decoherence rates by a different method (position space, ball distance distribution), at the GRW
            point, the Majorana 2022 bound and the XENONnT 2026 bound (PRL 136, 120201; arXiv:2506.05507).
  C [C]     energy budget of white-noise CSL heating over the age of the Universe vs an internal bath (CMB per baryon).
  D [C]     Diosi-Penrose E_G and tau for the proposed 1e-14 kg sphere; granular part at the XENONnT R0 bound.
  E [D, C]  GPP / Milburn total-H dephasing against the Bothwell et al. two-region ellipse observable (random
            second-pulse phase): algebra of the double commutator + Monte Carlo of common-mode vs independent jitter.
  F [D, C]  P6c: internal representatives of BOTH the never-switched and the instantly switched -10 H^2 mass
            (H is monotonic when Hdot < 0, so any switching time is a value of H).
  G [C]     spontaneous-CP criterion by brute-force grid (no BFGS) incl. near-boundary triples; the criterion at the
            two model points actually in the program (heavy spot point; 7737487 gauge point of gp_phase).
  H [C]     antiparticle budget of the tau = 20 orbit by a different formulation (polar ODE with conserved Q,
            Gauss-Legendre Fourier quadrature); definition dependence check.
  I [C]     release-patch count with own constants.
  J [C]     GPP exponents recomputed (Sr, molecule, proposed sphere).
  K [C]     rephasing count including a renormalizable inflaton coupling chi Sigma* Y^2.
  L [E]     missed postulate: unweighted no-boundary weight exp(8 pi^2 M^2/H^2) at the program's inflation scales.
"""
import sys
import json
import time
import numpy as np
from scipy import constants as sc
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.optimize import minimize

hbar, c, G, kB = sc.hbar, sc.c, sc.G, sc.k
u = sc.physical_constants['atomic mass constant'][0]
eV = sc.e
t_P = np.sqrt(hbar * G / c**5)
m0 = u


def conv(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, complex):
        return [o.real, o.imag]
    raise TypeError(type(o))


# ------------------------------------------------------------------ A
sx = np.array([[0, 1], [1, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)


def lindblad_exact(rho0, H, gamma, t):
    Id = np.eye(2)
    Ls = -1j * (np.kron(Id, H) - np.kron(H.T, Id)) + gamma * (np.kron(sz.T, sz) - np.kron(Id, Id))
    return (expm(Ls * t) @ rho0.reshape(-1, order='F')).reshape(2, 2, order='F')


def dephase(rho, f):
    return np.array([[rho[0, 0], f * rho[0, 1]], [f * rho[1, 0], rho[1, 1]]])


def collision_channel_run(rho0, H, gamma, tf, K, matched, strang):
    tau = tf / K
    theta = 0.5 * np.arccos(np.exp(-2 * gamma * tau)) if matched else np.sqrt(gamma * tau)
    f = np.cos(2 * theta)          # exact single-collision coherence factor (derived in part A_sym)
    if strang:
        Uh = expm(-1j * H * tau / 2)
        rho = rho0.copy()
        for _ in range(K):
            rho = Uh @ rho @ Uh.conj().T
            rho = dephase(rho, f)
            rho = Uh @ rho @ Uh.conj().T
    else:
        Uh = expm(-1j * H * tau)
        rho = rho0.copy()
        for _ in range(K):
            rho = Uh @ rho @ Uh.conj().T
            rho = dephase(rho, f)
    return rho


def full_dilation(psi0, H, gamma, tf, K, matched=True, strang=True):
    """Closed unitary whole: system qubit + K fresh ancillas, each used once. Returns reduced system state and the
    fidelity of the global inverse (applied as one exterior operation)."""
    tau = tf / K
    theta = 0.5 * np.arccos(np.exp(-2 * gamma * tau)) if matched else np.sqrt(gamma * tau)
    Uc = expm(-1j * theta * np.kron(sz, sx)).reshape(2, 2, 2, 2)
    Uh = expm(-1j * H * (tau / 2 if strang else tau))
    psi = np.zeros([2] * (K + 1), complex)
    idx0 = (0,) * (K + 1)
    idx1 = (1,) + (0,) * K
    psi[idx0] = psi0[0]
    psi[idx1] = psi0[1]

    def step(psi, k):
        psi = np.tensordot(Uh, psi, axes=([1], [0]))
        psi = np.tensordot(Uc, psi, axes=([2, 3], [0, k + 1]))
        psi = np.moveaxis(psi, [0, 1], [0, k + 1])
        if strang:
            psi = np.tensordot(Uh, psi, axes=([1], [0]))
        return psi
    for k in range(K):
        psi = step(psi, k)
    ax = list(range(1, K + 1))
    rho = np.tensordot(psi, psi.conj(), axes=(ax, ax))
    return rho


def part_A():
    import sympy as sp
    # symbolic single collision: system rho (x) |0><0|, U = exp(-i th sz (x) sx); partial trace over the ancilla
    th = sp.symbols('theta', real=True)
    r00, r11 = sp.symbols('r00 r11', real=True)
    r01 = sp.symbols('r01')
    Sz = sp.Matrix([[1, 0], [0, -1]])
    Sx = sp.Matrix([[0, 1], [1, 0]])
    M = sp.kronecker_product(Sz, Sx)
    U = sp.cos(th) * sp.eye(4) - sp.I * sp.sin(th) * M     # (sz (x) sx)^2 = 1
    rho = sp.Matrix([[r00, r01], [sp.conjugate(r01), r11]])
    anc = sp.Matrix([[1, 0], [0, 0]])
    big = U * sp.kronecker_product(rho, anc) * U.H
    red = sp.zeros(2, 2)
    for i in range(2):
        for j in range(2):
            red[i, j] = sp.simplify(big[2 * i + 0, 2 * j + 0] + big[2 * i + 1, 2 * j + 1])
    sym = dict(reduced_01=str(sp.simplify(sp.trigsimp(red[0, 1]))), reduced_00=str(sp.simplify(red[0, 0])))

    Delta, gamma, tf = 1.0, 0.5, 3.0
    H = 0.5 * Delta * sx
    out = dict(symbolic_single_collision=sym)
    psiP = np.array([1, 1], complex) / np.sqrt(2)
    out['track_initial_state_plus'] = dict(commutator_norm_H_rho0=float(np.abs(H @ np.outer(psiP, psiP) - np.outer(psiP, psiP) @ H).max()),
                                           exact_lindblad=lindblad_exact(np.outer(psiP, psiP.conj()), H, gamma, tf),
                                           pure_dephasing_prediction_x=np.exp(-2 * gamma * tf))
    for label, psi0 in (('plus_state_track', psiP), ('site_localised_0', np.array([1, 0], complex))):
        out[label] = _part_A_state(psi0, H, gamma, tf)
    out['qsd_vs_lindblad_site_localised'] = _qsd_check(np.array([1, 0], complex), H, gamma, tf)
    out['note'] = ('The track initial state |+> is an eigenstate of H = sigma_x/2, so H never acts on the ensemble and '
                   'the averaged dynamics is pure dephasing x(t) = exp(-2 gamma t). dev_vs_same_channel ~ 1e-16 only '
                   'checks the tensor contraction; dev_vs_exact_lindblad is the physical comparison.')
    return out


def _qsd_check(psi0, H, gamma, tf, ntraj=8000, dt=5e-4, seed=23):
    rng = np.random.default_rng(seed)
    lz = np.array([1.0, -1.0])
    A_ = np.tile(psi0, (ntraj, 1)).astype(complex)
    B_ = A_.copy()
    U0 = expm(-1j * H * dt)
    for _ in range(int(round(tf / dt))):
        dW = rng.normal(0, np.sqrt(dt), ntraj)
        # QSD (Ito, Euler-Maruyama, renormalised)
        Lm = np.sum(np.abs(A_)**2 * lz, axis=1)
        La = A_ * lz - Lm[:, None] * A_
        A_ = A_ + (-1j * (A_ @ H.T) * dt - 0.5 * gamma * (La * lz - Lm[:, None] * La) * dt + np.sqrt(gamma) * La * dW[:, None])
        A_ /= np.linalg.norm(A_, axis=1)[:, None]
        # random unitary (Stratonovich-exact phase kick, then free step)
        B_ = (B_ * np.exp(-1j * np.sqrt(gamma) * dW[:, None] * lz)) @ U0.T
    rq = np.einsum('ti,tj->ij', A_, A_.conj()) / ntraj
    rr = np.einsum('ti,tj->ij', B_, B_.conj()) / ntraj
    rl = lindblad_exact(np.outer(psi0, psi0.conj()), H, gamma, tf)
    return dict(ntraj=ntraj, dt=dt, dev_qsd=float(np.abs(rq - rl).max()), dev_random_unitary=float(np.abs(rr - rl).max()),
                sampling_scale=1 / np.sqrt(ntraj), rho_lindblad=rl)


def _part_A_state(psi0, H, gamma, tf):
    rho0 = np.outer(psi0, psi0.conj())
    rhoL = lindblad_exact(rho0, H, gamma, tf)
    rows = []
    for K in (6, 12, 24, 48, 96, 192, 384, 768, 1536):
        row = dict(K=K, tau=tf / K)
        for matched in (False, True):
            for strang in (False, True):
                r = collision_channel_run(rho0, H, gamma, tf, K, matched, strang)
                row[f"dev_{'matched' if matched else 'track_theta'}_{'strang' if strang else 'lie'}"] = float(np.abs(r - rhoL).max())
        rows.append(row)
    # fitted convergence orders
    Ks = np.array([r['K'] for r in rows], float)
    orders = {}
    for key in rows[0]:
        if key.startswith('dev_'):
            y = np.array([r[key] for r in rows])
            sl = np.polyfit(np.log(1 / Ks[3:]), np.log(y[3:]), 1)[0]
            orders[key] = float(sl)
    # one genuine closed-whole dilation with matched theta + Strang at K = 16 and K = 12 (track's K)
    dil = {}
    for K in (12, 16):
        rd = full_dilation(psi0, H, gamma, tf, K, matched=True, strang=True)
        rc = collision_channel_run(rho0, H, gamma, tf, K, True, True)
        dil[f'K={K}'] = dict(dev_vs_exact_lindblad=float(np.abs(rd - rhoL).max()),
                             dev_vs_same_channel=float(np.abs(rd - rc).max()))
    rd12 = full_dilation(psi0, H, gamma, tf, 12, matched=False, strang=False)
    dil['K=12 track construction (unmatched theta, Lie)'] = dict(dev_vs_exact_lindblad=float(np.abs(rd12 - rhoL).max()))
    return dict(convergence=rows, fitted_order_in_tau=orders, dilations=dil)


# ------------------------------------------------------------------ B
def ball_distance_pdf(s, R):
    """pdf of the distance between two independent uniform points in a ball of radius R (0 <= s <= 2R)."""
    x = s / R
    return np.where((s >= 0) & (s <= 2 * R), 3 * x**2 / R * (1 - 0.75 * x + x**3 / 16), 0.0)


def csl_rate_ball(lam, m, R, rC, d=None):
    """Position-space CSL decoherence rate (mass-proportional, white noise):
    Gamma = (lam/m0^2) int int rho(x) rho(y) [exp(-|x-y|^2/4rC^2) - exp(-|x-y-d|^2/4rC^2)];
    for d >> R + rC the cross term vanishes and Gamma = lam (m/m0)^2 <exp(-s^2/4rC^2)>_ball."""
    N = m / m0
    norm = quad(lambda s: ball_distance_pdf(s, R), 0, 2 * R, limit=200)[0]
    I = quad(lambda s: ball_distance_pdf(s, R) * np.exp(-s**2 / (4 * rC**2)), 0, 2 * R, limit=400,
             points=[min(2 * R, k * rC) for k in (1, 3, 10, 30)])[0]
    return lam * N**2 * I, norm


def part_B():
    rC = 1e-7
    lam = {'GRW': 1e-16, 'Adler': 4e-8, 'Majorana2022_nuclei_95': 0.494 * rC**2,
           'XENONnT2026_90': 3.0e-3 * rC**2, 'XENONnT2026_95': 3.7e-3 * rC**2}
    out = dict(lambda_at_rC_1em7=lam,
               XENONnT_quote='lambda/rC^2 < 3.0e-3 s^-1 m^-2 (90% CL), 3.7e-3 (95% CL); GRW (1e-16, 1e-7 m) excluded at '
                             '9.1 sigma; DP R0 > 4.9e-10 m (90%), 4.5e-10 m (95%); white (Markovian) noise; '
                             'PRL 136, 120201 (2026), arXiv:2506.05507 v2 [Q]',
               GRW_over_XENONnT95=1e-16 / lam['XENONnT2026_95'])
    # proposed sphere 1e-14 kg silica
    m = 1e-14
    R = (3 * m / (4 * np.pi * 2200.0))**(1 / 3)
    sph = dict(R=R)
    for k, l in lam.items():
        G_, norm = csl_rate_ball(l, m, R, rC)
        sph[k] = dict(rate=G_, time=1 / G_)
    sph['pdf_norm'] = norm
    sph['track_GRW_rate_kspace'] = 29724701.836339206
    sph['ratio_own_to_track_GRW'] = sph['GRW']['rate'] / sph['track_GRW_rate_kspace']
    out['sphere_1e-14kg'] = sph
    # Fein molecule: point-like (R ~ 1.5 nm << rC), separation 266 nm (grating period, [R] in track), flight 2 m
    mm = 25000 * u
    v = 2 * np.pi * hbar / (mm * 53e-15)
    T = 2.0 / v
    d = 266e-9
    mol = dict(v=v, T=T)
    for k, l in lam.items():
        mol[k] = l * (mm / m0)**2 * (1 - np.exp(-d**2 / (4 * rC**2))) * T
    mol['track_GRW_exponent'] = 3.4426131055520373e-10
    out['molecule_exponent'] = mol
    return out


# ------------------------------------------------------------------ C
def part_C():
    rC = 1e-7
    t_U = 13.8e9 * 365.25 * 86400
    res = {}
    for name, lam in (('GRW', 1e-16), ('XENONnT95', 3.7e-17), ('Majorana95', 0.494e-14)):
        P = 0.75 * lam * hbar**2 / (m0 * rC**2)          # W per nucleon (m = m0)
        res[name] = dict(W_per_nucleon=P, eV_per_nucleon_over_age=P * t_U / eV)
    # CMB energy per baryon today
    T0 = 2.7255
    rho_g = (np.pi**2 / 15) * (kB * T0)**4 / (hbar * c)**3       # J m^-3
    n_g = 2 * 1.2020569 / np.pi**2 * (kB * T0 / (hbar * c))**3    # m^-3
    eta_b = 6.1e-10                                              # Planck 2018 / BBN [R]
    res['CMB_eV_per_baryon'] = rho_g / (eta_b * n_g) / eV
    res['CMB_over_GRW_deposit'] = res['CMB_eV_per_baryon'] / res['GRW']['eV_per_nucleon_over_age']
    res['reading'] = ('white CSL at GRW deposits ~1e-7 eV per nucleon in 13.8 Gyr; any internal bath with more than '
                      'that capacity per nucleon reproduces it at every finite time, so "unbounded heating" has no '
                      'finite-time observable consequence.')
    return res


# ------------------------------------------------------------------ D
def part_D():
    m = 1e-14
    R = (3 * m / (4 * np.pi * 2200.0))**(1 / 3)
    d = 250e-6
    EG = G * m**2 * (6 / (5 * R) - 1 / d)
    mn = 20 * u
    out = dict(R=R, EG_J=EG, tau_DP_s=hbar / EG)
    for R0 in (2.54e-10, 4.5e-10):
        out[f'granular_EG_R0={R0:.2e}'] = (m / mn) * 1.2 * G * mn**2 / R0
    return out


# ------------------------------------------------------------------ E
def part_E(seed=11):
    rng = np.random.default_rng(seed)
    # (1) algebra: total-H dephasing L(rho) = -s[H,[H,rho]] leaves any O with [O,H]=0 invariant; for
    #     O = s+^(1) s-^(2), H = (w1 sz1 + w2 sz2)/2: Tr(O L(rho)) = -s (w1-w2)^2 Tr(O rho).
    w1, w2 = 1.3, 1.3 - 0.07
    Sp = np.array([[0, 1], [0, 0]], complex)
    I2 = np.eye(2)
    Hm = 0.5 * (w1 * np.kron(sz, I2) + w2 * np.kron(I2, sz))
    O = np.kron(Sp, Sp.conj().T)
    psi = np.kron(np.array([1, 1]) / np.sqrt(2), np.array([1, 1j]) / np.sqrt(2)).astype(complex)
    rho = np.outer(psi, psi.conj())
    L = -(Hm @ (Hm @ rho - rho @ Hm) - (Hm @ rho - rho @ Hm) @ Hm)
    algebra = dict(ratio=float(np.real(np.trace(O @ L) / np.trace(O @ rho))), expected=-(w1 - w2)**2,
                   single_atom_ratio=float(np.real(np.trace(np.kron(Sp, I2) @ L) / np.trace(np.kron(Sp, I2) @ rho))),
                   single_atom_expected=-(w1**2))
    # (2) Monte Carlo of the Bothwell two-region ellipse measurement (random second-pulse phase).
    #     GPP/Milburn map = average over a COMMON time shift dt of the whole sample (one clock), D = w^2 sd^2 / 2.
    Nat, shots = 2000, 4000
    D_single = 2.0                       # single-atom coherence factor exp(-2) = 0.135
    dphi = 0.6                           # differential phase between regions (magnetic gradient in the paper)

    def run(common):
        P1 = np.empty(shots)
        P2 = np.empty(shots)
        sd = np.sqrt(2 * D_single)       # phase jitter (w * dt) with variance 2 D
        for s in range(shots):
            theta = rng.uniform(0, 2 * np.pi)
            if common:
                j = rng.normal(0, sd)
                ph1 = np.full(Nat, j)
                ph2 = np.full(Nat, j + dphi)
            else:
                ph1 = rng.normal(0, sd, Nat)
                ph2 = rng.normal(0, sd, Nat) + dphi
            p1 = 0.5 * (1 + np.cos(ph1 - theta))
            p2 = 0.5 * (1 + np.cos(ph2 - theta))
            P1[s] = rng.binomial(1, p1).mean()
            P2[s] = rng.binomial(1, p2).mean()
        out = {}
        for nm, P in (('region1', P1), ('region2', P2)):
            noise = np.mean(P * (1 - P)) / Nat
            out[nm + '_contrast'] = float(np.sqrt(max(8 * (np.var(P) - noise), 0)))
        # relative phase from the ellipse (correlation)
        cc = np.cov(P1 - P1.mean(), P2 - P2.mean())
        out['cos_dphi_est'] = float(cc[0, 1] / np.sqrt(cc[0, 0] * cc[1, 1]))
        return out
    mc = dict(D_single_atom=D_single, expected_single_atom_coherence=np.exp(-D_single), true_cos_dphi=np.cos(dphi),
              common_mode_jitter_GPP_Milburn=run(True), independent_per_atom_jitter=run(False))
    # (3) what GPP predicts for the measured relative observable: D scales with (w1 - w2)^2.
    om = 2 * np.pi * 429228004229873.0
    T = 37.0
    D_track = 1.5 * t_P**(4 / 3) * T**(2 / 3) * om**2
    frac = 1e-15                          # generous bound on |w1 - w2| / w across the sample [E]
    pred = dict(D_track_single_atom_vs_clock=D_track, assumed_max_fractional_differential=frac,
                D_relative_observable_upper=D_track * frac**2,
                Milburn_tau_bound_track_s=2 * (T / 52.8) / (om**2 * T),
                reading='The 37 s contrast is measured on the relative phase of two sub-ensembles with a randomised '
                        'laser phase; a total-H (one-clock) dephasing is common-mode and leaves it unchanged, so this '
                        'measurement does not bound GPP or Milburn decoherence.')
    return dict(algebra=algebra, monte_carlo=mc, prediction=pred)


# ------------------------------------------------------------------ F
def part_F():
    # m^2 chi^2/2 inflaton, m = 1, reduced M_P = 1, slow-roll start chi0 = 15
    def rhs(t, y):
        chi, dchi = y
        rho = 0.5 * dchi**2 + 0.5 * chi**2
        return [dchi, -3 * np.sqrt(rho / 3) * dchi - chi]
    chi0 = 15.0
    y0 = [chi0, -np.sqrt(2.0 / 3.0)]       # slow-roll attractor dchi = -sqrt(2/3) m M_P
    sol = solve_ivp(rhs, (0, 300), y0, rtol=1e-11, atol=1e-13, dense_output=True, max_step=0.02)
    t = sol.t
    chi, dchi = sol.y
    K = 0.5 * dchi**2
    rho = K + 0.5 * chi**2
    H2 = rho / 3
    epsH = 3 * K / rho
    iend = int(np.argmax((epsH >= 1) & (t > 1)))
    H_end2 = H2[iend]
    # representative 1: never switched, c rho / M^2 = -10 H^2 identically (Friedmann)
    rep_never = -10.0 / 3.0 * rho / H2
    # representative 2: instant switch at end of inflation, written as a function of internal fields only
    rep_switch = -10.0 * np.where(H2 > H_end2, 1.0, 0.0)          # in units of H^2
    monotone = bool(np.all(np.diff(rho) <= 1e-15))
    after = t > t[iend]
    out = dict(t_end=float(t[iend]), chi_end=float(chi[iend]), H_end=float(np.sqrt(H_end2)),
               rho_monotone_nonincreasing=monotone,
               never_switched_rep=dict(inflation_mean=float(rep_never[~after].mean()), after_mean=float(rep_never[after].mean())),
               switched_rep=dict(inflation_value_last=float(rep_switch[iend - 1]), after_mean=float(rep_switch[after].mean()),
                                 after_max=float(rep_switch[after].max())),
               reading='rho (equivalently H) is a non-increasing internal quantity, so any switch time t_s equals '
                       'Theta(H - H(t_s)); both the never-switched and the instantly switched -10 H^2 histories have '
                       'internal representatives. P6c forbids neither.')
    return out


# ------------------------------------------------------------------ G
def V3(A, B, a, b, c):
    return -a * np.cos(A) - b * np.cos(B) + c * np.cos(A + B)


def canted_crit(a, b, c):
    ia, ib, ic = 1 / a, 1 / b, 1 / c
    return (ia < ib + ic) and (ib < ia + ic) and (ic < ia + ib)


def brute_canted(a, b, c, n=721):
    g = np.linspace(-np.pi, np.pi, n, endpoint=False)
    A, B = np.meshgrid(g, g, indexing='ij')
    VV = V3(A, B, a, b, c)
    i = np.unravel_index(np.argmin(VV), VV.shape)
    x0 = np.array([A[i], B[i]])
    res = minimize(lambda x: V3(x[0], x[1], a, b, c), x0, method='Nelder-Mead',
                   options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
    vmin = min(res.fun, VV.min())
    vcp = min(V3(p, q, a, b, c) for p in (0, np.pi) for q in (0, np.pi))
    scale = a + b + c
    return bool(vmin < vcp - 1e-10 * scale), float(vcp - vmin)


def heavy_min(cH=10.0):
    p = dict(lx=1e-3, ly=1.0, ls=1.0, kap=0.0034586847211042523, M=30.0, h=40.0, mu=7.159456159513793)

    def V(q):
        x, y, s = np.exp(q)
        return (-cH * x**2 + p['lx'] * x**4 + p['ly'] * y**4 + p['ls'] * (s**2 - p['h']**2 / 2)**2 + p['M']**2 * y**2
                - 2 * p['kap'] * y * x**3 - 2 * p['mu'] * s * y**2)
    best = None
    for x0 in (50, 80, 100):
        for y0 in (1, 3, 10):
            for s0 in (20, 28, 35):
                r = minimize(V, np.log([x0, y0, s0]), method='Nelder-Mead',
                             options=dict(xatol=1e-13, fatol=1e-13, maxiter=40000))
                if best is None or r.fun < best.fun:
                    best = r
    return np.exp(best.x), p


def part_G(seed=5):
    rng = np.random.default_rng(seed)
    n = agree = 0
    near = dict(n=0, agree=0)
    disagreements = []
    for _ in range(500):
        a, b, c = 10 ** rng.uniform(-1.5, 1.5, 3)
        cr = canted_crit(a, b, c)
        bf, gap = brute_canted(a, b, c)
        n += 1
        agree += (cr == bf)
        if cr != bf:
            disagreements.append([a, b, c, gap])
    # near-boundary triples: margin 0.5%-3% of the largest reciprocal, both sides
    for _ in range(200):
        b = 1.0
        c = 10 ** rng.uniform(-1, 1.5)
        ib, ic = 1 / b, 1 / c
        sgn = rng.choice([-1, 1])
        frac = rng.uniform(0.005, 0.03)
        # place 1/a near the boundary 1/a = ib + ic or |ib - ic|
        if rng.uniform() < 0.5:
            ia = (ib + ic) * (1 + sgn * frac)
        else:
            ia = abs(ib - ic) * (1 + sgn * frac)
        if ia <= 0:
            continue
        a = 1 / ia
        cr = canted_crit(a, b, c)
        bf, gap = brute_canted(a, b, c, n=1441)
        near['n'] += 1
        near['agree'] += (cr == bf)
    # model points
    (x, y, s), p = heavy_min()
    Mred_over_H = 2.435323e18 / 1.3335e11
    a_h = 2 * p['kap'] * y * x**3
    b_h = 2 * p['mu'] * s * y**2
    c_h = 2 * s * y * x**3 / Mred_over_H            # Sigma* Y X^3 / Lambda, |coefficient| = 1, Lambda = M_red
    # window in c at the heavy point's a/b
    ia, ib = 1 / a_h, 1 / b_h
    c_lo, c_hi = 1 / (ia + ib), 1 / abs(ia - ib)
    heavy = dict(min_fields_over_H=[x, y, s], hw_gates_quote=[80.20, 3.442, 28.311], a=a_h, b=b_h, c_natural=c_h,
                 a_over_b=a_h / b_h, c_over_b=c_h / b_h, criterion=canted_crit(a_h, b_h, c_h),
                 c_window=[c_lo, c_hi], Lambda_window_GeV=[2.435323e18 * c_h / c_hi, 2.435323e18 * c_h / c_lo])
    # gp_phase 7737487 gauge point (units H_inf = 4.7e13 GeV), quoted minimum
    Hinf = 4.7e13
    xg, yg, sg = 23.82417, 1.891386, 0.2669760
    kap_g, mu_g = 1e-3, 0.5e12 / Hinf
    a_g = 2 * kap_g * yg * xg**3
    b_g = 2 * mu_g * sg * yg**2
    c_g = 2 * sg * yg * xg**3 / (2.435323e18 / Hinf)
    gauge = dict(a=a_g, b=b_g, c_natural=c_g, a_over_b=a_g / b_g, c_over_b=c_g / b_g,
                 gp_phase_quoted_c_over_b=12.97178719966999, criterion=canted_crit(a_g, b_g, c_g))
    return dict(random=dict(n=n, agree=agree, disagreements=disagreements), near_boundary=near,
                heavy_spot_point=heavy, gp_phase_gauge_point=gauge)


# ------------------------------------------------------------------ H
def part_H():
    lam = 1e-3
    y0 = np.array([-57.31774205507421, 63.02971554304654, -351.9822243388189, 245.42938010389813])
    x, yy, vx, vy = y0
    Q = x * vy - yy * vx
    r0 = np.hypot(x, yy)
    pr0 = (x * vx + yy * vy) / r0
    E = 0.5 * (vx**2 + vy**2) + lam / 4 * r0**4

    def rhs(t, z):
        r, pr, th = z
        return [pr, Q**2 / r**3 - lam * r**3, Q / r**2]

    def ev(t, z):
        return z[1]
    ev.direction = 1
    th0 = np.arctan2(yy, x)
    sol = solve_ivp(rhs, (0, 10.0), [r0, pr0, th0], method='Radau', rtol=1e-12, atol=1e-12, events=ev,
                    dense_output=True)
    te = sol.t_events[0]
    t1, t2 = te[0], te[1]
    Tr = t2 - t1
    Om = (sol.sol(t2)[2] - sol.sol(t1)[2]) / Tr
    wr = 2 * np.pi / Tr
    # Gauss-Legendre quadrature of Fourier coefficients of f(t) = r e^{i(theta - Om (t - t1))} over one radial period
    xg, wg = np.polynomial.legendre.leggauss(400)
    nsub = 64
    edges = np.linspace(t1, t2, nsub + 1)
    ts = np.concatenate([0.5 * (b - a) * xg + 0.5 * (b + a) for a, b in zip(edges[:-1], edges[1:])])
    ws = np.concatenate([0.5 * (b - a) * wg for a, b in zip(edges[:-1], edges[1:])])
    Z = sol.sol(ts)
    f = Z[0] * np.exp(1j * (Z[2] - Om * (ts - t1)))
    ns = np.arange(-80, 81)
    cn = np.array([np.sum(ws * f * np.exp(-1j * k * wr * (ts - t1))) / Tr for k in ns])
    om_n = Om + ns * wr
    P = np.abs(cn)**2
    Ql = float(np.sum(om_n * P))
    nminus = float(-np.sum((om_n * P)[om_n < 0]))
    parseval = float(np.sum(P) / (np.sum(ws * np.abs(f)**2) / Tr))
    Emin = 0.75 * lam**(1 / 3) * abs(Q)**(4 / 3)
    # alternative particle definition (definition dependence): adiabatic number with the orbit-averaged
    # mean-field frequency w_eff^2 = lam <r^2>, n_- = <|w chi + i dchi|^2> / (4 w) style split of Q
    r_ = Z[0]
    w_eff = np.sqrt(lam * np.sum(ws * r_**2) / Tr)
    chi = r_ * np.exp(1j * Z[2])
    dr = Z[1]
    dchi = (dr + 1j * Q / r_) * np.exp(1j * Z[2])
    nm_ad = np.sum(ws * np.abs(w_eff * chi - 1j * dchi)**2) / Tr / (4 * w_eff)
    npl_ad = np.sum(ws * np.abs(w_eff * chi + 1j * dchi)**2) / Tr / (4 * w_eff)
    return dict(Q=Q, Q_lines=Ql, nbar_per_Q_lines=nminus / Q, parseval_ratio=parseval, E_over_Emin=E / Emin,
                track_value=1.7377883531119396, en_particles_value=1.7377883530390552,
                rel_diff_track=abs(nminus / Q - 1.7377883531119396) / 1.7377883531119396,
                adiabatic_meanfield=dict(
                    w_eff=float(w_eff),
                    note='Q = Im(chi* dchi) > 0 for e^{+iwt}: X number = <|w chi - i dchi|^2>/4w, Xbar = <|w chi + i dchi|^2>/4w',
                    X_over_Q=float(nm_ad / Q), Xbar_over_Q=float(npl_ad / Q),
                    check_X_minus_Xbar_over_Q=float((nm_ad - npl_ad) / Q)))


# ------------------------------------------------------------------ I
def part_I():
    Hi = 1.3335e11
    Mred = 2.435323e18
    H0 = 67.4 * 1e3 / 3.0856775814913673e22 * sc.hbar / sc.e * 1e-9      # GeV
    T0 = 2.7255 * sc.k / sc.e * 1e-9                                   # GeV
    gs0 = 2 + 7 / 8 * 2 * 3.044 * 4 / 11
    out = dict(H0_GeV=H0, T0_GeV=T0, gs0=gs0)
    for gstar in (106.75, 116.75):
        rho = 3 * Mred**2 * Hi**2
        TRH = (30 * rho / (np.pi**2 * gstar))**0.25
        a_ratio = (gs0 / gstar)**(1 / 3) * T0 / TRH
        N = (a_ratio * Hi / H0)**3
        out[f'g*={gstar}'] = dict(T_RH_instant=TRH, patches=N, residual=N**-0.5)
    return out


# ------------------------------------------------------------------ J
def part_J():
    def D(om, T):
        return 1.5 * t_P**(4 / 3) * T**(2 / 3) * om**2
    om_sr = 2 * np.pi * 429228004229873.0
    m_mol = 25000 * u
    v = 2 * np.pi * hbar / (m_mol * 53e-15)
    Tm = 2.0 / v
    om_mol = m_mol * 9.80665 * 266e-9 / hbar
    om_sph = 1e-14 * 9.80665 * 250e-6 / hbar
    req = np.sqrt(1e-3 / (1.5 * t_P**(4 / 3))) * hbar / eV
    return dict(t_P=t_P, Sr=D(om_sr, 37.0), molecule=D(om_mol, Tm), sphere=D(om_sph, 2.5), required_dE_eV_for_D_1em3_in_1s=req)


# ------------------------------------------------------------------ K
def part_K():
    import sympy as sp
    rows = {'kappa Y*X^3': (3, -1, 0), 'mu0 Sigma*Y^2': (0, 2, -1), 'c_chi chi Sigma*Y^2 (inflaton)': (0, 2, -1),
            'mu1 phi Sigma*Y^2': (0, 2, -1)}
    out = {}
    for name, keys in (('Z2_phi exact, no inflaton coupling', ['kappa Y*X^3', 'mu0 Sigma*Y^2']),
                       ('Z2_phi exact, + inflaton chi Sigma*Y^2', ['kappa Y*X^3', 'mu0 Sigma*Y^2', 'c_chi chi Sigma*Y^2 (inflaton)']),
                       ('Z2_phi broken, + inflaton', ['kappa Y*X^3', 'mu0 Sigma*Y^2', 'mu1 phi Sigma*Y^2', 'c_chi chi Sigma*Y^2 (inflaton)'])):
        Mx = sp.Matrix([rows[k] for k in keys])
        out[name] = dict(couplings=len(keys), rank=int(Mx.rank()), physical_phases=len(keys) - int(Mx.rank()))
    return out


# ------------------------------------------------------------------ L
def part_L():
    """[E] Missed postulate: no-boundary (Hartle-Hawking) state of the whole. Unweighted tree-level weight
    |Psi|^2 ~ exp(+24 pi^2 M^4 / V) = exp(8 pi^2 M^2 / H^2) (reduced Planck mass M) [R]: it prefers the lowest H."""
    Mred = 2.435323e18
    out = {}
    for name, H in (('heavy_window_H_i', 1.3335e11), ('gp_phase_H_inf', 4.7e13), ('TCC_like_1e8', 1e8)):
        out[name] = dict(H_GeV=H, exponent_8pi2M2_over_H2=8 * np.pi**2 * (Mred / H)**2)
    out['reading'] = ('relative weight of two supplied inflation scales is exp(8 pi^2 M^2 (1/H_1^2 - 1/H_2^2)); '
                      'without volume weighting the heavy window is suppressed by exp(-O(1e16)) against any lower '
                      'scale the potential allows. Measure-dependent; not evaluated further.')
    return out


def main(outpath):
    t0 = time.time()
    res = {}
    for name, fn in (('A_collision', part_A), ('B_csl', part_B), ('C_energy', part_C), ('D_dp', part_D),
                     ('E_gpp_ellipse', part_E), ('F_p6c', part_F), ('G_spontaneous_cp', part_G),
                     ('H_orbit', part_H), ('I_patches', part_I), ('J_gpp', part_J), ('K_phases', part_K), ('L_no_boundary', part_L)):
        t1 = time.time()
        res[name] = fn()
        print(f'{name} done in {time.time() - t1:.1f} s')
        print(json.dumps(res[name], indent=1, default=conv)[:2500])
    res['runtime_s'] = time.time() - t0
    with open(outpath, 'w') as f:
        json.dump(res, f, indent=1, default=conv)


if __name__ == '__main__':
    main(sys.argv[1])
