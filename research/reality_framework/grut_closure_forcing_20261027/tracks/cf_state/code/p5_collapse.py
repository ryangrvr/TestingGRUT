"""P5 (global purity + unitary whole): what it forbids and how the forbidden models are bounded (own code).

Part A [C, D]: four dynamics for a two-site 'position' qubit (H = (Delta/2) sigma_x tunnelling, localisation L = sigma_z)
  (1) Lindblad master equation  d rho/dt = -i[H, rho] + gamma (L rho L - rho);
  (2) nonlinear stochastic collapse unravelling (QSD / CSL type, Ito):
        d psi = [-iH dt - (gamma/2)(L - <L>)^2 dt + sqrt(gamma)(L - <L>) dW] psi ;
  (3) linear random-unitary unravelling (classical noise, no collapse): psi -> exp(-i sqrt(gamma) dW L) exp(-iH dt) psi;
  (4) a closed unitary whole: system + K fresh ancilla qubits (collision model), U_k = exp(-i theta sigma_z (x) sigma_x),
      theta^2 = gamma dt, so each collision multiplies the coherence by cos(2 theta).
  Claims checked: one-time statistics of (1)-(4) coincide (to sampling error / O(dt)); individual trajectories of (2)
  localise with Born frequencies while those of (3) do not; only (4) can be reversed (global inverse unitary returns
  the system to its pure initial state, fidelity 1).
Part B [C]: CSL and DP predictions for the molecule interferometer and the proposed 1e-14 kg superposition, at the
  GRW and Adler parameters and at the quoted current bounds; white-noise CSL energy non-conservation (heating rate).
Usage: python3 -B p5_collapse.py ../results/p5_collapse.json
"""
import sys
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
from common import hbar, G, u, m0_CSL, QUOTED, dump

sx = np.array([[0, 1], [1, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)


def lindblad_rho(rho0, H, L, gamma, t):
    d = 2
    Id = np.eye(d)
    # vec(A rho B) = (B^T kron A) vec(rho) with column stacking
    Lsup = -1j * (np.kron(Id, H) - np.kron(H.T, Id)) + gamma * (np.kron(L.T, L) - np.kron(Id, Id))
    v = expm(Lsup * t) @ rho0.reshape(-1, order='F')
    return v.reshape(d, d, order='F')


def part_A(seed=7):
    rng = np.random.default_rng(seed)
    Delta, gamma = 1.0, 0.5
    H = 0.5 * Delta * sx
    L = sz
    psi0 = np.array([1, 1], complex) / np.sqrt(2)
    rho0 = np.outer(psi0, psi0.conj())
    tf, dt = 3.0, 1e-3
    nsteps = int(round(tf / dt))
    U0 = expm(-1j * H * dt)
    ntraj = 20000
    # (2) nonlinear QSD collapse and (3) linear random unitary, same Wiener increments (vectorised over trajectories)
    lz = np.array([1.0, -1.0])
    A_ = np.tile(psi0, (ntraj, 1))
    B_ = np.tile(psi0, (ntraj, 1))
    for n in range(nsteps):
        dW = rng.normal(0, np.sqrt(dt), ntraj)
        Lm = np.real(np.sum(np.abs(A_)**2 * lz, axis=1))
        La = A_ * lz - Lm[:, None] * A_
        LLa = La * lz - Lm[:, None] * La
        HA = A_ @ H.T
        A_ = A_ + (-1j * HA * dt - 0.5 * gamma * LLa * dt + np.sqrt(gamma) * La * dW[:, None])
        A_ /= np.linalg.norm(A_, axis=1)[:, None]
        B_ = B_ @ U0.T
        ph = np.sqrt(gamma) * dW
        B_ = B_ * np.exp(-1j * ph[:, None] * lz)
    rho_qsd = np.einsum('ti,tj->ij', A_, A_.conj()) / ntraj
    rho_ru = np.einsum('ti,tj->ij', B_, B_.conj()) / ntraj
    final_Lexp_qsd = np.real(np.sum(np.abs(A_)**2 * lz, axis=1))
    final_Lexp_ru = np.real(np.sum(np.abs(B_)**2 * lz, axis=1))
    rho_L = lindblad_rho(rho0, H, L, gamma, tf)
    # (4) collision model: closed unitary whole, K ancillas, coarse step tau (theta^2 = gamma tau)
    K = 12
    tau = tf / K
    theta = np.sqrt(gamma * tau)
    Uc = expm(-1j * theta * np.kron(sz, sx))           # system (x) one ancilla
    Uh = expm(-1j * H * tau)
    dim = 2 ** (K + 1)
    state = np.zeros(dim, complex)
    # system is the most significant qubit; ancillas start in |0>
    state[0] = psi0[0]
    state[2 ** K] = psi0[1]
    psi = state.reshape([2] * (K + 1))
    ops = []
    for k in range(K):
        # free step on system
        psi = np.tensordot(Uh, psi, axes=([1], [0]))
        # collision with ancilla k (axis k+1)
        U4 = Uc.reshape(2, 2, 2, 2)
        psi = np.tensordot(U4, psi, axes=([2, 3], [0, k + 1]))
        psi = np.moveaxis(psi, [0, 1], [0, k + 1])
        ops.append(k)
    rho_sys = np.tensordot(psi, psi.conj(), axes=(list(range(1, K + 1)), list(range(1, K + 1))))
    # Lindblad comparison at the same coarse Trotterisation: per step coherence factor cos(2 theta) vs exp(-2 gamma tau)
    rho_L_coarse = rho0.copy()
    for k in range(K):
        rho_L_coarse = Uh @ rho_L_coarse @ Uh.conj().T
        rho_L_coarse = np.array([[rho_L_coarse[0, 0], np.cos(2 * theta) * rho_L_coarse[0, 1]],
                                 [np.cos(2 * theta) * rho_L_coarse[1, 0], rho_L_coarse[1, 1]]])
    # reversal of the closed whole: apply the inverse of every step in reverse order
    psi_r = psi.copy()
    for k in reversed(range(K)):
        U4i = Uc.conj().T.reshape(2, 2, 2, 2)
        psi_r = np.tensordot(U4i, psi_r, axes=([2, 3], [0, k + 1]))
        psi_r = np.moveaxis(psi_r, [0, 1], [0, k + 1])
        psi_r = np.tensordot(Uh.conj().T, psi_r, axes=([1], [0]))
    rho_back = np.tensordot(psi_r, psi_r.conj(), axes=(list(range(1, K + 1)), list(range(1, K + 1))))
    fid_back = np.real(np.vdot(psi0, rho_back @ psi0))
    # Born-rule check for the collapse trajectories (H = 0 run): fraction ending near +1 vs initial weight
    p_up = 0.3
    a0 = np.array([np.sqrt(p_up), np.sqrt(1 - p_up)], complex)
    nb = 4000
    A_ = np.tile(a0, (nb, 1))
    for n in range(4000):
        dW = rng.normal(0, np.sqrt(dt), nb)
        Lm = np.real(np.sum(np.abs(A_)**2 * lz, axis=1))
        La = A_ * lz - Lm[:, None] * A_
        LLa = La * lz - Lm[:, None] * La
        A_ = A_ + (-0.5 * gamma * LLa * dt + np.sqrt(gamma) * La * dW[:, None])
        A_ /= np.linalg.norm(A_, axis=1)[:, None]
    Lend = np.real(np.sum(np.abs(A_)**2 * lz, axis=1))
    ups = int(np.sum(Lend > 0))
    collapsed_fraction = float(np.mean(np.abs(Lend) > 0.99))
    stderr = 1.0 / np.sqrt(ntraj)
    out = dict(params=dict(Delta=Delta, gamma=gamma, t=tf, dt=dt, ntraj=ntraj, K_ancillas=K),
               rho_lindblad=rho_L, rho_qsd_avg=rho_qsd, rho_random_unitary_avg=rho_ru,
               dev_qsd=float(np.abs(rho_qsd - rho_L).max()), dev_ru=float(np.abs(rho_ru - rho_L).max()),
               sampling_error_scale=stderr,
               rho_collision_model=rho_sys, rho_lindblad_same_trotter=rho_L_coarse,
               dev_collision_vs_trotter_lindblad=float(np.abs(rho_sys - rho_L_coarse).max()),
               dev_collision_vs_exact_lindblad=float(np.abs(rho_sys - rho_L).max()),
               reversal_fidelity_closed_whole=fid_back,
               qsd_final_abs_L_mean=float(np.mean(np.abs(final_Lexp_qsd))),
               ru_final_abs_L_mean=float(np.mean(np.abs(final_Lexp_ru))),
               born_check=dict(p_up=p_up, fraction_collapsed_up=ups / nb, n=nb, stderr=np.sqrt(p_up * (1 - p_up) / nb), fraction_with_abs_L_gt_0p99=collapsed_fraction))
    return out


# ---------------- Part B: CSL and DP numbers ----------------
def csl_rate_point(lam, m, dx, rC):
    return lam * (m / m0_CSL)**2 * (1 - np.exp(-dx**2 / (4 * rC**2)))


def csl_rate_sphere(lam, m, R, dx, rC):
    """Gamma = lam rC^3/(pi^{3/2} m0^2) int d^3k e^{-k^2 rC^2} |mu(k)|^2 (1 - cos k.dx), homogeneous sphere.
    |mu(k)|^2 = m^2 [3 j1(kR)/(kR)]^2; angular average of (1 - cos k.dx) = 1 - sin(k dx)/(k dx)."""
    def f(k):
        x = k * R
        F = 3 * (np.sin(x) - x * np.cos(x)) / x**3 if x > 1e-4 else 1 - x**2 / 10
        y = k * dx
        ang = 1 - np.sin(y) / y if y > 1e-6 else y**2 / 6
        return 4 * np.pi * k**2 * np.exp(-k**2 * rC**2) * F**2 * ang
    kmax = 12 / rC
    pts = np.linspace(0, kmax, 4001)
    val = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        val += quad(f, a, b, limit=200)[0]
    return lam * rC**3 / (np.pi**1.5 * m0_CSL**2) * m**2 * val


def dp_EG_sphere(m, R, d):
    """Penrose E_G for a homogeneous sphere displaced by d >= 2R: G m^2 (6/(5R) - 1/d)."""
    return G * m**2 * (6 / (5 * R) - 1 / d)


def part_B():
    rC = 1e-7
    lam = {
        'GRW': QUOTED['CSL_GRW_lambda'][0],
        'Adler': QUOTED['CSL_Adler_lambda'][0],
        'bound_LPF': QUOTED['CSL_bound_LPF'][0],
        'bound_cantilever': QUOTED['CSL_bound_cantilever'][0],
        'bound_bulk_heating': QUOTED['CSL_bound_bulkheat'][0],
        'bound_Xray_GranSasso_white': QUOTED['CSL_bound_Xray_GranSasso'][0],
        'bound_Majorana_mp_nuclei_white': QUOTED['CSL_bound_Majorana_mp_nuclei'][0] * rC**2,
    }
    res = {'rC': rC, 'lambda_values': lam}
    # point-sphere check of the integral against the point formula (R << rC)
    m_test = 1e4 * u
    res['integral_check_point_limit'] = dict(
        integral=csl_rate_sphere(1.0, m_test, 1e-10, 3e-7, rC), point=csl_rate_point(1.0, m_test, 3e-7, rC))
    # Fein molecule (point-like on rC scale), separation = grating period
    m = QUOTED['Fein_mass_u'][0] * u
    d = QUOTED['Fein_grating_period'][0]
    v = 2 * np.pi * hbar / (m * QUOTED['Fein_lambda_dB'][0])
    T = QUOTED['Fein_length'][0] / v
    mol = {}
    for k, l in lam.items():
        G_ = csl_rate_sphere(l, m, 1.5e-9, d, rC)
        mol[k] = dict(rate=G_, exponent_over_flight=G_ * T)
    mol['flight_time'] = T
    mol['sensitivity_exponent'] = -np.log(QUOTED['Fein_visibility_fraction'][0])
    mol['lambda_reaching_sensitivity'] = mol['sensitivity_exponent'] / T / (csl_rate_sphere(1.0, m, 1.5e-9, d, rC))
    EG = G * m**2 * 6 / (5 * 1.5e-9)
    mol['DP_continuum_EG_J'] = EG
    mol['DP_time_s'] = hbar / EG
    res['molecule_Fein'] = mol
    # proposed 1e-14 kg sphere (silica density 2200), dx = 250 um
    mM = QUOTED['Bose_mass'][0]
    R = (3 * mM / (4 * np.pi * 2200.0))**(1 / 3)
    dx = QUOTED['Bose_dx'][0]
    mac = {'radius_m': R}
    for k, l in lam.items():
        G_ = csl_rate_sphere(l, mM, R, dx, rC)
        mac[k] = dict(rate=G_, time=1 / G_)
    EG = dp_EG_sphere(mM, R, dx)
    mac['DP_EG_J'] = EG
    mac['DP_time_s'] = hbar / EG
    mac['protocol_time_s'] = QUOTED['Bose_tau'][0]
    # granular (nuclear) DP part at R0 = bound: N_nuc * (6/5) G m_n^2 / R0, mean nucleus 20 u for SiO2
    mn = 20 * u
    R0 = QUOTED['DP_R0_bound_Majorana'][0]
    mac['DP_granular_EG_J_at_R0_bound'] = (mM / mn) * 1.2 * G * mn**2 / R0
    res['macroscopic_proposed'] = mac
    # white-noise mass-proportional CSL heating: dE/dt = 3 lam hbar^2 m/(4 m0^2 rC^2) per particle of mass m
    heat = {}
    for k, l in lam.items():
        heat[k] = dict(per_nucleon_W=3 * l * hbar**2 * u / (4 * m0_CSL**2 * rC**2),
                       earth_W=3 * l * hbar**2 * 5.97e24 / (4 * m0_CSL**2 * rC**2),
                       temperature_rate_K_per_s_gas=(2 / 3) * 3 * l * hbar**2 / (4 * m0_CSL * rC**2) / 1.380649e-23)
    res['white_CSL_heating'] = heat
    return res


def main(out):
    A = part_A()
    B = part_B()
    dump(dict(part_A=A, part_B=B, quoted={k: v for k, v in QUOTED.items() if k.startswith(('CSL', 'DP'))}), out)
    print('Part A: max |rho - rho_Lindblad|: QSD avg', A['dev_qsd'], ' random-unitary avg', A['dev_ru'],
          ' (1/sqrt(ntraj) =', A['sampling_error_scale'], ')')
    print('        collision model vs Trotterised Lindblad', A['dev_collision_vs_trotter_lindblad'],
          ' vs exact Lindblad', A['dev_collision_vs_exact_lindblad'])
    print('        reversal fidelity (closed whole):', A['reversal_fidelity_closed_whole'])
    print('        <|<L>|> at t_f: QSD', A['qsd_final_abs_L_mean'], ' random unitary', A['ru_final_abs_L_mean'])
    print('        Born check:', A['born_check'])
    print('Part B:')
    import json
    print(json.dumps(B, indent=1, default=float)[:6000])


if __name__ == '__main__':
    main(sys.argv[1])
