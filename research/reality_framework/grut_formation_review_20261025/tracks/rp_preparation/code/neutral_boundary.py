"""Track rp_preparation (Claude, issue #13): independent reproduction of the finite
gauge-neutral boundary example of Codex e160356 (BOUNDARY.md) and checks of the
symmetry argument behind "these conditions do not select the phase".

Own code. Method differs from Codex's boundary.py: ladder operators are built on the
full truncated three-mode Fock space as Kronecker products, the Hamiltonian is formed
from them, and everything is then restricted to the exactly neutral block
n_a + 3 n_b - 6 n_c = 0. Codex enumerates the neutral basis and writes the matrix
elements by hand. No code of Codex's or of our earlier tracks is reused.

Usage: python3 neutral_boundary.py OUT.json
"""
import sys, json, math, itertools
import numpy as np
import scipy.sparse as sp

CHARGES = np.array([1, 3, -6])
FREQ = np.array([1.0, 1.7, 2.3])
G1, G2 = 0.025, 0.04            # b^dag a^3 and c b^2 couplings
MEANS = np.array([6.0, 2.0, 2.0])  # |alpha|^2 of the pre-projection coherent state
N_CLOCK = 5
DT = 0.1


def ladder(cut):
    return sp.diags(np.sqrt(np.arange(1, cut + 1)), 1, format='csr')


def full_ops(cuts):
    dims = [c + 1 for c in cuts]
    ident = [sp.identity(d, format='csr') for d in dims]

    def embed(op, k):
        mats = list(ident)
        mats[k] = op
        return sp.kron(sp.kron(mats[0], mats[1]), mats[2], format='csr')

    A = [embed(ladder(c), k) for k, c in enumerate(cuts)]
    N = [embed(sp.diags(np.arange(d, dtype=float)), k) for k, d in enumerate(dims)]
    occ = np.array(list(itertools.product(*(range(d) for d in dims))))
    return A, N, occ, dims


def hamiltonian(A, N, chi=0.0):
    a, b, c = A
    H0 = FREQ[0] * N[0] + FREQ[1] * N[1] + FREQ[2] * N[2]
    T1 = b.conj().T @ (a @ a @ a)   # b^dag a^3  (charge +3-3 = 0)
    T2 = c @ (b @ b)                # c b^2      (charge -(-6) - 6 = 0)
    g2 = G2 * np.exp(1j * chi)
    H = H0 + G1 * (T1 + T1.conj().T) + g2 * T2 + np.conj(g2) * T2.conj().T
    if chi == 0.0:
        H = H.real
    return H.tocsr()


def coherent_full(occ, delta):
    alpha = np.sqrt(MEANS) * np.exp(1j * np.array([0.0, 0.0, delta]))
    logfact = np.array([math.lgamma(n + 1) for n in range(occ.max() + 1)])
    amp = np.ones(len(occ), dtype=complex)
    for k in range(3):
        amp *= alpha[k] ** occ[:, k] / np.exp(0.5 * logfact[occ[:, k]])
    return amp


def entropy(p):
    p = np.asarray(p, dtype=float)
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)))


class Model:
    def __init__(self, cuts, chi=0.0):
        self.cuts = tuple(cuts)
        A, N, occ, dims = full_ops(cuts)
        self.dims = dims
        self.occ = occ
        Hfull = hamiltonian(A, N, chi)
        q = occ @ CHARGES
        self.mask = q == 0
        idx = np.flatnonzero(self.mask)
        self.idx = idx
        self.d = len(idx)
        # charge conservation of the full H: no matrix element leaves the neutral block
        off = Hfull[np.flatnonzero(~self.mask)][:, idx]
        self.offblock_norm = float(abs(off).max()) if off.nnz else 0.0
        self.H = Hfull[idx][:, idx].toarray()
        self.nocc = occ[idx]
        self.L = np.diag((self.nocc[:, 0] + 3 * self.nocc[:, 1]).astype(float))
        self.Q = self.nocc @ CHARGES
        self.E, self.U = np.linalg.eigh(self.H) if chi == 0.0 else np.linalg.eigh(self.H)
        self.commHL = self.H @ self.L - self.L @ self.H
        self.is_real = chi == 0.0

    def state(self, delta):
        amp = coherent_full(self.occ, delta)
        trunc_norm2 = float(np.vdot(amp, amp).real)
        vec = amp[self.idx]
        norm2 = float(np.vdot(vec, vec).real)
        return vec / math.sqrt(norm2), norm2 / trunc_norm2, norm2 * math.exp(-MEANS.sum())

    def energy(self, psi):
        return float(np.vdot(psi, self.H @ psi).real)

    def torque(self, psi):
        # d<L>/dt = i <[H, L]>
        return float((1j * np.vdot(psi, self.commHL @ psi)).real)

    def torque_rho(self, rho):
        return float((1j * np.trace(rho @ self.commHL)).real)

    def marginals(self, psi):
        w = abs(psi) ** 2
        return [np.bincount(self.nocc[:, k], weights=w, minlength=self.dims[k]) for k in range(3)]

    def reduced_vn_entropies(self, psi):
        full = np.zeros(len(self.occ), dtype=complex)
        full[self.idx] = psi
        T = full.reshape(self.dims)
        out = []
        for k in range(3):
            Tk = np.moveaxis(T, k, 0).reshape(self.dims[k], -1)
            rho = Tk @ Tk.conj().T
            out.append(entropy(np.linalg.eigvalsh(rho)))
        return out

    def row(self, delta):
        psi, p_trunc, p_exact = self.state(delta)
        w = abs(psi) ** 2
        marg = self.marginals(psi)
        sm = np.diag([0.2, 0.8])
        rho_h = np.outer(psi, psi.conj())
        joint = np.kron(rho_h, sm)
        mi = entropy(np.linalg.eigvalsh(rho_h)) + entropy([0.2, 0.8]) - entropy(np.linalg.eigvalsh(joint))
        return dict(delta=delta,
                    neutral_projection_probability_truncated=p_trunc,
                    neutral_projection_probability_untruncated=p_exact,
                    mean_occupations=[float(w @ self.nocc[:, k]) for k in range(3)],
                    mean_gauge_charge=float(w @ self.Q),
                    energy=self.energy(psi),
                    light_charge=float(np.vdot(psi, self.L @ psi).real),
                    light_charge_torque=self.torque(psi),
                    mode_entropies_shannon_nats=[entropy(m) for m in marg],
                    mode_entropies_vonNeumann_nats=self.reduced_vn_entropies(psi),
                    hidden_SM_mutual_information=mi)

    # ---- internal clock --------------------------------------------------
    def clock(self, deltas, events):
        d, N = self.d, N_CLOCK
        step = (self.U * np.exp(-1j * DT * self.E)) @ self.U.conj().T
        Ws = [np.eye(d, dtype=complex)]
        for t in range(1, N):
            Ws.append(step @ Ws[-1])
        V = np.concatenate(Ws, axis=0) / math.sqrt(N)
        P = V @ V.conj().T
        C = np.zeros((N * d, N * d), dtype=complex)
        for t in range(N - 1):
            s0 = slice(t * d, (t + 1) * d)
            s1 = slice((t + 1) * d, (t + 2) * d)
            C[s0, s0] += 0.5 * np.eye(d)
            C[s1, s1] += 0.5 * np.eye(d)
            C[s1, s0] -= 0.5 * step
            C[s0, s1] -= 0.5 * step.conj().T
        evals = np.linalg.eigvalsh(C)
        kernel_dim = int(np.sum(evals < 1e-9))
        gap = float(np.sort(evals)[kernel_dim])
        # clock reversal R: |t> -> |N-1-t>; combined with complex conjugation K
        R = np.kron(np.eye(N)[::-1], np.eye(d))
        rk_invariance = float(np.linalg.norm(R @ C.conj() @ R - C))
        out = dict(isometry_error=float(np.linalg.norm(V.conj().T @ V - np.eye(d))),
                   constraint_on_isometry=float(np.linalg.norm(C @ V)),
                   constraint_spectrum_kernel_dim=kernel_dim,
                   constraint_first_gap=gap,
                   constraint_first_gap_analytic=1 - math.cos(math.pi / N),
                   reversal_conjugation_invariance_of_C=rk_invariance,
                   rows=[])
        for delta in deltas:
            psi, _, _ = self.state(delta)
            hist = V @ psi
            res = dict(delta=delta,
                       history_norm=float(np.vdot(hist, hist).real),
                       constraint_residual=float(np.linalg.norm(C @ hist)),
                       max_effect_commutator=0.0, max_born_mismatch=0.0, events=[])
            for name, proj in events:
                for t in range(N):
                    Pi = np.zeros((N * d, N * d))
                    sl = slice(t * d, (t + 1) * d)
                    Pi[sl, sl] = proj
                    eff = N * P @ Pi @ P
                    comm = float(np.linalg.norm(C @ eff - eff @ C))
                    born = float(np.vdot(hist, eff @ hist).real)
                    psi_t = Ws[t] @ psi
                    direct = float(np.vdot(psi_t, proj @ psi_t).real)
                    res['max_effect_commutator'] = max(res['max_effect_commutator'], comm)
                    res['max_born_mismatch'] = max(res['max_born_mismatch'], abs(born - direct))
                    if t == 2 and name == 'n_a>=6':
                        res['events'].append(dict(event=name, reading=t, physical=born, direct=direct))
            out['rows'].append(res)
        return out


def main(out_path):
    results = dict(scope=__doc__)
    deltas = [0.0, math.pi / 2, -math.pi / 2, math.pi]
    models = []
    for cuts in [(18, 8, 6), (24, 12, 9), (30, 15, 12), (36, 18, 15)]:
        m = Model(cuts)
        rows = [m.row(dl) for dl in deltas]
        plus, minus = rows[1], rows[2]
        # exact conjugation check: conj(psi_delta) == psi_{-delta}
        pp, _, _ = m.state(math.pi / 2)
        pm, _, _ = m.state(-math.pi / 2)
        conj_overlap = abs(np.vdot(pm, pp.conj()))
        # H real in the number basis
        h_real = float(abs(m.H.imag).max()) if np.iscomplexobj(m.H) else 0.0
        h_sym = float(abs(m.H - m.H.T).max())
        # delta scan: energy, torque and 6 dE/d delta
        grid = np.linspace(-math.pi, math.pi, 49)
        Es, Ts, dE6 = [], [], []
        h = 1e-5
        for dl in grid:
            psi, _, _ = m.state(dl)
            Es.append(m.energy(psi)); Ts.append(m.torque(psi))
            ep = m.energy(m.state(dl + h)[0]); em = m.energy(m.state(dl - h)[0])
            dE6.append(6 * (ep - em) / (2 * h))
        Es, Ts, dE6 = map(np.array, (Es, Ts, dE6))
        even_err = float(abs(Es - Es[::-1]).max())
        odd_err = float(abs(Ts + Ts[::-1]).max())
        grad_err = float(abs(Ts - dE6).max())
        imin = int(np.argmin(Es))
        # rotation identity exp(i theta L) psi_delta = psi_{delta + 6 theta}
        theta = 0.37
        rot = np.exp(1j * theta * np.diag(m.L)) * m.state(0.3)[0]
        tgt = m.state(0.3 + 6 * theta)[0]
        rot_overlap = abs(np.vdot(tgt, rot))
        # stationary states have zero torque
        eig_torques = [m.torque(m.U[:, k]) for k in range(m.d)]
        gibbs = []
        for beta in (0.1, 1.0, 10.0):
            w = np.exp(-beta * (m.E - m.E.min())); w /= w.sum()
            rho = (m.U * w) @ m.U.conj().T
            gibbs.append(dict(beta=beta, torque=m.torque_rho(rho)))
        # doubled conjugate example (explicit product state) for the two smaller cutoffs
        doubled = None
        if m.d <= 60:
            I = np.eye(m.d)
            Htot = np.kron(m.H, I) + np.kron(I, m.H)
            Lsg = np.kron(m.L, I) - np.kron(I, m.L)
            comm = Htot @ Lsg - Lsg @ Htot
            Psi = np.kron(pp, pm)
            Psi_rev = np.kron(pm, pp)
            Psi_sym = (Psi + Psi_rev); Psi_sym /= np.linalg.norm(Psi_sym)
            tq = lambda v: float((1j * np.vdot(v, comm @ v)).real)
            doubled = dict(initial_signed_light_charge=float(np.vdot(Psi, Lsg @ Psi).real),
                           signed_torque=tq(Psi), reversed_phase_signed_torque=tq(Psi_rev),
                           swap_symmetrized_signed_torque=tq(Psi_sym),
                           energy=float(np.vdot(Psi, Htot @ Psi).real),
                           twice_single_copy_torque=2 * plus['light_charge_torque'])
        clock = m.clock([0.0, math.pi / 2, -math.pi / 2, math.pi],
                        [('n_a>=6', np.diag((m.nocc[:, 0] >= 6).astype(float))),
                         ('n_b>=2', np.diag((m.nocc[:, 1] >= 2).astype(float))),
                         ('n_c>=2', np.diag((m.nocc[:, 2] >= 2).astype(float)))]) if m.d <= 60 else None
        models.append(dict(cutoffs=list(cuts), full_dimension=int(len(m.occ)), neutral_dimension=m.d,
                           H_offblock_norm=m.offblock_norm, H_imag_max=h_real, H_asymmetry=h_sym,
                           rows=rows,
                           conjugation_maps_plus_to_minus_overlap=float(conj_overlap),
                           rotation_identity_overlap=float(rot_overlap),
                           symmetry_checks=dict(energy_difference=abs(plus['energy'] - minus['energy']),
                                                torque_sum=abs(plus['light_charge_torque'] + minus['light_charge_torque']),
                                                occupation_difference=max(abs(a - b) for a, b in zip(plus['mean_occupations'], minus['mean_occupations'])),
                                                entropy_difference=max(abs(a - b) for a, b in zip(plus['mode_entropies_shannon_nats'], minus['mode_entropies_shannon_nats'])),
                                                shannon_vs_vonNeumann=max(abs(a - b) for r in rows for a, b in zip(r['mode_entropies_shannon_nats'], r['mode_entropies_vonNeumann_nats']))),
                           delta_scan=dict(grid=grid.tolist(), energy=Es.tolist(), torque=Ts.tolist(),
                                           six_dE_ddelta=dE6.tolist(), energy_even_error=even_err,
                                           torque_odd_error=odd_err, torque_minus_gradient_max=grad_err,
                                           energy_minimum_delta=float(grid[imin]), torque_at_energy_minimum=float(Ts[imin]),
                                           torque_amplitude=float(abs(Ts).max()),
                                           coherent_state_torque_amplitude=0.48 * 2 * math.sqrt(2)),
                           stationary_states=dict(max_eigenstate_torque=float(max(abs(t) for t in eig_torques)), gibbs=gibbs),
                           doubled_conjugate_example=doubled, clock=clock))
    results['models'] = models
    # cutoff convergence of the pi/2 row
    ref = models[-1]['rows'][1]
    results['cutoff_convergence_vs_largest'] = [
        dict(cutoffs=mm['cutoffs'], neutral_dimension=mm['neutral_dimension'],
             **{k: abs(mm['rows'][1][k] / ref[k] - 1) for k in ('energy', 'light_charge', 'light_charge_torque')})
        for mm in models[:-1]]

    # CP-violating coupling: complex c b^2 coefficient with phase chi
    cpv = []
    for chi in (0.3, 1.0):
        m = Model((18, 8, 6), chi=chi)
        grid = np.linspace(-math.pi, math.pi, 721)
        Es = np.array([m.energy(m.state(dl)[0]) for dl in grid])
        Ts = np.array([m.torque(m.state(dl)[0]) for dl in grid])
        imin = int(np.argmin(Es))
        pp = m.state(math.pi / 2)[0]; pm = m.state(-math.pi / 2)[0]
        cpv.append(dict(chi=chi, H_hermitian_error=float(abs(m.H - m.H.conj().T).max()),
                        energy_plus=m.energy(pp), energy_minus=m.energy(pm),
                        torque_plus=m.torque(pp), torque_minus=m.torque(pm),
                        energy_minimum_delta=float(grid[imin]), expected_pi_minus_chi=math.pi - chi,
                        torque_at_energy_minimum=float(Ts[imin]),
                        torque_amplitude=float(abs(Ts).max())))
    results['cp_violating_coupling'] = cpv

    # SM-correlated preparation: the SM bit records the sign of delta
    m = Model((18, 8, 6))
    pp = m.state(math.pi / 2)[0]; pm = m.state(-math.pi / 2)[0]
    corr = []
    for p in (0.5, 0.6, 0.9):
        joint = p * np.kron(np.outer(pp, pp.conj()), np.diag([1.0, 0.0])) + (1 - p) * np.kron(np.outer(pm, pm.conj()), np.diag([0.0, 1.0]))
        rho_h = p * np.outer(pp, pp.conj()) + (1 - p) * np.outer(pm, pm.conj())
        mi = entropy(np.linalg.eigvalsh(rho_h)) + entropy([p, 1 - p]) - entropy(np.linalg.eigvalsh(joint))
        corr.append(dict(p_plus=p, mutual_information_nats=mi, hidden_torque=m.torque_rho(rho_h),
                         sm_bit_asymmetry=2 * p - 1))
    results['sm_correlated_preparation'] = corr

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=1)
    for mm in models:
        r = mm['rows']
        print(mm['cutoffs'], 'd=', mm['neutral_dimension'], 'E(0,pi/2,pi)=', [round(x['energy'], 6) for x in (r[0], r[1], r[3])],
              'T(+pi/2)=', round(r[1]['light_charge_torque'], 8), 'S=', [round(s, 6) for s in r[1]['mode_entropies_shannon_nats']],
              'gradient id err=', mm['delta_scan']['torque_minus_gradient_max'], 'Emin at', round(mm['delta_scan']['energy_minimum_delta'], 4),
              'T there', mm['delta_scan']['torque_at_energy_minimum'])
        if mm['doubled_conjugate_example']:
            print('  doubled:', mm['doubled_conjugate_example'])
        if mm['clock']:
            c = mm['clock']
            print('  clock: iso', c['isometry_error'], 'CV', c['constraint_on_isometry'], 'kernel', c['constraint_spectrum_kernel_dim'],
                  'gap', c['constraint_first_gap'], c['constraint_first_gap_analytic'], 'RK', c['reversal_conjugation_invariance_of_C'],
                  'born', max(rw['max_born_mismatch'] for rw in c['rows']), 'comm', max(rw['max_effect_commutator'] for rw in c['rows']),
                  'event', c['rows'][1]['events'])
    print('cutoff convergence', results['cutoff_convergence_vs_largest'])
    print('CPV', cpv)
    print('SM corr', corr)


if __name__ == '__main__':
    main(sys.argv[1])
