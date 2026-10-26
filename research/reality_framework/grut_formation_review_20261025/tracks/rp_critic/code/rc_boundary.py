#!/usr/bin/env python3
"""rp_critic: own build of Codex's three-mode neutral boundary model (charges +1,+3,-6), to recompute the
load-bearing numbers (neutral dimension, E(delta), torque, projection probability) and check the two identities
rp_preparation states: T = 6 dE/ddelta and <c b^2>_proj = alpha_c alpha_b^2.  Independent of both codes:
the neutral block is enumerated directly and matrix elements are written from ladder-operator rules.
usage: rc_boundary.py OUT.json
"""
import sys, json, math, itertools
import numpy as np

def build(cut, delta, coup=(0.025, 0.04), freqs=(1.0, 1.7, 2.3), amps=(math.sqrt(6), math.sqrt(2), math.sqrt(2))):
    Na, Nb, Nc = cut
    basis = [(na, nb, nc) for na in range(Na + 1) for nb in range(Nb + 1) for nc in range(Nc + 1) if na + 3 * nb - 6 * nc == 0]
    idx = {s: k for k, s in enumerate(basis)}
    d = len(basis)
    H = np.zeros((d, d), dtype=complex)
    L = np.zeros(d)
    for s, k in idx.items():
        na, nb, nc = s
        H[k, k] = freqs[0] * na + freqs[1] * nb + freqs[2] * nc
        L[k] = na + 3 * nb
        # b^dagger a^3 : (na, nb, nc) -> (na-3, nb+1, nc), amplitude sqrt(na(na-1)(na-2)) sqrt(nb+1)
        if na >= 3 and nb + 1 <= Nb:
            t = (na - 3, nb + 1, nc)
            if t in idx:
                amp = coup[0] * math.sqrt(na * (na - 1) * (na - 2) * (nb + 1))
                H[idx[t], k] += amp; H[k, idx[t]] += amp
        # c b^2 : (na, nb, nc) -> (na, nb-2, nc-1), amplitude sqrt(nb(nb-1)) sqrt(nc)
        if nb >= 2 and nc >= 1:
            t = (na, nb - 2, nc - 1)
            if t in idx:
                amp = coup[1] * math.sqrt(nb * (nb - 1) * nc)
                H[idx[t], k] += amp; H[k, idx[t]] += amp
    # projected coherent state: amplitude prod alpha^n/sqrt(n!) with alpha_c = sqrt2 e^{i delta}
    ac = amps[2] * np.exp(1j * delta)
    psi = np.array([amps[0] ** na / math.sqrt(math.factorial(na)) * amps[1] ** nb / math.sqrt(math.factorial(nb))
                    * ac ** nc / math.sqrt(math.factorial(nc)) for (na, nb, nc) in basis], dtype=complex)
    norm_trunc = sum(abs(amps[0] ** na / math.sqrt(math.factorial(na)) * amps[1] ** nb / math.sqrt(math.factorial(nb))
                         * amps[2] ** nc / math.sqrt(math.factorial(nc))) ** 2
                     for na in range(Na + 1) for nb in range(Nb + 1) for nc in range(Nc + 1))
    pproj = float(np.vdot(psi, psi).real / norm_trunc)
    psi = psi / np.linalg.norm(psi)
    E = float(np.vdot(psi, H @ psi).real)
    Lop = np.diag(L)
    comm = H @ Lop - Lop @ H
    T = float((1j * np.vdot(psi, comm @ psi)).real)
    # <c b^2> in the projected state
    cb2 = np.zeros((d, d), dtype=complex)
    for s, k in idx.items():
        na, nb, nc = s
        if nb >= 2 and nc >= 1:
            t = (na, nb - 2, nc - 1)
            if t in idx:
                cb2[idx[t], k] = math.sqrt(nb * (nb - 1) * nc)
    cb2_exp = complex(np.vdot(psi, cb2 @ psi))
    occ = [float(sum(abs(psi[k]) ** 2 * basis[k][m] for k in range(d))) for m in range(3)]
    return dict(dim=d, pproj=pproj, E=E, T=T, cb2=cb2_exp, occ=occ, psi=psi, H=H, L=L)

if __name__ == '__main__':
    out = {}
    for cut in [(18, 8, 6), (24, 12, 9), (30, 15, 12)]:
        rows = {}
        for name, dl in [('0', 0.0), ('+pi/2', math.pi / 2), ('-pi/2', -math.pi / 2), ('pi', math.pi)]:
            r = build(cut, dl)
            rows[name] = dict(dim=r['dim'], pproj=r['pproj'], E=r['E'], T=r['T'], cb2=[r['cb2'].real, r['cb2'].imag], occ=r['occ'])
        # T = 6 dE/ddelta by central difference at +pi/2
        h = 1e-4
        Ep = build(cut, math.pi / 2 + h)['E']; Em = build(cut, math.pi / 2 - h)['E']
        dE = (Ep - Em) / (2 * h)
        # K-symmetry: E even, T odd; stationary-state torque
        r = build(cut, math.pi / 2)
        w, v = np.linalg.eigh(r['H'])
        Lop = np.diag(r['L']); comm = r['H'] @ Lop - Lop @ r['H']
        tmax = max(abs((1j * np.vdot(v[:, k], comm @ v[:, k])).real) for k in range(len(w)))
        out[str(cut)] = dict(rows=rows, six_dE_ddelta_at_plus_pi2=6 * dE, T_at_plus_pi2=rows['+pi/2']['T'],
                             exact_T_formula=-12 * 0.04 * 2 * math.sqrt(2), exact_E_amplitude=2 * 0.04 * 2 * math.sqrt(2),
                             E0_minus_Epi_over_2=(rows['0']['E'] - rows['pi']['E']) / 2,
                             E_even_check=rows['+pi/2']['E'] - rows['-pi/2']['E'], T_odd_check=rows['+pi/2']['T'] + rows['-pi/2']['T'],
                             max_eigenstate_torque=float(tmax), H_real=bool(np.max(np.abs(r['H'].imag)) == 0))
        print(cut, 'dim', rows['0']['dim'], 'pproj', rows['0']['pproj'], 'E', [rows[k]['E'] for k in rows], 'T(+pi/2)', rows['+pi/2']['T'],
              '6dE/dd', 6 * dE, '(E0-Epi)/2', out[str(cut)]['E0_minus_Epi_over_2'], '<cb2>(+pi/2)', rows['+pi/2']['cb2'],
              'occ', rows['+pi/2']['occ'], 'max eigenstate torque', tmax)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
