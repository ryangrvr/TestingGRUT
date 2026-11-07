"""GRUT detector histories in lattice QED (Schwinger model).

Instrument (GRUT P2 form): K_pm = (1 pm e^{i phi} W)/2, W = exp(i c F),
F = smeared electric field (Gaussian in physical units, centred, unit weight).
Two-instrument history at times 0 and T in the prepared vacuum.

Compare exact lattice probabilities with the quasi-free (two-point-only)
prediction built from the lattice's own <dF(t) dF(0)>. In the continuum
massless Schwinger model E = (g/sqrt(pi)) phi is a free field, so the
two-point summary is sufficient there (GRUT eq. 9-11 language: the retained
span is the Gaussian one).
"""
import itertools, json, sys
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
from schwinger import hamiltonian


def run(N, x, m, c=2.0, sigma_g=0.8, phi=0.0, Tmax=6.0, nT=31):
    H, L, st = hamiltonian(N, x, m)
    Hg = H / (2 * np.sqrt(x))                       # H/g, time in 1/g
    w, v = eigsh(Hg, k=1, which="SA", tol=1e-12)
    E0, om = w[0], v[:, 0].astype(complex)
    nl = N - 1
    pos = (np.arange(nl) - (nl - 1) / 2) / np.sqrt(x)   # link centres, units 1/g
    f = np.exp(-pos ** 2 / (2 * sigma_g ** 2)); f /= f.sum()
    F = L @ f                                      # diagonal smeared field
    Fm = float(np.real(np.vdot(om, F * om)))
    dF = F - Fm
    Ts = np.linspace(0, Tmax, nT)
    # two-point function G(T) = <dF e^{-i(H-E0)T} dF>
    u = dF * om
    evo = expm_multiply(-1j * Hg, u, start=0, stop=Tmax, num=nT, endpoint=True)
    G = np.array([np.vdot(u, e) * np.exp(1j * E0 * t) for e, t in zip(evo, Ts)])
    W = np.exp(1j * c * F)
    exact = {}
    for s1 in (1, -1):
        psi = 0.5 * (om + s1 * np.exp(1j * phi) * W * om)
        evo = expm_multiply(-1j * Hg, psi, start=0, stop=Tmax, num=nT, endpoint=True)
        for s2 in (1, -1):
            exact[(s1, s2)] = np.array([np.linalg.norm(0.5 * (e + s2 * np.exp(1j * phi) * W * e)) ** 2 for e in evo])
    # no first instrument: marginal of the second
    p2_alone = {s2: float(np.linalg.norm(0.5 * (om + s2 * np.exp(1j * phi) * W * om)) ** 2) for s2 in (1, -1)}

    # quasi-free prediction
    def qf(T_idx, s1, s2):
        g0, gT = G[0], G[T_idx]
        two = {(0, 0): g0, (1, 1): g0, (1, 0): gT, (0, 1): np.conj(gT)}   # <dF(t_a) dF(t_b)>, 1 = time T
        tot = 0
        A = [(1, None), (s1, (0, -1))]                  # 1 or s1 U1^dagger
        B = [(2, None), (s2, (1, +1)), (s2, (1, -1))]   # 2, s2 U2, s2 U2^dagger
        C = [(1, None), (s1, (0, +1))]
        for a, b, cc in itertools.product(A, B, C):
            coef = a[0] * b[0] * cc[0]
            ops = [o[1] for o in (a, b, cc) if o[1] is not None]
            val = 1.0 + 0j
            for (t, sg) in ops:
                val *= np.exp(1j * sg * (phi + c * Fm))
            ex = 0
            for i, (ti, si) in enumerate(ops):
                ex += -0.5 * (c * si) ** 2 * two[(ti, ti)]
                for (tj, sj) in ops[i + 1:]:
                    ex += -(c * si) * (c * sj) * two[(ti, tj)]
            tot += coef * val * np.exp(ex)
        return float(np.real(tot)) / 16
    dev = 0.0
    disturb = 0.0
    for k in range(nT):
        for s1 in (1, -1):
            for s2 in (1, -1):
                dev = max(dev, abs(exact[(s1, s2)][k] - qf(k, s1, s2)))
        for s2 in (1, -1):
            disturb = max(disturb, abs(exact[(1, s2)][k] + exact[(-1, s2)][k] - p2_alone[s2]))
    norm_err = max(abs(sum(exact[(a, b)][k] for a in (1, -1) for b in (1, -1)) - 1) for k in range(nT))
    return {"N": N, "x": x, "m_over_g": m, "L_g": N / np.sqrt(x), "links_in_2sigma": int(np.sum(np.abs(pos) < 2 * sigma_g)),
            "var_cF": float(c ** 2 * np.real(G[0])), "max_dev_exact_vs_twopoint": dev,
            "max_disturbance_by_first_instrument": disturb, "probability_norm_error": norm_err,
            "p_pp_series": [round(float(v), 6) for v in exact[(1, 1)]], "T": Ts.tolist()}


if __name__ == "__main__":
    out = []
    for (N, x) in [(10, 4.0), (14, 9.0), (20, 16.0)]:
        for m in (0.0, 0.5):
            r = run(N, x, m)
            print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k not in ("p_pp_series", "T")}, flush=True)
            out.append(r)
    json.dump(out, open("histories_results.json", "w"), indent=1)
