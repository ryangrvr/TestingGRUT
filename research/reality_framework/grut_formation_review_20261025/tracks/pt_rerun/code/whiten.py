"""Whitened-covariance test of a generated Fourier realization against C_k = F F^T / L^3 (own code)."""
import numpy as np
from scipy.fft import fftn
from scipy.stats import chi2 as CHI2


def whiten_test(z, d, wk, L, N):
    F = d[wk]; lookup = {int(x): i for i, x in enumerate(d['keys'])}
    c = fftn(z, axes=(-3, -2, -1), norm='forward')
    m = np.fft.fftfreq(N) * N; I, Jj, K = np.meshgrid(m, m, m, indexing='ij')
    first = np.where(I != 0, I, np.where(Jj != 0, Jj, K))
    key = (I ** 2 + Jj ** 2 + K ** 2).astype(int)
    half = (first > 0)
    stats = []; allT = 0.; alldof = 0; pseudo = []; outside = 0.
    kk = 2 * np.pi / L * np.sqrt(key)
    outside = float(np.max(abs(c[:, half & (kk > 8 + 1e-9)]))) if np.any(half & (kk > 8 + 1e-9)) else 0.
    diag_all = []
    for kv, i in lookup.items():
        sel = half & (key == kv)
        if not sel.any(): continue
        Ck = F[i] @ F[i].T / L ** 3
        Lc = np.linalg.cholesky(Ck)
        w = np.linalg.solve(Lc, c[:, sel])  # 6 x n
        n = w.shape[1]; Sm = (w @ w.conj().T) / n; Ps = (w @ w.T) / n
        diag_all.append(np.real(np.diag(Sm)))
        T = 2 * np.sum(abs(w) ** 2); allT += T; alldof += 12 * n
        stats.append(dict(k=float(2 * np.pi / L * np.sqrt(kv)), n=int(n), diag=np.real(np.diag(Sm)).tolist(),
                          offdiag_max=float(np.max(abs(Sm - np.diag(np.diag(Sm))))), pseudo_max=float(np.max(abs(Ps))),
                          chi2_p=float(CHI2.sf(T, 12 * n)) if T > 12 * n else float(CHI2.cdf(T, 12 * n))))
    D = np.concatenate(diag_all)
    # shells with >= 50 modes for a scatter statement
    big = [s for s in stats if s['n'] >= 50]
    Dbig = np.concatenate([s['diag'] for s in big])
    zero_mode = c[:, 0, 0, 0].real
    mean = np.r_[d['light_mean_q'], d['light_mean_p']]
    return dict(shells=len(stats), modes=int(sum(s['n'] for s in stats)), pooled_mean_abs2=float(allT / alldof),
                pooled_chi2_dof=alldof, pooled_chi2_two_sided_p=float(2 * min(CHI2.sf(allT, alldof), CHI2.cdf(allT, alldof))),
                diag_range_all=[float(D.min()), float(D.max())], diag_std_all=float(D.std()),
                diag_range_shells_ge50=[float(Dbig.min()), float(Dbig.max())], diag_std_shells_ge50=float(Dbig.std()),
                expected_std_shells_ge50=float(np.mean([1 / np.sqrt(s['n']) for s in big])),
                min_shell_p=float(min(s['chi2_p'] for s in stats)), n_shells_p_lt_0p01=int(sum(s['chi2_p'] < .01 for s in stats)),
                max_offdiag_shells_ge50=float(max(s['offdiag_max'] for s in big)), max_pseudo_shells_ge50=float(max(s['pseudo_max'] for s in big)),
                max_coeff_beyond_k8=outside, zero_mode_minus_mean_max=float(np.max(abs(zero_mode - mean))),
                low_k_shells=[dict(k=s['k'], n=s['n'], diag=[round(x, 3) for x in s['diag']]) for s in stats if s['k'] < 1.8])

