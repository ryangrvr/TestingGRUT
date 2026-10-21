"""Reproducible cross-checks quoted in NOTES section 2 (1D Fokker-Planck):
  * spectral-sum p_- vs the relative-precision delta-start identity (H = 300, 250 GeV);
  * tridiagonal eigenvalue vs MFPT Lambda_1, including where the eigenvalue fails (200, 150 GeV);
  * fresh-start lag N_lag = N - 2 p_-(N)/Lambda_1 at N = 50 (120-300 GeV), and the rejected
    MFPT-difference estimate of the same lag (it neglects fast-mode terms)."""
import numpy as np
from common import Spectator, V, dump
import fp1d


def main():
    rows = []
    for H in (300.0, 250.0, 200.0, 150.0, 120.0):
        sp = Spectator(H)
        yy, h = fp1d.make_grid(sp)
        m0 = int(np.argmin(np.abs(yy - V / H)))
        s, lam, vec = fp1d.odd_spectrum(sp, yy, h)
        L1, _ = fp1d.lambda1_mfpt(sp, yy, h)
        pd = fp1d.p_minus_delta(sp, yy, h, m0, [50.0])[50.0]["p_minus"]
        ps = fp1d.p_minus_from_spectrum(s, lam, vec, m0, 50.0)
        # rejected estimate: lag ~ T_m0 - E_pi[tau], from the birth-death sums (no huge terms)
        s_, up, dn, kill, k = fp1d.odd_chain(sp, yy, h)
        logpi = -s_
        rev = np.logaddexp.accumulate(logpi[::-1])[::-1]
        r = dn.copy()
        r[0] = kill[0]
        term = np.exp(rev - (logpi + np.log(r)))
        cs = np.cumsum(np.where(np.arange(len(s_)) >= 1, term, 0.0))
        pi = np.exp(-(s_ - s_.min()))
        pi /= pi.sum()
        lag_mfpt = float(np.sum(pi * (cs[m0] - cs)))
        rows.append({"H_GeV": H, "Lambda1_eig": float(lam[0]), "Lambda1_mfpt": L1,
                     "p_minus_N50_spectral": float(ps), "p_minus_N50_delta_identity": pd,
                     "N_lag_N50": 50.0 - 2 * pd / L1, "N_lag_mfpt_difference_estimate_REJECTED": lag_mfpt})
    dump("checks_1d.json", {"module": "inflation_spectator", "part": "1D cross-checks", "rows": rows})
    return rows


if __name__ == "__main__":
    for r in main():
        print(r)
