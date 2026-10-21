"""Parts 2+4 combined: opposite-vacuum fraction AFTER the post-inflation roll/oscillation, fresh start at +v.

End-of-inflation distribution P_end(phi) from the full (two-sided) SqRA Fokker-Planck chain, delta start at
phi = +v, N_* e-folds. Each end value phi_end is then evolved in an instantly starting radiation era
(post_inflation.oscillate: phi frozen at the end of inflation, Hubble friction included) to its trapped
vacuum. p_final = P(trapped in -v). Below |phi_end| ~ 2.5 v no end value crosses zero (checked), so the
correction only matters for H_inf >~ 500 GeV, where p_- is already > 1e-3.
Outside-action inputs: instant RD after inflation, phi_dot(end) = 0, no preheating/restoration (see NOTES).
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from common import Spectator, V, DIFF, PC_SITE_SC, dump
import fp1d
from post_inflation import oscillate


def full_chain_end(sp, h, ymax, y0, N):
    M = int(np.ceil(ymax / h))
    y = (np.arange(-M, M) + 0.5) * h
    s = sp.s(y)
    k = DIFF / h ** 2
    up = np.zeros_like(s)
    dn = np.zeros_like(s)
    up[:-1] = k * np.exp(-(s[1:] - s[:-1]) / 2)
    dn[1:] = k * np.exp(-(s[:-1] - s[1:]) / 2)
    # forward equation dp/dN = Q^T p, Q_{i,i+1} = up_i, Q_{i,i-1} = dn_i
    QT = diags([up[:-1], -(up + dn), dn[1:]], [-1, 0, 1], format="csr")
    p0 = np.zeros_like(s)
    p0[int(np.argmin(np.abs(y - y0)))] = 1.0
    p = expm_multiply(QT * N, p0)
    return y, np.clip(p, 0, None)


def sign_map(y, vH, ythr):
    out = np.sign(y).astype(float)
    big = np.abs(y) >= ythr
    cache = {}
    for i in np.where(big)[0]:
        a = round(abs(y[i]), 6)
        if a not in cache:
            cache[a] = oscillate(a, vH)[1]
        out[i] = np.sign(y[i]) * cache[a]
    return out


def main():
    rows = []
    for H in (500.0, 700.0, 1000.0, 1400.0, 2000.0, 3000.0, 5000.0, 10000.0):
        sp = Spectator(H)
        yy, h = fp1d.make_grid(sp, smax=40.0)
        vH = V / H
        row = {"H_GeV": H}
        for N in (50.0, 60.0):
            y, p = full_chain_end(sp, h, yy[-1], vH, N)
            sf = sign_map(y, vH, 2.4 * vH)
            p_end_neg = float(p[y < 0].sum())
            p_final = float(p[sf < 0].sum())
            row["N%d" % N] = {"p_minus_end_of_inflation": p_end_neg, "p_minus_after_trapping": p_final,
                              "prob_phi_end_above_2.4v": float(p[y > 2.4 * vH].sum())}
        rows.append(row)
        print(row, flush=True)
    out = {"module": "inflation_spectator", "part": "2+4 final vacuum fraction after post-inflation trapping",
           "label": "computed (1D FP end distribution x RD trapping map); instant RD, no preheating/restoration",
           "rows": rows}
    for N in ("N50", "N60"):
        Hs = np.array([r["H_GeV"] for r in rows])
        pf = np.array([r[N]["p_minus_after_trapping"] for r in rows])
        i = np.where((pf[:-1] < PC_SITE_SC) & (pf[1:] >= PC_SITE_SC))[0]
        out["H_perc_after_trapping_" + N] = (float(10 ** np.interp(np.log10(PC_SITE_SC), np.log10(pf[i[0]:i[0] + 2]),
                                                                   np.log10(Hs[i[0]:i[0] + 2]))) if i.size else None)
    dump("final_sign.json", out)
    return out


if __name__ == "__main__":
    o = main()
    print({k: v for k, v in o.items() if k.startswith("H_perc")})
