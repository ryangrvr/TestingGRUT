"""INDEPENDENT REPRODUCTION of S2-Sigma load-bearing claims (NR-Sigma-1..7). No S2-Sigma code imported.
Criteria (larger = better), computed with review routines (nr_sigmah.py: Moebius support weights, eigvalsh entropies):
  L = -k_eff (support-size mean, in blocks)          P = weight fraction on single-block supports (autonomy)
  M = -(dim of support-closed operator space)/4^n    G = -(fraction of block pairs that interact)  [NR-Sigma-2 only]
  Q = -mean_f [S_f(tau) - S_f(0)]                    R = max_f #{g != f : I(f:g)(tau) >= (1-delta) S_f(tau)}, S_f > 0.1
Candidates: 202 groupings x frames {id, 3 fresh H/S/CZ Cliffords, Haar, commutant e^{-iH 0.7}} = 1212."""
import itertools as it
import numpy as np
from scipy.linalg import expm
from nr_sigmah import (n, N, X, Z, op, ptrace_keep, support_weights, GROUPS, dtype, entropy_bits, rand_clifford, haar, mfi, prod_state)

SUBS = [tuple(s) for r in range(1, n + 1) for s in it.combinations(range(n), r)]
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

def h_only(w, g):
    blk = {q: i for i, b in enumerate(g) for q in b}
    present = set(); num = den = single = 0.0
    tol = 1e-9 * sum(w.values())          # relative: Moebius inclusion-exclusion leaves ~1e-10 cancellation noise
    for S, v in w.items():
        if v < tol: continue
        B = frozenset(blk[q] for q in S); present.add(B); num += len(B) * v; den += v
        if len(B) == 1: single += v
    closure = 0
    for S in SUBS:
        B = frozenset(blk[q] for q in S)
        if any(B <= P_ for P_ in present): closure += 3 ** len(S)
    pairs = {frozenset(p) for B in present for p in it.combinations(sorted(B), 2)}
    nb = len(g)
    return dict(L=-num / den, P=single / den, M=-closure / 4 ** n, G=-len(pairs) / (nb * (nb - 1) / 2))

def ent_table(rho):
    return {S: entropy_bits(ptrace_keep(rho, list(S))) for S in SUBS}

def state_crit(E0, Et, g, delta=0.1):
    blocks = [tuple(sorted(b)) for b in g]
    Q = -np.mean([Et[b] - E0[b] for b in blocks]); R = 0
    for f in blocks:
        Sf = Et[f]
        if Sf <= 0.1: continue
        cnt = sum(1 for h in blocks if h != f and Sf + Et[h] - Et[tuple(sorted(f + h))] >= (1 - delta) * Sf - 1e-9)
        R = max(R, cnt)
    return Q, R

def candidates(H, rho, frames, tau=1.0, delta=0.1):
    ev, V = np.linalg.eigh(H); Ut = (V * np.exp(-1j * ev * tau)) @ V.conj().T
    rows = []
    for fname, U in frames:
        w = support_weights(U.conj().T @ H @ U)
        r0 = U.conj().T @ rho @ U; rt = U.conj().T @ (Ut @ rho @ Ut.conj().T) @ U
        E0, Et = ent_table(r0), ent_table(rt)
        for g in GROUPS:
            row = h_only(w, g); row["Q"], row["R"] = state_crit(E0, Et, g, delta)
            row.update(frame=fname, g=g, d=dtype(g)); rows.append(row)
    return rows

KEYS5 = ["L", "Q", "R", "P", "M"]
def mat(rows, keys): return np.round(np.array([[r[k] for k in keys] for r in rows], float), 8)
def pareto_idx(V):
    out = []
    for i in range(len(V)):
        dom = ((V >= V[i]).all(1) & (V > V[i]).any(1)).any()
        if not dom: out.append(i)
    return out
def dominant_idx(V): mx = V.max(0); return [i for i in range(len(V)) if (V[i] == mx).all()]
def lab(r): return f"{r['frame']}:{r['d']}:" + "".join("(" + "".join(map(str, b)) + ")" for b in r["g"])

def norm(V, how):
    if how == "minmax": s = V.max(0) - V.min(0); s[s == 0] = 1; return (V - V.min(0)) / s
    if how == "zscore": s = V.std(0); s[s == 0] = 1; return (V - V.mean(0)) / s
    return np.argsort(np.argsort(V, 0), 0) / len(V)
TRANS = {"id": lambda x: x, "log1p": lambda x: np.log1p(x - x.min(0)), "sq": lambda x: (x - x.min(0)) ** 2, "sqrt": lambda x: np.sqrt(x - x.min(0))}

def make_frames(H, rng):
    return [("id", np.eye(N, dtype=complex))] + [(f"cl{i+1}", rand_clifford(rng)) for i in range(3)] + [("haar", haar(rng)), ("comm", expm(-1j * H * 0.7))]

def models(rng):
    gen = lambda d: (lambda A: (A + A.conj().T) / 2)(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    def pair(h4, a): return np.kron(np.kron(np.eye(2 ** a), h4), np.eye(2 ** (n - a - 2)))
    Hm = mfi()
    ring = sum(op({i: Z, (i + 1) % n: Z}) for i in range(n)) + sum(op({i: X}) for i in range(n))
    clus = sum(pair(gen(4), a) for a in (0, 2, 4)) + 0.05 * (op({1: X, 2: X}) + op({3: X, 4: X}))
    Wc = rand_clifford(rng, 30)
    conf = Hm + Wc @ Hm @ Wc.conj().T          # sum of a chain local in id and a chain local in a Clifford frame
    return {"MFI chain (positive control)": Hm, "TFIM ring (translation-symmetric)": ring, "clustered pairs + weak bridges": clus, "conflict: MFI + Clifford-rotated MFI": conf}

def states(H, rng):
    ev, V = np.linalg.eigh(H); pr = prod_state(1.4)
    hv = (lambda v: v / np.linalg.norm(v))(rng.normal(size=N) + 1j * rng.normal(size=N))
    rec = op({0: np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)}) @ np.eye(N, dtype=complex)[:, 0]
    g = (V * np.exp(-(ev - ev.min()))) @ V.conj().T; g /= np.trace(g)
    to_rho = lambda v: np.outer(v, v.conj())
    return {"product": to_rho(pr), "eigenstate": to_rho(V[:, N // 2].astype(complex)), "Haar": to_rho(hv), "record-forming": to_rho(rec), "Gibbs b=1": g}

if __name__ == "__main__":
    rng = np.random.default_rng(4242)
    MOD = models(rng)
    hdr("NR-Sigma-1 / 3 / 5: full criterion vector (L,Q,R,P,M), fresh frames; dominance, objective dependence, state dependence")
    store = {}
    for mname, H in MOD.items():
        frames = make_frames(H, rng); STS = states(H, rng); store[mname] = (H, frames, STS)
        for sname, rho in STS.items():
            rows = candidates(H, rho, frames); V = mat(rows, KEYS5)
            fr = pareto_idx(V); dom = dominant_idx(V)
            # objective conventions
            wins = set()
            for tn, T in TRANS.items():
                for nm in ("minmax", "zscore", "rank"):
                    wins.add(lab(rows[int(np.argmax(norm(T(V), nm).sum(1)))]))
            Wd = rng.dirichlet(np.ones(5), 3000); sw = np.argmax(norm(V, "minmax") @ Wd.T, 0); u, c = np.unique(sw, return_counts=True)
            lexw = set()
            for perm in it.permutations(range(5)):
                idx = np.arange(len(V))
                for k in perm:
                    col = V[idx, k]; idx = idx[col == col.max()]
                lexw.add(frozenset(lab(rows[i]) for i in idx))
            # monotone-invariance of the Pareto set
            inv = all(set(pareto_idx(np.round(T(V), 10))) == set(fr) for T in TRANS.values())
            dts = sorted({rows[i]["d"] for i in fr})
            print(f"[{mname[:30]:<30} | {sname:<14}] front {len(fr):4d} (dtypes {len(dts)}); dominant {'YES' if dom else 'none'}; "
                  f"scalar winners (4 transforms x 3 norms) {len(wins)}; simplex winners {len(u)} (top share {c.max()/c.sum():.2f}); "
                  f"lexicographic (120 orders) distinct winner sets {len(lexw)}; Pareto set monotone-invariant: {inv}")
            store[(mname, sname)] = rows

    hdr("NR-Sigma-2: local-dimension conflict (best dtype per criterion; clustered + MFI, product state)")
    for mname in ("clustered pairs + weak bridges", "MFI chain (positive control)"):
        rows = store[(mname, "product")]
        for k in ("L", "P", "G", "M", "Q", "R"):
            vals = np.round(np.array([r[k] for r in rows]), 8); best = vals.max(); bd = sorted({rows[i]["d"] for i in np.where(vals == best)[0]})
            print(f"  {mname[:30]:<30} criterion {k}: best value {best:+.4f} attained by dtypes {bd[:6]}{' ...' if len(bd) > 6 else ''}")

    hdr("NR-Sigma-4: scale dependence (clustered and MFI; record-forming state): tau, delta, GUE perturbation eps")
    for mname in ("clustered pairs + weak bridges", "MFI chain (positive control)"):
        H0, frames, STS = store[mname]; rho = STS["record-forming"]
        G = (lambda A: (A + A.conj().T) / 2)(rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))); G *= np.linalg.norm(H0) / np.linalg.norm(G)
        for eps in (0.0, 0.05, 0.3):
            for tau, delta in ((0.2, 0.1), (1.0, 0.1), (5.0, 0.1), (1.0, 0.3)):
                rows = candidates(H0 + eps * G, rho, frames, tau, delta); V = mat(rows, KEYS5)
                fr = pareto_idx(V); dom = dominant_idx(V)
                ew = lab(rows[int(np.argmax(norm(V, "minmax").sum(1)))])
                print(f"  {mname[:12]} eps={eps:<4} tau={tau:<3} delta={delta}: front {len(fr):4d}, dominant {'YES' if dom else 'no'}, equal-weight minmax winner {ew}")

    hdr("NR-Sigma-5b: H-only (L,P,M) vs H+psi winners across states, same H (MFI and clustered)")
    for mname in ("MFI chain (positive control)", "clustered pairs + weak bridges"):
        honly = None
        for sname in ("product", "eigenstate", "Haar", "Gibbs b=1"):
            rows = store[(mname, sname)]
            Vh = mat(rows, ["L", "P", "M"]); Vf = mat(rows, KEYS5)
            wh = lab(rows[int(np.argmax(norm(Vh, "minmax").sum(1)))]); wf = lab(rows[int(np.argmax(norm(Vf, "minmax").sum(1)))])
            print(f"  {mname[:12]} {sname:<11}: H-only equal-weight winner {wh:<34} | H+psi equal-weight winner {wf}")

    hdr("NR-Sigma-6: CPR positive control (d=2, n=6 qubit grouping only; H-only criteria; MFI)")
    H = MOD["MFI chain (positive control)"]; frames = store["MFI chain (positive control)"][1] + [("comm1.9", expm(-1j * H * 1.9))]
    qub = [[q] for q in range(n)]
    for fname, U in frames:
        w = support_weights(U.conj().T @ H @ U); v = h_only(w, qub)
        print(f"  frame {fname:<8}: L {v['L']:+.4f}  P {v['P']:.4f}  M {v['M']:+.5f}")
    print("  -> within the supplied (n=6, d=2, k) class the id frame and its commutant images are tied at the top (H-relative gauge); fresh Cliffords / Haar are dominated")

    hdr("NR-Sigma-7: commutant check (time-evolution subgroup) and an exploratory nonlinear-commutant probe (NOT adjudicated)")
    ev, V = np.linalg.eigh(H); psi = prod_state(1.4); rho = np.outer(psi, psi.conj())
    for s in (0.7, 1.9):
        W = expm(-1j * H * s)
        rw = candidates(H, rho, [("W", W)]); rshift = candidates(H, (V * np.exp(1j * ev * s)) @ V.conj().T @ rho @ ((V * np.exp(1j * ev * s)) @ V.conj().T).conj().T, [("id", np.eye(N, dtype=complex))])
        r0 = candidates(H, rho, [("id", np.eye(N, dtype=complex))])
        dh = max(abs(a[k] - b[k]) for a, b in zip(rw, r0) for k in ("L", "P", "M"))
        ds = max(abs(a[k] - b[k]) for a, b in zip(rw, rshift) for k in ("Q", "R"))
        d0 = max(abs(a[k] - b[k]) for a, b in zip(rw, r0) for k in ("Q", "R"))
        print(f"  W = exp(-iH*{s}): |H-only diff vs id| = {dh:.1e}; |state diff vs id at psi(-s)| = {ds:.1e}; |state diff vs id at psi| = {d0:.3f}")
    Wf = (V * np.exp(-1j * ev ** 2 * 0.05)) @ V.conj().T     # nonlinear f(E) = 0.05 E^2
    rwf = candidates(H, rho, [("Wf", Wf)]); dh = max(abs(a[k] - b[k]) for a, b in zip(rwf, r0) for k in ("L", "P", "M"))
    best = 9
    for s in np.arange(-6, 6.01, 0.25):
        Us = (V * np.exp(1j * ev * s)) @ V.conj().T
        rs = candidates(H, Us @ rho @ Us.conj().T, [("id", np.eye(N, dtype=complex))])
        best = min(best, max(abs(a["Q"] - b["Q"]) for a, b in zip(rwf, rs)))
    print(f"  W = exp(-i 0.05 H^2) (nonlinear commutant): |H-only diff vs id| = {dh:.1e} (gauge); min over time shifts s in [-6,6] of max|Q diff| = {best:.3f}"
          " -> not a time shift in this instance (exploratory; full-commutant question remains NOT adjudicated)")
