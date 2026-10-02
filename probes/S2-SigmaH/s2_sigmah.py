"""SCOUT-2 S2-SigmaH: joint factorization / boundary-condition selector.  Pre-registered (PROBE_CHARTERS, S2-SigmaH).
Reuses the S2-Sigma candidate machinery (6 qubit slots, 202 groupings) with frames:
id, 3 Cliffords, Haar, commutant e^{-iH 0.7}, and W_psi (Householder: psi(0) -> |0..0>)."""
import itertools as it
import sys
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
sys.path.insert(0, "../S2-Sigma")
import s2_sigma as SG                                    # n = 6, N = 64, GROUPINGS, pauli_coeffs, entropies, ...

n, N = SG.n, SG.N
rng = np.random.default_rng(1729)
X, Z, Yc, I2 = SG.X, SG.Z, SG.Y, SG.I2
op = SG.op
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
SUB = {s: i for i, s in enumerate(SG.SUBSETS)}
TGRID = np.arange(-10, 10.0001, 0.25); T0 = int(np.argmin(abs(TGRID)))

def mfi(hx=0.9045, hz=0.809):
    return sum(op({i: Z, i + 1: Z}) for i in range(n - 1)) + hx * sum(op({i: X}) for i in range(n)) + hz * sum(op({i: Z}) for i in range(n))
def prod_state(angles):
    out = np.ones(1, complex)
    for a in angles: out = np.kron(out, np.array([np.cos(a / 2), np.sin(a / 2)]))
    return out
def householder(psi):
    e = np.zeros(N, complex); e[0] = 1; ph = psi[0] / abs(psi[0]) if abs(psi[0]) > 1e-12 else 1
    v = psi - ph * e; nv = np.linalg.norm(v)
    if nv < 1e-12: return np.eye(N, dtype=complex)
    v /= nv; return np.eye(N) - 2 * np.outer(v, v.conj())
def kmax_of(c, g):
    T = np.stack([SG.NONID[:, b].any(1) for b in g], 1); size = T.sum(1)
    return int(size[(np.abs(c) > 1e-9) & (size > 0)].max())
def h_mdl(c, g, delta):
    w = c.copy(); w[0] = 0
    T = np.stack([SG.NONID[:, b].any(1) for b in g], 1)
    code = T.astype(np.int64) @ (1 << np.arange(len(g)))
    present = np.unique(code[(np.abs(w) > delta) & (code != 0)]); cov = np.zeros(len(code), bool)
    for S in present: cov |= (code & ~S) == 0
    return int((cov & (code != 0)).sum())

def state_tables(state, H, frames):
    """entropy tables E[frame][t] over the 63 subsets, plus global entropy (0 for pure)."""
    ev, V = np.linalg.eigh(H); out = {}
    pure = state.ndim == 1
    for fname, U in frames:
        rows = []
        for t in TGRID:
            if pure:
                st = V @ (np.exp(-1j * ev * t) * (V.conj().T @ state)); st = U.conj().T @ st
            else:
                Ut = (V * np.exp(-1j * ev * t)) @ V.conj().T; st = U.conj().T @ (Ut @ state @ Ut.conj().T) @ U
            e = SG.entropies(st); rows.append(np.array([e[s] for s in SG.SUBSETS]))
            if not pure and t == TGRID[0]:
                rows = [rows[0]] * len(TGRID); break      # Gibbs: stationary
        out[fname] = np.array(rows)
    Sg = 0.0 if pure else SG.ent(np.linalg.eigvalsh(state)) if hasattr(SG, "ent") else float(-(lambda p: (p[p > 1e-14] * np.log2(p[p > 1e-14])).sum())(np.linalg.eigvalsh(state)))
    return out, Sg

def candidates(H, state, frames, delta=1e-9):
    tabs, Sg = state_tables(state, H, frames)
    rows, meta = [], []
    for fname, U in frames:
        c = SG.pauli_coeffs(U.conj().T @ H @ U); E = tabs[fname]
        for g in SG.GROUPINGS:
            L, P, M = SG.h_criteria(c, g)
            Ct = sum(E[:, SUB[tuple(b)]] for b in g) - Sg
            imin = int(np.argmin(Ct))
            rows.append(dict(L=L, P=P, M=M, k=kmax_of(c, g), C0=Ct[T0], Cmin=Ct.min(), Cavg=Ct.mean(), Ctyp=float(np.median(Ct)),
                             tstar=TGRID[imin], Ct=Ct, E=E, c=c, g=g))
            meta.append((fname, SG.dtype_of(g), g))
    return rows, meta

def vec(rows, keys): return np.round(np.array([[r[k] for k in keys] for r in rows]), 9)
def front_report(rows, meta, Ckey):
    V = vec(rows, ["L"]) ; V = np.c_[V, -vec(rows, [Ckey])]
    fr = SG.pareto(V); dom = SG.dominant(V)
    frames = sorted({meta[i][0] for i in fr}); dts = sorted({meta[i][1] for i in fr})
    distinct = len({tuple(V[i]) for i in fr})
    return fr, dom, frames, dts, distinct

def lex(rows, order):
    idx = np.arange(len(rows))
    for key, sense in order:
        vals = np.round(np.array([rows[i][key] for i in idx]) * (1 if sense == "max" else -1), 9)
        idx = idx[vals == vals.max()]
    return idx
RULES = {"A (max L; min C)": [("L", "max"), ("C", "min")], "B (min C; max L)": [("C", "min"), ("L", "max")],
         "C (min k; min C; min MDL)": [("k", "min"), ("C", "min"), ("M", "max")], "D (min C; min k; max autonomy)": [("C", "min"), ("k", "min"), ("P", "max")]}
# note: M is stored as -dim/N^2 in S2-Sigma (larger = shorter description), so 'max M' = min description length

def analyse(name, H, state, extra_frames=(), show=True, Ckeys=("C0", "Cmin", "Cavg", "Ctyp")):
    psi0 = state if state.ndim == 1 else np.linalg.eigh(state)[1][:, -1]
    frames = [("id", np.eye(N)), ("cliff1", SG.CLIFFS[0]), ("cliff2", SG.CLIFFS[1]), ("cliff3", SG.CLIFFS[2]), ("haar", SG.HAAR),
              ("commutant", expm(-1j * H * 0.7)), ("W_psi", householder(psi0))] + list(extra_frames)
    rows, meta = candidates(H, state, frames)
    res = {"rows": rows, "meta": meta}
    if show: print(f"\n[{name}]  candidates {len(rows)}")
    for Ck in Ckeys:
        fr, dom, frs, dts, distinct = front_report(rows, meta, Ck)
        res[Ck] = (fr, dom)
        if show:
            print(f"   Pareto(L, -{Ck:<4}): front {len(fr):4d} (distinct vectors {distinct:3d}); dominant: "
                  f"{SG.label(meta, dom[0]) + (f' (+{len(dom)-1} ties)' if len(dom) > 1 else '') if len(dom) else 'NONE'}; frames {frs}; dtypes {dts}")
    return res

def lex_eps(res, Ck, show=True):
    rows, meta = res["rows"], res["meta"]
    for r in rows: r["C"] = r[Ck]
    out = {}
    for nm, order in RULES.items():
        w = lex(rows, order); out[nm] = w
        if show: print(f"   lexicographic {nm:<32} [{Ck}]: {len(w)} winner(s): " + "; ".join(SG.label(meta, i) for i in w[:3]) + (" ..." if len(w) > 3 else ""))
    Ls = np.round(np.array([r["L"] for r in rows]), 9); Cs = np.round(np.array([r["C"] for r in rows]), 9); wins = []
    for q in (0.5, 0.75, 0.9, 0.99):
        ok = np.where(Ls >= np.quantile(Ls, q) - 1e-9)[0]; best = Cs[ok].min()
        wins.append(frozenset(SG.label(meta, i) for i in ok[Cs[ok] == best]))
    for q in (0.01, 0.1, 0.25, 0.5):
        ok = np.where(Cs <= np.quantile(Cs, q) + 1e-9)[0]; best = Ls[ok].max()
        wins.append(frozenset(SG.label(meta, i) for i in ok[Ls[ok] == best]))
    distinct = sorted(set(wins), key=len)
    labs = distinct
    if show: print(f"   epsilon-constraint sweep [{Ck}] (8 tolerances, full tie sets): {len(distinct)} distinct winner SETS, sizes {[len(w) for w in distinct]};"
                   f" common to all: {len(frozenset.intersection(*wins))}; e.g. " + "; ".join(sorted(distinct[0])[:2]))
    return out, labs

if __name__ == "__main__":
    H = mfi(); ev, Vh = np.linalg.eigh(H)
    p_tilt = prod_state(np.full(n, 1.4))
    def energy(a): s = prod_state(a); return float(np.real(np.vdot(s, H @ s)))
    mf = prod_state(min((minimize(energy, rng.uniform(0, 2 * np.pi, n), method="BFGS") for _ in range(6)), key=lambda r: r.fun).x)
    gibbs = (Vh * np.exp(-(ev - ev.min()))) @ Vh.T; gibbs = gibbs / np.trace(gibbs)
    haar = (lambda v: v / np.linalg.norm(v))(rng.normal(size=N) + 1j * rng.normal(size=N))
    ring = sum(op({i: Z, (i + 1) % n: Z}) for i in range(n)) + sum(op({i: X}) for i in range(n))
    gk = rng.uniform(0.8, 1.2, n - 1)
    star = sum(gk[k - 1] * op({0: Z, k: Yc}) for k in range(1, n)) + sum(0.1 * op({q: Z}) for q in range(1, n))
    rec = op({0: np.array([[1, 1], [1, -1]]) / np.sqrt(2)}) @ np.eye(N)[:, 0].astype(complex)
    CASES = [("PC  positive control: MFI chain + product psi (same frame)", H, p_tilt),
             ("H1  conflict: MFI local in id, psi product in cliff1 frame", H, SG.CLIFFS[0] @ p_tilt),
             ("H3  translation ring (TFIM) + translation-invariant product psi", ring, p_tilt),
             ("H4  record-forming star + |+>|0..>", star, rec),
             ("H5a MFI mid-spectrum eigenstate", H, Vh[:, N // 2].astype(complex)),
             ("H5b MFI Gibbs beta = 1", H, gibbs.astype(complex)),
             ("H6  MFI Janus mean-field product state", H, mf),
             ("H7  MFI Haar state", H, haar)]
    hdr("SigmaH-1/2/5: Pareto fronts over (L_H, -C_psi) for four epoch prescriptions; dominance; local dims (H2 = W_psi frame present in every run)")
    RES = {}
    for nm, Hm, st in CASES:
        RES[nm] = analyse(nm, Hm, st)
    hdr("SigmaH-3/4: lexicographic priority rules and epsilon-constraint sweeps (C0 = supplied epoch; Cmin = epoch-free)")
    for nm, Hm, st in CASES:
        print(f"\n[{nm}]")
        for Ck in ("C0", "Cmin"):
            lex_eps(RES[nm], Ck)

    hdr("SigmaH-6: time translation psi -> e^{-iHs} psi, s = 2.25 (grid-aligned)")
    for nm, Hm, st in CASES[:2] + CASES[6:7]:
        s = 2.25; ev2, V2 = np.linalg.eigh(Hm); sts = V2 @ (np.exp(-1j * ev2 * s) * (V2.conj().T @ st))
        R2 = analyse(nm + " [shifted]", Hm, sts, show=False, Ckeys=("C0", "Cmin"))
        for Ck in ("C0", "Cmin"):
            a = RES[nm][Ck][1]; b = R2[Ck][1]
            la = {SG.label(RES[nm]["meta"], i).replace("W_psi", "W") for i in a}; lb = {SG.label(R2["meta"], i).replace("W_psi", "W") for i in b}
            ta = sorted({RES[nm]["rows"][i]["tstar"] for i in a}); tb = sorted({R2["rows"][i]["tstar"] for i in b})
            print(f"   {nm[:3]} {Ck:<4}: dominant set before {sorted(la)[:2] or 'NONE'} | after {sorted(lb)[:2] or 'NONE'} | same Sigma class: {la == lb and len(la) > 0}"
                  + (f" | t* before {ta[:3]} after {tb[:3]}" if Ck == "Cmin" else ""))

    hdr("SigmaH-7: time reversal (Theta = K, real H): psi vs conj(psi)")
    for nm, Hm, st in [CASES[0], CASES[6]]:
        Rr = analyse(nm + " [reversed]", Hm, np.conj(st), show=False, Ckeys=("Cmin",))
        a, b = RES[nm]["Cmin"][1], Rr["Cmin"][1]
        print(f"   {nm[:3]}: Cmin-dominant before {[SG.label(RES[nm]['meta'], i) for i in a][:2] or 'NONE'} t*={sorted({RES[nm]['rows'][i]['tstar'] for i in a})[:3]};"
              f" reversed {[SG.label(Rr['meta'], i) for i in b][:2] or 'NONE'} t*={sorted({Rr['rows'][i]['tstar'] for i in b})[:3]}")
        idr = [i for i, m in enumerate(RES[nm]["meta"]) if m[0] == "id" and len(m[2]) == n][0]
        print(f"        id-frame qubit TPS: C(t) at t = -4,-2,0,2,4: " + " ".join(f"{RES[nm]['rows'][idr]['Ct'][T0 + int(t/0.25)]:.3f}" for t in (-4, -2, 0, 2, 4))
              + "  -> symmetric about the special epoch: orientation not selected")

    hdr("SigmaH-8: add arrow quality R (forward relaxation from the Cmin epoch) as a third criterion")
    for nm, Hm, st in CASES:
        rows, meta = RES[nm]["rows"], RES[nm]["meta"]
        for h in (2, 5):
            for r in rows:
                i0 = int(round((r["tstar"] + 10) / 0.25)); seg = r["Ct"][i0: i0 + int(h / 0.25) + 1]
                r[f"R{h}"] = float(np.mean(np.diff(seg) >= -1e-9)) if len(seg) > 1 else 0.0
            V = np.c_[vec(rows, ["L"]), -vec(rows, ["Cmin"]), vec(rows, [f"R{h}"])]
            fr = SG.pareto(V); dom = SG.dominant(V)
            print(f"   {nm[:3]} horizon h={h}: Pareto(L, -Cmin, R) front {len(fr):4d}; dominant {'NONE' if not len(dom) else SG.label(meta, dom[0]) + f' (+{len(dom)-1})'}")

    hdr("SigmaH-9: records under ONE fragment rule (each factor = system; each other single factor = fragment) at t* + 2")
    for nm, Hm, st in [CASES[0], CASES[3], CASES[6]]:
        rows, meta = RES[nm]["rows"], RES[nm]["meta"]; recs = []
        for r in rows:
            i1 = min(len(TGRID) - 1, int(round((r["tstar"] + 2 + 10) / 0.25))); E = r["E"][i1]; best = 0
            for f in r["g"]:
                Sf = E[SUB[tuple(f)]]
                if Sf <= 0.1: continue
                cnt = sum(1 for h_ in r["g"] if h_ != f and Sf + E[SUB[tuple(h_)]] - E[SUB[tuple(sorted(f + h_))]] >= 0.9 * Sf)
                best = max(best, cnt)
            r["Rrec"] = best; recs.append(best)
        recs = np.array(recs); nf = np.array([len(m[2]) for m in meta])
        V = np.c_[vec(rows, ["L"]), -vec(rows, ["Cmin"]), recs]; fr = SG.pareto(V)
        print(f"   {nm[:3]}: max record count {recs.max()} attained by {np.sum(recs == recs.max())} candidates; mean record count by #factors: "
              + ", ".join(f"{k}:{recs[nf == k].mean():.2f}" for k in sorted(set(nf))) + f";  Pareto(L,-Cmin,Rrec) front {len(fr)}")

    hdr("SigmaH-10: MDL = #H support-closure params(delta) + c_psi * #amplitudes(delta) in a declared language")
    def amp_count(st, U, g, lang, delta):
        v = U.conj().T @ (st if st.ndim == 1 else np.linalg.eigh(st)[1][:, -1])
        if lang == "computational": w = v
        elif lang == "local Hadamard": w = op({q: np.array([[1, 1], [1, -1]]) / np.sqrt(2) for q in range(n)}) @ v
        else:   # product-aware: rotate each factor into the eigenbasis of its reduced state
            t = v.reshape([2] * n)
            for b in g:
                rest = [q for q in range(n) if q not in b]
                m = np.transpose(t, list(b) + rest).reshape(2 ** len(b), -1); rho = m @ m.conj().T
                _, Ub = np.linalg.eigh(rho); m = Ub.conj().T @ m
                t = np.transpose(m.reshape([2] * len(b) + [2] * len(rest)), np.argsort(list(b) + rest))
            w = t.reshape(-1)
        return int((np.abs(w) > delta).sum())
    for nm, Hm, st in [CASES[0], CASES[1], CASES[6], CASES[7]]:
        rows, meta = RES[nm]["rows"], RES[nm]["meta"]
        frames = dict([("id", np.eye(N)), ("cliff1", SG.CLIFFS[0]), ("cliff2", SG.CLIFFS[1]), ("cliff3", SG.CLIFFS[2]), ("haar", SG.HAAR),
                       ("commutant", expm(-1j * Hm * 0.7)), ("W_psi", householder(st if st.ndim == 1 else np.linalg.eigh(st)[1][:, -1]))])
        winners = set()
        sub = list(range(0, len(rows), 1))
        for lang in ("computational", "local Hadamard", "product-aware"):
            for delta in (1e-3, 1e-6, 1e-9):
                hc = np.array([h_mdl(rows[i]["c"], rows[i]["g"], delta) for i in sub])
                ac = np.array([amp_count(st, frames[meta[i][0]], rows[i]["g"], lang, delta) for i in sub])
                for cpsi in (0.25, 1, 4, 64, 1024):
                    tot = hc + cpsi * ac; i = sub[int(np.argmin(tot))]; winners.add((lang, delta, cpsi, SG.label(meta, i)))
        labs = sorted({w[3] for w in winners})
        print(f"   {nm[:3]}: 45 language/precision/cost settings -> {len(labs)} distinct MDL winners: " + "; ".join(labs[:5]) + (" ..." if len(labs) > 5 else ""))

    hdr("SigmaH-11: symmetry quotient on the translation ring (Sym(H,psi) contains translations)")
    rows, meta = RES[CASES[2][0]]["rows"], RES[CASES[2][0]]["meta"]
    fr, _ = RES[CASES[2][0]]["Cmin"]
    def translate(g, s): return sorted(sorted(((q + s) % n) for q in b) for b in g)
    def refl(g): return sorted(sorted(((-q) % n) for q in b) for b in g)
    cls = {}
    for i in fr:
        f, d, g = meta[i]
        key = (f, min(str(h_) for h_ in [translate(g, s) for s in range(n)] + [translate(refl(g), s) for s in range(n)])) if f == "id" else (f, str(g))
        cls.setdefault(key, []).append(i)
    print(f"   Cmin front: {len(fr)} bare TPS candidates -> {len(cls)} classes after quotienting id-frame groupings by translations+reflection")

    hdr("SigmaH-11b: quotient every Cmin Pareto front by local-unitary (LU) equivalence of frames for the same grouping")
    def lu_equiv(U1, U2, g):
        W = U1.conj().T @ U2
        for b in g:
            rest = [q for q in range(n) if q not in b]
            t = W.reshape([2] * (2 * n))
            m = np.transpose(t, list(b) + [x + n for x in b] + rest + [x + n for x in rest]).reshape(4 ** len(b), -1)
            if np.linalg.matrix_rank(m, tol=1e-8) != 1: return False
        return True
    for nm, Hm, st in CASES:
        rows, meta = RES[nm]["rows"], RES[nm]["meta"]; fr, dom = RES[nm]["Cmin"]
        psi0 = st if st.ndim == 1 else np.linalg.eigh(st)[1][:, -1]
        FR = dict([("id", np.eye(N)), ("cliff1", SG.CLIFFS[0]), ("cliff2", SG.CLIFFS[1]), ("cliff3", SG.CLIFFS[2]), ("haar", SG.HAAR),
                   ("commutant", expm(-1j * Hm * 0.7)), ("W_psi", householder(psi0))])
        classes = []
        for i in fr:
            f, d, g = meta[i]
            for cl in classes:
                j = cl[0]
                if meta[j][2] == g and lu_equiv(FR[meta[j][0]], FR[f], g): cl.append(i); break
            else: classes.append([i])
        reps = [SG.label(meta, cl[0]) for cl in classes]
        fr_frames = sorted({meta[cl[0]][0] for cl in classes}); fr_d = sorted({meta[cl[0]][1] for cl in classes})
        print(f"   {nm[:3]}: Cmin front {len(fr)} bare candidates -> {len(classes)} LU classes; frames {fr_frames}; dtypes {fr_d}")
