#!/usr/bin/env python3
"""
V0-1 blind reproduction of P-17 acceptance tests T1-T8 (exact rational Taylor checks).

Independent design (written from the spec's model definition only):

  * System A (Hamiltonian, full phase space): power-series recursion for
        q' = p,  p' = -V'(q) + sum_j c_j (x_j - c_j q / w_j^2),
        x_j' = p_j,  p_j' = -w_j^2 x_j + c_j q,
    with q(0)=a, p(0)=0, x_j(0) = c_j a / w_j^2 + xi_j, p_j(0) = eta_j   (shifted preps)
    or x_j(0) = X_j, p_j(0) = eta_j                                       (product prep).
    Coefficients are exact polynomials in QQ[a, beta, xi_j, eta_j] (sympy sparse rings).
    E[.] is then taken monomial-by-monomial with exact moments (Gaussian or random-phase).

  * GLE comparators B, B_R, D, B+slip: power-series recursion DIRECTLY for
        q'' = -V'(q) - int_0^t gamma(t-s) q'(s) ds + F(t) + S(t),
    using int_0^t (t-s)^k s^l ds = k! l!/(k+l+1)! t^(k+l+1).  No bath variables are used.
      - B   : F symbolic through its derivatives f_k = F^(k)(0); E[.] by Wick/Isserlis with
              Cov(f_k, f_l) = T * (-1)^l * gamma^(k+l)(0)  (stationary covariance T*gamma).
      - B_R : F(t) = sum_j c_j sqrt(2T)/w_j cos(w_j t + phi_j), phi_j iid uniform, exogenous.
      - D   : F = 0.
      - slip: S(t) = -gamma(t) a (deterministic).

Statistic: E[q^(n)(0)] = n! E[q_n] as polynomial in (beta, T, a).
Usage: python3 v0_1_taylor_checks.py <nmodes: 1|2> <Nmax>
"""
import sys, time
from functools import lru_cache
from math import factorial
from sympy import QQ, Rational, symbols, factor, expand, sympify
from sympy.polys.rings import ring

NM = int(sys.argv[1]) if len(sys.argv) > 1 else 2
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else (8 if NM == 2 else 12)

K = QQ(23, 10)                       # V = K q^2/2 + beta q^4  ->  V' = K q + 4 beta q^3
ALL_MODES = [(QQ(1), QQ(1, 2)), (QQ(2), QQ(3, 4))]   # (omega_j, c_j)
MODES = ALL_MODES[:NM]

# Result ring QQ[a, beta, T]
Rres, ra, rb, rT = ring("a,beta,T", QQ)
Ssym_a, Ssym_b, Ssym_T = symbols("a beta T")


def to_expr(p):
    return factor(p.as_expr(Ssym_a, Ssym_b, Ssym_T))


def dfact2(m):  # (m-1)!! for even m, with (-1)!! = 1
    r = 1
    for k in range(m - 1, 0, -2):
        r *= k
    return r


def cs_moment(i, k):
    """E[cos^i(phi) sin^k(phi)], phi ~ U(0,2pi): (i-1)!!(k-1)!!/(i+k)!! if i,k even, else 0."""
    if i % 2 or k % 2:
        return QQ(0)
    num = dfact2(i) * dfact2(k)
    den = 1
    for m in range(i + k, 0, -2):
        den *= m
    return QQ(num, den)


def cube_coeff(c, m):
    """m-th coefficient of (sum c_k t^k)^3 using available c[0..m]."""
    s = 0
    for i in range(m + 1):
        for j in range(m + 1 - i):
            s += c[i] * c[j] * c[m - i - j]
    return s


# -------------------------------------------------------------------- System A (Hamiltonian)
def hamiltonian_series(prep_kind):
    """Taylor coefficients q_0..q_NMAX of q(t) in QQ[a,beta,xi_j,eta_j].
    prep_kind: 'shifted' -> x_j(0)=c_j a/w_j^2 + xi_j ; 'product' -> x_j(0)=xi_j."""
    names = ["a", "beta"] + [f"xi{j}" for j in range(NM)] + [f"eta{j}" for j in range(NM)]
    R, *g = ring(",".join(names), QQ)
    a, b = g[0], g[1]
    xi = g[2:2 + NM]
    eta = g[2 + NM:2 + 2 * NM]
    q = [a]; p = [R(0)]
    x = []; pj = []
    for j, (w, c) in enumerate(MODES):
        x0 = (QQ(c) / QQ(w) ** 2) * a + xi[j] if prep_kind == "shifted" else xi[j]
        x.append([x0]); pj.append([eta[j]])
    for m in range(NMAX):
        # derivative coefficients at order m -> new coefficient m+1
        rhs_p = -QQ(K) * q[m] - 4 * b * cube_coeff(q, m)
        for j, (w, c) in enumerate(MODES):
            rhs_p += QQ(c) * (x[j][m] - (QQ(c) / QQ(w) ** 2) * q[m])
        qn = p[m] * QQ(1, m + 1)
        pn = rhs_p * QQ(1, m + 1)
        for j, (w, c) in enumerate(MODES):
            xn = pj[j][m] * QQ(1, m + 1)
            pjn = (-QQ(w) ** 2 * x[j][m] + QQ(c) * q[m]) * QQ(1, m + 1)
            x[j].append(xn); pj[j].append(pjn)
        q.append(qn); p.append(pn)
    return R, q


def expect_bath(poly, law):
    """E over (xi_j, eta_j). law: 'gauss' (xi~N(0,T/w^2), eta~N(0,T)) or
    'phase' (xi=sqrt(2T)/w cos phi, eta=-sqrt(2T) sin phi)."""
    out = Rres(0)
    for monom, coeff in poly.terms():
        ea, eb = monom[0], monom[1]
        ex = monom[2:2 + NM]; ee = monom[2 + NM:2 + 2 * NM]
        val = QQ(coeff); Tpow = 0; zero = False
        for j, (w, c) in enumerate(MODES):
            i, k = ex[j], ee[j]
            if i % 2 or k % 2:
                zero = True; break
            if law == "gauss":
                val *= QQ(dfact2(i) * dfact2(k)) / QQ(w) ** i
                Tpow += (i + k) // 2
            else:
                # (sqrt(2T)/w)^i (-sqrt(2T))^k E[cos^i sin^k]
                val *= QQ(2) ** ((i + k) // 2) / QQ(w) ** i * QQ(cs_moment(i, k))
                Tpow += (i + k) // 2
        if zero or val == 0:
            continue
        out += val * ra ** ea * rb ** eb * rT ** Tpow
    return out


# -------------------------------------------------------------------- GLE comparators
def gamma_coeffs(N):
    """Taylor coefficients g_k of gamma(t) = sum c^2/w^2 cos(w t)."""
    g = []
    for k in range(N + 1):
        if k % 2:
            g.append(QQ(0))
        else:
            s = sum(QQ(c) ** 2 / QQ(w) ** 2 * QQ(w) ** k for (w, c) in MODES)
            g.append(s * QQ((-1) ** (k // 2), factorial(k)))
    return g


def gle_series(R, a, b, Fc, Sc):
    """q_0..q_NMAX from q''=-Kq-4bq^3 - int gamma(t-s)q'(s)ds + F + S; Fc,Sc Taylor coeff lists."""
    g = gamma_coeffs(NMAX)
    q = [a, R(0)]
    for m in range(NMAX - 1):
        # memory term coefficient of t^m: sum_{k+l+1=m} g_k d_l k! l!/m!,  d_l=(l+1)q_{l+1}
        mem = R(0)
        for l in range(m):
            k = m - 1 - l
            d_l = (l + 1) * q[l + 1]
            mem += g[k] * QQ(factorial(k) * factorial(l), factorial(m)) * d_l
        rhs = -QQ(K) * q[m] - 4 * b * cube_coeff(q, m) - mem + Fc[m] + Sc[m]
        q.append(rhs * QQ(1, (m + 1) * (m + 2)))
    return q


def gle_gauss_B(with_F=True, slip=False):
    nf = NMAX - 1  # f_0..f_{NMAX-2}
    names = ["a", "beta"] + [f"f{k}" for k in range(nf)]
    R, *gg = ring(",".join(names), QQ)
    a, b = gg[0], gg[1]
    f = gg[2:]
    Fc = [(f[k] * QQ(1, factorial(k)) if with_F else R(0)) for k in range(nf)] + [R(0)] * 2
    g = gamma_coeffs(NMAX)
    Sc = [(-g[k] * a if slip else R(0)) for k in range(NMAX + 1)]
    q = gle_series(R, a, b, Fc, Sc)
    return R, q, nf


def gamma_deriv0(m):
    """gamma^(m)(0)."""
    if m % 2:
        return QQ(0)
    return sum(QQ(c) ** 2 / QQ(w) ** 2 * QQ(w) ** m for (w, c) in MODES) * QQ((-1) ** (m // 2))


def make_wick(nf):
    C = [[(QQ((-1) ** l) * gamma_deriv0(k + l)) for l in range(nf)] for k in range(nf)]  # /T

    @lru_cache(maxsize=None)
    def wick(idx):  # idx: sorted tuple of f-indices; returns coefficient of T^(len/2)
        if not idx:
            return QQ(1)
        if len(idx) % 2:
            return QQ(0)
        i0 = idx[0]; rest = idx[1:]
        s = QQ(0)
        seen = set()
        for pos, j in enumerate(rest):
            # pairing i0 with element at pos; identical indices give identical sub-results
            cij = C[i0][j]
            if cij == 0:
                continue
            sub = rest[:pos] + rest[pos + 1:]
            s += cij * wick(sub)
        return s
    return wick


def expect_wick(poly, nf, wick):
    out = Rres(0)
    for monom, coeff in poly.terms():
        ea, eb = monom[0], monom[1]
        fm = monom[2:]
        idx = []
        for k, e in enumerate(fm):
            idx += [k] * e
        if len(idx) % 2:
            continue
        v = wick(tuple(idx))
        if v == 0:
            continue
        out += QQ(coeff) * v * ra ** ea * rb ** eb * rT ** (len(idx) // 2)
    return out


def gle_phase_BR():
    names = ["a", "beta"] + [f"u{j}" for j in range(NM)] + [f"v{j}" for j in range(NM)]
    R, *gg = ring(",".join(names), QQ)
    a, b = gg[0], gg[1]
    u = gg[2:2 + NM]; v = gg[2 + NM:]
    # F(t)=sum c/w (u cos wt - v sin wt), u=sqrt(2T)cos phi, v=sqrt(2T) sin phi
    Fc = []
    for k in range(NMAX + 1):
        s = R(0)
        for j, (w, c) in enumerate(MODES):
            cw = QQ(c) / QQ(w) * QQ(w) ** k / QQ(factorial(k))
            if k % 4 == 0: s += cw * u[j]
            elif k % 4 == 1: s += -cw * v[j]
            elif k % 4 == 2: s += -cw * u[j]
            else: s += cw * v[j]
        Fc.append(s)
    Sc = [R(0)] * (NMAX + 1)
    q = gle_series(R, a, b, Fc, Sc)
    return R, q


def expect_uv(poly):
    out = Rres(0)
    for monom, coeff in poly.terms():
        ea, eb = monom[0], monom[1]
        eu = monom[2:2 + NM]; ev = monom[2 + NM:]
        val = QQ(coeff); Tp = 0; zero = False
        for j in range(NM):
            i, k = eu[j], ev[j]
            mom = cs_moment(i, k)
            if mom == 0:
                zero = True; break
            val *= QQ(2) ** ((i + k) // 2) * QQ(mom)
            Tp += (i + k) // 2
        if zero:
            continue
        out += val * ra ** ea * rb ** eb * rT ** Tp
    return out


def deriv_list(qs):
    return [qs[n] * factorial(n) for n in range(len(qs))]


def main():
    t0 = time.time()
    print(f"=== V0-1 P-17 Taylor checks: {NM} mode(s) {MODES}, n <= {NMAX}; K=23/10, M=1 ===")
    gam0 = gamma_deriv0(0)
    print(f"gamma(0) = {gam0}")

    RA, qA = hamiltonian_series("shifted")
    EA = deriv_list([expect_bath(c, "gauss") for c in qA])
    ER = deriv_list([expect_bath(c, "phase") for c in qA])
    print(f"[A,R done {time.time()-t0:.1f}s]")
    RP, qP = hamiltonian_series("product")
    EP = deriv_list([expect_bath(c, "gauss") for c in qP])
    print(f"[product done {time.time()-t0:.1f}s]")

    RB, qB, nf = gle_gauss_B(True, False)
    wick = make_wick(nf)
    EB = deriv_list([expect_wick(c, nf, wick) for c in qB])
    print(f"[B done {time.time()-t0:.1f}s]")
    RBs, qBs, nf2 = gle_gauss_B(True, True)
    EBs = deriv_list([expect_wick(c, nf2, wick) for c in qBs])
    print(f"[B+slip done {time.time()-t0:.1f}s]")
    RD, qD, _ = gle_gauss_B(False, False)
    ED = deriv_list([Rres(0) + expect_wick(c, nf, wick) for c in qD])
    RBR, qBR = gle_phase_BR()
    EBR = deriv_list([expect_uv(c) for c in qBR])
    print(f"[D, B_R done {time.time()-t0:.1f}s]")

    def show(label, lst):
        print(f"\n--- {label} ---")
        for n in range(NMAX + 1):
            d = lst[n]
            print(f"  n={n:2d}: {'0' if d == 0 else to_expr(d)}")

    show("T1: A - B", [EA[n] - EB[n] for n in range(NMAX + 1)])
    AD = [EA[n] - ED[n] for n in range(NMAX + 1)]
    show("T2: A - D", AD)
    lead = next(n for n in range(NMAX + 1) if AD[n] != 0)
    pred = Rres(factorial(6)) * (-rb * ra * rT * gam0 * QQ(1, 10))
    print(f"\n--- T3: leading nonzero A-D at n={lead}: {to_expr(AD[lead])};"
          f" 6!*[-beta a T gamma(0)/10] = {to_expr(pred)}; equal: {AD[lead] == pred and lead == 6}")
    show("T4: A - R", [EA[n] - ER[n] for n in range(NMAX + 1)])
    show("T5: R - B_R", [ER[n] - EBR[n] for n in range(NMAX + 1)])
    show("T6: B_R - B", [EBR[n] - EB[n] for n in range(NMAX + 1)])
    show("T7a: product - B", [EP[n] - EB[n] for n in range(NMAX + 1)])
    show("T7b: product - (B + slip -gamma(t) a)", [EP[n] - EBs[n] for n in range(NMAX + 1)])

    # T8 controls: substitute beta=0 or T=0
    print("\n--- T8: controls ---")
    ok_b = all(AD[n].as_expr(Ssym_a, Ssym_b, Ssym_T).subs(Ssym_b, 0) == 0 for n in range(NMAX + 1))
    ok_T = all(AD[n].as_expr(Ssym_a, Ssym_b, Ssym_T).subs(Ssym_T, 0) == 0 for n in range(NMAX + 1))
    print(f"  beta=0 => A-D=0 at all n: {ok_b}")
    print(f"  T=0    => A-D=0 at all n: {ok_T}")
    # also directly: beta=0 Hamiltonian is linear in bath data -> zero-mean deviation
    print("\n--- raw E_A[q^(n)(0)] (thermal) for reference ---")
    for n in range(NMAX + 1):
        print(f"  n={n:2d}: {to_expr(EA[n])}")
    print(f"\n[total time {time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
