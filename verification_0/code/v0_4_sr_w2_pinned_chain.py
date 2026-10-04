#!/usr/bin/env python3
"""V0-4 second reader, witness W2 (E-1, E-14, E-19: the pinned conservative chain).

Fresh code. Stand-in model built from the record's stated data only:
  K = c I - T on a half-line chain, c = 2 + pin, T = nearest-neighbour hopping (off-diagonal 1),
  so spec(K) = [pin, pin + 4] (record: [0.3, 4.3] at pin = 0.3; K_11 = 2.3 = omega_s^2).

Part (a)  vary the spectral bottom (pin = supplied S-1 datum) at fixed earned predicates:
   r = (K_inf^{-1})_11 = (c - sqrt(c^2 - 4))/2       (record: r = (2.3 - sqrt(1.29))/2)
   X_J(inf) = (T_s - T_b) + T_b r^2 / 2               (record LS-2 closed form)
   zero of X_J at T_s/T_b = 1 - r^2/2                 (the visible SIGN of X_J is pin-dependent)
   kappa = Im F(omega_s^2) / (2 omega_s) = 1/(2 sqrt(c))   (record: kappa = 1/(2 sqrt 2.3))
   relaxational memory k(tau) = [exp(-K tau)]_11:  k_pin = exp(-pin tau) k_0  (E-1 factorization);
   late-time decay rate = spectral bottom = pin.
Part (b)  vary the EDGE EXPONENT at fixed band edges and fixed locality (nearest-neighbour, conservative):
   Gegenbauer-graded Jacobi chain: end-site spectral measure ~ (1 - x^2)^(alpha - 1/2) on the band;
   alpha = 1 is the uniform chain (record: square-root edges, t^{-3/2} tails); alpha = 2 gives
   (1 - x^2)^{3/2} edges. Native oscillator response phi(t) = [cos(sqrt(K) t)]_11 tail ~ t^{-(alpha+1/2)}.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, expm


def r_closed(pin):
    c = 2 + pin
    return (c - np.sqrt(c*c - 4))/2


def r_numeric(pin, N=4000):
    c = 2 + pin
    d = np.full(N, c)
    off = np.full(N - 1, -1.0)
    # solve K x = e1 by Thomas algorithm
    from scipy.linalg import solve_banded
    ab = np.zeros((3, N))
    ab[0, 1:] = off
    ab[1, :] = d
    ab[2, :-1] = off
    e1 = np.zeros(N); e1[0] = 1
    return solve_banded((1, 1), ab, e1)[0]


def memory_kernel(pin, taus, N=600):
    c = 2 + pin
    w, V = eigh_tridiagonal(np.full(N, c), np.full(N - 1, -1.0))
    wt = V[0, :]**2
    return np.array([np.sum(wt*np.exp(-w*t)) for t in taus]), w.min()


def gegenbauer_K(alpha, N, pin=0.3):
    n = np.arange(1, N)
    beta = n*(n + 2*alpha - 1)/(4.0*(n + alpha)*(n + alpha - 1))
    off = 2*np.sqrt(beta)          # -> 1 as n -> infinity
    c = 2 + pin
    return np.full(N, c), -off


def tail_exponent(alpha, N=6000, pin=0.3):
    d, off = gegenbauer_K(alpha, N, pin)
    lam, V = eigh_tridiagonal(d, off)
    wt = V[0, :]**2
    om = np.sqrt(lam)
    t = np.linspace(40, 400, 7201)
    phi = np.array([np.sum(wt*np.cos(om*tt)) for tt in t])
    a = np.abs(phi)
    # envelope: local maxima of |phi|
    pk = [i for i in range(1, len(a) - 1) if a[i] >= a[i-1] and a[i] >= a[i+1]]
    tp, ap = t[pk], a[pk]
    # upper envelope: keep running maxima in log bins
    bins = np.logspace(np.log10(40), np.log10(400), 25)
    et, ea = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (tp >= lo) & (tp < hi)
        if m.any():
            j = np.argmax(ap[m]); et.append(tp[m][j]); ea.append(ap[m][j])
    s = np.polyfit(np.log(et), np.log(ea), 1)[0]
    edge_lo, edge_hi = lam.min(), lam.max()
    return -s, edge_lo, edge_hi


def main():
    out = []
    P = out.append
    P("W2  pinned conservative chain: E-1 / E-14 / E-19 (fresh stand-in code)")
    P("=" * 100)
    P("(a) vary the spectral bottom pin (supplied S-1) -- predicates unchanged, edge-read quantities move")
    pairs = [(2, 1), (10, 1), (0.5, 1), (0.1, 1), (0.9, 1)]
    for pin in (0.1, 0.3, 0.6, 1.0):
        rc, rn = r_closed(pin), r_numeric(pin)
        kap = 1/(2*np.sqrt(2 + pin))
        xs = [(Ts - Tb) + 0.5*Tb*rc**2 for Ts, Tb in pairs]
        P(f"  pin={pin:4.2f}  band=[{pin:.2f},{pin+4:.2f}]  r={rc:.6f} (numeric {rn:.6f})  kappa={kap:.6f}"
          f"  X_J zero at T_s/T_b={1-rc**2/2:.5f}")
        P("            X_J(inf) at (T_s,T_b)=" + ", ".join(f"{p}:{x:+.5f}" for p, x in zip(pairs, xs)))
    P("  record check (pin=0.3): r = %.6f, |X_J| (2,1),(10,1),(1/2,1),(1/10,1) = %s ; kappa = %.6f = 1/(2 sqrt 2.3)"
      % (r_closed(0.3), [round(float(abs((a - b) + 0.5*b*r_closed(0.3)**2)), 5) for a, b in pairs[:4]], 1/(2*np.sqrt(2.3))))
    P("  -> the predicate NET-ARROW (X_J, X_sigma of fixed forward sign on the declared pairs) holds at every pin;")
    P("     the numbers, and the sign of X_J for a near-equilibrium pair such as (0.9, 1), are functions of the")
    P("     band bottom: at (0.9,1) X_J changes sign where r^2/2 = 0.1, i.e. it is an edge READ, not a FIX.")
    taus = np.array([5.0, 10.0, 20.0, 40.0])
    k0, l0 = memory_kernel(0.0, taus)
    for pin in (0.3, 0.6):
        kp, lp = memory_kernel(pin, taus)
        ratio = kp/(k0*np.exp(-pin*taus))
        rate = -np.log(kp[-1]/kp[-2])/(taus[-1] - taus[-2])
        P(f"  E-1 memory: pin={pin}: k_pin/(e^(-pin tau) k_0) = {np.round(ratio, 12).tolist()} ;"
          f" lambda_min = {lp:.4f} ; late log-rate(20->40) = {rate:.4f}")
    P("  -> both pins are EXPONENTIAL-grade (P_memory holds); the decay rate is the supplied spectral bottom.")
    P("")
    P("(b) vary the edge exponent at fixed band edges, fixed locality (nearest-neighbour), fixed conservativity")
    for alpha in (1.0, 2.0):
        ex, lo, hi = tail_exponent(alpha)
        P(f"  alpha={alpha}: end-site measure ~ (1-x^2)^{alpha-0.5:.1f}; spectrum [{lo:.4f},{hi:.4f}];"
          f" native phi(t) envelope exponent = {ex:.3f} (theory alpha+1/2 = {alpha+0.5})")
    P("  -> 'native g=1 -> non-Markovian dissipation only, branch-cut tail' holds for both; the tail exponent")
    P("     (t^-3/2 in the record) is a read-out of the supplied chain's edge exponent, not fixed by S5-1.")
    print("\n".join(out))


if __name__ == "__main__":
    main()
