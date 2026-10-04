#!/usr/bin/env python3
"""V0-4 second reader, witness W3 (E-21: sector-selection S-1, per-branch counting + min-dominance).

Fresh code. Two influence spectral densities (low-frequency / spectral-bottom edge exponent is the
datum) that obey the SAME per-branch counting law:
    branch 1 : exponent s        (base class)
    branch 2 : exponent s + 4    (order-2 symmetry cancellation, the record's '+4')
  scalar channel  J(w) = a w^s + b w^(s+4)
  matrix channel  J(w) = [[a w^s, c w^(s+2)], [c w^(s+2), b w^(s+4)]]  (observable = smaller eigenvalue
                  and = larger eigenvalue; both reported)
Only the AMPLITUDES (a, b, c) differ: the record says amplitudes are sector-supplied (S-7/S-11; PC-D
"lambda J admissible for every lambda > 0"). We report the true w -> 0 exponent and the effective
two-point slope on windows [w0, 2 w0].
"""
import numpy as np


def slope(f, w0):
    w1, w2 = w0, 2*w0
    return (np.log(f(w2)) - np.log(f(w1)))/np.log(2)


def main():
    s = 3.0
    out = []
    P = out.append
    P("W3  E-21 branch counting with amplitude-dependent masking (fresh code)")
    P("=" * 100)
    P(f"per-branch exponents fixed by the counting law: s = {s}, s+4 = {s+4}")
    cases = [(1.0, 1.0), (1e-2, 1.0), (1e-4, 1.0), (0.0, 1.0)]
    windows = [1e-3, 1e-2, 1e-1, 0.45, 0.9]
    P("scalar J = a w^s + b w^(s+4)")
    for a, b in cases:
        f = lambda w, a=a, b=b: a*w**s + b*w**(s+4)
        true = s if a > 0 else s + 4
        sl = [slope(f, w0) for w0 in windows]
        P(f"  a={a:7.0e} b={b}: true w->0 exponent = {true:.0f};  window slopes " +
          ", ".join(f"[{w0:g},{2*w0:g}]:{x:.3f}" for w0, x in zip(windows, sl)))
    P("matrix channel (eigen-slopes), a=1e-4, b=1, c varied:")
    for c in (0.0, 1e-3, 5e-3):  # matrix positivity (cone) needs c^2 < a*b = 1e-4
        def ev(w, c=c, k=0):
            M = np.array([[1e-4*w**s, c*w**(s+2)], [c*w**(s+2), w**(s+4)]])
            return np.sort(np.abs(np.linalg.eigvalsh(M)))[k]
        sl_lo = [slope(lambda w: ev(w, k=0), w0) for w0 in windows]
        sl_hi = [slope(lambda w: ev(w, k=1), w0) for w0 in windows]
        P(f"  c={c:6.0e}: small-eig slopes " + ", ".join(f"{x:.3f}" for x in sl_lo) +
          " | large-eig slopes " + ", ".join(f"{x:.3f}" for x in sl_hi))
    P("READING")
    P("  * Every model obeys the identical per-branch counting law (s and s+4).")
    P("  * The observable spectral-edge exponent (w -> 0) is s if the base branch has non-zero amplitude and")
    P("    s+4 if it vanishes; on any finite window the effective exponent interpolates (non-integer values")
    P("    such as the record's +2.216 shift are of this kind). The amplitude ratio is supplied, so the")
    P("    counting law SELECTS AMONG classes (CON) and does not FIX the edge exponent.")
    print("\n".join(out))


if __name__ == "__main__":
    main()
