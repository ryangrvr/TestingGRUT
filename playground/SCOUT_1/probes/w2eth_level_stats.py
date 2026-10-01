"""SCOUT-1 W2-ETH addendum: does genericity alone fix a universal weight-0 number?
Level-spacing ratio <r> = < min(s_n, s_{n+1}) / max(s_n, s_{n+1}) > over the middle half of the spectrum.
  GOE (generic, time-reversal real H): <r> = 0.5307;  GUE: 0.5996;  Poisson (integrable): 2 ln 2 - 1 = 0.3863.
Symmetries removed by weak random site fields with OBC (no translation, no reflection).
  generic:     H = sum Z Z + sum g_i X_i + sum h_i Z_i               (real => GOE)
  integrable:  H = sum Z Z + sum g_i X_i                              (free fermions, random g_i => Poisson-like)
  generic + complex: add sum b_i Y_i  (breaks the antiunitary => GUE)
"""
import numpy as np

rng = np.random.default_rng(17)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0 + 0j, -1.0]); I2 = np.eye(2)


def op(o, i, L):
    out = np.array([[1.0 + 0j]])
    for k in range(L):
        out = np.kron(out, o if k == i else I2)
    return out


def ham(L, g, h, b):
    Hm = np.zeros((2 ** L, 2 ** L), complex)
    for i in range(L - 1):
        Hm += op(Z, i, L) @ op(Z, i + 1, L)
    for i in range(L):
        Hm += g[i] * op(X, i, L) + h[i] * op(Z, i, L) + b[i] * op(Y, i, L)
    return Hm


def rbar(E):
    E = np.sort(E)
    n = len(E)
    E = E[n // 4: 3 * n // 4]
    s = np.diff(E)
    s = s[s > 1e-12]
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return r.mean()


L = 11
print(f"L = {L} OBC, middle half of the spectrum; references GOE 0.5307, GUE 0.5996, Poisson 0.3863")
for label, gen in (
    ("generic real (GOE expected)", lambda: (0.9045 + 0.2 * rng.standard_normal(L), 0.809 + 0.2 * rng.standard_normal(L), np.zeros(L))),
    ("integrable random-TFIM (Poisson expected)", lambda: (0.9045 + 0.3 * rng.standard_normal(L), np.zeros(L), np.zeros(L))),
    ("generic complex (GUE expected)", lambda: (0.9045 + 0.2 * rng.standard_normal(L), 0.809 + 0.2 * rng.standard_normal(L), 0.5 + 0.2 * rng.standard_normal(L))),
):
    vals = []
    for rep in range(3):
        g, h, b = gen()
        E = np.linalg.eigvalsh(ham(L, g, h, b))
        vals.append(rbar(E))
    print(f"  {label:44s}: <r> = {np.round(vals, 4)}  mean {np.mean(vals):.4f}")
print("  coupling magnitudes vary across realizations; <r> does not (universality); its value is set by the")
print("  supplied antiunitary-symmetry class (Dyson's threefold way): genericity fixes the value WITHIN a class.")
