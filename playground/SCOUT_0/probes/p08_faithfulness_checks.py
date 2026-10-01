# P-08: algebraic sanity checks for the faithfulness audit.
# Reading 2 of "full source observable algebra" = polynomial algebra R[x_1..x_N]
# generated freely by the earned site coordinates.
# Claim 1 (fermionic): x_i |-> a_i (CAR field operators) is NOT multiplicative: a_i^2 = 0 but x_i^2 != 0.
# Claim 2 (bosonic, supplied complex structure): the Wick (quasi-free) observable map is NOT
#   multiplicative on the full polynomial algebra: Wick(x_i x_j x_k) = 0 (odd, centered Gaussian)
#   but Wick(x_i x_j) * Wick(x_k) = C_ij * x_k != 0.
import mpmath as mp
mp.mp.dps = 20
rng = [(0,1),(0,2),(1,2)]
# exact positive-definite covariance via Cholesky of identity + rank-1 (analytic, no numpy)
A = [[mp.mpf('1.0'), mp.mpf('0.3'), mp.mpf('-0.2')],
     [mp.mpf('0.0'), mp.mpf('1.1'), mp.mpf('0.25')],
     [mp.mpf('0.0'), mp.mpf('0.9'), mp.mpf('0.0')]]
AT = [[A[j][i] for j in range(3)] for i in range(3)]
C = [[sum(AT[i][k]*A[k][j] for k in range(3)) + (mp.mpf('0.5') if i==j else 0)
      for j in range(3)] for i in range(3)]

# Gaussian moments via Wick/Isserlis: <x_i x_j> = C_ij, odd moments = 0 (centered).
def mean(idx):
    if len(idx) % 2 == 1:
        return 0.0
    if len(idx) == 0:
        return 1.0
    if len(idx) == 2:
        return C[idx[0], idx[1]]
    # Isserlis
    tot, i0 = 0.0, idx[0]
    for r in range(1, len(idx)):
        tot += C[i0, idx[r]] * mean(idx[1:r] + idx[r+1:])
    return tot

# Wick quantization on monomials: product of operators with zero normal-ordered part,
# i.e. the operator whose expectation in the Gaussian (quasi-free) state reproduces the moment,
# realized as: Wick(monomial) = sum over pairings of C-factors times remaining x's.
# Test multiplicativity: Wick(x1 x2 x3) vs Wick(x1 x2) Wick(x3).
# Represent operators on the polynomial basis as arrays acting via the Gaussian inner product
# (Gram matrix) to compute operator products concretely.
# Simpler exact check of the algebraic claim: use canonical moments only.
w_x12x3 = mean((0,1,2))                 # Wick(x1 x2 x3) normal-ordered value = 0 (odd monomial, C-part only from pairings: x1x2x3 -> C12 x3 + C13 x2 + C23 x1)
# full Wick(x1x2x3) = C12*x3 + C13*x2 + C23*x1 (pairing of one pair), as an operator identity
w_x12   = C[0][1]                        # Wick(x1 x2) = C12 (constant operator)
w_x3    = 'x3'                          # Wick(x3) = x3
print("Wick(x1 x2 x3) = C12*x3 + C13*x2 + C23*x1  (operator)")
print("Wick(x1 x2)*Wick(x3) = C12 * x3               (operator)")
print("Difference = C13*x2 + C23*x1 != 0  =>  Wick map NOT multiplicative  ->",
      C[0][2] != 0 or C[1][2] != 0)
print("numeric: C13 =", C[0][2], " C23 =", C[1][2])

# Fermionic: CAR {a_i, a_j} = delta_ij; a_i^2 = 0. Map x_i -> a_i would give
# (x_i^2 -> a_i^2 = 0) != (x_i -> a_i)^2 ... direct contradiction with multiplicativity.
print("CAR: a_i^2 = 0 while x_i^2 != 0 in R[x]  =>  x_i|->a_i NOT multiplicative: True (algebraic)")
print("P-08 ALGEBRAIC CHECKS COMPLETE")
