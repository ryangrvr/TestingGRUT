# Numerical arithmetic scope

Inherit the frozen S8 ARITHMETIC.md contract: binary64 round-to-nearest,
correctly rounded square root, standard absolute BLAS reduction bounds,
positive long-double norm accumulation with an explicit conversion fee,
underflow pads, checked SVD/eigenvector Gram and residual bounds, full gate
audits, finite-receiver cut and time-grid accounts. No library factorization
output is treated as exact. Dependency hashes identify every imported source.

The new production fee is explicit in certificate.py. For stored X,K,Z,
R_hat=fl(X-KZ), the true residual differs by at most

    gamma_(16r+32) (||X||_F+||K||_F||Z||_F) + underflow,

where r=K's column count. For w_hat=fl(K^dagger R_hat), include

    gamma_(16n+32)||K||_F||R_hat||_F + sqrt(1+eta_K)*radius_R

plus its underflow pad; n=K's row count. These overcount ordinary complex
dot products. Positive norms enclose the stored arrays separately. Unknown
physical prefix/suffix Grams are handled by SWEEP_CERTIFICATE.md. SVD
reconstruction and existing gate/input arithmetic remain separate fees.

Positive squared-loss sums are rounded upward within each eligible batch;
their square-root increments enter the cumulative ledger. In-range defects,
input/gate errors and normalization costs are always linear. Every reversal
starts a new batch. No arithmetic or projection fee is retroactively lowered
in S8. A stopped partial sweep has no target-time probability certificate.

Analytic proof constants in formation_controls.constants use exact rational
arithmetic. Small-matrix, spectral quadrature and reference-walk diagnostics
use ordinary double precision and are not outward-rounded all-N proofs.
This is not formal verification of BLAS/compiler/hardware or external review.
