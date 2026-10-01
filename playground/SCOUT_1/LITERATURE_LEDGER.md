# SCOUT-1 LITERATURE LEDGER (imports with source-access grade)

**Grades:**
- **PRIMARY-READ:** the primary text was read in this campaign.
- **PRIMARY-VERIFIED-EXTERNAL:** an external auditor read it.
- **STANDARD-TEXTBOOK:** classical result; statement reproduced from standard knowledge and checked by explicit
  computation where marked ✓.
- **SECONDARY:** search excerpt only.

The sandbox blocks arXiv and several publishers.

| Import | Statement used | Grade | Used in |
|---|---|---|---|
| Curie's principle / Buckingham-π | symmetric premises cannot yield a conclusion that breaks the symmetry; dimensionful quantities fixed only relative to a reference scale | STANDARD-TEXTBOOK ✓ (fixed-point lemma, sympy) | W1-C |
| Spectral theorem, cyclic vector | `(K, e_r)` with `e_r` cyclic is determined up to unitary equivalence fixing `e_r` by the spectral measure `μ_r` | STANDARD-TEXTBOOK ✓ (random hidden U, moments to 2e-15) | W1-C, W1-I |
| Jacobi/Lanczos uniqueness | the Jacobi matrix generated from `(K, e_r)` is a function of `μ_r` alone | STANDARD-TEXTBOOK ✓ (Lanczos coefficients agree to 2e-15) | W1-C, W1-I |
| Lieb–Schultz–Mattis / Oshikawa / Yamanaka–Oshikawa–Affleck | U(1) × translation at non-integer filling ν: no unique gapped ground state; low-energy state at momentum 2πν | STANDARD-TEXTBOOK ✓ (momentum-resolved ED, 10 models) | W1-A |
| Stieltjes/Jacobi inverse spectral theorem; Krein string; Hankel moment problem | end-site measure of a Jacobi matrix determines it up to signs; moments m_0..m_2k fix k layers; Hankel matrices are exponentially ill-conditioned | STANDARD-TEXTBOOK ✓ (reconstruction 4e-12; depth test; cond up to 2e34) | W1-I |
| Pusz–Woronowicz 1978; Lenard 1978 | completely passive ⟺ Gibbs at a single β ∈ [0, ∞] | STANDARD-TEXTBOOK ✓ (exact ergotropy of N copies; dense-bath single-copy window) | W1-P |
