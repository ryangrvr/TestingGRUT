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
| Friedan–Qiu–Shenker 1984; GKO coset | unitary Virasoro reps with c < 1 occur only at c = 1 − 6/(m(m+1)) with Kac weights | STANDARD-TEXTBOOK ✓ (table; lattice c = 0.501 at non-integrable point) | W1-R |
| Calabrese–Cardy entanglement scaling | S(l) = (c/6) log[(2N/π) sin(πl/N)] (OBC), (c/3) log(L/π) (PBC half chain) | STANDARD-TEXTBOOK ✓ | W1-R |
| XXZ Bethe ansatz (v, x₁) | v = π√(1−Δ²)/(2 arccos Δ); x₁ = (π − arccos Δ)/(2π) | STANDARD-TEXTBOOK ✓ (ED L ≤ 20) | W1-R |
| U(1) current algebra ⇒ c ≥ 1 | gapless 1D U(1)-conserving CFT contains a U(1) Kac–Moody algebra (c ≥ 1) | STANDARD-TEXTBOOK (not separately checked) | W1-R |
| van Hove / image method / threshold resonance (Levinson) / Bessel processes | γ = d/2 − 1 + #Dirichlet on hypercubic lattices; point defects relevant (d=1), marginal (d=2), irrelevant unless resonant (d=3); radial dimension D gives return probability t^{−D/2} | STANDARD-TEXTBOOK ✓ (Bessel exact; G₀ quadrature; expm_multiply n = 20000) | W1-S |
| Gleason 1957; Busch 2003; Gisin 1990 / Polchinski 1991 | noncontextual frame functions (d ≥ 3) / additive effect probabilities (d ≥ 2) are Born; nonlinear outcome weights signal under composition | STANDARD-TEXTBOOK ✓ (d=2 non-Born family; d=3/POVM failures; signalling spread up to 3.5e-2) | W1-G |
| Hardy 2001; Chiribella–D'Ariano–Perinotti 2011; Masanes–Müller 2011 (local tomography) | K_AB = K_A K_B holds for complex, fails for real QM | STANDARD-TEXTBOOK ✓ (counting; rebit witness) | W1-L |
| quaternionic composites not locally tomographic | — | SECONDARY (not checked) | W1-L |
| Renou et al., Nature 600, 625 (2021) | network experiment excludes real QM under independent sources | SECONDARY (not readable from sandbox) | W1-L |
| Gibbons–Hawking 1977; Bunch–Davies; Unruh–DeWitt detector | dS Bunch–Davies state is KMS at T = H/2π along geodesics | STANDARD-TEXTBOOK ✓ (detailed balance to 7 digits) | W1-F |
| α-vacua non-thermal (Mottola; Allen) | dS-invariant non-Hadamard states are not KMS at H/2π | SECONDARY (not computed) | W1-F |
| ETH (Deutsch 1991; Srednicki 1994; Rigol–Dunjko–Olshanii 2008); GGE for integrable chains | non-integrable closed systems relax locally to Gibbs at β fixed by energy; integrable ones to a GGE | STANDARD-TEXTBOOK (ED-indicated L ≤ 12) | W2-ETH |
