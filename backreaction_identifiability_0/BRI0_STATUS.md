# BRI0 STATUS

**State:** **BRI1-X1-THEOREM PROVED — BRI-E2+O EARNED AT CLASS-THEOREM LEVEL FOR X1** (X1 ∉ E₂± for all
sufficiently large finite N_B; the witness vanishes O(N_B⁻¹) in the reservoir limit, which is E₁). The frozen-τ grade stays
**PF4Q-I / X1-PF-INDETERMINATE**. INTERNALLY PROVED / NOT EXTERNALLY REVIEWED. AWAITING OWNER REVIEW.

**Branch:** `grut-backreaction-identifiability-0`.

**Frozen parent:** `scout-0 @ ab2da47407bb670d94e2a52c87599fa13fd8ab99` (unmodified).

**Programme:** BRI0 — BACK-REACTION IDENTIFIABILITY 0. Lane 3; a **new-premise campaign** downstream of P-17 / P-18 and
the SCOUT-0 saturation map. It does not reopen or weaken the SCOUT-0 verdict. It is not part of DA0.

**Boundaries:** Charter 0 `89b8236`; Scope Repair 01 `8df4b88`; Scope Repair 02 and Candidate 1 `6220a9e` (accepted);
PF4Q declaration `6d25e7e`; PF4Q amendment `35d6985`; PF4Q run `d5a0bdb` (accepted); analytic escape theorem = this
commit. See `BRI0_LEDGER.md` and `BRI1_PF4Q.md`.

**Accepted (owner ruling, Scope Repair 01, `BRI0_CHARTER.md` §R):**
- **Clamp family 𝒫 = 𝒳 accepted** as the primary interventional protocol family. The object is the interventional
  environment-force law. The at-rest clamp is the default reference only at times where its force variance is non-zero.
- **E₂± accepted as primary** (signed affine class: G ∈ ℝ \ {0}). The positive-G class E₂ is retired. The ladder is
  **E₁ ⊊ E₂± ⊊ E_univ**, with strictness re-proved reflection-safe.
- **Zero-amplitude firewall:** E₂± membership requires a common degeneracy set. A protocol-dependent zero is graded
  **BRI-DEG** and never earns BRI-E2+ alone.
- **BRI-C1 and BRI-C2± accepted** as exact at the fdd level (a path-law upgrade only under continuous / càdlàg
  regularity):
  - E₂± membership holds exactly when there are **causal** signs s_q with s_q·Z_q having a common law.
  - **Correction recorded:** the per-tuple reflection-orbit condition is necessary but **not** sufficient (explicit
    counterexample). The primary witness is a reflection-orbit violation. "Causal-sign inconsistency" (**BRI-CSI**) is a
    distinct escape mode whose grade needs an owner ruling.
- **Witnesses made reflection-safe:** even standardised cumulants and |odd| ones, |correlations|, and dependence modulo
  reflections. Raw skew-sign or copula-sign changes never count.
- **BRI-UPPER earned at its precise scope** (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED): the continuous-time
  deterministic causal parent class lies in E_univ, back-reaction included. The discrete-time randomisation result for
  stochastic kernels is kept separately. **No continuous-time general-stochastic claim.** Quantum is out of scope.

**Scope Repair 02 (§R2):**
- E₂± membership requires **one shared causal sign functional** S_t[q_[0,t]].
- Positives: **BRI-E2+O** (orbit / shape) and **BRI-E2+C** (causal-affine coherence failure). Both are genuine BRI-E2+.
- Finite CSI certificate registered (BRI-O11).

**BRI1 analytic escape theorem (`BRI1_ANALYTIC_ESCAPE_THEOREM.md`):**
- **T1:** c_{P1}(t, t; t) < 0 on some (0, δ). This uses an exact Taylor expansion through t⁷ with an integrable
  remainder.
- **T2:** κ₃(F_{P1,N}(t*)) ≠ 0 for all N_B ≥ N₀.
- **T3:** the P0 skewness is exactly 0, so P1 shows a reflection-safe orbit violation and X1 ∉ E₂±.
- **T4:** the reservoir limit is the protocol-independent Gaussian, i.e. E₁.
- Existential δ and N₀; no new computation.
- Earned: the identifiability theorem relative to E₂± only. **Not** GRUT physics, primitive randomness, ontology, or an
  escape from E_univ.

**BRI1-PF4Q (`BRI1_PF4Q.md`), accepted; certified-τ pipeline not to be built now:**
- **LEMMA BRI1-R1 proved:** κ₃(F) = K/N_B + O(N_B⁻²).
- **V1** certified (arb); **V3** passes; **V2** led to a pre-run Method-A amendment.
- **All 20 frozen K_q[a,b,c] are non-zero as evidence:** two independent methods agree to about 10⁻¹¹, with |K| from
  5.65·10⁻³ to 0.25.
- **No certified enclosure:** the validated-ODE / certified-quadrature machinery is unavailable, so this is **PF4Q-I**.
- The PF-A path (a certified pipeline) is an owner decision.

**Candidate 1 (`BRI1_CANDIDATE_CHARTER.md`), analytic preflight only:**
- **C1:** the A = q + q³/3 harmonic control lies **in E₂± exactly** (proof plus an exact series check). Structurally it
  is BRI-E1.
- **X1:** the finite Duffing bath, ε = N_B^{−1/2}, is well posed and reciprocal.
  - The exact i.i.d. structure under clamp gives κ_n(F) = N_B^{1−n/2}κ_n(X^ε).
  - The first non-affine departure is O(N_B^{−1}), and X1 → E₁ in the reservoir limit.
  - The O(N_B^{−1}) odd-cumulant coefficient is proved non-zero as a functional, but **not certified at the frozen τ**.
- **CSI** is structurally impossible at τ.
- **Terminal: X1-PF-INDETERMINATE.** Resolution proposed: a certified deterministic quadrature, pending owner
  authorisation.
- **Frozen:** P0 / P1 / P2, τ = (π, 3π/2, 2π), D_orb, N_B ∈ {1, 2, 4, 8, 16, 32, 64}. No sample count.

No Monte Carlo, no finite-N_B simulation, no numerical D_orb. Deterministic quadrature only, for PF4Q. No PR. No merge.
