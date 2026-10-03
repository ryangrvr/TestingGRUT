# RA0 · G5 — METASTABLE BLIND-RECOVERY THEOREM (result)

**Central question.** Under the G4-T hypotheses with fixed rank K, does η_N/g_N → 0 make the cubic moment tensor of the
canonical slow space near-odeco, with primitive components converging to the block-indicator idempotents? And does that
make blind partition recovery theorem-grade?

**Answer: yes, in the sense of vanishing misclassified π-mass.**
- Exact eventual recovery is **not** proved. It needs a new pointwise margin hypothesis (H4).
- All proofs below are written out here (labelled PROVED HERE), and **not externally reviewed**.
- Numerical audit: `g5/g5_recovery.py` + `g5/g5_recovery.log` (independent code path, not independent reviewer).

**Notation** (all on L²(π)):
- A = A_𝔅 = span{1_{B_i}}, the hidden partition algebra, with p_i = π(B_i), Λ = p_min^{−1/2} and λ_i = p_i^{−1/2} ∈ [1, Λ].
- V = V_N is the slow space at the certified rank K.
- s = ‖E_V − E_A‖, and C_V = sup_{f ∈ V, ‖f‖ = 1} ‖f‖_∞.
- ‖·‖ on 3-tensors is the injective norm. For symmetric tensors it equals max_{‖a‖ = 1} |T(a, a, a)| (Banach's theorem).

## G5.1 Hypotheses (PRICED)

| | hypothesis | role |
|---|---|---|
| H1 | reversible, irreducible P_N | π derived; spectral projectors π-orthogonal |
| H2 | hidden proof partition 𝔅_N with fixed K; η_N/g_N → 0 | gives s ≤ η/(g − η) → 0 (G4-T, Davis–Kahan); certifies the rank-K cut |
| H3′ | p_min ≥ p_* > 0, and **s_N² · C_{V_N} → 0** | a **weakening** of G4-T's H3 (C_V bounded). G5 never uses Δ_alg. A crude sufficient condition is s² · π_min^{−1/2} → 0 |
| (H4) | pointwise margin, needed **only** for exact recovery: max_i ‖ẽ_i − 1_{B_i}‖_∞ < 1/2 eventually | not assumed by the main theorem |

**The family is supplied.** 𝔅_N appears **only in proofs**. The recovery map R below never sees it.

## G5.2 LEMMA G5-1 (exact tensor; PROVED HERE, elementary)

For u_i = 1_{B_i}/√p_i (a π-orthonormal basis of A):

  T_A(a, b, c) := E_π[abc] = Σ_i λ_i (a·u_i)(b·u_i)(c·u_i), with **λ_i = p_i^{−1/2}**.

This is orthogonally decomposable (odeco).

**Proof.** E_π[u_i u_j u_l] = δ_{ijl}·p_i·p_i^{−3/2}. ∎

**The compressed product.** Because A is an algebra, the compressed product a∗b := E_A(ab) equals ab.
- In u-coordinates, (x∗y)_i = λ_i x_i y_i.
- Its idempotents are the solutions of λ_i x_i² = x_i, i.e. x_i ∈ {0, 1/λ_i}. There are exactly **2^K** of them:
  x_S = Σ_{i∈S} u_i/λ_i = 1_{∪_{i∈S} B_i}.
- ‖x_S‖² = Σ_{i∈S} p_i ≤ 1.

**Primitive idempotents.** These are x_{{i}} = 1_{B_i} = u_i/λ_i, i.e. tensor component divided by its coefficient.

**The multiplication operator.** L_{x_S} := T_A(x_S, ·, ·) = diag(λ_i x_i) is the orthogonal projection onto
span{u_i : i ∈ S}. Its eigenvalues are 1 (with multiplicity |S|) and 0.

**Uniqueness of the odeco form.**
- The eigenvectors of T_A (T_A(v, v) = μv, ‖v‖ = 1, μ > 0) are exactly v_S ∝ x_S (Robeva's formula, rederived from
  λ_i c_i² = μ c_i).
- The rank-one components are the |S| = 1 eigenvectors. Hence the decomposition is unique up to permutation.
- For odd order, the sign is fixed by requiring λ_i > 0.

## G5.3 LEMMA G5-2 (canonical transport and tensor perturbation; PROVED HERE)

Let s < 1. Let W: A → V be the polar factor of E_V E_A|_A, i.e. W = E_V E_A (E_A E_V E_A)^{−1/2}. W is canonical (fixed by
the two projectors) and unitary from A onto V.

1. **‖(W − I)a‖ ≤ √2·s‖a‖.** Use principal vectors a_j ∈ A, v_j ∈ V with W a_j = v_j and E_A v_j = cos θ_j a_j. The
   vectors v_j − a_j are mutually orthogonal, and ‖v_j − a_j‖² = 2 − 2cos θ_j ≤ 2 sin² θ_j ≤ 2s².
2. **‖E_A(W − I)a‖ ≤ s²‖a‖.** E_A(W − I)a_j = (cos θ_j − 1)a_j, and 1 − cos θ ≤ sin² θ.
3. **Transported tensor.** Define T̃(a, b, c) := T_V(Wa, Wb, Wc) = E_π[(Wa)(Wb)(Wc)]. W intertwines the compressed products:
   W(a ∗̃ b) = (Wa) ∗_V (Wb), where x ∗_V y = E_V(xy). So idempotents, primitivity and the spectra of multiplication
   operators correspond exactly under W.
4. **Bound: ‖T̃ − T_A‖ ≤ ε_bd := s²·(11Λ + 2C_V).**

**Proof of 4.** Write δa = (W − I)a for unit a, b, c ∈ A, and use ‖a‖_∞ ≤ Λ on A. Expand
E[(a + δa)(b + δb)(c + δc)] − E[abc] and bound each group of terms:
- **First-order terms** (3 of them, e.g. E[δa·bc]): bc ∈ A **because A is an algebra**, so E[δa·bc] = ⟨E_A δa, bc⟩, with
  absolute value ≤ s²·‖b‖_∞‖c‖ ≤ Λs². The first-order terms therefore contribute only at order s².
- **Second-order terms** (3 of them): each ≤ ‖δa‖‖δb‖‖c‖_∞ ≤ 2s²Λ.
- **Third-order term:** ≤ ‖δa‖‖δb‖‖δc‖_∞ ≤ 2s²(C_V + Λ).

Summing: 3Λs² + 6Λs² + 2s²(C_V + Λ) = s²(11Λ + 2C_V). ∎

**Remarks.**
- This is the step where the **algebra structure of A** is used: it makes the perturbation **second order** in s.
- C_V enters only through the third-order term, which is why H3′ needs only s²·C_V → 0.

**Numerics (F1).** ‖W − I‖ ≈ s, with ≤ √2·s holding at every M. The true ε scales ∝ s², e.g. 1.27e-2 → 2.3e-4 for
M = 20 … 160.

## G5.4 Source audit

| source | what it actually proves | used? |
|---|---|---|
| Anandkumar, Ge, Hsu, Kakade & Telgarsky, *JMLR* 15 (2014) | an **algorithmic, randomized** guarantee: the robust tensor power method with random restarts and deflation returns, with high probability, components close to those of a perturbed odeco tensor (a Wedin-type analogue) | **not used.** It is about the output of a randomized algorithm, not an identifiability statement about the fixed-point set; per the ruling it is not upgraded |
| Mu, Hsu & Goldfarb, *SIAM J. Matrix Anal. Appl.* 36 (2015) 1638–1659 | **deterministic** robustness of successive rank-one approximation for nearly odeco tensors | not used. It concerns SROA outputs (successive best rank-one approximations), a different object from our idempotent set. It corroborates stability |
| Auddy & Yuan, *Information and Inference* 12 (2023) 1044–1072 | **deterministic** perturbation bounds for singular values and vectors of odeco tensors, with no gap condition | not used (full text not read). Its gap-free character is consistent with Lemma G5-3, whose constants do not depend on spacing between the λ_i |
| Robeva, *SIAM J. Matrix Anal. Appl.* 37 (2016) 86–102 | the eigenvectors of odeco tensors | rederived in Lemma G5-1 |
| Newton–Kantorovich theorem (standard; Kantorovich 1948, Ortega–Rheinboldt) | existence, uniqueness and error bound for a zero near an approximate solution | **used** in Lemma G5-3 |
| Banach's theorem on symmetric multilinear forms (norm attained on the diagonal) | injective norm of a symmetric tensor = max |T(a, a, a)| | used only in the numerics |

**Verdict.** G5 does not rely on an external tensor-recovery theorem. Identifiability is proved directly (Lemma G5-3), so
the terminal SOURCE THEOREM MISMATCH does not apply. Full texts were not re-read in this environment; citations rest on
metadata and abstracts.

## G5.5 The recovery map R (deterministic; no random starts in the definition)

Input: V (any k-dimensional unital slow space) and π.

1. **Tensor.** T_V(a, b, c) = E_π[abc] on V. Equivalently, define the compressed product x ∗_V y = E_V(xy).
2. **Idempotents.** ℐ(V) = {y ∈ V : E_V(y²) = y}, the zero set of a polynomial system fixed by (V, π).
3. **Primitive idempotents.** 𝒫(V) = {y ∈ ℐ(V) : the operator z ↦ E_V(yz) has **exactly one** eigenvalue > 1/2}.
   - The operator is symmetric, so its spectrum is real.
   - 1/2 is the midpoint of the spectrum {0, 1} of an exact projection; it is a mathematical constant, not a tuned
     threshold.
   - In an exact algebra this criterion is equivalent to "rank-one multiplication operator", i.e. primitive.
4. **Rounding.** If ℐ(V) is finite and |𝒫(V)| = k = dim V, set Π(V) := the argmax partition x ↦ argmax_i ẽ_i(x) over
   𝒫(V) = {ẽ_1, …, ẽ_k}. This is defined off the **tie set** (points where the maximum is attained twice), which is
   reported. **Otherwise R refuses.**

**What R uses and does not use.**
- k = dim V comes from the certified cut.
- No K, label, coordinate, ε or objective enters.
- The numerical implementation finds ℐ(V) by Newton's method from many seeded starts. That is a solver for a fixed
  polynomial system, and the theory does not depend on it. The audit checks that the count found equals the theorem's
  2^k.

## G5.6 LEMMA G5-3 (idempotent stability for near-odeco products; PROVED HERE)

**Setting.** Let T = T_A + E on ℝ^K, where T_A is the odeco tensor of Lemma G5-1 with λ_i ∈ [1, Λ], E is symmetric, and
‖E‖ ≤ ε.

**Conclusions.**
- **(a) Local.** If ε ≤ 1/(8Λ), each exact idempotent x_S has a perturbed idempotent x̃_S with ‖x̃_S − x_S‖ ≤ (8/3)ε. It
  is the unique idempotent in the open ball of radius 1/(3Λ) around x_S.
- **(b) Global.** If ε < ε_glob := 1/(36·K^{3/2}·Λ²), then ℐ(T) = {x̃_S : S ⊂ {1, …, K}}: exactly 2^K idempotents and no
  others.
- **(c) Primitivity.** Under (b), T(x̃_S, ·, ·) has exactly |S| eigenvalues > 1/2. So 𝒫(T) = {x̃_{{i}}}: exactly K
  primitive idempotents.

**Proof of (a): Newton–Kantorovich for F(x) = T(x, x) − x.**
- F′(x_S) = D_S + 2E(x_S, ·), where D_S = diag(±1). Since ‖x_S‖ ≤ 1, β := ‖F′(x_S)^{−1}‖ ≤ 1/(1 − 2ε).
- ‖F(x_S)‖ = ‖E(x_S, x_S)‖ ≤ ε, so η₀ ≤ ε/(1 − 2ε).
- F′ is Lipschitz with constant L = 2‖T‖ ≤ 2(Λ + ε).
- With ε ≤ 1/(8Λ) and Λ ≥ 1: h = βLη₀ ≤ 2(Λ + ε)ε/(1 − 2ε)² ≤ 1/2.
- Hence there is a zero within r₋ ≤ 2η₀ ≤ (8/3)ε. It is unique within r₊ ≥ 1/(βL) ≥ (3/4)/(2·(9/8)Λ) = 1/(3Λ).

**Proof of (b): every idempotent lies near some x_S.**
- Let x = T(x, x), x ≠ 0, and write x̂ = x/‖x‖. Then ‖T_A(x̂, x̂)‖ = (Σ λ_i² x̂_i⁴)^{1/2} ≥ 1/√K. So
  1/‖x‖ = ‖T(x̂, x̂)‖ ≥ 1/√K − ε, giving ‖x‖ ≤ 2√K (since ε ≤ 1/(2√K)).
- Coordinatewise, x_i(λ_i x_i − 1) = −E(x, x)_i =: −r_i, with ‖r‖ ≤ ε‖x‖² ≤ 4Kε.
- Let t_i = min(|x_i|, |x_i − 1/λ_i|). Then t_i² ≤ |r_i|/λ_i ≤ |r_i|.
- Choosing the nearer point coordinate by coordinate gives an S with
  ‖x − x_S‖² ≤ Σ|r_i| ≤ √K‖r‖ ≤ 4K^{3/2}ε, i.e. ‖x − x_S‖ ≤ 2K^{3/4}√ε < 1/(3Λ).
- By the uniqueness in (a), x = x̃_S. The 2^K points x̃_S are distinct: the x_S are separated by at least 1/Λ, and each
  moves by at most (8/3)ε.

**Proof of (c): Weyl.**
- ‖T(x̃_S, ·, ·) − T_A(x_S, ·, ·)‖ ≤ Λ‖x̃_S − x_S‖ + ε‖x̃_S‖ ≤ (8Λ/3)ε + 2ε < 1/2 for ε < ε_glob.
- T_A(x_S, ·, ·) is a projection of rank |S|. ∎

### THEOREM G5 (METASTABLE BLIND RECOVERY; PROVED HERE, not externally reviewed)

Assume H1, H2 and H3′ with fixed K, and let s_N ≤ η_N/(g_N − η_N) → 0. Then for all N with
ε_N := s_N²(11Λ + 2C_{V_N}) < ε_glob:

1. **Rank.** The rank-K cut diverges (G4-T), and V_N is the slow space.
2. **Map defined.** The compressed product on V_N has exactly 2^K idempotents, exactly K of them primitive, so R is
   defined. This follows from Lemma G5-3 on the transported tensor, using Lemma G5-2(3).
3. **Idempotent convergence.** After relabelling, the primitive idempotents ẽ_i = W x̃_{{i}} satisfy
   ‖ẽ_i − 1_{B_i}‖ ≤ (8/3)ε_N + √2·s_N·√p_i → 0.
4. **Vanishing misclassification.** The misclassified-or-tied π-mass of Π_N = R(V_N) satisfies
   m_N ≤ (256/9)·K·ε_N² + 8·s_N² → 0.
5. **Projector convergence.** gap(A_{𝔅_N}, A_{Π_N}) ≤ 2Λ√m_N, and so gap(V_N, A_{Π_N}) ≤ s_N + 2Λ√m_N → 0.

**Proof of 4 (LEMMA G5-4, rounding).**
- If x ∈ B_j is misassigned or tied, some i ≠ j has ẽ_i(x) ≥ ẽ_j(x).
- With d_l = ẽ_l − 1_{B_l}, this gives d_i(x) − d_j(x) ≥ 1, so Σ_l d_l(x)² ≥ 1/2.
- Chebyshev: m ≤ 2Σ‖d_l‖², and (a + b)² ≤ 2a² + 2b² finishes the bound.

**Proof of 5.**
- Each recovered block has mass ≥ p_i − m > 0, so the two partition algebras have equal dimension.
- For a unit a = Σγ_i 1_{B_i}/√p_i, compare â = Σγ_i 1_{B̂_i}/√p_i. Then |a − â| ≤ 2Λ on the misassigned set and 0
  elsewhere.
- The final bound follows from the triangle inequality for the gap metric. ∎

**Sense of recovery earned:**
- projector convergence;
- L² convergence of the canonical primitive idempotents to the block indicators;
- **vanishing misclassified π-mass**, with the tie set included.

**Not earned:** eventual exact block equality.

## G5.7 Margin issue

- **Vanishing π-mass needs no margin.** The rounding lemma (G5-4) converts L² control into a mass bound by Chebyshev,
  using only the fixed gap of 1 between the indicator values 0 and 1.
- **Exact eventual recovery** needs the pointwise statement (H4): max_i ‖ẽ_i − 1_{B_i}‖_∞ < 1/2.
  - The term ‖W(x̃ − x)‖_∞ ≤ C_V·(8/3)ε is controlled under bounded C_V.
  - The term ‖(W − I)1_{B_i}‖_∞, uniform convergence of slow eigenfunctions to block constants, is **not** implied by
    H1–H3′. L²-small errors can sit on small-mass states.
  - H4 is therefore a **new supplied pointwise / margin hypothesis**. It may be provable family-by-family (e.g. from the
    eigen-equation in mean-field blocks), but that is not done here.
- **F1 diagnostic** (not a proof): the sup error is 0.178 → 0.084 → 0.040 → 0.020 for M = 20 … 160, i.e. < 1/2. This is
  consistent with the observed exact recovery.

## G5.8 Rank-2 control

**PROP G5-R2 (PROVED HERE, holds for every 2-dimensional unital V; no asymptotics).** Let V = span{1, f}, with f of
π-mean 0 and variance 1, and s = π(f³). Then ℐ(V) = {0, 1, e₊, e₋} with

  e_± = (1 ∓ s/√(s² + 4))/2 ± f/√(s² + 4).

Both e_± are primitive, and argmax(e₊, e₋) = {f > s/2}.

**Proof.**
- Write y = α + βf. Using E_V(f²) = 1 + s f, the condition E_V(y²) = y becomes α² + β² = α and 2αβ + sβ² = β.
- For β ≠ 0: α = (1 − βs)/2 and β²(s² + 4) = 1.
- **Primitivity:** on the basis (1, e), z ↦ E_V(ez) has matrix [[0, 0], [1, 1]], with eigenvalues 0 and 1. For y = 1 the
  operator is I, with eigenvalues 1 and 1.
- **Rounding:** e₊ > e₋ ⇔ e₊ > 1/2 ⇔ f > s/2. ∎

**Consequence.** For K = 2, the G5 recovery map **is identical** to PROP G4-P's rounding, and G4-P supplies a dimension-free
bound with no class hypotheses.

**Numerics (F2 coarse level, M = 5, 20, 80):** 4 idempotents; the partition equals the G4-P rounding; the closed form
matches to ≤ 4e-14.

## G5.9 Rank-3 numerical audit (F1)

Bounds are reported as computed, including the loose ones. Proof-side quantities use the hidden labels, after R has
returned.

| M | idempotents / primitive (blind) | s (DK bound) | ε true | ε_bd = s²(11Λ + 2C_V) | ε_glob | L² idempotent error (bound) | sup error | misclassified (bound with true ε / with ε_bd) |
|---|---|---|---|---|---|---|---|---|
| 5 | 8 / 3 | 0.997 (none) | 3.70 | 29.4 | 1.4e-3 | 0.60 (10.8) | 1.03 | 0.31 (vacuous) |
| 10 | 8 / 3 | 0.992 (none) | 5.32 | 33.6 | 1.2e-3 | 0.62 (15.1) | 1.09 | 0.33 (vacuous) |
| 20 | 8 / 3 | 0.145 (13.0) | 1.27e-2 | 0.575 | 1.2e-3 | 0.066 (0.170) | 0.178 | 0 (0.18 / 28) |
| 40 | 8 / 3 | 0.068 (0.86) | 3.10e-3 | 0.129 | 1.2e-3 | 0.031 (0.072) | 0.084 | 0 (0.038 / 1.45) |
| 80 | 8 / 3 | 0.033 (0.30) | 9.1e-4 | 3.0e-2 | 1.2e-3 | 0.015 (0.034) | 0.040 | 0 (8.9e-3 / 0.088) |
| 160 | 8 / 3 | 0.016 (0.13) | 2.3e-4 | 7.5e-3 | 1.2e-3 | 0.0075 (0.016) | 0.020 | 0 (2.2e-3 / 6.9e-3) |

**Reading.**
- **Every proved inequality holds.** The idempotent count found is 2^k = 8 and the primitive count is k = 3 at every M.
  Lemma G5-3 guarantees this only once ‖E‖ < ε_glob; at smaller M it holds without a guarantee.
- **The theorem's hypothesis is met late.**
  - Lemma G5-3, applied with the *true* ε, holds from M = 80 (9.1e-4 < 1.2e-3).
  - Theorem G5's stated sufficient condition (ε_bd < ε_glob) is not yet met at M = 160 (7.5e-3); extrapolating
    ε_bd ∝ s² ∝ M^{−2}, it would be met around M ≈ 400.
  - With the fully analytic s ≤ η/(g − η) (DK, about 8× looser), the guarantee starts far later.
- **The bounds are asymptotically correct in rate but loose in constants.** This is preserved as found.

## G5.10 General-conjecture firewall

- CONJECTURE G4-C (every fixed-rank almost-subalgebra is near a partition algebra) is **not used and remains OPEN**.
- G5 bypasses it: it uses s → 0 (Davis–Kahan) and the algebra structure of the **hidden** A, inside the proof only. It
  never uses Δ_alg.
- **Caveat.** R is not, by itself, a metastability detector.
  - Forced onto the diffusive ladder (β = 0, rank 3), R returns a partition (8 idempotents, 3 primitive; gap 0.577).
  - On SSEP it refuses (1515 idempotents found, 1084 "primitive" ≠ 3).
  - So refusal of non-metastable families still rests on the **certified-cut rule** of G4, not on R.

## Terminal

**METASTABLE BLIND RECOVERY PROVED.** The proof is written here and not externally reviewed.

**Sense:** vanishing misclassified π-mass, L² convergence of the canonical primitive idempotents, and projector
convergence. **Exact eventual recovery is not proved; it needs the supplied pointwise margin H4.**

**Grade: CONDITIONAL DYNAMICAL PARTITION DERIVATION — THEOREM-GRADE IN G4-T CLASS** (in the vanishing-mass sense).

**Still recorded:**
- the family is supplied;
- reversibility and irreducibility are supplied;
- the regularity hypotheses (p_min, s²C_V) are supplied;
- the divergence certificate is a proof about the family;
- no frozen GRUT residual is eliminated globally;
- the partition is a macrostate partition, so no subsystem split is derived;
- no consciousness claim is made.
