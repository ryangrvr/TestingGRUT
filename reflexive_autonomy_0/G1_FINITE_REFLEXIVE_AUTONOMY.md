# RA0 · G1 — FINITE REFLEXIVE AUTONOMY (result)

> **Repaired by RA0 REPAIR 01** (RA1-01 … 06; `RA0_CORRECTION_LEDGER.md`). Where wording differs, the ledger takes
> precedence. Logs are kept as emitted.

**Evidence:**
- the propositions below (proofs given);
- the exhaustive finite checks in `g1/g1_finite.py` + `g1/g1_finite.log` (independent code path, not independent
  reviewer; illustration only).

**Status labels:**
- **KNOWN** = established Markov-chain theory;
- **PROP** = proposition proved here (elementary; not externally reviewed);
- **NUM** = numerical illustration.

**Setup.**
- Finite X, |X| = n ≥ 2.
- Row-stochastic P.
- Partition Π with blocks B_1 … B_k.
- Y_t = Π(X_t).
- **Order:** Π ⪯ Π′ iff Π refines Π′; ⊥ = discrete, ⊤ = indiscrete.
- **Versions:** F_P = Version A (present-inclusive); F⁺_P = Version B (strict-future). Note F_P(Π) = Π ∧ F⁺_P(Π).

## G1.1 Exact fixed-point characterization

**Definition (KNOWN, Kemeny–Snell 1960).** Π is **strongly lumpable** for P iff for all blocks B_i, B_j and all
x, x′ ∈ B_i:

    P(x, B_j) = P(x′, B_j).

The lumped kernel is then Q(i, j) := P(x, B_j) for any x ∈ B_i.

**PROP 1 (Version A).** Π = F_P(Π) iff Π is strongly lumpable.

*Proof.*
- (⇒) If Π is fixed, x ∼ x′ for all x, x′ in one block. In particular the laws of Y_{t+1} agree, so P(x, B_j) = P(x′, B_j).
- (⇐) Assume lumpability. By induction on k, for x ∈ B_i:

      P(Y_{t+1} = j₁, …, Y_{t+k} = j_k | X_t = x) = Q(i, j₁) Q(j₁, j₂) ⋯ Q(j_{k−1}, j_k).

  For k = 1 this is the definition. For the step, sum over y ∈ B_{j₁} of P(x, y) · P(rest | X_{t+1} = y). The inductive
  hypothesis makes the second factor depend only on j₁, giving Q(i, j₁) · (product).
- So the law of the whole future block process depends only on the block of x. With the present included (Y_t = i), x ∼ x′
  holds exactly when they are in the same block. Hence F_P(Π) = Π. ∎

**PROP 2 (Version B).** Π = F⁺_P(Π) iff (i) Π is strongly lumpable **and** (ii) the lumped kernel Q has pairwise
distinct rows.

*Proof.*
- **Fixed ⇒ (i).** Π ⪯ F⁺_P(Π) means same-block states have equal Y_{t+1}-laws, which is lumpability as in PROP 1.
- **Fixed ⇒ (ii).** F⁺_P(Π) ⪯ Π means states in different blocks have different future laws. Under (i), the future law from
  x ∈ B_i is the law of the Q-chain path (Z₁, Z₂, …) started at Z₀ = i. That law is determined by Q(i, ·), since the path
  law is Q(i, ·) followed by Q-transitions, and conversely Q(i, ·) is its first marginal. So distinct future laws ⇔
  distinct rows of Q.
- **(i) + (ii) ⇒ fixed**, by the same two facts. ∎

So Version B fixed points are the lumpable partitions whose lumped chain has no two "predictively identical" macrostates.
Fix(F⁺_P) ⊆ Fix(F_P).

**NUM.** Over all 52 partitions of 5 states, there were **0 mismatches** with PROP 1 / PROP 2. Families checked:
20 random Dirichlet chains, 20 sparse chains, all 120 bijections, the i.i.d. chain, the identity.

**KNOWN (identification).** For Markov chains, probabilistic bisimulation is the same notion as lumpability (Larsen–Skou,
Inf. Comput. 94 (1991)). Version-A fixed points are exactly the **unlabelled bisimulation partitions** of P.
- In bisimulation theory a *labelling* (an observation map) is normally supplied, and the canonical object is the coarsest
  bisimulation refining it.
- **Without a supplied labelling, the coarsest bisimulation is ⊤.** That is the formal shadow of the selector-screen
  finding: the canonical object needs the observable channel that this construction tried to avoid.

## G1.2 Lattice properties

**PROP 3 (monotonicity).** F⁺_P and F_P are order-preserving.

*Proof.*
- If Π ⪯ Π′ then Π′ = g ∘ Π for some map g. The Π′-future process is a function of the Π-future process, so equal
  Π-future laws imply equal Π′-future laws. Hence F⁺_P(Π) ⪯ F⁺_P(Π′).
- F_P = id ∧ F⁺_P is a meet of order-preserving maps, so it is order-preserving. ∎

**NUM:** 358 refining pairs × 3 chains, 0 violations for either version.

**PROP 4 (universal endpoints).**
- (a) ⊤ is fixed by both versions for every P.
- (b) ⊥ is fixed by F_P for every P.
- (c) ⊥ is fixed by F⁺_P iff P has pairwise distinct rows.

*Proof.*
- (a) The ⊤-process is constant.
- (b) Lumpability of ⊥ is vacuous.
- (c) Apply PROP 2 with Q = P. ∎

**COROLLARY 5 (Tarski gives existence, never selection).** Fix(F_P) and Fix(F⁺_P) are complete lattices (Knaster–Tarski).
- **Version A:** μF_P = ⊥ ≠ ⊤ = νF_P for **every** P on |X| ≥ 2. The least and greatest fixed points never coincide.
- **Version B:** νF⁺_P = ⊤ always. μF⁺_P = νF⁺_P happens only when ⊤ is the **only** fixed point. That is the trivial
  partition, e.g. the i.i.d. chain (NUM: |Fix_B| = 1).
- So the coincidence μ = ν can **never** single out a nontrivial partition.

This supersedes the weaker scratch wording ("μF = νF only when there is no nontrivial solution").

**NUM (degenerate cases).**
- i.i.d. chain: *every* partition is Version-A fixed (52 / 52), and Version B has only ⊤.
- A chain with two equal rows loses ⊥ in Version B.

## G1.3 Generic no-go

**PROP 6 — KNOWN / REDERIVED IN RA0** [RA1-01]. The codimension (k − 1)(n − k) matches the classical
lumpability-testing literature: the asymptotic chi-square test of a fixed k-block lumping has (k − 1)(n − k) degrees of
freedom (owner-cited; Statistics & Probability Letters, ScienceDirect S0167715203001263; not re-read here). The internal
proof is retained, and the measure-zero corollary stands. **No novelty is claimed.** Let 𝒮_n be the row-stochastic polytope, of dimension n(n − 1) inside the affine space 𝒜_n of real matrices
with unit row sums. For a partition Π with 1 < k < n blocks, the lumpable set L_Π is contained in an affine subspace of
𝒜_n of **codimension exactly (k − 1)(n − k) ≥ 1**.

Hence the set of chains with *any* nontrivial lumpable partition, ∪_{1<k<n} L_Π, is a finite union of proper affine
slices. It has **Lebesgue measure zero** in 𝒮_n.

*Proof.*
- Fix one representative r_i in each block B_i. Lumpability is equivalent to the linear equations

      ⟨P_x − P_{r_i}, 1_{B_j}⟩ = 0    for x ∈ B_i \ {r_i}, j = 1 … k − 1

  (the j = k equation follows from the row sums).
- There are Σ_i (|B_i| − 1)(k − 1) = (n − k)(k − 1) such functionals.
- They are linearly independent modulo the row-sum constraints:
  - functionals for different x involve the distinct row P_x;
  - for fixed x, the vectors 1_{B_1}, …, 1_{B_{k−1}} are independent and not in span{1} (the row-sum direction), because
    k ≥ 2 and every block is proper.
- Each functional is non-constant on 𝒜_n: move mass in row x from outside B_j into B_j.
- The intersection of a proper affine subspace with 𝒮_n has measure zero. There are finitely many partitions. ∎

**COROLLARY 7 (generic finite no-go).**
- Statement: for any probability law on 𝒮_n absolutely continuous with respect to Lebesgue measure (e.g. i.i.d.
  Dirichlet rows), almost surely Fix(F_P) = {⊥, ⊤}. Since distinct rows also hold a.s., Fix(F⁺_P) = {⊥, ⊤} too.
- So exact reflexive autonomy is **nongeneric**.
- NUM: 0 / 1000 random Dirichlet(1) chains had a nontrivial fixed point. The codimension formula was verified by rank for
  every nontrivial partition of 5 states.

**Qualifications (stated explicitly).**
- The theorem is about the full-dimensional interior. Zero entries (sparse / local dynamics), deterministic maps,
  reducible chains (closed classes are lumpable) and symmetric chains all lie on lower-dimensional strata, where nontrivial
  fixed points can exist.
- E.g. bijections on 5 states have 4 – 10 Version-A fixed points (scratch run).
- **Physical models are typically nongeneric in exactly these ways** (locality ⇒ zeros, symmetry ⇒ equivariance), so the
  no-go says *exact* autonomy needs special structure in D. It does **not** say physical systems lack it.

## G1.4 Symmetry

**PROP 8.** Let a group G act on X by permutations with P(gx, gy) = P(x, y) (equivariance). For every subgroup H ≤ G, the
H-orbit partition X/H is strongly lumpable (so it is Version-A fixed). It is Version-B fixed iff its lumped rows are
distinct.

*Proof.* For an H-orbit O and h ∈ H:

    P(hx, O) = Σ_{y∈O} P(hx, y) = Σ_{y∈O} P(x, h⁻¹y) = P(x, O),

because h⁻¹O = O. ∎

**Consequence.** Richer symmetry → more subgroups → **several, generally incomparable**, nontrivial fixed points.

**NUM.** A circulant (Z₆-equivariant) chain on 6 states has exactly two nontrivial Version-B fixed points: the Z₂-orbit
partition (0, 1, 2, 0, 1, 2) and the Z₃-orbit partition (0, 1, 0, 1, 0, 1). They are **incomparable**, so there is no
canonical choice between them.

**Lesson: symmetry → structure, not unique selection.** And the structure is the symmetry already present in D.

## G1.5 Orientation

**PROP 9 (bijective deterministic dynamics).** For a bijection T on finite X, Fix(F_{P_T}) = Fix(F_{P_{T⁻¹}}), and the
same holds for F⁺.

*Proof.*
- Lumpability of P_T means each block maps *into* a single block: T(B_i) ⊆ B_{σ(i)}.
- Since T is surjective, every block B_j = ∪_{σ(i)=j} T(B_i). So every j has a preimage under σ, and σ is a surjection,
  hence a bijection of the finite block set.
- Then |T(B_i)| = |B_i| forces T(B_i) = B_{σ(i)}. So T⁻¹(B_j) = B_{σ⁻¹(j)}, and Π is lumpable for T⁻¹. By symmetry the
  converse holds.
- The lumped kernels are the permutation matrices of σ and σ⁻¹, whose rows are distinct (relevant for Version B). ∎

**NUM:** all 120 bijections on 5 states, both versions, sets identical.

**Reversible chains.** With detailed balance π(x)P(x, y) = π(y)P(y, x), the time reversal P*(x, y) = π(y)P(y, x)/π(x)
equals P, so the fixed sets are trivially identical (NUM: P* = P).

**Non-reversible chains.**
- P* requires π. For an **irreducible** chain π is unique (Perron–Frobenius), so it is *derived from P*, not supplied.
  Otherwise P* is measure-priced.
- **NUM (one explicit example):** Π = {0,1 | 2,3,4} is lumpable for P but **not** for P*. Fix_A(P) = {⊤, Π, ⊥} while
  Fix_A(P*) = {⊤, ⊥}. A coarse-graining can be autonomous in one time direction only.
- **Reading.** This asymmetry is a property of the **supplied irreversible P**. It is the frozen "class-relative orientation
  carrier" in D (Bridge B4-CUC-2). The construction is orientation-covariant: it reports which coarse-grainings are
  autonomous in each direction and **prefers neither**. Calling the autonomous direction "forward" would be an
  interpretive identification (S1).
- **No orientation is selected.** Related literature on lumpability and reversibility (Marin & Rossi, Acta Informatica
  2017) is SOURCE LOCATED only.

## G1.6 Frozen-twin control

**PROP 10.** F_P is a function of P alone. Any structure that is not reflected in P cannot be distinguished, e.g. a
different labelling of the same P, or a different net / frame attached to the same generator.

F_P **can** rank candidate decompositions, but only through exact autonomy of their factor partitions. A factor partition
is lumpable iff that factor's dynamics receive no back-action from the rest.

**NUM (4 states = {0,1}², decompositions D1 = (a, b) and D2 = (a ⊕ b, b)):**

| dynamics | D1 factor fixed | D2 factor fixed |
|---|---|---|
| factor autonomous in D1 (here realised by a product kernel) | **True** | False |
| generic interacting | False | False |

**Reading.**
- The scratch claim "both twins necessarily receive the same answer" is **refined**:
  - the operator is identical;
  - it selects a decomposition **only when that decomposition has an exact autonomous factor** [RA1-02]. That means
    no back-action onto the coarse variable: P(A_{t+1} | A_t, B_t) = P(A_{t+1} | A_t). The rest may still depend on A,
    so one-way / triangular interacting dynamics qualify. Full factorization P_A ⊗ P_B is **not** required;
  - for interacting dynamics it accepts neither.
- Every frozen Σ witness involves interacting nets (W-Σ1: coupled isospectral oscillator nets; W-Σ2 / W-Σ3: Pareto /
  dimension ties with interaction).
- So on the frozen witnesses **finite reflexive autonomy does not select Σ**. Where it does select (autonomous factors), it
  recovers factors with no back-action. **Finite reflexive autonomy detects autonomous factors, not necessarily
  independent factors** [RA1-02]. This is a dynamics-conditional derivation (cf. frozen E-S2 / CPR
  class), not a selection among interacting twins.

## G1 terminal

> **FINITE REFLEXIVE AUTONOMY NO-GO**

1. **Fixed points reduce to known theory.** Version A = strong lumpability = unlabelled probabilistic bisimulation.
   Version B = lumpability + predictively distinct lumped states.
2. **Universal endpoints** (⊥ and ⊤ in A; ⊤ always in B) block any Tarski selection. μ = ν never yields a nontrivial
   partition.
3. **Generic finite dynamics** (measure one in the interior of the stochastic polytope) have **no** nontrivial exact
   autonomous coarse-graining. The codimension is exactly (k − 1)(n − k).
4. **Structured dynamics** (symmetry, sparsity, determinism) can have **several**, incomparable, nontrivial fixed points.
5. **The structure found reflects symmetry / conservation / autonomous (no-back-action) factors already present in D:**
   RELOCATION into D [RA1-02].
6. **Orientation is not selected.** Bijections and reversible chains are exactly symmetric. Non-reversible asymmetry is
   supplied by P.
7. **Σ twins** that are interacting, or encoded outside P, are not distinguished.

**Consciousness firewall.** Finite reflexive predictive autonomy is **insufficient** as a mathematical model of conscious
access:
- any conserved quantity or symmetry orbit qualifies;
- a memory clause M ⊆ Π either imports an internal / external split (S5) or is downstream of Π.

**No consciousness claim.**
