# RA0 FINAL HANDOFF — GRUT REFLEXIVE AUTONOMY CONSTRUCTION 0

**Branch:** `grut-reflexive-autonomy-0`, based on `grut-selector-screen-1-frozen @ a119454`. To be frozen as
`grut-reflexive-autonomy-0-frozen`.

**Status of results.** New theory construction. It is non-canonical, it is not a GRUT prediction, and it does not update
canonical GRUT.
- All proofs: **INTERNALLY PROVED / NOT EXTERNALLY REVIEWED**.
- All numerics: independent code path, not independent reviewer, illustration only.
- `RA0_CORRECTION_LEDGER.md` (RA1, RA2, RA3) takes precedence over earlier wording. Logs are kept as emitted.

## The construction, gate by gate

### G1 — FINITE REFLEXIVE AUTONOMY NO-GO (`G1_FINITE_REFLEXIVE_AUTONOMY.md`)

- **Known theory.** Fixed points of the reflexive access operator are the strongly lumpable partitions (Kemeny–Snell),
  i.e. probabilistic bisimulations (Larsen–Skou). The future-only version adds the condition that the lumped kernel has
  distinct rows.
- **Universal endpoints** (the trivial partition ⊤ and the discrete partition ⊥) block any Tarski-style selection.
- **Generic no-go.** Nontrivial exact autonomy has measure zero: codimension (k − 1)(n − k), KNOWN / REDERIVED [RA1-01].
- **Non-uniqueness.** Symmetry creates multiple fixed points (orbit partitions).
- **Autonomy means no back-action** onto the coarse variable, not independence [RA1-02].
- **No orientation selector.** Orientation is covariant: Fix(T) = Fix(T⁻¹), and reversible chains satisfy P* = P.
- **No frozen Σ selector.**

### G2 — CONDITIONAL MACROVARIABLE DERIVATION (`G2_ASYMPTOTIC_AUTONOMY.md`)

- **Emergence.** Hydrodynamic-type structure emerges, e.g. SSEP density modes and an asymptotically autonomous local
  density.
- **No ε is needed.** Every criterion is a limit.
- **Still supplied:** the macro criterion and the family / scaling.
- **No uniqueness.** There is no unique asymptotic slow subspace [RA1-05]. What is canonical is the finite-size lowest
  eigenspace [RA1-04].

### G3 — CONDITIONAL DYNAMICAL HIERARCHY DERIVATION (`G3_SPECTRAL_HIERARCHY.md`)

- **Projectors.** Diverging adjacent spectral ratios yield basis-free spectral projectors. Nested separations yield
  filtrations.
- **Canonicality.** At finite N it is earned. Across N it is conditional on rank sequences being spectrally determined
  [RA2-02].
- **SSEP** is a multiscale continuum (adjacent ratios ≤ 4 in the slow window; counting profile 2⌊√x⌋ [RA2-03]), so it
  has no canonical cut.
- **Family / scaling are priced.** The positive controls detect supplied separation; they do not create it [RA2-04].

### G4 — CONDITIONAL DYNAMICAL PARTITION DERIVATION (`G4_METASTABLE_ALGEBRA.md`)

- **Finite algebra ↔ partition.** A finite unital commutative function algebra is the same thing as a partition (Gelfand,
  finite case).
- **Nearly-uncoupled class G4-T** (Weyl + Davis–Kahan):
  - a certified diverging cut;
  - slow space → partition algebra;
  - asymptotic product closure under the regularity hypothesis H3.
- **Rank 2 proved directly** (PROP G4-P, dimension-free).
- **False positives** (SSEP, independent particles, ladder β ≤ 2) are refused by the certified-cut rule.
  - Algebraization ≠ diverging cut: the ladder at β = 2 algebraizes, but is refused.
- **CONJECTURE G4-C** (a general almost-algebra lies near a partition algebra) **remains OPEN**.

### G5 — METASTABLE BLIND RECOVERY PROVED (`G5_METASTABLE_RECOVERY_THEOREM.md`)

**The chain**, under the priced G5 hypotheses:

  supplied metastable family
  → canonical diverging cut
  → canonical slow spectral space
  → near-odeco cubic product (‖T̃ − T_A‖ ≤ s²(11Λ + 2C_V))
  → canonical primitive idempotents (Lemma G5-3)
  → macro-partition (argmax).

**Earned:**
- projector convergence;
- L² convergence of the idempotents;
- **vanishing misclassified π-mass**;
- **asymptotic uniqueness within the G5 class** (Corollary G5-U);
- rank 2 reduces exactly to G4-P.

**Not earned:**
- **no exact-recovery claim** without the pointwise margin H4;
- recovery is **blind**, but family **certification is not** [RA3-03].

**Finite-M audit** [RA3-01, RA3-02]:
- ε_num ≤ ε_true ≤ ε_bd.
- The rigorous condition is not reached through M = 160; M ≈ 400 is an extrapolation.
- The numerical idempotent enumeration is evidence, not proof of completeness.

## Information accounting

**EARNED.** Within the certified G5 metastable class, a macrostate partition can be recovered from dynamics alone,
without supplying that partition to the recovery map. That partition is asymptotically unique within the class.

**STILL SUPPLIED:**
- the system family;
- reversibility;
- irreducibility;
- the metastable scale-separation hypothesis and its certificate;
- the fixed-rank asymptotic class;
- p_min regularity;
- s²C_V regularity.

**NOT EARNED:**
- generic partition selection;
- subsystem factorization Σ;
- frozen A_partition, globally;
- orientation;
- A_resolution, globally;
- consciousness, awareness or selfhood;
- definite quantum outcomes;
- Born weights.

**TRUE COMPRESSION relative to frozen GRUT: 0.** No frozen witness pair or canonical residual item has been eliminated
globally.

**New constructive result:** CONDITIONAL DYNAMICAL PARTITION DERIVATION IN A CERTIFIED METASTABLE CLASS. This is **not** a
GRUT prediction.

## Consciousness entry ruling

**MATHEMATICALLY ELIGIBLE AS A NEW HYPOTHESIS SECTOR — NOT DERIVED.**

**Reason.** In the certified metastable class, RA0 supplies an endogenous macro-algebra and macro-partition **before any
observer is introduced**.

**What a future hypothesis (opened elsewhere, never on RA0) may ask.** Whether memory, self-reference, internal predictive
access or recursive modelling can be defined entirely inside such an endogenous algebra, without supplying a subject or
system boundary.

**Limits on that eligibility:**
- It is **not** evidence that the construction describes consciousness.
- Consciousness must not be used to choose the partition. That is A_partition's former job, and in this class dynamics
  now does it.

## Constraints honoured throughout

- No modification of any frozen parent branch or of canonical GRUT (`ryangrvr/GRUT-RAI @ b935099…`).
- No PR, no merge, no canonical update.
- No quantum, collapse, Born-rule or consciousness work on RA0.
