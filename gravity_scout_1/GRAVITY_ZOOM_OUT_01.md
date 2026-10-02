# GRAVITY-SCOUT-1 ZOOM-OUT 01 (after G1 d_min) — HARD STOP for owner review

**Branch:** `gravity-scout-1` (from `qft-scout-1-frozen @ 5b07573`). All frozen inputs untouched. **Not frozen.**

> **Repaired by GRAVITY REPAIR 01** (`GRAVITY_CORRECTION_LEDGER.md`). Where this file and the ledger differ, the ledger
> takes precedence. G1 is accepted provisionally after the repair.

## 1. Does gravity force d_min > 0?

**At the inclusion level: no positive universal minimum splitting length** [GR1-01, GR1-07]:

| regime | result |
|---|---|
| fixed-background LCQFT (Fewster) | d_min^incl = 0 under the split assumptions |
| semiclassical backreaction | no tested argument derives a positive universal inclusion-level scale (no stronger claim) |
| perturbative gravity (Donnelly–Giddings, O(κ), arbitrary U_ε) | **NO POSITIVE COLLAR SCALE DERIVED; SUBSYSTEM NOTION REPLACED / REPRICED** |
| fine-grained gravity (Raju, scoped examples) | **ORDINARY d_min^incl MAY BE NOT APPLICABLE / BLOCKED** |

The stronger issue is that **the ordinary split structure may be replaced or may fail** in dynamical gravity, rather than
acquire a finite d_min.

**At the preparation level: a preparation-class constraint, conditionally** [GR1-02, GR1-03, GR1-05]:
- Preparing H_cross = 0 across a collar costs energy. The exact finite-lattice minimum is α₃ N ħc A/d³, with α₃ ≈ 0.0068
  per massless scalar.
- With universal coupling, the spherical trapped-sphere criterion and premise **L_all**, horizon-free preparation needs
  **d ≳ d\*(R) ~ (N ℓ_P² R)^{1/3}**.
  - The scaling is the robust part.
  - The prefactor (8π α₃)^{1/3} is model-dependent.
- d ≥ d\* is **not excluded by this test**. That is necessary, not sufficient.
- Terminal: **PREPARATION-CLASS CONSTRAINT — CONSTRAINED-NONUNIQUE — CONDITIONAL ON L_all — HEURISTIC.**

## 2. Was the Planck scale inserted by hand?

**No.**
- d\* comes from a consistency condition, not from dimensional analysis.
- Its candidate scaling ℓ_P^{2/3} R^{1/3} gives 10⁶ – 10²⁰ ℓ_P for 10⁻¹⁵ m ≤ R ≤ Hubble (prefactor-dependent).
- **Every scale in it is supplied:** ħ, G, c, R and N. So it is a preparation-class narrowing, not a compression.
- **No novelty claim for the scaling** [GR1-09]: the (ℓ_P² R)^{1/3} form has Károlyházy / Ng–van Dam analogues.
- The routes that do return ℓ_P:
  - dimensional analysis and the Jacobson cutoff reading are RELOCATION;
  - the Bousso-type step is a **HEURISTIC APPLICATION OF A CONJECTURE** [GR1-08]. It is sub-Planckian (≈ 0.1 √N ℓ_P)
    and supplies no G1 evidence.

## 3. Did gravity narrow "the simplest new degree of freedom"?

**Only as a preparation-class constraint. The deeper effect is on the subsystem notion itself.**
- **As a property of the net** (is 𝒜(S) ⊂ 𝒜(S_d) split?), gravity derived no positive scale. In dynamical gravity the
  ordinary split structure is replaced (perturbatively) or may fail (fine-grained, scoped examples) [GR1-01].
- **As a property of preparations** (can a zero-correlation state across the collar exist without a horizon?), the
  collapse test excludes d < d\*(R), conditional on L_all. d ≥ d\* is merely not excluded [GR1-03].
- The preparation bound is ordinary semiclassical GR + QFT reasoning: **not GRUT-distinctive and not observable.**

## 4. Controls

| control | result |
|---|---|
| flat / AQFT, inclusion level | d_min^QFT = 0 (owner-stated; Fewster) |
| flat, preparation level (Part A) | E_prod(D = 1) converges (0.0528 / 0.0540 / 0.0548); sharp cut diverges (≈ 0.11/a) ⇒ d_min^prep,flat = 0 |
| G → 0 | d* → 0 |
| QNEC | non-gravitational; not used |

## 5. Scorecard

| measure | value |
|---|---|
| d_min^incl | fixed background: 0; perturbative gravity: no positive scale, subsystem notion replaced; fine-grained (scoped): ordinary split may fail / not applicable [GR1-01] |
| (H_cross = 0, d)\|_horizon-free | **PREPARATION-CLASS CONSTRAINT — CONSTRAINED-NONUNIQUE — CONDITIONAL ON L_all — HEURISTIC** [GR1-02] |
| CONDITIONALLY SELECTED | no |
| TRUE COMPRESSION | **0** |
| Distinctive predictions | **ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS** |
| Status vs canonical GRUT | auxiliary (QG-1 – QG-7, QP-5 supplied) |

## 6. Questions for the owner (answered by owner ruling; see `GRAVITY_CORRECTION_LEDGER.md`, decisions 1 – 3)

1. **Does a preparation-level narrowing count as narrowing d?** The owner's G1 question allowed "otherwise constrain the
   admissible collar scale". The verdict reports both levels and keeps them separate (GF-7).
2. **Premise L.** Should a later step try to upgrade L, the localization of product-state energy, using quantum energy
   inequalities? Or should GS-P1 stay heuristic-grade?
3. **G2 framing.** QG-7 suggests that in dynamical gravity the split itself becomes "split modulo gravitational dressing /
   charges". G2 (𝒩) would then ask whether the intermediate type-I factor survives gravitational dressing and what it
   costs. This is a scale-free question, not a d question.
4. **Larger campaign?** Gravity does narrow the collar without a hand-inserted Planck scale, but only through supplied
   scales and only at the state level. Whether that justifies a larger gravity campaign is the owner's call.

## 7. Sources

| source | grade |
|---|---|
| Fewster arXiv:1501.02682; Jacobson arXiv:1505.04753; QNEC arXiv:1509.02542; Bousso hep-th/9905177; Susskind–Uglum PRD 50, 2700 | owner-scoped (as in the charter) |
| Donnelly–Giddings PRD 98, 086006 (arXiv:1805.11095) | PRIMARY / SOURCE-TEXT VERIFIED (owner) [GR1-10] |
| Raju arXiv:2110.05470 | PRIMARY / SOURCE-TEXT VERIFIED (owner), scoped examples [GR1-10] |
| Witten arXiv:2112.12828; Chandrasekaran–Longo–Penington–Witten arXiv:2206.10780 | PRIMARY ARXIV ABSTRACT VERIFIED [GR1-10] |
| Braunstein–Pirandola–Życzkowski arXiv:0907.1190 | PRIMARY ARXIV RECORD VERIFIED; conceptual precedent only [GR1-09] |
| Hayward PRD 53, 1938 (1996) (spherical Misner–Sharp criterion) | PRIMARY JOURNAL ABSTRACT VERIFIED [GR1-10] |
| Ng–van Dam gr-qc/9906003 | SOURCE LOCATED (novelty comparison only) [GR1-09] |

Numerics are an **independent code path, not an independent reviewer.**

**HARD STOP.**
