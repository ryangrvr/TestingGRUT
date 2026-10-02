# GRAVITY-SCOUT-1 ZOOM-OUT 01 (after G1 d_min) — HARD STOP for owner review

**Branch:** `gravity-scout-1` (from `qft-scout-1-frozen @ 5b07573`). All frozen inputs untouched. **Not frozen.**

## 1. Does gravity force d_min > 0?

**At the inclusion level: no.**
- QFT on curved backgrounds (Fewster) and semiclassical backreaction leave the splitting-distance infimum at 0.
- Perturbative quantum gravity (Donnelly–Giddings dressing) reprices the *split notion itself* at every d. It does not
  introduce a scale.
- **d_min^incl: NO CONSTRAINT.**

**At the preparation level: yes, conditionally.**
- Preparing H_cross = 0 across a collar costs energy. The exact lattice minimum is α₃ N ħc A/d³, with α₃ ≈ 0.0068 per
  massless scalar.
- With universal coupling and the spherical no-trapped-surface criterion, horizon-free preparation requires

      d ≥ d*(R) = (8π α₃ N ℓ_P² R)^{1/3}.

- **d_min^prep: CONSTRAINED-NONUNIQUE (conditional, heuristic grade: premise L, spherical / hoop).**

## 2. Was the Planck scale inserted by hand?

**No.**
- d* comes from a consistency condition, not from dimensional analysis.
- It scales as ℓ_P^{2/3} R^{1/3}, which is 10⁶ – 10²⁰ ℓ_P for 10⁻¹⁵ m ≤ R ≤ Hubble.
- **But every scale in it is supplied:** ħ, G, c, R and N. So it is a narrowing, not a compression.
- The routes that do return ℓ_P (dimensional analysis, the Jacobson cutoff reading, the Bousso entropy bound) are
  RELOCATION or CONJECTURE-CONDITIONED. The Bousso route is even sub-Planckian (≈ 0.1 √N ℓ_P), so it gives no admissible
  semiclassical constraint.

## 3. Did gravity narrow "the simplest new degree of freedom"?

**Partly, and only on the state side.**
- The new degree of freedom (d) has two faces:
  - as a property of the **net** (is 𝒜(S) ⊂ 𝒜(S_d) split?), gravity does not touch it;
  - as a property of **preparations** (can a zero-correlation state across the collar exist without a horizon?), gravity
    cuts off d < d*(R).
- This is the first gravitational narrowing in the program. It is ordinary semiclassical GR + QFT physics: **not
  GRUT-distinctive and not observable.**

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
| d_min^incl | **NO CONSTRAINT** |
| d_min^prep | **CONSTRAINED-NONUNIQUE (conditional / heuristic)** |
| CONDITIONALLY SELECTED | no |
| TRUE COMPRESSION | **0** |
| Distinctive predictions | **ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS** |
| Status vs canonical GRUT | auxiliary (QG-1 – QG-7, QP-5 supplied) |

## 6. Questions for the owner

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
| Donnelly–Giddings PRD 98, 086006 (arXiv:1805.11095) | SOURCE LOCATED, TEXT NOT RE-READ |
| Witten arXiv:2112.12828; Chandrasekaran–Longo–Penington–Witten arXiv:2206.10780 | SOURCE LOCATED, TEXT NOT RE-READ |
| Braunstein–Pirandola–Życzkowski PRL 110, 101301 | SOURCE LOCATED, TEXT NOT RE-READ (novelty comparison only) |
| Misner–Sharp spherical trapped-surface criterion | CITED FROM STANDARD LITERATURE |

Numerics are an **independent code path, not an independent reviewer.**

**HARD STOP.**
