# GRAVITY-SCOUT-1 — GRAVITY / CAUSAL-STRUCTURE ENVELOPE (charter)

**Status:** non-canonical playground campaign, authorized by owner ruling ("AUTHORIZE GRAVITY / CAUSAL-STRUCTURE CAMPAIGN",
this thread, after the QFT-SCOUT-1 freeze).

| field | value |
|---|---|
| Branch | `gravity-scout-1`, created from `qft-scout-1-frozen @ 5b07573bf45650d332df14258554cc8428a0706c` |
| Frozen inputs (read-only) | `qft-scout-1-frozen @ 5b07573`; `grut-bridge-1-frozen @ 753b90a`; `scout-2-reviewed @ 1f08ae1`; canonical GRUT `GRUT-RAI/master-w25bu9 @ b935099` |
| Not modified | scout-0, scout-1, scout-2, scout-2-review, scout-2-reviewed, grut-bridge-1, grut-bridge-1-frozen, qft-scout-1, qft-scout-1-frozen, GRUT-RAI |

## 0. What this campaign opens

**This campaign explicitly opens the gravity envelope.**
- Bridge-1 ran under a "do not reopen gravity" rule, and QFT-SCOUT-1 inherited it.
- The canonical record holds **GRAVITY_UNDERDETERMINED**. The gravity quadruple (coupling, spin-2, universal reach, causal
  cone) is **supplied** in canonical Ledger A-8.
- Opening it here means that **gravitational inputs are consumed as priced premises (QG-1 …)**. It does **not** mean
  canonical GRUT has selected its gravitational layer.
- Every result stays **auxiliary to canonical GRUT** unless the gravitational layer itself is selected. The QFT lift
  (QP-5) is also still supplied.

## 1. Central question and order

The QFT-SCOUT-1 handoff left **H_cross CONSTRAINED-NONUNIQUE**: a product (zero-correlation) boundary is admissible only
with a **split buffer (d, 𝒩)**, never across a sharp non-type-I cut in the normal class. The buffer data (d, 𝒩) is the
simplest new degree of freedom that QFT-SCOUT-1 created.

| gate | target | question | stop |
|---|---|---|---|
| **G1** | **d_min** | **Does gravity force a nonzero minimum splitting distance d_min, or otherwise constrain the admissible collar scale?** | **HARD STOP → `GRAVITY_ZOOM_OUT_01.md`** |
| **G2** | **subsystem structure** (𝒩 and its replacement) | **What replaces QFT split / type-I subsystem structure when gravitational gauge constraints are imposed, and does the replacement remove any frozen residual information or only relocate it into charges, dressing, boundary observables, or resolution?** (owner reframing, with GRAVITY REPAIR 01) | **HARD STOP → `GRAVITY_ZOOM_OUT_02.md`** |
| G3 | orientation | — | not authorized (owner review after G2) |
| G4 | modular pattern / dimension | — | not authorized |

Owner ruling: "hard stop after d_min. If gravity cannot narrow that simplest new degree of freedom without inserting the
Planck scale by hand, we should know before building a larger gravity campaign."

## 2. Definitions (fixed before any G1 numbers were read)

**Splitting distance** (Fewster, locally covariant split analysis):
- An inclusion 𝒜(S) ⊂ 𝒜(S_d) is **split** if there is a type-I factor 𝒩 with 𝒜(S) ⊂ 𝒩 ⊂ 𝒜(S_d).
- d(S) := inf { d > 0 : 𝒜(S) ⊂ 𝒜(S_d) is split }, where S_d is S thickened by a collar of width d.

**Two levels of "minimum collar"** are distinguished. Conflating them is a firewall violation (GF-7):

| level | symbol | meaning |
|---|---|---|
| inclusion level | **d_min^incl** | the infimum of collar widths for which the **inclusion** is split. A property of the net; no state is prepared |
| preparation level | **d_min^prep** | the infimum of collar widths for which a **normal product state across the collar** (H_cross = 0) can be **physically prepared** under the stated gravitational premises (e.g. without forming a trapped surface) |

Split ⇔ a normal product state exists (Buchholz–Summers / Summers, QFT-SCOUT-1 G1). Existence is not preparability:
gravity can act on the energy that a product state carries even if it cannot act on the inclusion.

**Flat / AQFT control** (owner-stated):
**d_min^QFT = 0 at the level of the split-distance infimum, under the stated split-property assumptions.**
- This does **not** mean "QFT requires a finite buffer scale".
- The sharp boundary d = 0 admits no normal product state (QFT-SCOUT-1 G1-P1).
- Arbitrarily small d > 0 may admit split inclusions.
- Under timeslice and local quasi-equivalence, sufficiently strong distal splitting gives a vanishing splitting distance
  for balls, and uniform distal splitting gives the full split property (Fewster).

## 3. Priced gravitational premises (QG)

| ID | premise | status |
|---|---|---|
| **QG-1** | Newton coupling G | supplied (canonical Ledger A-8) |
| **QG-2** | dynamical spin-2 metric; Einstein equations | supplied (A-8) |
| **QG-3** | universal coupling: all energy gravitates | supplied (A-8) |
| **QG-4** | Lorentzian causal cone; globally hyperbolic backgrounds where used | supplied (A-8); not derived (SCOUT did not derive Lorentz structure) |
| **QG-5** | semiclassical coupling G_μν = 8πG⟨T_μν⟩_ω (state-dependent backreaction), with its validity regime (curvature scales ≫ ℓ_P) | supplied hypothesis |
| **QG-6** | an admissibility criterion for preparations: no trapped surface. Spherical symmetry: Misner–Sharp 2Gm(r)/(c²r) < 1. General shapes: the hoop conjecture | spherical case: classical-GR result, imported; general case: **conjecture** |
| **QG-7** | perturbative quantum-gravity gauge constraints (gravitational dressing, Gauss law at O(G)) | supplied hypothesis; sources located |

ħ enters through QP-5 (the supplied quantum lift) and c through QG-4. **ℓ_P = √(ħG/c³) is therefore assembled entirely from
supplied inputs.**

## 4. Success ladder and G1 terminals

The ladder is unchanged: FORBIDDEN < CONSTRAINED-NONUNIQUE < CONDITIONALLY SELECTED < TRUE COMPRESSION.

Preregistered G1 terminals (owner):

| terminal | meaning |
|---|---|
| **NO CONSTRAINT** | arbitrarily small positive collars remain admissible |
| **RELOCATION** | d is identified with ℓ_P or with another already supplied scale |
| **CONSTRAINED-NONUNIQUE** | gravity proves d ≥ d_min > 0, or rules out part of the collar family |
| **CONDITIONALLY SELECTED** | a unique d follows from explicit gravitational and state premises |
| **TRUE COMPRESSION** | the split scale is fixed without importing equivalent scale information |

A terminal is reported **separately for each level** (incl / prep).

## 5. Smuggling firewalls

| ID | firewall |
|---|---|
| **GF-1** | **Planck insertion.** d = ℓ_P = √(ħG/c³) obtained by dimensional analysis alone is **RELOCATION / SCALE INSERTION**: ħ, G and c are all supplied |
| **GF-2** | **Susskind–Uglum.** The horizon-entropy UV divergence renormalizes G (PRD 50, 2700 (1994)). A UV divergence is therefore not automatically evidence of a literal cutoff d ~ ℓ_P |
| **GF-3** | **Bousso** (hep-th/9905177) is a **conjecture**, not theorem-grade input. A route that uses it is graded **CONJECTURE-CONDITIONED** |
| **GF-4** | **QNEC** (arXiv:1509.02542) does not involve gravity. It is a **non-gravitational control** only and can never count as gravity selecting d |
| **GF-5** | **Jacobson entanglement equilibrium** (arXiv:1505.04753) assumes small-ball geometry and the equilibrium condition. It relates gravity and entanglement but does **not** select a split buffer |
| **GF-6** | **Fewster** (arXiv:1501.02682; White Rose eprint 84764) transports and preserves split structure. It supplies **no** universal gravitational minimum collar |
| **GF-7** | **Level separation.** An inclusion-level statement and a preparation-level statement are never substituted for each other |
| **GF-8** | **Real progress** requires a bound of the form d_min ≥ F(G, ħ, state, curvature, …) from a **physical consistency condition**, or a theorem excluding split preparations below some scale. Even then it is **CONSTRAINED-NONUNIQUE**, not compression |
| **GF-9** | **Non-gravitational control at every gate.** The flat AQFT result (QFT-SCOUT-1) and the G → 0 limit are run alongside every gravitational claim |
| **GF-10** | **Auxiliary status.** Results are auxiliary to canonical GRUT unless the gravitational layer is selected |
| **GF-11** | **Lattice controls** are illustrations of continuum statements. Numerics are an independent code path, not an independent reviewer |

## 6. Source grading

The grades are as in QFT-SCOUT-1:
- **PRIMARY / SOURCE-TEXT VERIFIED** (by the owner);
- **SOURCE LOCATED, TEXT NOT RE-READ** (direct fetches from this environment are blocked by the egress proxy);
- **CITED FROM STANDARD LITERATURE**;
- **GRAVITY-SCOUT PROPOSITION**, proved or sketched here with its review status stated;
- **HEURISTIC**, an estimate whose premises are listed and whose grade is never upgraded silently.

## 7. Payoff bar

The bar is carried unchanged: forced · residual-input invariant · convention invariant · non-trivial · distinctive vs the
ordinary comparison class (here: **semiclassical GR + QFT with the same premises**) · observable · falsifiable.

Baseline: **ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS.**

## 8. Rules carried

- No PR unless asked. No merge. No canonical edit.
- Commit trailers as before. No model identifiers in commits or artifacts.
- Hard stop at each zoom-out for owner review.
