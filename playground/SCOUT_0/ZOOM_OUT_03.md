# SCOUT_0 ZOOM-OUT REVIEW 03 — after P-02, P-01, P-24

CANDIDATE structure only.

## What was banked

1. **P-02:** the SF-1 split collapses for free bosons. The carrier is **exclusion** (hard-core bosons
   reproduce SF-1 exactly). Theorem: the IR class is a function of the occupation-edge velocity
   `v_b = ε′(k_b^∞)`.
2. **P-01:** both E-1 floor edges factor through the retained local spectral measure `μ_r`.
   - Gap → memory ports only to amenable baths; non-amenable baths give pin-free exponential memory.
   - Passivity → positivity: necessity needs every mode to be visible; asymmetry breaks it.
3. **P-24:** `Σ₀ = μ₀/2` is census-confirmed as the record's only distinctive (conditional) lock.

## 1. Repeated structure — the edge dictionary

P-02 and P-01 found the same thing from opposite ends: **the law-level invariants of quadratic
parents are local spectral data at an edge.**

| Probe | Edge | Datum | Law invariant |
|---|---|---|---|
| P-02 (sectors) | occupation edge `k_b^∞` | edge velocity `ε′(k_b^∞)` | `z = 1` vs `2`; soft-point count |
| P-01 (memory) | bottom of `supp μ_r` | edge location `λ₀` | exponential vs algebraic memory |
| L0-1a / S5-1 (record) | bottom edge at 0 | edge exponent `γ` (`μ_r ~ λ^γ`) | tail `k(τ) ~ τ^{−(1+γ)}` (Karamata) |
| P-01 (positivity) | visible spectrum | sign of visible weights | CM vs not |

**Check (`probes/z03_edge_dictionary.py`):** on the pin-free chain, `μ_r([0, x]) ∝ x^{3/2}`
(`γ = 1/2`), and the local slope of `ln k / ln τ` is `−1.493, −1.497, −1.498` on `τ = 20–160`.
This is the record's S5-1 native `t^{−3/2}` and L0-1a's algebraic tail, from the same edge datum.

> **Edge dictionary (standard: van Hove / Tauberian).** For quadratic parents, the IR law class
> seen by a retained observable is fixed by (edge location, edge exponent, edge velocity, visible
> sign) of the relevant local spectral measure. Everything else is fibre data.

This is the **quotient principle** (ZOOM_OUT_02) made concrete: the quotient is the **edge data of the
local spectral measure**.

## 2. Assumptions that keep reappearing

- **Quadratic (free) parents:** every law-class result in P-01/P-02 and SF-1 is for free models.
- **One retained site / one readout:** the local measure is site-relative (E-7/C-7).
- **A supplier of the edge:** pin vs non-amenability (memory); exclusion vs interaction (sector).
  The supplier is never earned.

## 3. Route families killed (for reconnaissance)

- **"Conservation alone forms law classes."** It needs an edge-moving mechanism (P-02).
- **"E-1's edges are properties of gap and passivity as such."** They are properties of the chain's
  spectral measure (P-01).
- **Hunting further record-wide locks by text search** is saturated (P-24). New locks would have to be
  derived, not found.

## 4. Candidate class theorems

- **Edge dictionary** (§1): standard as mathematics. In GRUT terms it would unify E-1, E-14 (S5-1
  `t^{−3/2}`), E-13 (SF-1) and P-01/P-02 as one statement about edge data. Formalizing it would turn
  several separate floor and sector certificates into corollaries.
- **Restriction principle (conjecture):** a GRUT mechanism restricts the effective law **iff** it
  fixes edge data that the supplied layers would otherwise leave free. No earned mechanism currently
  fixes edge data. The pin, `K_b`, statistics and the sector are all supplied.

## 5–6. Parameters that cancel; empirical invariants

Unchanged from ZOOM_OUT_02: K1 `Σ₀ = μ₀/2` (conditional) is the only distinctive lock; K5 (no
crossing) is on exposure. Edge exponents are amplitude-free (`τ^{−3/2}` does not care about
coupling strength) but standard (van Hove), so they are not GRUT-distinctive.

## 7. Where to search next

- **P-02b (interacting bosons; Bogoliubov `c = √(gn)`):** tests whether "density → IR velocity"
  rather than exclusion is the general carrier. It is the cleanest test of the edge dictionary's
  sector half beyond free parents. Requires lifting the interaction fence (owner/auditor call).
- **Does any earned structure fix edge data?** This is a read-only audit: do the earned predicates
  (E-1…E-22) pin any edge location, exponent or velocity independently of supplied inputs? If none
  does, the restriction-principle conjecture predicts the P-23 non-predictivity statement, and the two
  would merge into one structural theorem.
- **P-15** (nonlinear branching → Born), last in the queue as before.

## Decision

**Recommended next:** the read-only **edge-data audit** (cheap, decisive for merging P-23 with the
edge dictionary). Then P-02b if the interaction fence is lifted, then P-15. Owner/auditor to choose.
