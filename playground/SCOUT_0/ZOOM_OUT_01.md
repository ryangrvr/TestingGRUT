# SCOUT_0 ZOOM-OUT REVIEW 01 — after P-06, P-06b, P-08, P-09

Standing rule (auditor): every 3–4 completed probes, step back. If a local question does not
materially change the broader map, bank it and move on. Everything here is CANDIDATE structure,
never GRUT physics.

## What was banked from this neighbourhood

1. **P-06b (family theorem):** `f_p` CM ⟺ `0 < p ≤ 1` in the Lorentzian-exponent family.
2. **P-08 (structural negative):** faithfulness is not an earned selector. Making it precise exposes
   supplied choices: image class and algebra type. Λ-F fails the conditional `ℝ[x]` criterion at
   theorem grade.
3. **P-09 (structural negative, same mechanism):** earned locality (the declared commuting site net
   and `G(K_b)`) is passed by all four lifts. Every pruning notion needs a new choice: bath
   factorization, site-local `J`, or tensor vs graded locality.

## 1. Repeated mathematical structure

- **Descent invariance (lift interface).** In P-08 and P-09 every earned criterion turned out to be a
  function of data descended from the one-particle contraction `e^{−Kt}`, which every lift
  reproduces. So it is constant on the lift set. Discrimination always needs structure on the lift
  side, and the record lists all of that as supplied.
- **Positive-measure (Bernstein/Stieltjes) representations.** These appear in P-06b's proof. They
  are also the hidden content of the baseline's "positivity/passivity/accretivity" seed (accretive
  `K` is needed for Sz.-Nagy; CM is the Level-0 generator's signature).

## 2. Assumptions that keep reappearing

These are representational choices with no earned source: readout identification, image class,
algebra type, complex structure, statistics/grading, bath factorization. All six live in one place,
the step from the classical substrate to a lift. That step is S-3, already SUPPLIED.

## 3. Route family that can now be killed (for reconnaissance)

**"Find one more earned structural criterion that selects the lift."** If a candidate criterion is
built only from `K`-descended data, descent invariance says it cannot select. Any other candidate
needs a supplied choice, which the R-0 fence forbids. This is a reconnaissance kill, not a theorem.
The lift neighbourhood is closed to further scout probes unless a premise changes.

## 4. Candidate class theorems

- **CM ⟺ Lorentzian-mixture spectral density.** This is the general form of P-06b. A stationary
  kernel `f(t) = ∫cos(ωt) S(ω) dω` (even, `S ≥ 0`, finite mass) is CM on `t > 0` iff
  `S(ω) = (1/π)∫ s/(s²+ω²) μ(ds)` with `μ ≥ 0`. In words: `S` is a positive superposition of
  zero-centred Lorentzians. Proof sketch: Bernstein plus uniqueness of the Fourier transform. For
  divergent members (`p ≤ 1/2`) `μ` is only locally finite, and that case needs care. The math is
  standard (REDISCOVERED-KNOWN); it is new in GRUT.
  P-06b is the special case. Its weight `(s²−κ²)^{−p}` on `s > κ` reproduces `(ω²+κ²)^{−p}` up to
  the constant `π/(2 sin πp)`. This was checked numerically for `p ∈ {0.25, 0.5, 0.8, 0.95}`, with
  relative error ≤ 6e-12. The class theorem supersedes P-06c (two-scale mixtures are trivially
  inside it), so **P-06c is retired**.
- **Descent-invariance lemma (lift interface).** Any predicate that factors through
  lift ↦ (one-particle semigroup, site labels) is constant on the lift set. The statement is
  trivial. Its content is the classification of which recorded criteria factor this way: D-1, the
  D-6 fact, P-08 Reading 1, and P-09 L-E. It is banked as a structural observation.

## 5. Parameters that might cancel between observables

None identified yet. The CM test is scale-free (`κ` drops out under `t ↦ κt`). Any CM-based
observable built from shape alone, not amplitude, would therefore be free of the bath width. This
should be carried into P-23 as a design rule.

## 6. Empirical invariant that might survive

This is a hypothesis, to be tested in P-17/P-23. The class theorem turns "the environment supports a
Level-0 CM generator" into a line-shape statement about measurable noise spectra: the environment's
spectral density must be a zero-centred Lorentzian mixture. A resonance peak away from `ω = 0`, or
a tail decaying faster than `ω⁻²` (the `p > 1` side), rules it out. Unlike an amplitude, this test
is free of nuisance parameters.

## 7. Areas of physics to search next

- **Open-system spectral densities (P-17).** Look at Drude–Lorentz, Ohmic with exponential cutoff,
  and structured baths. A Drude-cutoff Ohmic bath has `J(ω)/ω ∝ 1/(ω²+κ²)`, which is exactly
  `p = 1`: the edge member of P-06b. The standard regularized Caldeira–Leggett bath therefore sits
  on the CM boundary. P-17 should check this first.
- **1/f noise (McWhorter).** 1/f spectra arise as Lorentzian superpositions, so they lie inside the
  CM class. This is a second physical environment family.
- **Hydrodynamic long-time tails.** S5-1's native `t^{−3/2}` branch cut is the d = 3
  Alder–Wainwright tail. That tail is a non-CM-generator mechanism with a large literature, and a
  candidate for P-18's homogenization null.

## Scope correction (auditor, applied in P-17)

The Drude/P-06b correspondence is to the **friction/memory kernel** `γ(t) = (2/π)∫ J(ω)/ω cos ωt dω`,
not automatically to `J(ω)`, to the finite-temperature quantum noise correlation, or to a measured
noise spectrum. Four objects are kept separate: `J(ω)`, `γ(t)`, `C_F(t)`, and the reduced stochastic
representation. Classical thermal `C_F = Tγ` inherits CM; the quantum symmetrized noise kernel of the
same Drude bath is **not** CM at low temperature (negative Matsubara weights; `P17_RESULT.md` §1). §6's
line-shape hypothesis is therefore restricted to the classical friction kernel.

## Decision

The lift interface is banked and closed. Next, in priority order: **P-17**, then P-18, P-23,
P-01/P-02, P-15. P-17 starts with the Drude edge observation above, which costs almost nothing and
connects P-06b to a physical bath.
