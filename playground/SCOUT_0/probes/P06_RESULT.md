# SCOUT_0 W1 P-06 RESULT — the CM/non-CM boundary of the S5-1 memory kernel

**Charter:** `PROBE_CHARTERS.md` P-06 (frozen; criteria unchanged after result).
**Script:** `probes/p06_cm_boundary.py`; raw outputs `probes/p06_result_raw.json`,
`probes/p06_persistence.json`. All numbers mpmath dps=22, CM test = Bernstein
derivative-sign test (-1)^n f^(n)(t) ≥ 0, n ≤ 6, t-grid 0.08–12 log-spaced.

## Setup (recorded objects only)

- Spectral density of the S5-1 scaled-limit kernel: Lorentzian κ/(x²+κ²), κ = 1/(2√2.3) ≈ 0.32969 (`S5_OWNER_RULING_03.md:49`).
- Controlled deformation family: Lorentzian exponent p, S_p(x) ∝ (x²+κ²)^{-p}; inverse transform f_p(t) ∝ t^{p-1/2} K_{p-1/2}(κt). p=1 is the recorded kernel.

## Self-tests (all pass; tester live)

- S-1: p=1 reduces to e^{-κt} to 1.3e-23 — PASS CM (the recorded S5-1 scaled kernel is CM).
- S-2 (mutation): planted non-CM function e^{-t}cos(t) — tester FLAGGED it (evidence: 4th derivative violation at t=0.08). Tester is live; result not void.
- S-3: p=1/2 is K_0(κt), which has the exact positive-measure Laplace representation K_0(t)=∫e^{-t cosh s}ds — PASS CM (analytic sanity agreeing with numerics).

## Family sweep

| p | CM? | evidence |
|---|-----|----------|
| 0.05–1.0 (9 values) | **YES** | no violation, n≤6, t∈[0.08,12] |
| 1.05, 1.1, 1.2, 1.3 | NO **at small t only** | 6th-derivative violation at t=0.08; **no violation found at t ∈ {0.5, 1, 2}** |
| 1.5 | NO at small t | same pattern |
| 2.0 | **NO — proven** | analytic: f(t)=e^{-κt}(t+1/κ) ⇒ f''(t)=e^{-κt}(κ²t−κ) < 0 on (0, 1/κ)≈(0,3.03). Independent of any small-t artifact. |

## Verdict (under the frozen success criteria)

**Success criterion (a) MET — with a sharpened structure:**

1. **The recorded S5-1 scaled-limit kernel (p=1) is completely monotone.** This is a check the canonical record never states; it is consistent with (and does not alter) S5-1's terminal, whose failure modes were the *rate* (κg²→0 in physical time) and the *class* (underdamped, not first-order CM Level-0), not the scaled kernel's sign structure.
2. **CM is not isolated — it is an open half-neighborhood on the heavy-tail side.** The whole sub-Lorentzian family 0 < p ≤ 1 is CM at tested resolution. The CM property survives deforming the spectrum to decay MORE slowly.
3. **The boundary is crossed on the light-tail side, with a proven counterexample at p=2** (f'' < 0 on (0,3.03), exact). At finite tested resolution the last fully-CM member is p=1.0 and non-CM appears by p=1.05 — but only in the small-t region, so the honest bracket is:
   - **CM for all t>0: p ∈ (0, 1] at tested resolution (n≤6), p=2 provably non-CM.**
   - **Small-t non-CM: p ∈ [1.05, 1.5] (observed), proven mechanism plausibly the t^{p-1/2} UV singularity (t^α with α∈(0,1) fails CM at the origin by direct derivative computation — flagged, not yet fully proven for the full kernel).**

**Physical reading.** Within the Lorentzian-exponent family, the S5-1 parent sits exactly at the *terminal fully-CM member*: any sharpening of the UV tail (p>1) destroys CM, first at small t, and by p=2 destroys it on an open interval. De-sharpening (p<1, heavier spectral tails, stronger memory) keeps CM. **Serves the unselected S-7 directly:** the d_mono boundary question ("what sets the monotone boundary") now has a concrete, computed instance — in this family the monotone boundary is the Lorentzian exponent itself, and the recorded parent sits on it.

## Hostile attack on this result (immediate, per standing loop)

- **H1 (normalization):** CM is invariant under positive rescaling; normalization dropped. Sound.
- **H2 (small-t artifact):** the p∈[1.05,1.5] rows rest on one grid point (t=0.08) at order 6. Mitigation: p=2's analytic proof is artifact-free, so the *existence* of a crossing is secure; only its *location* between 1.0 and 2.0 is numerically bracketed. Recorded as such — not overstated.
- **H3 (grid/precision):** dps=22, n≤6. A violation at higher order or between grid points is not excluded for p≤1. Mitigation: the analytic anchors (S-1 exact exponential, S-3 exact K_0 representation) bound the trust region; both pass.
- **H4 (standard physics):** the CM of e^{-κ|t|} and of K_0 is textbook (Bernstein/Laplace-representation theory). The *new* content is only the family statement (open half-neighborhood + proven counterexample + the recorded parent being the terminal CM member). Classified for the literature ledger: **KNOWN-BUT-NEW-IN-GRUT** (the family boundary applied to the S5-1 kernel is not in the canonical record).
- **H5 (does this rescue S5-0?):** NO. CM of the scaled kernel does not derive the Level-0 generator — the canonical failures (rate κg²→0; OTHER-CLASS underdamped structure) are untouched. Charter forbids that reinterpretation explicitly.

## What changes in the dependency map (scout-level, not canonical)

- The supplied generator layer (A-3/S5-0) gains a neighboring fact: the parent's scaled kernel is CM and terminally so in its natural deformation family. If a future canonical campaign wants S-7 (d_mono), the Lorentzian exponent is the first computable order parameter for that boundary.

## Follow-ups spawned (for the Wave-2 queue)

- P-06b: prove the small-t mechanism (t^α singularity kills CM for α∈(0,1)) for the full kernel — closes the [1.0, 2.0] bracket to a theorem.
- P-06c: does the same boundary structure hold for a *different* deformation direction (e.g., two-Lorentzian mixtures with separated scales)? Tests whether "Lorentzian exponent = d_mono order parameter" is family-specific.

**Status: P-06 COMPLETE (success criterion (a) met, hostile-checked). Next: P-08.**
