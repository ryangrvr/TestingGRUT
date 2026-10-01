# SCOUT-1 ZOOM-OUT 03 — information accounting (after AUDIT REPAIR 01, W2-ETH, W2-FP, W2-TC1, W2-DT, W2-QE)

**Central question:** has SCOUT-1 found a selector that **removes** input information, rather than exchanging one
supplied primitive for another?

## 1. The accounting table

Categories:
- **VF** value fixed;
- **RF** relation fixed;
- **CD** class discretized;
- **PF** preparation forgotten dynamically;
- **EX** supplied primitive exchanged;
- **TR** tuning relocated;
- **NR** true net reduction of input information.

| Probe | Category | What enters | What leaves | Net |
|---|---|---|---|---|
| W1-C | — (no-go) | — | — | proves that G-moved values cannot be selected by the earned layer |
| W1-A | **RF** (Δq ≡ 2πν mod G*) | U(1), translation, ν | the momentum shift as an independent datum | 0: the output is a function of the supplied ν |
| W1-I | gauge | frame + boundary readout | K beyond `μ_r` | 0: gauge fixing |
| W1-P | **EX** | complete-passivity preparation | the `T(ω)` profile (a function) | trades a function for a *principle*. The principle is supplied |
| W1-R | **CD** + TR | criticality (tuned), c < 1, symmetry | continuum c → discrete set | relevant direction tuned; irrelevant couplings forgotten (see §2) |
| W1-S | **CD** | translation, integer d, boundary class | γ continuum → `{d/2 − 1 + k}` | 0: the inputs re-encode γ |
| W1-G | **EX** | Hilbert space + noncontextuality / composition | `h(p)` | exchange |
| W1-L | **EX** (within the declared class) | composition + local tomography | ℝ vs ℂ | exchange; the three-for-one collapse is a hypothesis |
| W1-F | **RF** (`T/H = 1/2π`) | dS background + Hadamard | T as an independent datum | 0: T is in units of the supplied H |
| W2-ETH | **PF** (local) | non-integrability (genericity) | preparation detail beyond the energy | **partial NR**: generic Hamiltonians carry less information than a prepared KMS state; fine-grained information conserved; β supplied |
| W2-QE | **PF** (coarse) | pilot-wave law, mixing ψ, coarse scale, no-microstructure assumption | coarse disequilibrium | **not NR**: the disequilibrium moves below the cell scale (fails refinement) |
| W2-FP | **TR** | exact conservation + infinite time-scale separation | `h − h_c` tuning | relocation to G-fixed zeros |
| W2-DT | **NR** | an interacting, asymptotically free sector | **the dimensionless coupling** (consumed by the anomaly) | **TRUE NR: −1 dimensionless parameter**; the remaining scale is only the unit (TC-1-unfixable anyway) |

## 2. Answer to the central question

> **Yes, once — and in a precise sense.**
>
> **The one clean net reduction (W2-DT).** Dimensional transmutation **consumes** a dimensionless parameter. The
> quantum theory has strictly fewer free dimensionless inputs than the classical one, and the leftover scale is
> exactly the datum TC-1 says no theory can fix (the unit).
>
> **The closest runner-up (W2-ETH).** It is a **partial** reduction. Generic dynamics makes preparation detail
> locally irrelevant, leaving only the conserved energy, but fine-grained information persists.
>
> **Everything else** is exchange, relation, discretization by supplied classes, or relocation.

**Unifying observation.** Both genuine reductions are **RG / coarse-graining phenomena of interacting or
non-integrable dynamics.**
- RG universality (W1-R) also forgets infinitely many irrelevant couplings, at the cost of tuning the finitely
  many relevant ones. Net reduction for IR observables = (irrelevant forgotten) − (relevant tuned).
- Dimensional transmutation is the extreme case: zero relevant dimensionless directions remain.

**And all of these mechanisms are absent from GRUT's earned floor, by construction.**
- The floor is quadratic: E-3, linearity ⇒ exact reduction.
- Quadratic (Gaussian) theories:
  - have no running and no anomaly, so no DT;
  - do not thermalize to KMS (they go to GGE, W2-ETH);
  - keep every spectral datum marginal (EDA-01: `μ_r` free; W1-S: a continuum γ).

This is the **E-3 tension**, now seen three times (W2-ETH, W2-DT, and EDA-01 in hindsight):

> **The linearity that makes GRUT's reductions exact is the same property that prevents every
> information-reducing mechanism found so far.**

## 3. The eight standing questions (brief)

1. **Reducing freedom?** Yes. Values, relations and classes are fixed. There is one true NR (W2-DT) and one
   partial (W2-ETH).
2. **Selectors = new supplied assumptions?** Mostly, yes: EX / RF / CD. The exceptions are of the genericity /
   anomaly type.
3. **Same obstruction?** Yes: the scale (TC-1, now robust to anomalies), and the E-3 linearity wall.
4. **Layers constraining each other?** W1-A × W1-R exclusion (scoped). E-3 × all NR mechanisms.
5. **New quotient?** The **RG-covariant family**: for interacting theories, the quotient becomes
   "fixed point + relevant data". The Gaussian floor has no such collapse.
6. **Prediction path?** Not yet distinctive. A GRUT-specific path would require a GRUT-specific interacting
   sector whose RG fixed point is determined by GRUT structure.
7. **Kill?** "Born from dynamics" as an *earned* route (W2-QE fails refinement). "SOC as untuned selection"
   (relocates).
8. **What would most change the view?** An **earned-compatible nonlinearity** inside the record's declared
   classes that flows to a fixed point with **no** free dimensionless data. That would be a GRUT-internal NR.

## 4. Theorem status

- **TC-1 — earned-layer non-selection under G:**
  - CANONICAL-CANDIDATE;
  - survived the record census (W2-TC1) and the anomaly hostile (W2-DT), with G read RG-covariantly;
  - external review still required.
- **TC-2 (price classification):** still a classification, not a theorem.
- **New candidate TC-3 (E-3 obstruction, conjecture).** Within the quadratic / Gaussian class:
  - no mechanism reduces the number of free quotient data;
  - every quotient datum is exactly marginal under coarse-graining.

  Grade: CONJECTURE, supported by EDA-01, W1-S, W2-ETH and W2-DT. Proof sketch target: Gaussian fixed points
  have no running of quadratic data. Spectral measures are preserved by decimation up to the Schur-complement
  map, which is a bijection on the relevant class. To be tested in Wave 3.

## 5. Wave 3 (spawned)

| ID | Question | Gate |
|---|---|---|
| **W3-NL** | Inside the record's **declared** nonlinear class (S2 C-B nonlinear drift, E-15), extended to a chain / conserved field: does any nonzero nonlinearity flow to a universal fixed point whose exponents are fixed with **no** free dimensionless data, turning the nonlinearity strength into a mere crossover scale (a classical analogue of DT)? Test: symmetric vs asymmetric exclusion (EW z = 2 → KPZ z = 3/2; width exponent β = 1/4 → 1/3), several asymmetries. | none (declared class) |
| **W3-TC3** | Test / prove TC-3: decimation of quadratic chains (Schur complement) preserves the free data of `μ_r`, and no Gaussian coarse-graining reduces it. | none |
| W3-PIN | The pin decides universality classes in several places (W1-S, W2-FP's zero, momentum conservation). Is the pin = 0 fork (E-2 Outcome A, held outside the floor) the unique entry point for any RG reduction in GRUT? | owner-level fork (analysis only) |

**Order:** W3-NL → W3-TC3 → ZOOM_OUT_04.
