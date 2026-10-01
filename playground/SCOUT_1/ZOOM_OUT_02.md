# SCOUT-1 ZOOM-OUT 02 — Wave-1 synthesis (after W1-R, W1-S, W1-G, W1-L, W1-F)

## The Wave-1 ledger

| Probe | Target (Q) | Weight | Outcome | Principle Π | Supplied structure it needs |
|---|---|---|---|---|---|
| W1-C | all G-moved components | ≠ 0 | **NON-SELECTION THEOREM** | — | — |
| W1-A | Q4 soft momenta | 0 | **FIXED** `2πν·ℤ` | anomaly / LSM | U(1) sector, lattice translation, filling |
| W1-I | Q1 beyond `μ_r` | — | NON-SELECTION (gauge fixing) | — | frame + boundary readout |
| W1-P | Q2 profile `T(ω)` | 0 (shape) | **FIXED** to a single β | complete passivity (Kelvin) | bath Hamiltonian, preparation class |
| W1-R | Q4 central charge | 0 | **DISCRETE** (Kac) | unitarity + conformal | criticality (tuned), c < 1, symmetry class |
| W1-S | Q3 edge exponent γ | 0 | **DISCRETE** `d/2 − 1 + k` | van Hove | translation, integer d, boundary class |
| W1-G | Q5 `h(p)` | 0 | **FIXED** `h = p` | noncontextuality / no-signalling | Hilbert space, d ≥ 3 or POVMs, tensor composition |
| W1-L | Q1 number field | 0 | **FIXED** ℂ | local tomography | tensor composition |
| W1-F | Q7 `T/H` | 0 | **FIXED** `1/2π` (retired) | KMS / Hadamard | dS background, state class |

## The eight questions

**1. Are we reducing freedom?**
**Yes — six weight-0 components were fixed or discretized.** This is the first time in the program that
quotient values were made to *choose*. Every fix is conditional on supplied structure, and **none** is
produced by the earned layer E-1…E-22.

**2. Are the selectors just new supplied assumptions?**
Each selector is **(supplied structure) × (principle) × (imported theorem)**. The principles fall into **five
premise classes**, which coincide with the record's supplied-primitive categories:

| Premise class | Probes | Record primitive it prices |
|---|---|---|
| symmetry / sector layer | W1-A, W1-S | S-1 graph/lattice, S-7 sector |
| preparation / state layer | W1-P, W1-F | S-8 noise, S-11 preparation, S-13 state |
| composition layer | W1-G, W1-L | A-6 quantum class, lift prices, `p = |α|²` |
| fixed-point tuning | W1-R (and the W1-C corollary) | pin = 0 / criticality (held outside the floor) |
| gauge / frame | W1-I | S-10 readout, frame (not a selector) |

**This is the main structural output of Wave 1:** a **price table** saying which supplied layer buys which
quotient value.

**3. Is the same obstruction recurring?**
Yes, twice:
- **(a) The scale.** G-moved values are never fixed (W1-C; TC-1).
- **(b) The supplied-premise wall.** Every fixed weight-0 value needs at least one premise from the
  five classes.

**4. Have layers begun constraining one another?**
Yes, in three ways:
- **W1-A:** two supplied layers jointly fix what neither fixes alone.
- **W1-R × W1-A — exclusion:**
  - the U(1) layer that fixes soft momenta forces `c ≥ 1`, which removes unitarity rigidity;
  - so on one 1D parent the two selectors cannot both act;
  - **a cross-layer constraint that *removes* selection power**.
- **W1-G + W1-L — sharing:** one composition premise fixes three components (number field, `h(p)`, and the
  `p = |α|²` chain of P-15).

**5. New quotient objects?**
- the affine shape class of `μ_r`;
- access depth (2k moments ⟺ k layers);
- the Jacobi tail (γ is encoded in the `1/n²` tail of the Lanczos couplings, an infinitely deep datum);
- the β-profile shape.

**6. Prediction path worth prioritizing?**
**None distinctive.** Every fixed value is universal standard physics (LSM, Kac, van Hove, Born, ℂ, KMS) and is
shared by every theory with the same supplied layers. The empirical standard fails at "a standard alternative
differs" in every case. The only route that touches GRUT-specific content (Q6, gravity) is the
**positivity-bounds** candidate (C17/C18). It is literature-gated (below).

**7. Kill which route families?**
- **KILL:**
  - scale selection (W1-C);
  - topology/access selection (W1-I);
  - the KMS/FDT prediction route (W1-F);
  - "earned layer alone fixes a weight-0 IR datum" (W1-S's continuum).
- **KEEP:** the composition subfamily, as an owner option; fixed-point selection (Wave 2).

**8. What would most change our view?**
A **dynamical** principle that replaces one of the supplied *prices* — tuning or preparation — with an
attractor:
- **criticality without tuning** (self-organized criticality / RG attractor): removes W1-R's "tuned" price;
- **Born without supplied equilibrium** (relaxation to quantum equilibrium): removes the preparation price
  in de Broglie–Bohm-type parents of W1-G;
- **KMS without supplied passivity** (thermalization / ETH): removes W1-P's preparation price.

If dynamics can remove the price, GRUT-style dynamical structure could *earn* a value. If the attractor itself
needs supplied structure (a slow-drive separation, a coarse-graining, an integrability-breaking condition), the
wall is confirmed at the dynamical level too.

## Theorem recognition

**TC-1** (ZOOM_OUT_01; G-moved components) still stands. Hostile search W2-TC1 is queued.

**TC-2 candidate — the price theorem (empirical classification, not a proof).** Across nine probes, every
selector that fixes a weight-0 quotient component takes the form

    (supplied structure from one of the five premise classes) × (principle) × (imported theorem).

No selector draws only on the earned layer. The **formal** version — "any selector satisfying P1…Pn is
equivalent to supplying Q" — is **proved only for G-moved Q (TC-1)**. For weight-0 Q it is a classification
over the tested selectors and is **not** promoted. Promotion would need either a proof that the earned layer
cannot distinguish members of each premise class, or the Wave-2 dynamical attack failing.

## Wave 2 (spawned; highest-information first)

| ID | Question | Why now | Gate |
|---|---|---|---|
| **W2-QE** | Does de Broglie–Bohm relaxation drive a non-equilibrium `ρ ≠ |ψ|²` to Born **without** supplied equilibrium? What does it need (mode number, coarse-graining)? | attacks the preparation price at Q5 directly; computable (2D box, many modes) | none |
| **W2-FP** | Can a slowly driven, locally relaxing dynamics (sandpile / SOC) reach criticality with **no tuned parameter**? Is a separation of time scales then the hidden price? | attacks the tuning price (W1-R, W1-C corollary) | none |
| **W2-ETH** | Does a closed non-integrable chain thermalize a subsystem to KMS (removing W1-P's preparation price)? What does it need (non-integrability, energy-shell choice)? | attacks the preparation price at Q2 | none |
| W2-TC1 | hostile search of E-11/E-22/gravity inputs for an earned absolute reference that would break TC-1 | theorem hygiene | none |
| W2-POS | positivity/causality bounds on Horndeski α-functions: do they exclude one K1 side (Q6 sign)? | the only Q6-touching route | **literature-gated** (primary texts unreachable); could be attempted from first principles later |
| C29 | TKNN / Chern quantization | expected TP-1 (U(1) + gap + topology) | low priority |

**Order:** W2-QE → W2-ETH → W2-FP → ZOOM_OUT_03 → W2-TC1. W2-POS stays parked unless the owner provides access
or the texts.
