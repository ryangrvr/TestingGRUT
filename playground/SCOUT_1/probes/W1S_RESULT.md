# SCOUT-1 W1-S RESULT — does spectral dimension fix the edge exponent? (C20)

**Charter:** `PROBE_CHARTERS.md` §W1-S. Preregistered outcomes: CROSS-LAYER RESTRICTION (`γ ∈` discrete set;
dimension supplied) / NOT FIXED.

**Files:** `w1s_spectral_dimension.py`, with log `w1s_spectral_dimension.log`.

**Method:**
- exact return probabilities via Bessel functions (image method for Dirichlet/Neumann boundaries);
- threshold Green's functions for point defects;
- `expm_multiply` on radial chains (n = 20000, τ up to 6400).

**Labels:**
- CROSS-LAYER RESTRICTION (supplied layers) **and** NON-SELECTION (earned layer alone);
- KNOWN RESULT IMPORT: van Hove; image method; threshold resonances / Levinson; Bessel processes — STANDARD-TEXTBOOK ✓;
- NEW-IN-GRUT: the E-14 boundary-class reading.

## 0. Verdict

> **Both, and the split is exact.**
>
> 1. **Translation invariance + integer dimension d + boundary class ⇒ γ is discrete:**
>    `γ = d/2 − 1 + k`, where k is the number of Dirichlet-constrained directions at the retained site.
>    - This holds to 4 digits for d = 1, 2, 3 (bulk, Neumann face, Dirichlet face, Dirichlet corner) and is
>      **independent of coupling strength and anisotropy** (`J_y/J_x = 0.2`; `1 : 0.5 : 0.1`).
>    - **Point defects:**
>      - d = 1: relevant — any V ≠ 0 moves γ from −½ to +½;
>      - d = 2: marginal — γ = 0 with `1/log²` corrections (consistent with SCOUT-0 EDA W-A's drifting −1.16);
>      - d = 3: irrelevant unless tuned to the threshold resonance (γ → −½).
>
>      So the discrete set is `{d/2 − 1 + k}`, plus the tuned resonance values and logarithmic marginal cases.
> 2. **The earned predicates alone (local, passive, bounded couplings, gapped with a pin) allow a continuum.**
>    - Radial chains of real dimension D realize `γ = D/2 − 1` for **any real D ≥ 1**:
>
>      | D | 1.5 | 2.5 | 3.0 | 3.7 |
>      |---|---|---|---|---|
>      | measured γ (converging) | −0.248 | 0.261 | 0.516 | 0.875 |
>      | target `D/2 − 1` | −0.25 | 0.25 | 0.50 | 0.85 |
>
>    - Every member is nearest-neighbour, symmetric positive-semidefinite, with couplings in [0.39, 1].
>    - This upgrades SCOUT-0 EDA-01's two-point witness to an **explicit continuum family**.
>
> **The discreteness comes entirely from supplied layers.** Translation invariance makes the dimension an
> integer, and the boundary class sets k.

## 1. Results

**Hypercubic sites.** Local slopes of `k₀(τ)` at τ = 50–3200, against the prediction `−(γ+1)`:

| Site | Slopes | γ |
|---|---|---|
| d = 1 bulk | −0.5007, −0.5002, −0.5000 | −½ |
| d = 1 Dirichlet end | −1.4980 → −1.4999 | +½ |
| d = 1 Neumann end | −0.4993 → −0.5000 | −½ |
| d = 2 bulk (iso / aniso 0.2) | −1.0001 / −1.0003 | 0 |
| d = 2 Neumann edge | −1.0000 | 0 |
| d = 2 Dirichlet edge | −1.9999 | 1 |
| d = 2 Dirichlet corner | −2.9997 | 2 |
| d = 3 bulk (iso / aniso 1 : 0.5 : 0.1) | −1.5001 / −1.5006 | ½ |
| d = 3 Dirichlet face | −2.5000 | 3/2 |
| d = 3 Dirichlet corner | −4.4996 | 7/2 |

**Threshold Green's function `|G₀|`** at s = 10⁻² … 10⁻⁵ below the edge:

| d | `|G₀|` | Behaviour |
|---|---|---|
| 1 | 4.99, 15.8, 50.0, 158.1 | `∝ s^{−1/2}` |
| 2 | 0.64, 0.83, 1.01, 1.19 | `∝ log` |
| 3 | 0.245 → 0.2525 | finite: the Watson value `G_c = 0.2527` |

**Radial chains** (τ = 100 → 6400): the slopes converge monotonically to `−D/2` for D = 1, 1.5, 2, 2.5, 3,
3.7. Corrections are slow (log / power) for D > 2, as expected near a non-integer-dimension Bessel process.

## 2. Ten-point hostile test of "spectral dimension ⇒ γ"

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | With P = "local, passive, gapped": **yes**, a continuum. With P = "+ translation invariance + d + boundary class": no. |
| 2 | P silently contains Q? | **Yes, largely.** The boundary class (Dirichlet count) and d are γ, re-encoded (`γ = d/2 − 1 + k`). |
| 3 | Representation-dependent? | No. γ is a property of `μ_r` (W1-C eligible, G-invariant). |
| 4 | Physical or gauge? | Physical (memory-tail exponent). |
| 5 | Standard? | Yes. |
| 6 | Parent variation? | Coupling strength and anisotropy do not move γ. Boundary and defect class do. |
| 7 | Composition? | Product geometries add Dirichlet counts and dimensions (that is how the table is built). |
| 8 | Coarse-graining? | γ is an IR/RG-fixed datum. Marginal defects in d = 2 produce logs, i.e. the coarse-graining flow is visible. |
| 9 | Unique or stationary? | Unique per (d, k) class. Resonant cases need tuning. |
| 10 | Boundary condition selecting? | **Yes, literally:** a Dirichlet vs Neumann end moves γ by 1. |

## 3. What this means

- **Q3's eligible target (γ) is not earned-fixed.** It is discretized only by supplied translation invariance +
  dimension + boundary class.
- **Synthesis with W1-I.** In Lanczos form, every radial chain has couplings `b_n → 1` with a `1/n²` tail
  whose coefficient carries D. So γ is encoded in the **asymptotic tail** of the Jacobi parameters: an
  infinitely deep datum (standard: Jacobi-parameter tails control edge behaviour).
  - By the W1-I access-depth theorem, **no finite-order short-time criterion can constrain γ**.
  - This is a structural reason for EDA-01's moment-reader/edge-reader split.
- **E-14 reading (NEW-IN-GRUT placement).**
  - S5-1's native `t^{−3/2}` branch cut is the **d = 1 Dirichlet-end** value (γ = ½).
  - The record's parent has end diagonal `2 + pin`, the Dirichlet type, with `K₁₁ = 2.3`.
  - A Neumann end (diagonal `1 + pin`, `K₁₁ = 1.3`) gives the same chain and the same earned predicates, but
    a `t^{−1/2}` tail and a different κ.
  - So the exponent 3/2 in E-14 is a **supplied boundary-class datum**, not an earned one. This is
    consistent with EDA-01's READ classification and makes its upstream trace explicit.
- **TP-1 instance 5.** Another case of "selection = supplied symmetry (translation) × supplied class
  (dimension, boundary) × imported theorem (van Hove)". The rigidity (discreteness) is real, and it is
  bought entirely with supplied layers.

**Status: W1-S COMPLETE — CROSS-LAYER RESTRICTION (γ = d/2 − 1 + k; all inputs supplied) + NON-SELECTION
(earned predicates allow a continuum γ = D/2 − 1, explicit radial-chain family). E-14's 3/2 traced to the
supplied Dirichlet end.**
