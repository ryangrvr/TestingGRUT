# SCOUT-1 W1-I RESULT — does access + topology fix the one-particle contraction? (C08, C09, C33)

**Charter:** `PROBE_CHARTERS.md` §W1-I. Preregistered outcomes: CROSS-LAYER SELECTION (topology + access) /
NON-SELECTION (counterexample). Both may hold at different scopes.

**Files:** `w1i_inverse_spectral.py`, with log `w1i_inverse_spectral.log`.

**Labels:**
- KNOWN RESULT IMPORT: Jacobi/Stieltjes inverse spectral theorem; Krein string (continuum analogue); Hankel
  moment problem — STANDARD-TEXTBOOK ✓;
- NON-SELECTION THEOREM (for values);
- identifiability statement (for the chain gauge slice);
- NEW-IN-GRUT: the access-depth reading of E-15.

## 0. Verdict

> **Both outcomes, at different scopes. In substance: NON-SELECTION.**
>
> 1. **End-site readout + declared chain topology ⇒ K is fixed by `μ_r`**, up to the sign gauge `b_i → ±b_i`.
>    Reconstruction error is `4·10⁻¹²` (n = 30, random couplings).
> 2. **But "chain" is not a physical selector. It is the Lanczos gauge slice of the hidden-unitary orbit.**
>    Every (graph, root) has a chain representative with the identical `μ_r`:
>    - a star tree (uniform couplings) ≡ a chain with `b₁ = √2`;
>    - a 2D square-lattice corner ≡ a chain with `b = √2, √3, 1.826, …`.
>
>    Declaring the topology is declaring the frame. W1-C already showed that U moves everything in K except
>    `μ_r`. So identifiability here = **gauge fixing**, not selection.
> 3. **Interior readout breaks even the chain-slice identifiability.**
>    - In the finite case, 1716 distinct (non-mirror) chains share the centre `μ_r` (13 atoms split 6/7), and
>      8190 ordered splits exist in total.
>    - In the infinite case there is a continuous family (split σ as `fσ + (1−f)σ`).
>    - The readout site is a supplied primitive (S-10). This confirms E-7/E-8's site-relativity.
> 4. **Finite access depth theorem (exact).**
>    - Moments `m₀ … m_{2k}` fix exactly the first k Jacobi layers. This was verified to `10⁻¹⁴` for k = 3, 6, 10.
>    - Two chains equal on layers < 5 have `m₀ … m₁₀` equal to `10⁻¹⁵`, and they first differ at `m₁₁`.
>    - The moment route is exponentially ill-conditioned: `cond(H_k)` = 9·10⁸, 1.5·10¹⁸, 6·10²⁶, 2·10³⁴ for
>      k = 4, 8, 12, 16.
>    - The atom/weight route is stable: about ×100 amplification at depth 28.
>
> **Nothing here fixes a value of Q1 or Q3.** The positive statement is that Q3 (`μ_r`) is a **complete
> invariant** of Q1 modulo the supplied frame. That sharpens the quotient but does not make GRUT choose.

## 1. Results in detail

| Test | Result |
|---|---|
| chain, end site, random `a_i ∈ 2.3 ± 0.5`, `|b_i| ∈ [0.5, 1.5]`, n = 30 | Stieltjes reconstruction: max error 3.9·10⁻¹² (a), 1.5·10⁻¹² (|b|) |
| sign gauge `S K S`, `S₀₀ = 1` | same `μ_r` (S is a frame-preserving hidden U) |
| centre readout, chain A (unit centre couplings) vs chain B (re-split atoms) | identical centre `μ_r` and identical spectrum. `‖K_A − K_B‖_F = 3.22`; B is not the mirror of A. B's couplings include 0.012 and 0.254 (still local, passive, gapped) |
| star tree `0–1–{2,3}`, uniform couplings | 3-atom root measure = the chain `(a = 2.3³, b = 1, √2)` |
| 12×12 square-lattice corner | Lanczos chain `b₀…b₆ = 1.414, 1.732, 1.826, 1.889, 1.916, 1.939, 1.950` |
| finite moments | `m₀…m_{2k}` ⟺ layers `0…k−1` (exact) |
| noisy weights/atoms ε | median `|δb_j|`: ε = 10⁻⁸ → 9·10⁻⁹ (j = 0), 1.5·10⁻⁶ (j = 28); ε = 10⁻⁴ → 10⁻⁴ … 1.4·10⁻² |

## 2. Ten-point hostile test of the candidate selector "end readout + chain topology ⇒ K"

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | Not on the chain slice (modulo signs). Yes once the topology is not declared (§1, star / 2D corner). |
| 2 | P silently contains Q? | **Yes.** The topology declaration is the frame in which K is written. It is the gauge choice that U-invariance (W1-C) leaves free. |
| 3 | Representation-dependent? | K beyond `μ_r` is representation (frame)-dependent. `μ_r` is not. |
| 4 | Physical or gauge? | For a retained-site observer the K-entries beyond `μ_r` are **gauge**. They become physical only if other sites are also accessed (more readout = supplied access). |
| 5 | Standard? | Yes: Stieltjes/Jacobi, Krein, Gesztesy–Simon local Borg–Marchenko. |
| 6 | Parent variation? | Random chains, star tree, 2D corner, interior site: the pattern holds. |
| 7 | Composition? | Gluing two half-chains at a readout site loses identifiability (§1 centre test). |
| 8 | Coarse-graining? | Coarse access (low moments) fixes only the top layers. Access depth = half the moment order. |
| 9 | Unique or stationary? | Unique on the slice (modulo signs). |
| 10 | Boundary condition selecting? | **Yes:** "end site" is a boundary declaration. It is the readout choice that makes the slice unique. |

## 3. What this means

- **Q1 → Q3 collapse explained.** SCOUT-0 found that every lift shares K and that observables read `μ_r`.
  W1-C + W1-I give the exact reason:
  - `μ_r` is the complete invariant of `(K, e_r)` under frame changes;
  - the chain is merely its canonical representative;
  - the "one-particle contraction" quotient **is** `μ_r` for a single-site observer.
- **TP-1 extends.** The would-be selector (topology) is a supplied frame. Its "fixing" is gauge fixing. This
  is the second instance of "selection = composing supplied layers", here with the twist that the composed
  layer is a **gauge choice** rather than a physical symmetry.
- **E-15 reinterpreted (NEW-IN-GRUT placement):**
  - The C-B discriminator reads `K₁₁ = a₀`, the **depth-0** Jacobi layer.
  - The short-time expansion to order `t^{2k}` sees only k layers of the environment chain, and inferring
    deeper layers from data is exponentially ill-conditioned.
  - So any earned *short-time* criterion structurally cannot constrain the IR edge data (Q3's eligible
    target, γ). This explains, from first principles, EDA-01's split of the earned layer into
    "moment-readers" and "edge-readers".
- **Empirical consequence:** none distinctive. Inverse-spectral identifiability is a universal property of
  Jacobi operators.

**Status: W1-I COMPLETE — NON-SELECTION (values). Chain-slice identifiability = gauge fixing (supplied frame +
boundary readout). Finite-access depth theorem (exact, known). Interior-readout counterexample family.**
