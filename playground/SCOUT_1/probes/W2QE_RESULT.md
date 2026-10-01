# SCOUT-1 W2-QE RESULT — Born as a dynamical attractor? (de Broglie–Bohm relaxation, 2D box)

**Charter:** `PROBE_CHARTERS.md` §W2-QE, plus the **fine-grained firewall**. The firewall was pre-registered in
`8edb077`, before any output was read.

**Files:**
- `w2qe_firewall.py`: exact backward-trajectory evaluation of `f = ρ/|ψ|²` on a 96² grid.
- Logs `w2qe_firewall_{modes_1_2, modes_4_9, modes_16_25, micro, dt}.log`.
- `w2qe_bohm_relaxation.py` / `.log`: the independent forward-particle cross-check, appended below when complete.

**Labels:**
- **A. COARSE-GRAINED RELAXATION — PREMISE-PRICED**;
- KNOWN RESULT IMPORT: Valentini 1991 coarse H-theorem; Valentini–Westman 2005 — STANDARD, reproduced ✓;
- NON-DISTINCTIVE.

## 0. Verdict

> **A. COARSE-GRAINED RELAXATION — PREMISE-PRICED. Not B. Not "Born derived from dynamics without preparation".**
>
> - **Fine-grained:** `H_fine` is conserved within trajectory/grid error in every run (e.g. 16 modes: 0.778 →
>   0.787 → 0.762; 9 modes: 0.621 → 0.626 → 0.635). It is exactly conserved for the stationary single mode.
> - **Coarse-grained:** H̄ decays strongly for many-mode ψ. With 16 modes at C = 4: 0.397 → 0.057 → 0.009 at
>   t = 0, π, 4π.
> - **The decay is a coarse-graining effect.** At 16 modes and t = 4π:
>
>   | cells C | 4 | 8 | 16 | 32 |
>   |---|---|---|---|---|
>   | H̄ | 0.009 | 0.026 | 0.062 | 0.176 |
>   | fraction of H_fine retained | 1 % | 3 % | 8 % | 23 % |
>
>   Refinement recovers the disequilibrium, which has moved to small scales. As cell size → 0, H̄ → `H_fine`, a
>   conserved quantity. **Apparent relaxation does not survive refinement** (firewall item 8).

## 1. The eight firewall reports

| # | Report | Finding |
|---|---|---|
| 1 | fine-grained | `H_fine` is constant to the numerical error of the grid evaluation (norm drift ≤ 5 % at t = 4π from unresolved fine structure). 1 mode: exactly constant (1.9422). |
| 2 | coarse H̄ | Decays for ≥ 2 modes. Faster with more modes. |
| 3 | cell-size dependence | **Strong.** At t = 4π, 25 modes, C = 4 → 32: 0.007, 0.027, 0.069, 0.201. |
| 4 | mode count | At C = 8, t = 4π: 1 mode 1.087 (no change); 2 modes 0.604 → 0.277; 4 modes 0.459 → 0.115; 9 modes 0.473 → 0.059; 16 modes 0.521 → 0.026; 25 modes 0.642 → 0.027. |
| 5 | initial microstructure | A smooth start and a start with k = 24 microstructure relax alike at the coarse level (C = 8: 0.027 vs 0.026). |
| 5′ | the hidden-disequilibrium control | An `|ψ₀|²` start carrying only sub-cell microstructure (`f₀ = 1 + 0.9 sin 24x sin 24y`) has `H_fine = 0.112` (conserved: 0.108, 0.112) but **H̄ ≈ 0 at C ≤ 16 from t = 0**. Coarse graining **cannot see** fine disequilibrium. Coarse "equilibrium" therefore does not certify Born at the fine level. |
| 6 | recurrence / failure | ψ recurs exactly at t = 4π (half-integer energies), but the trajectory map does not. No coarse recurrence is seen by t = 4π. The single-mode control never relaxes (velocity field zero). |
| 7 | low-mode control | 2 modes give incomplete relaxation (C = 4: 0.43 → 0.14). |
| 8 | refinement | Fails to survive (item 3). |

**Numerical convergence:** dt = 0.001 vs 0.002 at t = 4π, 16 modes: H̄(C = 4…32) = 0.0055 / 0.026 / 0.064 /
0.168 vs 0.0087 / 0.026 / 0.062 / 0.176. Agreement is to ≤ 0.01.

## 2. The ledger (what remains supplied)

| Price | Status |
|---|---|
| the pilot-wave guidance law | supplied dynamics |
| the wavefunction (mixing: many modes, generic phases) | supplied. The mechanism fails for 1 mode and is weak for 2 |
| the coarse-graining scale | supplied. The result scales with it and vanishes under refinement |
| initial "no microstructure" assumption | supplied. Sub-cell disequilibrium is invisible (item 5′) and fine disequilibrium is never destroyed |
| the measure `|ψ|²` itself as the reference | it is the equivariant measure of the flow. The *target* is structural, and that is what the coarse H-theorem uses |

## 3. What this means

- **Same pattern as W2-ETH, one level sharper.** In both, preparation information becomes **locally/coarsely
  irrelevant** under generic (mixing / non-integrable) dynamics, while the fine-grained information is conserved
  (unitarity there, Liouville-type transport here).
  - In W2-ETH the "coarse graining" is the restriction to a subsystem.
  - In W2-QE it is the cell size, and refinement explicitly recovers the hidden disequilibrium.
- **Information accounting:** **PREPARATION FORGOTTEN DYNAMICALLY — at the coarse level only.** The supplied
  preparation is not removed; it is moved below the observation scale. The new supplied data are the mixing ψ
  and the coarse-graining scale.
- **For Q5 / P-15:** Born frequencies can be *effectively* reached at coarse resolution from non-Born starts.
  The outcome law is not *derived*. W1-G's premise-priced selection (noncontextuality / composition) remains
  the theorem-level route.

**Status: W2-QE COMPLETE — A. COARSE-GRAINED RELAXATION, PREMISE-PRICED (pilot-wave law, mixing ψ, coarse
graining, no-microstructure assumption). Fine-grained disequilibrium conserved; relaxation fails refinement.**

## Addendum — forward-particle cross-check (`w2qe_bohm_relaxation.log`) and the method cross-check (`w2qe_crosscheck.py`, `.log`)

**Forward particles** (N = 20000, C = 12, M = 4 box):

| Run | H̄ |
|---|---|
| 16 modes, t = 0 → 11 | 0.802 → 0.432 → 0.288 → 0.196 → 0.150 → 0.112 → 0.065 → 0.044 |
| same, dt = 0.001 | agrees to ≤ 0.002 |
| equilibrium start | flat at the sampling bias (0.003–0.004) |
| 16 modes, random amplitudes | 1.086 → 0.133 |
| 1 mode | 1.3608, constant |
| 2 modes (1,2)+(2,1) | **0.346 → 0.336**: essentially **no** relaxation |

The forward method reproduces the firewall pattern: many modes relax at the coarse level; the controls do not.

**The two-mode discrepancy, resolved.** In the firewall run, the 2-mode case relaxed partially (C = 8: 0.60 → 0.28).
Running forward particles on the **firewall's own box** gives, at t = π:

| C | forward | backward (dt 0.002) | backward (dt 0.0005) |
|---|---|---|---|
| 8 | 0.380 | 0.407 | 0.412 |
| 4 | 0.231 | 0.249 | 0.250 |

- **The methods agree** to within the backward evaluation's grid error. That error shows up as norm drift (~4 %)
  and as the slight upward drift of `H_fine` on an under-resolved grid.
- **The difference is physical.** Equal-energy modes give a **static** velocity field. The density then spreads only
  along closed streamlines, and how much coarse spreading occurs depends on the relative phase (1.37 vs 0.45 rad)
  and on the cell geometry.
- So low-mode relaxation is **not robust**: it ranges from none to partial depending on supplied ψ details.
  This **strengthens** the mixing-ψ price on the ledger.

**Adjudication unchanged: A. COARSE-GRAINED RELAXATION — PREMISE-PRICED.**
