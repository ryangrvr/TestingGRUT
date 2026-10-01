> **AUDIT REPAIR 01:** the selected object is sharpened. The theorem-level datum is the **LSM/Oshikawa momentum
> (anomaly) shift `Δq = 2πν mod reciprocal lattice`**, together with the **obstruction to a unique symmetric gapped
> ground state** at fractional cell filling.
> - "All `q ∈ 2πν·ℤ` are soft" is **withdrawn** as a universal claim.
> - The obstruction can be discharged by gaplessness, by ground-state degeneracy / translation breaking, or (in
>   d > 1) by other anomaly realizations.
> - Only in the declared **1D gapless branch** (YOA; the tested ED parents) does it appear as a low-energy state at
>   `2k_F = 2πν`.
>
> Classification: **CROSS-LAYER SELECTION OF MOMENTUM SHIFT / ANOMALY DATUM**. The ED evidence stands for the tested
> parents.

# SCOUT-1 W1-A RESULT — LSM/Oshikawa: symmetry + sector fix the soft-point momenta (C01; C02/C03 hostile)

**Charter:** `PROBE_CHARTERS.md` §W1-A. Preregistered outcomes: CROSS-LAYER SELECTION / NOT FIXED / UNRESOLVED.

**Files:**
- `w1a_lsm_soft_points.py` — momentum-resolved ED on periodic rings; `.log`;
- `w1a_summary.py` — the diagnostic table below plus the U(1)-breaking hostile; `.log`.

**Sanity:** the sector-resolved ground energies equal the full-Fock ED to 10 digits in four cases: F, HCB, BH and
staggered F + V.

**Labels:**
- CROSS-LAYER SELECTION (supplied layers);
- KNOWN RESULT IMPORT: Lieb–Schultz–Mattis 1961; Oshikawa 2000; Yamanaka–Oshikawa–Affleck 1997 — STANDARD-TEXTBOOK, ED-confirmed ✓;
- NON-DISTINCTIVE.

## 0. Verdict

> **CROSS-LAYER SELECTION — a dimensionless Q4 component IS fixed.**
>
> - **The fix.** Given a U(1) charge × lattice translation and filling ν, the soft-point momenta are fixed at
>   `q* ∈ 2πν·ℤ (mod 2π)`. This holds independent of statistics (fermions, hard-core bosons, soft-core bosons)
>   and of interaction (V = −1, 0, 1, 3; U = 1, 4).
> - **What stays free.** The velocity varies from 1.03 to 2.47 across the same models (Q4's law-relevant
>   `v_IR`).
> - **The status of the selector.** Every input (the U(1) sector, the lattice translation group, ν) is
>   **supplied** in the GRUT record (S-1 graph, S-7 sector).
> - **The status of the fixed object.** What is fixed is a **relation between two supplied/quotient
>   components** (soft momentum ≡ 2π × filling), not a value produced from earned structure.
> - **Consistency with W1-C.** The fixed component is weight-0 and G-invariant, as W1-C requires of anything
>   fixable.

## 1. Diagnostic

Let `g(q)` be the lowest excitation energy at momentum q relative to the ground state. A **soft** q has `g·L`
bounded as L grows; a **hard** q has `g·L ∝ L`. The low-q sound mode (`q = 2π/L`) is soft in every gapless
case, so the test is whether `q = 2πν` is soft while the midpoint `q = πν` is hard.

| Model | ν | `g·L` at `q = πν` (midpoint) | `g·L` at `q = 2πν` | v (L = 12) |
|---|---|---|---|---|
| free F | 1/2 | 16 → 24 → 32 → 40 (hard) | 0 (shell degeneracy) | 1.91 |
| F, V = +1 | 1/2 | 20.7 → 51.9 (hard) | 0 | 2.47 |
| HCB | 1/2 | 20.9 → 45.8 (hard) | 12.25 → 12.52 (soft) | 1.98 |
| BH U = 1 (nmax 3) | 1/2 | 12.3 → 17.7 (growing) | 18.4 → 19.0 (soft) | 1.03 |
| BH U = 4 (nmax 3) | 1/2 | 16.3 → 25.1 (hard) | 16.5 → 17.3 (soft) | 1.61 |
| free F | 1/4 | grows (shell oscillation) | 0 / 8.8 alternating with N parity (soft) | — |
| HCB | 1/4 | 8.7 → 18.4 (hard) | 8.66 → 8.86 (soft) | — |
| BH U = 4 | 1/4 | ~5 (dilute: L below healing length) | 9.37 → 9.74 (soft) | — |
| HCB V = 1 | 1/3 | 14.9 → 26.7 (hard) | 10.66 → 10.84 (soft) | — |
| F V = 3 | 1/3 | 19.7 → 38.2 (hard) | 0 / 10.2 alternating (soft) | — |

L ranges: 8–20 (F, HCB at ν = 1/2), 8–14 (BH), 8–24 (ν = 1/4), 9–21 (ν = 1/3).

Fixed soft points: q = 2πν in every row. The BH U = 1 midpoint grows more slowly (the dilute crossover). It is
still well above the soft line at the largest L, and it is reported as is.

## 2. Hostile tests

| Hostile | Result | Reading |
|---|---|---|
| Staggered potential δ = 0.5 at ν = 1/2 (cell = 2, **integer** cell filling) | every `g·L` grows ∝ L (F: gap = 1.000 = 2δ; HCB similar) | the selection **disappears** when the supplied translation group is coarsened to one where the filling is an integer |
| Staggered δ = 0.5 at ν = 1/4 (cell filling 1/2) | `q = 2πν` stays soft (0 / 8.3 alternating) | survives. The fixed datum lives **mod the reciprocal lattice of the supplied translation group** |
| U(1) broken by pairing Δ = 0.3 (BdG) | gapped at generic μ (min E = 0.51–0.60) | without the charge layer no soft momentum is fixed |
| | gapless only at fine-tuned μ = ±2, with soft k ∈ {0, π} | the only survivors are the **fixed points of k → −k** — the W1-C fixed-point corollary again |
| Vary ν (grand canonical, Δ = 0) | `k_F = πν` moves continuously with μ | the output is fixed **only relative to** the supplied sector |

## 3. Ten-point hostile test of the selector `{U(1), T, ν} ⇒ q* = 2πν`

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | No: ten models, three statistics, interactions from −1 to 4. |
| 2 | P silently contains Q? | **Partly.** ν is an input and the output is 2πν. The selector is a fixed **map** sector → momentum. Without supplied ν nothing is fixed. |
| 3 | Representation-dependent? | Momenta are defined mod the supplied translation group (staggered test). |
| 4 | Physical or gauge? | Physical: structure-factor singularities. |
| 5 | Standard? | Yes: LSM/Oshikawa/YOA. |
| 6 | Parent variation? | Passes (§1). |
| 7 | Composition? | Decoupled chains give the union of soft sets. The anomaly is additive (standard; not run). |
| 8 | Coarse-graining? | **Fragile.** Blocking by b sites leaves a constraint only while bν ∉ ℤ. The staggered test is the explicit loss. |
| 9 | Unique or stationary? | The theorem forces gaplessness **or** degeneracy, with a low state at 2πν. It guarantees the set `2πν·ℤ` is soft; it does not exclude extra soft points. |
| 10 | Boundary condition selecting? | The Oshikawa flux argument uses a ring. The bulk statement is BC-independent (standard). |

## 4. What this means for the campaign

- **First explicit `P ⇒ Q = Q*` in SCOUT-1.** It is a genuine value fixing of a dimensionless Q4 component, and
  it survives interaction and statistics. That is stronger than anything in SCOUT-0, where statistics moved
  the law class (P-02) and interaction restored the split (P-02b).
- **But the selector composes supplied layers** (sector × translation × filling). It is a **fixed relation
  between supplied data**, structurally like K1 (a fixed relation among quotient components that the theory
  does not itself fix). It does not fix the law-relevant components `v_IR` and Luttinger K.
- **Pattern note for the theorem search.** W1-C and W1-A both say the same thing:
  - an invariant selector can only fix G-invariant data;
  - the G-invariant data actually fixed are fixed by **composing** supplied symmetry/sector layers.

  Tentative pattern **TP-1**: *"selection = anomaly/symmetry matching between supplied layers; the output is
  a function of the supplied layer data and lives modulo the supplied group."* Kept open until W1-I/W1-P.
- **GRUT relevance.** E-13 (SF-1) earned "sector ⇒ law class". W1-A adds that sector + translation fix the
  soft momenta universally, while the velocity that P-02 showed decides the law class stays free. This is
  a **NEW-IN-GRUT placement** of a known theorem.

**Status: W1-A COMPLETE — CROSS-LAYER SELECTION (supplied layers), KNOWN RESULT IMPORT, NON-DISTINCTIVE.**
