# SCOUT-1 QUOTIENT MAP

Each inherited quotient is split into **dimensionful** and **dimensionless** components. The distinction is
load-bearing for W1-C: the earned layer has no intrinsic unit.

| # | Quotient object | Components (dimensionful / dimensionless) | Upstream supplied variables | Earned constraints known | Known degeneracies | What would count as FIXING it | Candidate fixing principles | Easiest counterexample to fixing | Empirical consequence if fixed |
|---|---|---|---|---|---|---|---|---|---|
| Q1 | one-particle contraction `e^{−Kt}` / lift-equivalence class | rates, couplings / topology, number field, statistics | `K`, graph, pin, lift, complex structure | passivity, locality, linearity (E-1–E-3); descent invariance (SCOUT-0 P-08/P-09) | all lifts share `K` (P-08/P-09); cospectral vertices share `μ_r` | `K` determined from earned data; or the lift's field/statistics forced | inverse spectral theory + topology (C08); reconstruction axioms / local tomography (C32); spin-statistics (C15) | scaling `K → λK`; cospectral graphs | which `K` / which lift |
| Q2 | effective forcing law (memory kernel + force law) | `J(ω)` amplitude and width / FDT ratio, cumulant ratios | bath, preparation, `T(ω)`, coupling | none earned (FDT borrowed, NO_GO 7) | P-17: ontology invisible; mode-temperature freedom (P-17 converse) | `T(ω)` forced constant (KMS), or cumulants forced | complete passivity ⇒ KMS (C12); local detailed balance (C22); Gallavotti–Cohen (C23) | non-equilibrium Gaussian bath | FDT ratio; fluctuation relations |
| Q3 | local spectral data `μ_r` (edge + moments) | edge location `λ₀`, moments / edge exponent `γ`, visible sign | `K`, pin, graph, readout | E-1 implications; EDA-01: FIX = 0 | chain vs half-plane, same predicates (EDA W-A) | `γ` or `λ₀` forced | spectral dimension (C20); sum rules (C19); Krein inverse problem (C33) | half-plane vs chain witness | memory-tail exponent |
| Q4 | many-body IR data | `v_IR` / Luttinger `K`, central charge `c`, soft-point momenta (lattice units), fixed-point type | sector, statistics/interaction, dimension, symmetry | none earned (SFG-0: no GRUT formation variable) | P-02b: exclusion ↔ interaction interchangeable | soft points / `c` / `K` forced | LSM/Oshikawa (C01); Haldane (C02); FQS unitarity (C04); anomaly matching (C16); TKNN (C29) | translation-breaking defect | structure-factor soft points; `c` |
| Q5 | outcome branch generator | rate `σ` / hitting law `h(p)` (scale function), absorption | stochastic law, `p = |α|²`, Hilbert space | none earned (Born BORROWED, NO_GO 7) | many parents share `2b/σ²` (P-15) | `h` forced | decomposition independence (P-15 corollary, C11); Gleason/Busch (C10) | nonlinear drift (CE-06) | Born frequencies |
| Q6 | gravity effective relation | `α_M, α_B` amplitudes / ratio `(Σ−1)/(μ−1)`, signs | `p_tt`, slip, `x` kernel, background | none earned (K1 conditional; K1-HS degenerate) | luminal Horndeski fills both K1 sides | relation or sign forced | positivity bounds (C17); GW-speed + stability (C18); dispersion/causality sign (C35) | K1-HS EFT-stable families | μ, Σ |
| Q7 | thermal scale | `T` / `T/H` | background state | KMS `T = H/2π` given dS and Hadamard (canonical record: "forced uniquely" given background) | background supplied | `T` forced without the background | KMS + Hadamard (C13) | non-dS background | noise/dissipation ratio |
| Q8 | anomaly ratio `α = a/c` | — / dimensionless `1/3` | IR carrier (conformal mode) | CONDITIONAL theorem (rung9a: IF conformal-mode carrier THEN 1/3) | carrier supplied | carrier forced | anomaly matching (C16) | another IR carrier | μ family amplitude (via `x`) |

**Reconnaissance observation (a priori, to be tested in W1-C):**
- Every earned predicate inspected so far is **scale-free**: passivity, locality, linearity, CM, gap > 0,
  positivity, and the cone's existence.
- So only **dimensionless** quotient components are even eligible to be fixed by earned structure.
- The dimensionless components are: exponents, soft momenta in lattice units, `c`, Luttinger `K`, `h(p)`,
  ratios and signs, `T/H`, and `a/c`.
