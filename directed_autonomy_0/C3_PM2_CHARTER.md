# DA0 · C3 — PRIMARY MODEL 2 CHARTER: ADAPTIVE CONSERVED-FLOW NETWORK

**Status:** CHARTER ONLY — **NOT RUN.** No PM2 calculation or simulation exists.

> **PM2 CHARTER REPAIR 01** (PM2R-01 … 07) was applied **before any computation**.
> - **Owner ruling:** charter CONDITIONAL PASS; frozen law ACCEPTED. After this repair is committed, **deterministic
>   Stage A is approved to run** without another stop.
> - **Repaired text (§R below) supersedes** the original sections wherever they differ.

**Owner ruling:** **M2-D SELECTED as C3 Primary Model 2. M2-B not selected.**

**Base:** `grut-directed-autonomy-0 @ b71e9307`.

**Preserved:**
- the PM1 terminal, unchanged;
- preflight `dde3b949…`;
- Repair 02 `31b5714e…`;
- C2 boundary `99428ff…`;
- RA0-frozen `ab4fd86…`.

**Constraints:** C3-B closed; no PR; no merge.

**Labels (as in DA0):** KNOWN / REDERIVED, PROVED HERE (internally), NUMERICAL ILLUSTRATION, CONJECTURE.

## Purpose

Test whether a bounded local adaptive transport law can generate **growing collective architecture** from endogenous
competition for conserved flow, without supplying:
- modules;
- hierarchy;
- branch labels;
- a patterned drive.

This is a **construction test**. Known transport-network results (Hu–Cai; Bohn–Magnasco; optimal channel networks) are
**precedent and comparators, not discoveries** of this programme.

## PM2-F0 — The single frozen law

**Substrate.** A fixed microscopic graph G_N with edge lengths L_e > 0 and conductances C_e ≥ 0.

**Fast, quasi-static flow (X).** Given C, the node potentials p and edge flows Q solve:

  Q_e = C_e (p_i − p_j)/L_e  for e = (i, j),  Σ_{j∼i} Q_ij = h_i  (Kirchhoff)

with the supplied source field h (PM2-F1). Equivalently, Q minimises the dissipation Σ_e L_e Q_e²/C_e under Kirchhoff.
Edges with C_e = 0 carry no flow. Q is unique whenever the support graph connects every source to the outlet.

**Declared functional** (Hu–Cai class: dissipation + concave material cost):

  **E(C) = Σ_e L_e [ Q_e(C)²/C_e + ν·C_e^γ ]**, with γ ∈ (0, 1) (concave cost).

**Derivative.** By the envelope theorem for the Kirchhoff minimisation:

  ∂E/∂C_e = L_e [ −Q_e²/C_e² + νγ·C_e^{γ−1} ].

**Adaptation law (W) — the exact law, frozen.** It is the gradient flow of E in the positive-orthant metric
g_e = L_e/(κC_e):

  **dC_e/dt = κ [ Q_e²/C_e − νγ·C_e^γ ]**.

- **Local:** each edge uses only its own Q_e, C_e and L_e.
- **Positivity:** C_e > 0 is preserved, and C_e = 0 is invariant (pruned edges stay pruned).
- **Lyapunov:** dE/dt = −Σ_e (κC_e/L_e)(∂E/∂C_e)² ≤ 0, so E is a Lyapunov function.
- **Stationary edges** satisfy Q_e² = νγ·C_e^{γ+1}, i.e. C_e ∝ |Q_e|^{2/(γ+1)} (Murray-type), or C_e = 0.
- **Attribution.** This is a gradient form of the Hu–Cai dissipation-plus-cost adaptation class. It was written here from
  the declared functional; **an exact match to the published Hu–Cai equation is not claimed** (full text not re-read).

**Regimes:**
- γ < 1: concave cost, **tree-favouring** (loops are suppressed in the optimisation picture; Bohn–Magnasco).
- γ > 1: convex cost, distributed / loop-favouring.

**Fixed constants, independent of N:** γ = 1/2 (primary); ν = 1; κ = 1; L_e = 1. **No constant is varied after results.**

## PM2-F1 — Primary geometry and drive (fixed before results)

**Substrate.** The 2D square lattice on an L × L grid of nodes with **open boundary**; N = L² nodes, bounded degree ≤ 4.
This is an existence / construction test, **not** a claim that fundamental space is 2D. The lattice and its embedding are
**supplied (priced)**.

**Drive (single convention).**
- **One outlet** at the corner node (0, 0).
- **Homogeneous injection** at every other node: h_i = 1/(N − 1) for i ≠ outlet, and h_outlet = −1, so the total flux is
  fixed at 1.

**Forbidden:** patterned source fields, basin templates, target hierarchies, hand-chosen terminals, food / demand maps,
and spatially varying κ, ν or γ. The outlet and outer boundary are **supplied O(1) information**, recorded in the
accounting.

**Initial ensemble.**
- **(i)** Exactly uniform C_e = 1 (deterministic).
- **(ii)** A preregistered perturbed ensemble C_e = 1 + δξ_e, with ξ_e i.i.d. uniform on [−1, 1], δ = 0.01, and seeds
  1 … 20 (`numpy.random.default_rng(seed)`). The perturbation law is supplied and priced. Statements about the ensemble
  are ensemble statements (C2-F7 discipline).

**Sizes** (computation stage, when authorised): L ∈ {16, 24, 32, 48, 64}.

## PM2-F2 — Stage A is deterministic

There are no fluctuating sinks, no conductance noise, and no stochastic sources in Stage A. Noise or fluctuations may be
chartered later (**PM2-BIS**) only if Stage A earns a reason for them. They are not part of any Stage A positive.

## PM2-F3 — Topology itself is not enough

**Things that do not count:**
- the number of active edges;
- graph depth growing with diameter;
- the number of possible spanning trees;
- "every spanning tree is a local minimum".

None of these satisfies C3-A1. A positive must **exceed generic tree combinatorics** (PM2-F4).

### Canonical collective units

**Branches.** For a tree-like stationary network rooted at the outlet, take an active edge e and remove it. **B_e** is the
component not containing the outlet, and the branch mass is **M_e = Σ_{i∈B_e} h_i**. It is derived from the generated
topology and the supplied outlet, not supplied as a module.

**Loopy stationary networks** (if any arise, which is unexpected at γ < 1). The preregistered generalisation is the
**flow-weighted upstream set**: for each active edge e, M_e = |Q_e|, which equals the branch mass on trees. No other
decomposition may be chosen post hoc.

### Elementary architecture transitions (canonical)

**Landscape nodes.** Symmetry-inequivalent (PM2-F5) **stable stationary networks** (local minima of E).

**Landscape edges.**
- **Adjacency:** minimum-barrier / minimum-action connections on E between minima.
- **Computational proxy (registered as an upper bound):** for tree minima differing by one **fundamental swap** (add a
  non-tree edge f, remove an edge of the cycle created), take the maximum of E along the straight conductance path
  between the two minima. A mountain-pass computation (e.g. nudged elastic band, a declared numerical method) may refine
  this. The proxy never **defines** adjacency by clustering.

**Rerouted mass** (the PM2 collectivity measure). For adjacent tree minima T_a and T_b:
- **Absolute:** m_reroute(a, b) = the total source mass whose outlet path changes between T_a and T_b.
- **Relative:** R(a, b) = m_reroute / (total source mass).

**Grading:**
- **Collective:** the absolute rerouted mass of minimum-barrier elementary transitions **grows with N**, judged by the
  sign of a fitted scaling exponent across sizes (a diagnostic; no post-hoc threshold).
- **System-scale:** R = O(1), reported separately.
- **Microscopic:** swaps rerouting O(1) mass are microscopic and fail.

## PM2-F4 — Random-tree firewall (hard comparators)

| ID | comparator | purpose |
|---|---|---|
| **R-TREE** | the uniform spanning-tree ensemble on the same graph, rooted at the same outlet (Wilson's algorithm; seeds 1 … 20) | generic random tree combinatorics |
| **SP-TREE** | a deterministic shortest-path (BFS) tree on the same graph and outlet, with ties broken lexicographically by node index | simple geometric routing |

**Primary invariants** (chosen now; no post-hoc statistic selection):

| ID | invariant | positive requires |
|---|---|---|
| **P1** | **branch-mass tail exponent** τ, from P(M_e ≥ m) ~ m^{−τ} over active edges (fit window: the middle decade of m at the largest L, declared now) | the PM2 τ differs from **both** R-TREE and SP-TREE, with non-overlapping ensemble confidence intervals at the two largest sizes |
| **P2** | **Hack-type exponent** η_H, from the mean outlet-path length of B_e's root vs M_e: ℓ ~ M^{η_H} | the same comparison rule as P1 |
| **P3** | **rerouting-mass exponent** for minimum-barrier elementary transitions: m_reroute ~ N^{ρ} | **ρ > 0** (collective) |

**Secondary (reported, not used for the terminal unless the primaries pass):**
- Strahler order of the outlet (canonical for rooted trees), vs both controls;
- the barrier scaling exponent for elementary transitions.

## PM2-F5 — Symmetry quotient

**Quotient group.** The symmetries preserving the substrate, the homogeneous drive and the corner outlet: the
**reflection across the diagonal through (0, 0)** (order 2). Translations are broken by the outlet.

**Counting rule.** Architectures are counted modulo this group. Symmetry copies are not distinct.

**Uniform initial condition (i).** It is symmetric, so the deterministic flow preserves the symmetry. A symmetric
stationary state that is a **saddle** of the symmetry-broken dynamics is reported as such, and is never counted as an
architecture.

## PM2-F6 — Landscape / multiplicity firewall

**KNOWN comparator.** In the concave regime, spanning trees can be local minima. **Rediscovering multiplicity earns
nothing.** PM2 must establish at least one **additional growing collective property** generated by the adaptive law: P1,
P2 or P3 exceeding the comparators, plus growing barriers for a full positive.

## Stage A — preregistered questions (exactly these)

1. Does deterministic local adaptation converge to a sparse / tree-like architecture from homogeneous initial
   conductance (ensembles (i) and (ii)), without a supplied hierarchy?
2. Does the branch / source-mass structure have an asymptotic collective invariant (P1 / P2) that differs from R-TREE
   **and** SP-TREE?
3. Does the number or depth of collective, symmetry-inequivalent structures grow with N? Count = distinct minima reached
   by ensemble (ii), a sampling lower bound; depth = outlet Strahler order.
4. Do elementary adaptive reroutings involve source mass growing with N (P3)?
5. Do barriers between collective alternatives grow with N? This may initially be analytic / scaling-only, or use the
   upper-bound proxy. **C3-A1 is not weakened if barrier scaling is unavailable**; a partial terminal is reported
   instead.

## Controls

| ID | control |
|---|---|
| **A0** | frozen conductances C ≡ 1 (the full lattice; no topology selection) |
| **A1** | adaptation decoupled from flow feedback: Q fixed at its uniform-network value Q⁰ throughout, i.e. dC_e/dt = κ[(Q_e⁰)²/C_e − νγC_e^γ] |
| **A2** | an explicitly supplied hierarchical tree (an H-tree-type construction on the same lattice), to validate the detectors. Graded ARCHITECTURE SUPPLIED |
| **PM1** | the adaptive Hebbian comparator (conceptual: adaptation without conserved-flow competition) |
| **K1 – K4** | permanent DA0 comparators |
| **R-TREE, SP-TREE** | PM2-F4 |

No further controls are added without a stated purpose.

## Preregistered PM2 terminals

| ID | terminal | requirement |
|---|---|---|
| **PM2-A** | **ENDOGENOUS GROWING COLLECTIVE FLOW ARCHITECTURE — CONSTRUCTED** | a bounded fixed local law; homogeneous fixed drive; no supplied hierarchy; P1 or P2 beyond R-TREE / SP-TREE; P3 with ρ > 0; growing architecture count or depth; growing barrier / stability scale. **Satisfies C3-A1 for PM2** |
| **PM2-B** | **GENERATED TREE / NETWORK — GENERIC TOPOLOGICAL HIERARCHY ONLY** | |
| **PM2-C** | **LOCAL CHANNEL SELECTION / MICROSCOPIC BOOKKEEPING** | |
| **PM2-D** | **KNOWN OPTIMAL-NETWORK STRUCTURE RECOVERED — NO NEW C3 COMPRESSION** | |
| **PM2-INDETERMINATE** | stated technical obstruction only | |

**Partial terminal.** If P1 / P2 / P3 pass but barriers are unavailable, report **PM2-A PARTIAL — BARRIERS
UNESTABLISHED**. It is not counted as C3-A1.

## Information accounting

**Supplied (priced):**
- the microscopic graph and embedding;
- the outlet and outer boundary;
- the homogeneous source law;
- the conductance state W = C and the flow state X = (p, Q);
- the adaptation equation;
- γ, ν, κ, L_e;
- the initial ensemble (uniform, plus δ-perturbation with preregistered seeds);
- the power supplied by the drive, i.e. the dissipation Σ L Q²/C, and the material-maintenance cost.

**Derived only if earned:** the active topology; the branch hierarchy; collective rerouting structure; the metastable
landscape; scaling exponents.

**Not earned automatically:** subsystem factorisation; consciousness; C3-B; a generic architecture selector; TRUE
COMPRESSION.

## Interpretation, recorded in advance

**What a PM2-A positive would not mean:** that GRUT discovered transport-network hierarchy.

**The narrower construction result it would establish:**

> A fixed, physically realised local law with conserved-flow competition can generate growing collective architecture
> without that architecture being supplied explicitly.

That would show the C2 architecture primitive is **replaceable in principle** by a bounded adaptive physical mechanism:
local adaptation + conservation + resource cost. Whether that mechanism belongs in the final working GRUT law remains a
later question.

## Hard stop

Charter only. No PM2 calculation or simulation.

**Sources** (metadata / abstracts; full texts not re-read):
- Hu & Cai, *PRL* 111, 138701 (2013);
- Bohn & Magnasco, *PRL* 98, 088702 (2007) (cond-mat/0607819);
- Rigon, Rinaldo, Rodríguez-Iturbe et al., *Water Resour. Res.* 29 (1993), optimal channel networks.


---

## §R — PM2 CHARTER REPAIR 01 (supersedes conflicting text above)

### PM2R-01 — Fixed source density, and exact scaling covariance

**Drive (replaces PM2-F1's drive):**
- h_i = +1 at every non-outlet node; h_outlet = −(N − 1), fixed by Kirchhoff conservation.
- Same corner outlet, same homogeneous bulk rule.
- On a tree: **S_e = M_e = |B_e|**.
- Rerouted mass = the number of unit-source nodes whose outlet path changes; R = m_reroute/(N − 1).
- Initial condition (unchanged): C(0) = 1 (uniform) or 1 + δξ.
- (The old total-flux normalisation made m_reroute ≤ 1 at every N, so ρ > 0 was impossible.)

**PROP PM2-S (exact scaling covariance; PROVED HERE).** Two facts:
- **Linearity:** for fixed C, Q(C, h) is linear in h.
- **Ratio invariance:** Q(λC, h) = Q(C, h), because flows depend only on conductance ratios.

Let C* be stationary for h, so Q*_e² = νγ C*_e^{γ+1} on active edges and C*_e = 0 otherwise. Under h → ah, put
**C′ = a^{2/(γ+1)} C***. Then:
- Q′ = aQ*, and Q′² = a²Q*² = νγ (a^{2/(γ+1)}C*)^{γ+1}. So C′ is stationary.
- The support (topology) is unchanged.
- **E(C′; ah) = a^{2γ/(γ+1)} E(C*; h)**, since a²/a^{2/(γ+1)} = a^{2γ/(γ+1)} = (a^{2/(γ+1)})^γ.
- For γ = 1/2: **C → a^{4/3}C and E → a^{2/3}E.** ∎

**Caveat on trajectories.** Trajectories map onto trajectories (with a time rescaling) **only if the initial condition is
rescaled too**. The set of stationary points is covariant; *which* one is reached from C(0) = 1 can depend on the initial
scale. The fixed-density drive with C(0) = 1 is therefore what is frozen here.

**Reduced energy.** Ē = E / H_N^{2γ/(γ+1)}, with H_N = Σ_{h_i>0} h_i = N − 1, i.e. Ē = E/(N − 1)^{2/3}. Any barrier
diagnostic reports **both raw and reduced** values. Trivial drive scaling is not a barrier result.

### PM2R-02 — Pruning is absorbing: tree-to-tree swaps are not Stage-A dynamics

**Boundary convention (lower-semicontinuous):**
- Q_e²/C_e = 0 when Q_e = 0 and C_e = 0;
- Q_e²/C_e = +∞ when C_e = 0 but nonzero flow would be required.

**Terminology.** A fundamental swap is a **static landscape adjacency / rerouting comparison**. Straight-path or
mountain-pass quantities are **static energy-barrier diagnostics**. Static barrier paths must keep every source connected
to the outlet (the straight path between two trees keeps the union of their supports positive for t ∈ (0, 1)). **Neither
is a physical elementary transition of the deterministic law.**

**Stage-A ceiling (maximum possible grade):** **PM2-A PARTIAL — ENDOGENOUS GROWING COLLECTIVE FLOW ARCHITECTURE
CONSTRUCTED; DYNAMICAL BARRIERS / REVERSIBLE TRANSITIONS UNESTABLISHED.** This is **not C3-A1**. If reached, it may justify
a later, separately chartered PM2-BIS (a physical fluctuation / edge-revival mechanism), which is **not added now**.

### PM2R-03 — Initial-condition information firewall

The i.i.d. perturbation carries O(N) random numbers, so:
- **A seed-specific final topology is not derived architecture.** The seed selects an instance, as in spontaneous
  symmetry breaking.
- **The number of distinct trees across seeds is a sampling diagnostic only** (bounded by 20). It is not a growth
  criterion.
- **What can be earned** is an **ensemble-stable structural law** of the fixed dynamics: branch scaling, path / basin
  scaling, rerouting scaling, and canonical depth beyond the controls.

**Stage-A question 3 (replaced):** "Does an ensemble-stable canonical depth / hierarchy measure (outlet Strahler order)
grow with N beyond R-TREE and SP-TREE, robustly across seeds?"

### PM2R-04 — P1, repaired

**Per-tree fit.**
- S_e = |B_e|. For each tree, the CCDF P(S_e ≥ s) is evaluated at s_k = unique integer rounds of 20 log-spaced points in
  **[√(N/10), √(10N)]** (one decade centred at √N, fixed algebraically).
- τ = −(OLS slope of log P vs log s).
- Edges within one tree are not independent samples.

**Uncertainty:**
- **PM2 perturbed** and **R-TREE**: one τ per seed (20 each); report the mean and 95% t-interval.
- **SP-TREE**: a single τ_SP, with no interval.
- **The uniform-start PM2 tree** is reported separately.

**Pass at L = 48 and 64:**
- the PM2 interval does not overlap the R-TREE interval; **and**
- τ_SP lies outside the PM2 interval.

### PM2R-05 — P2, repaired

**Branch length.** For edge e = (u → v), with v downstream:

  L_branch(e) = max_{w ∈ B_e} d_T(w, v) = 1 + height(u).

This is the longest upstream path draining through e.

**Fit.** L_branch ~ S_e^{η_H}, by an OLS fit of log L_branch vs log S_e over edges with S_e in the P1 window. One exponent
per tree; the same uncertainty and comparison rules as P1.

### PM2R-06 — P3, repaired (structural diagnostic only)

**Swap set (canonical, complete; no selection by rerouted mass).** Every non-tree lattice edge f, paired with every tree
edge e on the cycle that f closes. Every spanning tree defines a stationary point of the law (on a tree the flows are
fixed, Q_e = S_e, and C_e = (S_e²/νγ)^{1/(γ+1)}), so **every swap is a valid stationary-tree candidate**.

**Exact identity (PROVED HERE).** Removing e re-attaches B_e through f. Exactly the nodes of B_e change their outlet path,
so **n_reroute(f, e) = S_e**.

**Statistic.** Per tree, the **median** of n_reroute over the full swap set (the mean is also reported). Fit
median ~ N^ρ across L; R = n/(N − 1). The same statistic is reported for R-TREE and SP-TREE as structural context.

**Scope.** ρ > 0 is a **structural** collectivity diagnostic. It establishes neither metastability nor transitions.

**Static-barrier diagnostic (reported; it cannot upgrade Stage A).**
- 10 swaps per PM2 tree, drawn uniformly from the full swap set with `default_rng(1000 + seed)`.
- For each swap, take the maximum of E along the straight conductance path (11 points, t = 0, 0.1, …, 1) between the two
  stationary trees, minus E(T_a). Report it raw and reduced.

### PM2R-07 — A1 scope

A1 freezes Q at the uniform-network flow Q⁰, which already reflects the outlet, boundary and geometry. **A1 is not an
architecture-free null.** It tests only whether conductance → flow feedback matters. **A0 is the true no-adaptation
null.**

### Revised Stage-A positive

A Stage-A positive requires **all** of the following:
1. sparse / tree architecture from the fixed homogeneous law;
2. P1 or P2 distinct from **both** R-TREE and SP-TREE (PM2R-04 / 05);
3. P3 with ρ > 0;
4. canonical depth (Strahler) growing beyond the controls;
5. ensemble robustness.

**Barrier diagnostics never upgrade Stage A.** The maximum grade is the PM2R-02 ceiling.

**Implementation parameters** (numerical only; not model parameters): explicit Euler with adaptive step (max relative
change 5% per step); pruning floor C < 10⁻¹² · max C → 0; convergence when the relative stationarity residual on active
edges is < 10⁻⁶ and the topology has been unchanged for 2000 steps (or a step cap, reported).

**Stage-A grid (fixed):**
- L ∈ {16, 24, 32, 48, 64};
- PM2 uniform start, plus perturbed seeds 1 … 20 (δ = 0.01);
- R-TREE seeds 1 … 20; SP-TREE;
- A0, A1, A2. **A2** = a deterministic recursive-bisection spanning tree (supplied hierarchy). Each rectangle is split
  along its longer side; the far half attaches across the cut at the mid-cut cell; recursion runs from the outlet.

**No changes** to γ, ν, κ, δ, geometry or drive after results. No noise. No C3-B.
