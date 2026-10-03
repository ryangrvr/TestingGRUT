# DA0 · C3 — PRIMARY MODEL 2 CHARTER: ADAPTIVE CONSERVED-FLOW NETWORK

**Status:** CHARTER ONLY — **NOT RUN.** No PM2 calculation or simulation exists.

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
