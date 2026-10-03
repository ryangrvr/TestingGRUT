# DA0 · C3 — PRIMARY MODEL 2: CANDIDATE AUDIT (proposal only; NOT RUN; no model selected without owner review)

**Mechanism class (owner).** **Local adaptation with endogenous competition / resource constraint.**
- The fixed microscopic law contains a bounded local resource or homeostatic constraint, so that adaptive interactions
  **compete** rather than each simply reinforcing its observed correlation.
- **Firewall:** frustration or heterogeneity may arise as a consequence. It must **not** be encoded as a target pattern,
  module, hierarchy, sign schedule, basin or macrostructure.
- Any conservation / normalisation variable or resource field is a **new priced primitive**.
- Physical realisability outside an algorithmic / neural metaphor is audited for each candidate.

**Sources:** metadata and abstracts only. Literature shows mechanism and precedent; it does not substitute for a C3
result. No novelty is claimed.

## 0. A structural lesson carried from PM1 (PROVED HERE, one line)

**Envelope persistence under sign-blind constraints.** Suppose the competition constraint depends on J only through |J|.
Examples: a per-site budget Σ_{j∼i} J_ij² = κ, or Σ_j |J_ij| = κ.
- Then the feasible set S is invariant under J → |J|.
- PF5-E gives 𝓕_T(J) ≥ 𝓕_T(|J|) pointwise. Hence **the global minimisers over S are still unfrustrated.**

**Consequence.** Any candidate that hopes to *beat the unfrustrated envelope* needs **sign-sensitive** competition, or a
different substrate. Sign-blind normalisation can only reshape magnitudes; local minima are not covered by this
statement.

## 1. Candidates

| ID | candidate (fixed local law) | priced primitive(s) | generated or encoded? | physical realisation outside the neural metaphor | main risk |
|---|---|---|---|---|---|
| **M2-A** | **site-normalised Hebbian Ising:** PM1 dynamics projected onto a per-site coupling budget (Oja-type: Σ_{j∼i} J_ij² = κ) | per-site budget κ; the projection law | no pattern is encoded | **weak:** synaptic scaling is neural. A material analogue (conserved bond "stiffness" per node) is plausible but not established | **§0: sign-blind, so the envelope persists.** Likely winner-take-all bond selection (dimer / chain covers), which is microscopic bookkeeping (C3-F1), or near-uniform order |
| **M2-B** | **balanced local field (zero-sum) homeostasis:** Hebbian drift projected onto Σ_{j∼i} J_ij = 0 per site | per-site zero-sum constraint; it **breaks the Ising gauge symmetry** (priced) | no pattern is encoded; mixed signs are forced | **weak to moderate:** E/I balance is neural. A charge-neutrality-type constraint is a loose physical analogue | **sign-sensitive, so it escapes the envelope.** But zero-sum is satisfiable **without** frustration: on the cubic lattice J_x = +2a, J_y = J_z = −a has every plaquette unfrustrated. So frustration is not forced. It is also borderline "encoding" (it forces sign heterogeneity by fiat) |
| **M2-C** | **conserved, diffusing coupling resource:** Σ J locally conserved, with transfer between neighbouring bonds (model-B-like) | conservation law; transfer kinetics | none encoded | moderate (conserved local material) | conserved dynamics typically **coarsens**: domain count decreases, so architecture shrinks |
| **M2-D** | **adaptive transport / flow network with flux conservation and local material cost.** Conductances C_e (W) adapt to their own flow Q_e (X, fast: Kirchhoff / Darcy potential flow under a homogeneous source field and an O(1) outlet), e.g. dC_e/dt = a·Q_e^{2γ}C_e^{…} − b·C_e (Hu–Cai type), or concave cost Σ C_e^γ (γ < 1; Bohn–Magnasco) | local adaptation law; cost exponent γ; homogeneous injection; outlet / boundary; optional source fluctuations | **competition is endogenous:** conserved flux makes branches compete for flow. The hierarchy is **not** in the rule | **strong:** river networks (erosion: optimal channel networks), vasculature remodelling, *Physarum*, leaf venation | (i) the architecture lives in W itself (topology), so it must pass C3-F1 with *collective* units (sub-basins / branches with growing participation), not per-edge labels; (ii) the drive must be homogeneous, and the outlet is O(1) supplied boundary information; (iii) deterministic gradient flows make local minima absorbing, so metastability and barriers need priced fluctuations; (iv) much is KNOWN, so C3 adds only the audit against C3 criteria |

**Literature precedent (M2-D):**
- Hu & Cai, *PRL* 111, 138701 (2013): local adaptation dynamics minimises global dissipation. It produces hierarchical
  loop structures under flow fluctuations, and a tree / loop transition.
- Bohn & Magnasco, *PRL* 98, 088702 (2007): concave cost gives trees with many local optima, and a tree-to-loops
  transition.
- Rigon, Rinaldo, Rodríguez-Iturbe et al., *Water Resour. Res.* (1993), optimal channel networks: local optima of energy
  dissipation reproduce river-network hierarchy (Horton-type laws).

## 2. Assessment against the C3-A entry requirements

| requirement | M2-A | M2-B | M2-C | M2-D |
|---|---|---|---|---|
| fixed local rule, bounded description, parameters fixed in N | yes | yes | yes | yes. Kirchhoff is a local conservation law; flows are determined globally, like ⟨s_b⟩ in PM1 |
| no hierarchy / basin / module labels | yes | yes | yes | yes (the outlet is an O(1) boundary datum, priced) |
| Markov augmented state | yes | yes | yes | yes. Deterministic or stochastic adaptation; X is quasi-static |
| thermodynamic / energetic cost stated | budget maintenance | constraint maintenance | transport cost | **explicit:** dissipation plus material cost; the driving flux is supplied |
| plausible mechanism for **growing collective** architecture | low (§0; bookkeeping risk) | uncertain (frustration not forced) | low (coarsening) | **highest:** hierarchical branching with sub-basin sizes at all scales, and many local optima under concave cost (precedent) |
| scalable method (no 2^N) | mean-field / landscape | landscape | PDE / coarsening theory | **yes:** the X-solve is linear (sparse Laplacian); the landscape is the dissipation-plus-cost functional |
| continuity with PM1 / Ising lineage | high | high | medium | **low:** a new substrate |

## 3. Proposal (subject to owner review; nothing selected or run)

**Proposed Primary Model 2: M2-D, an adaptive flow network with conserved-flux competition and local material cost.**

Why:
1. **Competition is endogenous and physical.** Flux conservation makes channels compete for a conserved resource. No
   pattern, schedule or hierarchy is written into the law.
2. It is the only candidate with **strong physical realisations outside the neural metaphor** (rivers, vasculature,
   *Physarum*).
3. It has **known precedent for multiplicity** (many local optima under concave cost) **and hierarchy** (Horton-like
   branching). These are exactly the C3-A ingredients PM1 lacked.
4. Its landscape (dissipation + cost) is a **physically derived functional**. The C3R-05 detector transfers directly:
   symmetry-inequivalent local minima, canonical minimum-barrier connections, diverging barriers.

**Registered risks / hostile outcomes, if chartered:**
- **C3-F1 hazard.** Architecture lives in W (topology).
  - Required: collective macro-units (sub-basins) with participation growing with N, and elementary rerouting
    transitions with growing PR_edge and contrast.
  - Per-edge topology labels do not count.
- **Metastability needs fluctuations.** Without noise, local minima are absorbing. Fluctuating sources (Hu–Cai /
  Katifori-type) or explicit noise are priced primitives, and barriers must diverge.
- **Drive and boundary.** The source field must be homogeneous (no spatial pattern). The outlet and boundary are priced
  O(1) data.
- **Known theory.** Hierarchy and local-optima multiplicity are KNOWN. C3 can earn only an audited classification
  against C3-A1 (generated vs supplied; collective; growing barriers). It cannot claim discovery.
- **C3-B** stays closed regardless.

**Fallback (if the owner prefers Ising lineage): M2-B.** It is sign-sensitive and so escapes the envelope, but it is
weakly physical, does not force frustration, and is borderline on encoding.

**Rejected as primary:**
- M2-A, by the §0 envelope persistence plus the bookkeeping risk;
- M2-C, by the coarsening risk.

## 4. Next step (owner decision)

The owner chooses the Primary Model 2 candidate, or none. If chosen, a C3 PM2 charter addendum would preregister:
- the exact law, constants, drive, boundary and noise;
- the symmetry quotient;
- the collective units and edge definition;
- the barrier criterion;
- controls: A0 (frozen W), A1 (decoupled adaptation), A2 (supplied hierarchy), **PM1** (the comparator), and K1 – K4.

This happens **before** any computation.
