# 05 — The road to completion: what is actually open, posed sharply

GRUT III has five postulates and a list of data. A theory of everything would derive the postulates from fewer
and the data from none. Here is the honest distance, item by item, with the next computation for each.
Nothing below is a "ladder"; each item is a question with a definite form of answer.

## O1 — The area-law coefficient: where G comes from (the central open problem)

**Pose.** P2 uses S = ηA with η = 1/(4Għ). Area-proportionality is a theorem of P0 (entanglement entropy of a
region's complement); the *coefficient* is UV-divergent in QFT, ∝ Σ_species Λ_UV². Newton's constant renormalises
with the same species sum (Susskind–Uglum 1994; Jacobson 1994; Cooperman–Luty 2014 for the non-minimal
subtleties). So "G from the medium" is Sakharov's induced gravity made exact: the finite value of η requires the
UV completion in which S = A/4G is exact rather than leading-order.
**Next computation (in reach).** Using the species table of D5, verify that the one-loop matter contribution to
δ(1/G) and the entanglement-entropy area coefficient agree species-by-species in the same regularisation
(scalars with the ξ = 1/6 contact term; Weyl fermions; vectors with the Kabat contact term). This is a consistency
theorem for P2, not a derivation of G, and it is worth having in `calc/`.
**Form of a real answer.** A UV-finite theory in which (i) horizon thermality is exact, (ii) the area law has a
finite coefficient, (iii) P0 is recovered below the cutoff. AdS/CFT (Ryu–Takayanagi) satisfies all three in AdS;
GRUT III is agnostic about the completion and demands only these three properties of it.

## O2 — The value of Λ

**Pose.** D3 leaves Λ as an integration constant: the equation of state is silent on it, as unimodular gravity is.
**Dead end, recorded to prevent re-derivation.** The identity F = E − TS = 0 for the static patch (D4) looks like
a selection principle ("the medium relaxes to zero free energy") but holds identically for every H; it selects
nothing.
**Position.** GRUT III treats Λ as a datum and makes the sharp prediction F1 (w ≡ −1). If DESI DR3 establishes
evolution at ≥5σ, the equation-of-state reading of GR is incomplete at zeroth order and P2 must be extended (a
second horizon scale, i.e. a relaxor, re-enters with a price). That is the one place GRUT III can lose soon.

## O3 — The outcome rule: can P5 be derived from P0 + P1?

**Pose.** P5 (non-contextuality) buys Born via Gleason. In GRUT III an observer's statistics are functions of
the state restricted to their wedge (P1), and Einstein causality (P0) forbids signalling from the causal
complement. Any decomposition of the observer's reduced state that is *induced by an operation in the complement*
therefore cannot change the observer's outcome frequencies — which is decomposition-independence for that class
of decompositions, i.e. SCOUT_0 P-15's martingale condition, i.e. non-contextuality restricted to remotely
induced decompositions.
**Next computation (theoretical, high value).** Determine whether the remotely-inducible decompositions of a
reduced state exhaust the projector lattice needed by Gleason (dimension ≥ 3). If yes, P5 follows from P0 + P1
and the ledger falls to +4; if no, the residue is the precise content of the measurement postulate. Known
partial results: Gisin 1990 (nonlinear evolution + collapse ⇒ signalling), Zurek's envariance argument, and
the Hughston–Jozsa–Wootters theorem (every decomposition of ρ_W is remotely inducible given a purification) —
the last suggests the answer is yes, and should be checked carefully for its assumptions.

## O4 — The arrow (P4)

Unchanged and honest: the direction is the choice of the future horizon as the traced-out one, equivalently the
Hadamard condition with positive frequencies toward the future. The physics behind it is the cosmological
initial state. GRUT III locates the input; it does not remove it. No computation proposed.

## O5 — The memory dictionary (replaces the finite-memory axiom)

What "memory" means in GRUT III, so that the word is never again used for four different things:

| object | memory | status |
|---|---|---|
| response of a conformal medium at a point | none (odd polynomial in ω) | THEOREM (D2) |
| noise of operator of dimension Δ at a point, dS geodesic | exponential, rate (Δ−k)H after polynomial cancellations; 4H for T̂ | DERIVED (D2) |
| response of a homogeneous/extended mode | power law, light cone to the horizon (48/πΔ⁵ …) | record T3 chain |
| matter thresholds | power law × oscillation at 2mᵢ, exponents 7/2, 5/2, 3/2 | DERIVED (D5) |
| gravitational (minimally coupled) sector in dS | non-stationary, secular IR growth | record (tt_worldline); the open quantum-gravity IR problem |

**Next computation (in reach).** Generalise D1–D2 to any FRW a(t): for conformal matter along a comoving worldline
G = 1/(4π² a a'(η−η')²) exactly; compute the (non-stationary) noise correlation time in the radiation and matter
eras and show it tracks 1/H(t). This replaces "finite memory" by a derived, epoch-resolved function.

## O6 — The one risk taken: a₀ (see `06_CONJECTURE_A0.md`)

Next computation: Jacobson's Clausius construction with the de Sitter–Unruh temperature √(a² + H²)/2π and a
Rindler patch of finite size inside the cosmological horizon. Either an acceleration scale with a coefficient
emerges, or the construction shows that no local-equilibrium correction appears at a ∼ cH₀ — both outcomes are
informative and both are cheap relative to another audit wave.

## O7 — Owed independent verifications of this proposal

1. The closed forms of D2 by a second method (explicit residue sum at the order-8 poles).
2. The species kernel of D5 against Han–Lykken–Zhang (1999) and Giudice–Rattazzi–Wells (1999) primary texts,
   including the massive-vector threshold polynomial, which is **not computed here**.
3. The Clausius bookkeeping of D3 against Jacobson (1995) and Eling–Guedens–Jacobson (2006) with the
   non-equilibrium (shear) term carried explicitly.
4. The record correction in D1 (0.0339 vs 0.0341) by re-running `calc/worldline_reduction.py` with the
   one-sided/symmetrised distinction made explicit.
