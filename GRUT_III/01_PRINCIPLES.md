# 01 — Principles of GRUT III

Every statement carries a tier: **THEOREM**, **DERIVED** (by a script in `calc/`), **POSTULATE** (counted),
**DATUM** (measured constant, not counted), **CONJECTURE**. "Replaces" names the old register node(s) retired.

## Data (taken from nature; not postulates)

ħ, c, G, k_B; the cosmological constant Λ; the particle spectrum (masses, spins, multiplicities, gauge group).
Rationale: the charter's own ruling (2026-08-04) that "an observed value is a datum every candidate theory must
reproduce, not a discretionary posit" — applied consistently. The record applied it to Λ and not to G or ħ; that
inconsistency inflated its ledger by two.

## P0 — Covariant quantum vacuum (POSTULATE, +1)

Physics is a local relativistic quantum field theory on a Lorentzian manifold: a net of local operator algebras
with Einstein causality, a Poincaré-covariant (locally) vacuum satisfying the spectrum condition, and states
of Hadamard form. *Replaces* the record's rung-1 package (system/bath split, Gaussian truncation, Lorentzian
causal structure, 4d-covariant gauge-orbit zero: +4) and rung 6's "quantization" (+1): the Hilbert-space
structure, the causal structure, and the gauge-orbit statement are parts of one framework assumption; the
Gaussian truncation is not assumed at all (linear response is a derived approximation).

## P1 — Horizon relationality (IDENTIFICATION, 0)

For any observer, the **bath** is the algebra of the causal complement of the region they can access, and the
**responsive medium** is the restriction of the global state to their accessible wedge. The system/bath split is
the causal structure; it is observer-relative; nothing is "integrated out" by hand.

*Consequences (THEOREM):*
- Wedge in the Minkowski vacuum: the restricted state is KMS at β = 2π/a with respect to the boost flow
  (Bisognano–Wichmann 1975/76). This is the Unruh temperature as a theorem of P0.
- Bifurcate Killing horizon: at most one Hadamard state invariant under the horizon Killing flow; it is KMS at
  β = 2π/κ (Kay–Wald 1991).
- De Sitter, geodesic observer: the Bunch–Davies state is thermal at T = H/2π in proper time
  (Gibbons–Hawking 1977); verified in closed form in `v1` (§D1–D2).
- The KMS/FDT lock (the record's rung 2, "−1 derived") is inherited: the noise kernel is coth(βω/2) × the
  dissipative kernel, with β fixed by the horizon, not supplied.

*Replaces* `rung1_ontology_finite_memory` (+1) and `background_time_translation_flow` (+1): the "time-translation
flow" is the horizon's Killing flow; "finite memory" is no longer an axiom (see D2 for what is actually true).

## P2 — Local equilibrium of the medium (POSTULATE, +1)

Across every local causal horizon the medium's entropy is proportional to area, S = η A, and energy flux through
the horizon obeys Clausius, δQ = T δS, with T the KMS temperature of P1.

*Status of the parts:* area-proportionality of the entanglement entropy of a region's complement is a theorem of
quantum field theory (Bombelli–Koul–Lee–Sorkin 1986; Srednicki 1993) with a UV-sensitive coefficient; the
coefficient η = 1/(4Għ) is the **datum** G. The Clausius relation — that the medium is in local equilibrium —
is the one genuine structural postulate.

*Consequence (DERIVED, `v2`):* G_ab + Λ g_ab = 8πG T_ab, with Λ an **integration constant** (D3).
*Replaces* `rung5_gr_limit` (+2 assumed) → +1 postulate + 1 datum; closes NO_GO_LEDGER entry "GR borrowed".

## P3 — Universal coupling (POSTULATE, +1)

The medium couples to matter only through the stress tensor: the operational equivalence principle. In the
fluctuation sector this fixes the vertex (κ/2) h_μν T^μν, κ² = 32πG, so that the medium's response kernels are
the retarded and Hadamard ⟨T T⟩ functions and nothing else.

*Consequence (DERIVED, `v4`):* Im Π(s) = c_eff(s) G s², the species-resolved Weyl-anomaly coefficient; all
"spectral edge data" fixed by the spectrum (D5). *Replaces* `p_tt_ansatz` (+1): the tensor/scalar split of the
kernel is computed from ⟨TT⟩, not chosen; for conformal matter the spin-0 part vanishes identically.

## P4 — Time orientation (POSTULATE, +1)

The Hadamard (positive-frequency) condition is imposed toward the future: the horizon traced out is the future
horizon, β > 0, dissipation damps toward the future. This is the record's Past Hypothesis relocated to the one
place it enters. It is still an input; GRUT III does not derive the arrow's direction and says so (D8).
*Replaces* `arrow_of_time` (+1) at equal price, with the location now exact.

## P5 — Non-contextual outcome statistics (POSTULATE, +1)

The probability of a measurement outcome depends only on the projector (effect) measured, not on which other
compatible observables are co-measured. *Consequence (THEOREM):* for Hilbert spaces of dimension ≥ 3 the
probability measure is Tr(ρP) — the Born rule (Gleason 1957; Busch 2003 for effects, all dimensions).
*Replaces* the supplied `p = |α|²` of rung 6 (+1) by a principle at equal price, and discharges SCOUT_0's
P-15 requirement for "an earned reason for a martingale branch coordinate": non-contextuality is exactly
decomposition-independence of outcome statistics, which P-15 proved equivalent to the martingale property.

## Ledger

| | count |
|---|---|
| GRUT register net (2026-09-30) | **+16** |
| GRUT III structural postulates P0, P2, P3, P4, P5 | **+5** |
| Data (not counted) | ħ, c, G, k_B, Λ, spectrum |

The derivation of the reduction, node by node, is `04_LEDGER_MIGRATION.md`. The reduction is not cosmetic:
eleven units disappear because kernels, temperatures, memory times, the GR limit, the KMS lock, the TT
structure and the two-scale dark-energy ansatz are **computed** from P0–P3 in `calc/`, and because two units
were data all along.
