# Independent review of Codex e160356: physical curvature reheating, constrained preparation, portal sources

26 October 2026. Claude, issue #13. Tracks: `tracks/rp_reheating`, `rp_preparation`, `rp_portals`, `rp_critic`. The verdict
tables are in [tracks/rp_critic/CONSOLIDATED.md](tracks/rp_critic/CONSOLIDATED.md).

All of this is our own code, run by same-session AI agents. It is not outside replication. It covers one supplied
tree-level model at the heavy-mediator spot point. Nothing here is observational evidence, a validation, an exclusion
or a GRUT prediction. Candidate 1 is unchanged; the vacuum density is underived. The separate kinetic calculation
(`tracks/kin_*`) was not touched.

## Result

**Every number in e160356 that was checked reproduces**: the six reheating runs by three independent engines (to the
printed digits), the boundary model by three independent builds, the portal sources by two momentum-resolved solvers.
No numerical error was found. Four corrections in substance or scope:

1. **The yield deficit (18–78% of Y_Q) is a property of Codex's grid, not of the model.** The charge is fixed early
   (90% by τ ≲ 3, while radiation is ≤ 13% of the density) and then diluted by the entropy made afterwards, by a factor
   D = 2.3–13.4 on Codex's slow-reheating grid (Γ_I/m_I ≤ 0.15). With a matter era of about one Hubble time
   (Γ_I ≈ 0.2–0.3 m_I) the same model, with the same hidden couplings, scale and phase, gives 1.01–1.05 Y_Q for
   m_I = 3–50 H_i, at T_RH ≈ 2–5×10¹⁴ GeV. "No retuning performed" is accurate; "none reaches Y_Q" needs the grid label.
   Overreach check on ourselves: the yield is monotonic in Γ_I, so a crossing exists for every m_I ≥ 3; the content is
   that the model's perturbative range spans 0.18–1.15 Y_Q, not a prediction of Y_Q.
2. **The m_I = 10 Y excitation is a resonance:** 2m_I ≈ ω_YL = 22.25 H_i, peaking at N_YL = 2.22 for m_I = 12. Under our
   ledger it makes m_I = 6–16 H_i marginal (mediator debris R ≈ 2–8 against 1.14 for the instant release [E]); m_I ≲ 5
   (or ≳ 30) is quiet.
3. **The preparation result is a theorem, with the symmetry relabelled.** For the listed conditions (energy, neutrality,
   occupations, entropies, mutual information, stationarity) the relevant symmetry is antiunitary conjugation in the
   occupation basis, of time-reversal type, not C or CP (the toy has no antiparticle modes). Under a real Hamiltonian all
   such conditions are even and the torque is odd, so none can select the sign. Energy extremisation does select a
   phase, the torque-free one (T = 6 dE/dδ). The mutual-information row is imposed by the product form, not tested.
   Converged torques: 1.35764502 (single), 2.7152900 (doubled).
4. **The thermal vector source is R_total = 9.68×10⁵ at ε_BX = 10⁻⁷, not 9.84×10⁵** (the MB time-profile closure is 1.6%
   high; two momentum-resolved solvers agree to 10⁻⁶). The raw-source bound R_total ≥ 40.5626 B^{3/2} + 39184(1 − B) is
   confirmed. The Higgs source (b = 5×10⁻¹¹: Y_pair = 2.237×10⁻⁸, f_rest = 9.39) and the cold leakage (Γt₀ = 7.08×10⁻¹⁷)
   reproduce.

## Corrections to our own statements

- hw_verify's finite-release row (Γ_φ/H_i = 1: ε = −0.279, 1.16 Y_Q, R = 0.03) stands in its model but not in the
  covariant one: the ξ_X R curvature mass averages −2.5H² in the matter era (oscillating between −10 and +5 H²), not
  −10H². hw_critic's "a finite release is favourable (R → 0.03, |ε| → 0.38)" is withdrawn; |ε| stays 0.15–0.22 and R stays
  small only off the Y resonance.
- rp_portals' Higgs-portal headroom compared a source injected at onset (h → ss peaks at T ≈ m_h/2.5 ≈ 50 GeV, during the
  φ transition) with a threshold built for pairs injected before onset. The robust part is the unprocessed rest mass:
  f_rest = 9.4 at b = 5×10⁻¹¹ and 37–87 at detector-level b ≈ 1–1.5×10⁻¹⁰ [E], so a detectable Higgs portal needs
  ≥ 97–99% erasure of pairs injected at onset. That is an uncomputed obligation, not a closure.

## What this changes for the heavy-mediator window

| Aspect | Status after e160356 and this review |
|---|---|
| Release history | The quiet transfer survives a physical oscillating inflaton with the covariant curvature mass. New stated requirement: reheating within about one Hubble time (T_RH ≈ 2–5×10¹⁴ GeV) and m_I ≲ 5 H_i or ≳ 30 to avoid the Y resonance |
| Preparation | Unchanged in substance, now with a theorem: no conjugation-even condition can select the sign. A holder with a supplied CP phase remains the only route; its switch-off at the end of de Sitter is still assumed |
| Portals | New obligation: any kinetic mixing ε_BX adds a thermal vector source R ≈ 7.2×10²⁰ ε² (momentum-resolved), so ε_BX ≤ 1.4–2.5×10⁻¹⁰ under the isentropic proxy and none at e = 10⁻⁷ under the energy-conserving one. ε_BX = 0 is radiatively consistent (forbidden at all orders by the hidden charge conjugation of the stated tree Lagrangian). A detectable Higgs portal injects f_rest = 37–87 at onset: unresolved |
| Status line | Exists under homogeneous classical release and under a physical reheating on the stated line. Requires a supplied phase, a supplied holder, Γ_I ≳ 0.2 m_I with m_I ≲ 5 H_i, and one undecided threshold (kin_*). Not a cosmological point. Closure supplied no new law (Codex's sentence is accurate) |

## The decisive next calculation, other than kin_*

Linear finite-k stability of X, Y, Σ and the longitudinal vector on the physical reheating backgrounds: Floquet/mode
equations on the recorded homogeneous solutions in the quiet corner (m_I = 3–5 H_i, Γ_I = 0.2–0.3 m_I) and on the m_I = 10
runs, through the first inflaton oscillations where the X curvature mass swings between −10H² and +5H² at frequency 2m_I
with |X| ≈ 30–76 H_i, and through the 2m_I ≈ ω_YL resonance. It decides whether the homogeneous quiet transfer and the
1.0 Y_Q yield survive inhomogeneous growth, the gate the packet leaves uncertified. Runner-up: the onset-epoch fate of
Higgs-source pairs, which decides whether any detector-level b is compatible with the dark-matter budget.
