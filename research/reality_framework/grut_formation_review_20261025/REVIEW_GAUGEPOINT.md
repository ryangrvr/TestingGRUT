# Independent review of Codex b5174e5: the six requirements at the 7737487 Z6 gauge point

25 October 2026. Claude, issue #13.

**Sources.** Five tracks and a critic took part: `tracks/gp_matching`, `gp_phase`, `gp_stability`, `gp_walls`,
`gp_leftovers` and `gp_critic`. The verdict table is in [tracks/gp_critic/CONSOLIDATED.md](tracks/gp_critic/CONSOLIDATED.md).

**Scope.** All of this is our own code, written by same-session AI agents, and it is not outside replication. It
covers one supplied tree-level model with homogeneous calculations plus estimates. Nothing here is a validation, an
exclusion or a GRUT prediction.

## Result

On homogeneous calculations plus estimates, under the supplied release (H_inf = 4.7e13 GeV, c_H = 1), the stated point is
not a viable genesis point. Two requirements pass, two fail, and two remain open or partial.

| Requirement | Codex | Review |
|---|---|---|
| Potential stability | PASS (tree) | **PASS** at tree level and in a one-loop estimate. λ_Y reaches 4π near 1.4e15 GeV, close to the inflationary \|X\|, so perturbative control is marginal. A matching scale must be stated |
| Fixed initial phase | NOT SUPPLIED | **FAIL at this release.** The aligned zero-charge solution is exact but unstable: perturbations grow by e^29–e^30, which leaves local \|ε\| ≈ 0.1–0.2 with chaotic sign. Inflation does not hold the phase at its minimum (⟨cos ψ₂⟩ = 0.26). Neither a supplied DRT phase nor CP-violating higher operators produces a resolved net charge (8,832 runs, \|⟨ε⟩\| ≤ 0.034) |
| Bias to remove walls | ABSENT (variant passes) | **OPEN.** The tilt algebra is exact; it is Saikawa eq. 2.23 plus Δ/2. With a consistent g* the minimum T_ann is 6.6 MeV (central range 14–24 MeV), and the pass holds only under the standard scaling law. The cores attract walls by about 1.4 TeV, and wall–object capture has not been computed |
| Symmetric leftovers | OPEN | **FAIL on present estimates.** The net ε is ≤ 0.03, while ≥ 0.990 is needed (≥ 0.19 in the best equilibrium case). Two new sources are each sufficient on their own: longitudinal-vector production during inflation (about 1060 pairs per net charge) and the Y condensate (0.5 pairs per net charge if it never thermalizes) |
| Gauged objects | PARTIAL | **PARTIAL**, confirmed. The screened bound is 6.56e-20 GeV; the new unscreened bound is 3.4e-7 GeV, still 2.4e-10 of the margin. The fluctuation spectrum is open |
| Corrected survival | PASS | **PASS**: c₆ = 2.1068e5 and Γt₀ = 2.0e-18, robust to dimension 5–7 operators |

## Corrections for Codex

1. **Phase.** "Own alignment generates zero charge" holds for the mean only. The real subspace is unstable. The phase
   item is a fail at this release, not merely "not supplied".
2. **Isocurvature.** The phase masses do supply an isocurvature prediction: the light inflationary mode is the phase
   ψ₂ that sources the X charge. Since the release is chaotic, the operative failure is the absence of a correlated
   net charge.
3. **Phase count.** There are 0 physical CP phases at dimension ≤ 4 with Z2_φ, and 8 / 16 / 53 at dimension ≤ 5 / 6 /
   7. The wall tilt breaks Z2_φ, which makes φΣ*Y² renormalizable and adds a phase. The kinetic-mixing and CP audits
   must be redone.
4. **Mod 6.** The invariant is (Q_X + 3Q_Y) mod 6 on the whole closed slice. Drop the three-core grouping.
5. **Strings.** The Z6 strings have Gμ ≈ 3.1e-15, conditional on their forming, just below the LVK O3 limit of about
   4e-15.
6. **Walls.** Use 6.6 MeV, with a consistent g*. List wall–object capture as open.
7. **Leftovers.** The ledger contains no V₀. The debris, binding-pair and remainder energies are alternatives, not
   additive. The thresholds are 0.990 (non-relativistic count) and 0.988 (adiabatic). Add the longitudinal-vector and
   Y-condensate sources.
8. **Matching.** The radial companion is compulsory (|β₆| ≥ 2|α₆|). No angular-only AD window, including our earlier
   one, can be imported.

## The decisive next calculation

Scan the full three-field homogeneous release over H_inf, from the EFT-reliable range (about 3e11/√c_H GeV) up to
4.7e13 GeV, at c_H = 1 and 10. Record for each value:
- the net ε;
- its change under a 1e-9 rad shift, as a chaos indicator;
- the yield against Y_Q;
- the Y-condensate energy;
- the eaten-mode amplitude.

The construction survives only if some H_inf gives a release that is regular and productive (Y ≥ Y_Q) and safe for
symmetric leftovers (ε ≳ 0.19, or ≳ 0.01 with a demonstrated erasure). If none does, the one-mediator
"creates and protects" link fails at these couplings.

So far the reliable-EFT releases reach ε = 0.010–0.024 with Y = 0.013–0.42 Y_Q, so no computed regime meets both
conditions.
