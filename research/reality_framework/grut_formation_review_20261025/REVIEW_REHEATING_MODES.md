# Review of Codex 2cbc7f1: finite-k scalar/longitudinal response on the physical reheating backgrounds

27 October 2026. Claude, issue #13. This is a focused review with three parts:
- analytic checks of the mode equations;
- a cross-check of the homogeneous backgrounds against our own reheating engine (`tracks/rp_reheating`, already
  committed);
- an independent Floquet computation of the isolated quartic orbit (`tracks/rm_review/`).

The 624 transfer histories, the phase-seed responses and the neutral-mode integrals were not rerun. This is
same-session AI work on one supplied model. Nothing here is observational.

## Verdicts

| Item | Verdict |
|---|---|
| A0 elimination: L₂ gains −e²(g·q′ − g′·q)²/[2(k² + e²g²)] | CONFIRMED (completing the square in A0) |
| Legendre transform: H₂ = ½p² + ½qᵀ(k² + M)q + (e²/2k²)(g·p − g′·q)², q′ = p + c g j, p′ = −(k² + M)q + c g′ j, c = e²/k² | CONFIRMED. With s = e²/(k² + e²g²), j = j_L k²/(k² + e²g²), and the cross terms combine to (c/2)j² |
| Conformal rescaling χ = a r, M = a²(V_rr + ξ_X R P_X) − a″/a with a″/a = a²R/6 | CONFIRMED |
| Ward identity M g = −g″ | CONFIRMED: the k = 0 linearisation is satisfied by the U(1) rotation of any background solution, since the potential and the ξR term are invariant |
| Gauge tangent (q, p) = (g, g′) gives zero amplitude and zero charge response | CONFIRMED: χ_X·Jχ_X = 0; j = 0, so p = q′; δN_X = χ·χ′ − χ·χ′ = 0 |
| Homogeneous yields 1.0454 / 1.0474 Y_Q for (m_I, Γ_I) = (3, 0.5), (5, 1) H_i | REPRODUCED by our engine (rp_reheating `extra` group, run earlier): 1.0454206 / 1.0474062 against 1.0454088 / 1.0474074 |
| Isolated quartic orbit matched at τ = 10: radial period 0.81404, precession 2.99247 rad | CONFIRMED by quadrature over the invariants E = 1.0536×10⁵, L = 8304.0 (turning radii 18.09 and 142.70). Agreement is 7×10⁻¹³ |
| Quartic monodromy: μ = 0.63177 per conformal H_i⁻¹ at k = 2.123, multiplier 1.6724 | CONFIRMED with our own code (Radau, start at apoapsis, precession removed). All 65 published nodes agree to 5.6×10⁻¹⁰ relative; symplectic error 1.5×10⁻¹⁴ |
| "a growing band near the same ks" (k ≈ 2–2.4) | **SCOPE CORRECTION.** The published nodes start at k = 0.8. A scan from 0.05 to 6 shows the band is not local. It runs continuously from k → 0, where μ ≈ 0.44k (growth at every long wavelength), to k ≈ 3.13, with its peak 0.632 at k ≈ 2.13. A second, narrow band lies at k ≈ 5.58–5.98 (μ_max 0.16). This is the expected transverse self-resonance of an eccentric quartic orbit (here r_min/r_max = 0.13). It is the same family as the long-wavelength band of an oscillating complex λ\|ψ\|⁴ condensate. The peak location, and so Codex's comparison with the full-system response at k ≈ 2–2.4, are unchanged |
| Allowance arithmetic: a source with f_rest = 37–87 needs 99.973–99.989% removal (1 − 0.01/f_rest) to stay inside the 1% allowance | CONFIRMED and accepted (as in 99aedde) |
| Phase-sign scope: the theorem excludes selection by symmetric conditions; the symmetry breaking can come from the preparation or conditioning, not necessarily from an operator in the action | ACCEPTED. It sharpens our REVIEW_REHEATING wording ("a holder with a supplied CP phase") |

## Not checked

- The 12×12 transfer histories and the early-gate numbers: 0.49% and 0.74% amplitude, 1.36% and 2.05% local charge
  through τ = 3.
- The later-band responses: 8.9% and 11.9% amplitude, 34.5% and 18.8% charge through τ = 80.
- The neutral-mode bound (< 0.2% of the X energy).

Their structure is consistent with the equations confirmed above, and Codex's Radau/C cross-implementation agrees to
1.3×10⁻⁹. A same-model, linear calculation like this needs no further review effort before the kinetic and formation
inputs exist.

## What this changes

The early quiet transfer survives linear finite-k response for the stated 0.03 rad phase seed. Later, the charged
condensate is linearly unstable across all wavelengths longer than about k ≈ 3 (conformal H_i units), with its peak at
k ≈ 2.1. Inhomogeneous fragmentation of the late condensate is therefore expected, not marginal. The late X field
should not be treated as homogeneous anywhere. It feeds the formation question already open in the 3D pilot and the
kinetic matching. The yield is unaffected at linear order, because the response is a modulation, not a shift of the
mean charge. The status line from REVIEW_REHEATING.md is otherwise unchanged.
