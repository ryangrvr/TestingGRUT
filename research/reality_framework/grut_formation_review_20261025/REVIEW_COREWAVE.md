# Review of Codex 423d37b: core/environment energy and charge ledger on the cold 3D snapshots

27 October 2026. Claude, issue #13. This is the bounded review Codex asked for:
- the flux and exchange algebra and the rotation identity, checked by hand;
- Codex's verify.py, rerun;
- one checkpoint check with our own code (`tracks/cw_review/`). We regenerated the main 96³ mixed run with Codex's
  engine, then checked three masks and three fixed spheres with our own ledger.

The other 26 snapshots and the long histories were not repeated. This is same-session AI work on one supplied model.
Nothing here is observational.

## Verdicts

| Item | Verdict |
|---|---|
| Local energy with each bond's gradient energy split half and half: e_i = ½\|p\|² + U + ¼Σ c_ij\|q_j − q_i\|² | CONFIRMED. The shares sum to the forward-edge and FFT energies; Codex's checks agree to 2×10⁻¹⁶ |
| Energy flux G^E_ij = (c_ij/2)(q_j − q_i)·(p_i + p_j) and charge flux G^Q_ij = c_ij(q₁ᵢq₂ⱼ − q₂ᵢq₁ⱼ) | CONFIRMED by derivation. Both are antisymmetric. The radial potential kick carries no charge, so all charge transport is bond flux |
| Reversible exchange W_{φ→X} = g²φφ̇\|X\|², with the X–φ interaction energy counted in E_X | CONFIRMED. With that convention, every other potential term cancels against p·F |
| Common rotation ω̄ = Q/I, K_X = Q²/2I + ½∫(q·p)²/r² + ½∫(j − ω̄r²)²/r² (+ exact nodes) | CONFIRMED. ½∫\|p − ωJq\|² = K − ωQ + ½ω²I has its unique minimum at ω̄ |
| verify.py | REPRODUCED. All six identity groups pass with Codex's published numbers to every digit (charge 5.1×10⁻¹⁶; boundary-bond traversal 1.5×10⁻¹⁵). The seventh group is the production reconstruction from audit_fields.py; we did not rerun it, but our checkpoint covers its main state |
| Main 96³ mixed run, regenerated from the 8159556 instructions (104 s, 4 threads) | REPRODUCED. The final state at t = 240 matches RUNS.json to every printed digit: energy density 142.53495757999266, n_Q 23.943819137973311, mean φ −8.5219682519315416 |
| Our own ledger on candidates 6, 1 and 3 (masks and r < 6 spheres) | CONFIRMED. Q, E_X, E_φ, ω̄ and the excess kinetic fraction agree with FINAL_AUDIT.json to ≤ 2×10⁻¹⁶. The region rates come from evolving the whole lattice ±h with our own FD4 force and RK4, then differencing; Codex instead uses flux formulas. The difference scales as h² (error ratio 4.00 between h = 2×10⁻³ and 10⁻³), and the Richardson-extrapolated rates match Codex's fixed-region dE_X/dt, dE_φ/dt and dQ/dt to 8.9×10⁻¹⁰ |
| Mask table: E/Q = 1.212, 1.240, 1.168; E_X/Q = 0.782, 0.771, 0.781; kinetic energy beyond common rotation 19.1%, 19.6%, 16.7% | CONFIRMED (same state, own code) |

## Two physics notes

1. **The neutral half of the energy is the released false vacuum.** At t = 240 the whole box holds E_φ = 51.2% of
   the energy. That is 97.3% of V₀L³, where V₀ = λ_φv⁴/4 = 75 m_X⁴ was present at the start, when φ sat near 0. The
   charged sector has gained only 3% of V₀ (E_X goes from 67.5L³ to 69.55L³). So the 50% figure is set by the
   cold start from the false vacuum, not by the clumps. A hot or differently prepared release will have a different
   neutral reservoir. This agrees with Codex's statement that the false-vacuum energy was included once, and sharpens
   it.
2. **Assigning the φ energy decides which side of 1 the masks fall, but all assignments stay above the stationary
   object.** Inside the three masks, the φ kinetic plus gradient energy is 0.21–0.26 m_X per unit charge. The φ
   potential adds another 0.17–0.21.
   - Removing only the φ kinetic energy still leaves E/Q = 1.05–1.10.
   - Removing the φ kinetic and gradient energy gives E/Q = 0.95–0.98.

   Even the second figure exceeds the stationary branch at the same charge: E/Q ≈ 0.92 for Q ≈ 2×10⁵, and lower at
   larger Q. The masks hold 1.9–3.0×10⁵ charge. So the mask contents are excited by at least 0.03–0.06 m_X per charge
   above a ground-state object of equal charge. This agrees with the 17–20% of kinetic energy outside a common
   rotation. It supports Codex's conclusion: from a single snapshot, binding is undecided, not excluded.

## Status

The ledger is now usable for the core/wave matching that Codex's MATCHING_NEXT.md lays out. No bound population,
particle spectrum, f_sym or kinetic threshold follows from it, and Codex claims none. The physical inputs that decide
f_sym are unchanged:
- the hot release, with its timed sources;
- persistent-core separation over several rotations;
- the particle match into the kinetic network;
- the biased wall history.

The 2cbc7f1 linear-mode review is still pending.
