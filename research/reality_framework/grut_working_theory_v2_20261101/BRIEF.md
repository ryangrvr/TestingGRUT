# Brief: repair GRUT-M1 where it can be repaired, and compute the corrected physical picture (1 November 2026)

Owner request: "fix it, work out the working theory and correct the physical picture." This pass computes the repairs
and the corrected chain. A synthesis (WORKING_THEORY.md in this directory) is written afterwards from your results.
Codex's packets are read through git only; never edit them.

## What is established (read these first; do not redo them)

- research/reality_framework/grut_formation_review_20261025/REVIEW_M1.md: the M1 review (one vacuum, the Q = 2×10⁵
  object, calibration, holder, transport, frozen mass test 0.0667805(1−η) + 2.13356p ≤ 0.01; no computed population
  passes yet; holder EFT scale M_Y/ξ ≈ H_i; parsec dip unobservable; gas–core scattering omitted).
- REVIEW_TRANSPORT.md (same directory): the 1.75 X̄ budget is a prethermal plateau; φ channels relax it before onset
  (10–10⁴ e-folds); equilibrium X̄ per charge 0.069–0.096 at the pilot's energy, 0.19–0.26 at R = 1.14, 0.37–0.50 at
  R = 4.05, 2.3–2.9 at R = 40.6–48.5; the φ transition does not return J_r (untilted: random split; M1 tilt:
  resonance capture, no slow limit, sensitive to ε within ±5%).
- REVIEW_KINETIC.md, HANDOFF_KINETIC.md: case-L gas-only annihilation from a forced chemical-equilibrium start gives
  f = 0.0046 at total heat R = 4.97 (the e = 10⁻⁶ band), 0.0048 at R = 5, 0.029 at R = 7, 1.49 at R = 41.7 (e = 10⁻⁷);
  p = f/2.02855.
- research/reality_framework/grut_closure_forcing_20261027/SYNTHESIS.md: closure supplies no law; λ_Y = λ_Σ = 1 hits a
  Landau pole at 5.7×10¹⁵ GeV and the band λ_Y = λ_Σ ∈ [0.24, 0.57] keeps the window; Higgs survival needs ξ_H ≈
  0.05–1 (M1 supplies 0.1).
- Codex's M1: d4747cb (grut_matched_candidate_20261027), ed73791 (grut_m1_dynamics_20261027), bf8ab3b, c481530,
  94f660d, f8186fa (grut_m1_stabilized_transport_20261029: Option 1A holder (ξR/M_Y)(e^{iα}Σ*Y² + h.c.),
  ξ = 30.74672874, α = 2.19214407, μ_eff(12H_i²) = 10i H_i).

## Our engines and codes (own code; reuse freely, copy before modifying)

- `research/reality_framework/grut_formation_review_20261025/tracks/rp_reheating/code/rp.py` + `rpreheat.c`:
  homogeneous Jordan-frame reheating with an oscillating inflaton (mode 0; m_I = b H_i, Γ_I = gam H_i), the X, Y, Σ
  potential with κ and μ, ξ_X R; prepare() solves the conditional minimum at held phases; run() integrates. Its runs
  reproduce Codex (yields 1.0408 Y_Q,old at τ = 800 for (3, 0.5)).
- `tracks/m1_holder/code/rpholder.c`: the same engine with μ_eff(R) (Option 1A); reproduces +1.806458%.
- `tracks/kin_source`, `kin_rates`, `kin_thermo`, `kin_solve`: the case-L kinetic network (sources, rates including
  capture kernels and object branch E(Q), Bose state functions, moments solver with formation rules).
- `tracks/pt_m1tilt/code` and `skeptic/code`: homogeneous X + φ transition under the M1 potential.
- `tracks/m1_vacuum/code`: M1 Q-ball BVP solver (E(Q), Ω(Q), profiles), bulk state, dispersion.
- `tracks/m1_transport/code`, `skeptic/code`: cross sections, X–core three-channel scattering, free streaming.

## Units and key numbers

m_X = g v = 0.0924089637 GeV, g = 0.1, v = 10 m_X, λ_X = 10⁻³, λ_φ = 0.03, ε = 0.002; object Q = 2×10⁵ with
E/Q = 0.9374 m_X; Q_min = 5.67×10⁴, binding above 7.13×10⁴; bulk density 449.31 m_X³, spinodal 233.2 m_X³.
H_i = 1.342842735×10¹¹ GeV (M1), M_Y = 30 H_i, h = 40 H_i, μ = 7.159456 H_i, κ = 3.4587×10⁻³; inflaton
(m_I, Γ_I)/H_i = (3, 0.5); Y_Q = 5.07941423×10⁻⁹; D = 4.4×10⁻¹⁰ GeV per unit entropy.

## Rules

- Write only in `research/reality_framework/grut_working_theory_v2_20261101/tracks/<key>/` (NOTES.md, code/,
  results/*.json; nothing over 2 MB); scratch in /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/<key>/.
- Git is read-only. No commits.
- At most 2 threads or worker processes; other jobs share the 4 cores. No 3D lattice runs above 48³. Run in the
  foreground (Bash timeout up to 600000 ms per call; chunk long runs with restart files); never end with a
  background job running.
- Claim labels [D] [C] [Q] [E] [A] [R]. A repair is a model choice, not a consequence of closure; say so. Never claim
  validation, exclusion of the real world or a new theory. Data comparisons are retrospective.
- Return a digest of at most about 700 words.
