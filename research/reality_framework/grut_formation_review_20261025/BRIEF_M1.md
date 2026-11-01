# Brief: focused review of Codex's GRUT-M1 chain (27–30 October 2026), issue #13

1 November 2026. Claude was unavailable 27–31 October. Codex (with a third AI, "Argon"/Gemini, supplying notes)
built a model variant it calls **GRUT-M1**. The packets are on separate branches; read them through git only
(`git show <sha>:<path>`, `git archive`), never check them out or edit them. Run `git branch -r` and `git log --all
--since=2026-10-30` first: a newer isocurvature packet (adiabatic residual 4×10⁻¹² for a unit adiabatic mode; false
residual 3.02 without curvature forcing; hidden-scalar charge power 3.435×10⁻⁷⁰; 38/38 checks) may have been pushed
after f8186fa.

| Commit | Packet | Main claims |
|---|---|---|
| d4747cb | research/reality_framework/grut_matched_candidate_20261027/ (CLAUDE_HANDOFF.md, REPORT.md) | V_φ = λ_φ/4(φ²−v²)² − εv(φ−v)³, ε = 0.002, λ_φ = 0.03: V(v) = V′(v) = 0, vacuum curvature unchanged, derivative (u−1)[λu² + (λ−3ε)u + 3ε] with negative discriminant ⇒ +v the only stationary point ⇒ no φ domain walls. New Q = 200000 object: mass 17324.83235 GeV, E/(Qm_X) = 0.93739999, Ω = 0.88172049, central S = 0.12577; dE/dQ = ω to 2.9×10⁻⁷, dΩ/dQ < 0; fixed-Q radial Hessian positive (0.2233–0.2244); exact cold-fluid dispersion still grows for dilute densities (F1 peak Γ = 0.11428 m at k = 0.42 m). Calibration: H_i = 1.342842735×10¹¹ GeV (×1.00700905 from the old heavy point), Y_Q = 5.07941423×10⁻⁹, M_Y = 4.0285×10¹², h = 5.3714×10¹², μ = 9.614×10¹¹, κ = 3.4587×10⁻³. Charge loss Γt₀ = 6.27×10⁻¹⁷. Frozen mass test **0.0667805(1−η) + 2.13356p ≤ 0.01** (η = fraction of charge in reference objects, p = leftover pairs per charge): needs η ≥ 0.85026 at p = 0 or p ≤ 0.004687 at η = 1. ξ_H = 0.1 supplied; cutoff 10¹⁵ GeV; zero SM portals in the conservative point |
| ed73791 | grut_m1_dynamics_20261027/ | 13 periodic 3D runs: supplied object stays bound (E/Q 0.943; core 0.938); cold fields form merging compact regions, dilute masks hold 30–32% of charge; bulk Ω = 0.77586471, E/Q → 0.77586 at infinite size; the hot elliptic control (shape action 1.75, E/Qm = 5.79) is not spatially converged, 0.6% of charge in dense masks. 9 open spherical runs with an absorber: emitted negative charge 0.002/0.008/0.030/0.095 of Q for 5/10/20/40% pulses |
| bf8ab3b, c481530, 4fe963f, 12b5696, 8f67a83, 6dc9556, 6c504de | grut_m1_scattering_20261028/, grut_m1_kinetics_20261028/, grut_m1_expansion_clock_20261028/, later M1 packets | Γ(φ→XX̄) = 1.875659×10⁻⁴ m_X; three-channel radial scattering of X̄ on the core; conditional cold-completion mass excess 0.0275/0.1127/0.4903/2.482% for the four pulses; 10% pulse kinetics: N̄ 0.5589 from 1695, Q_core/Q₀ = 0.9419, cold excess 0.7266%; expanding-core ledger total 0.894%; a synchronous local source reaches ≤ 0.51% of volume |
| 94f660d | grut_m1_synthesis_20261029/REPORT.md | H/m_X = 1.93×10⁻²¹ at the 20 MeV clock; refit H_i −0.568% per 0.858% mass; sextic curvature holder ξ₆ = 4×10⁸, Λ = 10¹⁵ GeV gives δ = −2.10 (not −π/2) and +9.93% charge; free X: 6.595% of charge, v_rms = 0.0684c, 35.5 comoving pc to equality if collisionless, but Γ_coll/H = 15.7 at 2 MeV; CLASS did not converge |
| f8186fa | grut_m1_stabilized_transport_20261029/ (REPORT.md, DERIVATIONS.md, FOLLOWUP_NOTES.md, STABILIZED_CANDIDATE.json, TRANSFER_COMPARISON.json) | Option 1A holder: V ⊃ (ξR/M_Y)(e^{iα}Σ*Y² + h.c.), ξ = 30.74672874, α = 2.19214407, i.e. μ_eff(R)/H_i = μ/H_i + (10i − μ/H_i)R/(12H_i²); at R = 12H_i² the stationary point has amplitudes (87.19556, 5.75425, 28.38718)H_i, phase arg Σ − 6 arg X = −π/2 exactly, physical Hessian (13.0914, 14.5892, 538.147, 1546.66, 3239.89)H_i²; V bounded below at fixed R (≥ −285370 H_i⁴); late charge −8267.437 vs −8120.740 (+1.80645%) ⇒ H_i → 1.32691046×10¹¹ (−1.186%). Spectator charge power 1.236×10⁻¹⁰ (4 e-folds), 7.59×10⁻¹⁶ (8), 1.44×10⁻⁷⁰ (50). Gas: σ = 2.16627561×10⁻⁷/m_X², 7.91173 collisions per particle over a factor 10 in a, viscosity proxy crosses H at 0.62156 MeV, kurtosis 1.917 → 1.673; BGK hierarchy with 5.87674% gas: P_mix/P_cold = 0.997879 (628 pc), 0.979835 (209 pc), 0.893210 (62.8 pc), 0.849806 (20.9 pc), 0.764169 (6.3 pc). Codex: "a concrete conditional forecast that can now be challenged. It is not yet an earned, distinctive GRUT prediction" |

Codex's bounded request (d4747cb): "refute the one-vacuum argument and scalar BVP/Hessian first; then check the new
mass-normalization chain and included charge-loss normalization." The later packets were handed to Argon while
Claude was unavailable; no independent review of them exists.

## Our own results that bear on M1

- REVIEW_EARLY.md, REVIEW_TIMING.md and the 5f5513b review tracks (`tracks/pt_*`): the early classical X̄ budget is
  about 1.75 per net charge, set by the orbit's eccentricity; `tracks/pt_relax` finds it is a prethermal plateau that
  relaxes (classically, on the lattice) toward the equilibrium value of 0.07–0.11 per charge at the pilot's own energy,
  or 1.9–2.9 if the e = 10⁻⁷ vector heat thermalises; `tracks/pt_matching` finds the φ transition turns the coherent
  part into a random split (0.70 X̄ + 0.47 φ quanta per charge for the τ-20 orbit) rather than preserving the radial
  action; `tracks/pt_m1tilt` (running) tests this under M1's tilted potential.
- The kinetic case-L work (`tracks/kin_*`, REVIEW_KINETIC.md): gas-only annihilation erases pairs to f = 0.0046 at
  total heat R = 5 and 0.029 at R = 7; e = 10⁻⁷ (R ≈ 41) fails by 150×.
- The closure pass (research/reality_framework/grut_closure_forcing_20261027/SYNTHESIS.md): the heavy window needs
  ξ_H ≈ 0.05–1 for Higgs survival (Codex now supplies ξ_H = 0.1); λ_Y = λ_Σ = 1 hits a Landau pole at 5.7×10¹⁵ GeV;
  TCC would exclude the window; the stated Lagrangian has no CP phase, so the asymmetry's sign must be supplied (M1
  now supplies it through the complex coefficient e^{iα}).
- Our reheating engine: `tracks/rp_reheating/code/rp.py` (own code, reproduces Codex's e160356 and (3, 0.5)/(5, 1)
  runs); it contains the Y, Σ and μ terms and the Ricci scalar, so the Option 1A operator can be added as
  μ → μ_eff(R).
- Our Q-ball/object code: `tracks/kin_thermo` (object branch E(Q)), `tracks/gp_stability`, `tracks/independent_solver`.

## Rules

- Write only in `research/reality_framework/grut_formation_review_20261025/tracks/<key>/` (NOTES.md, code/,
  results/*.json; nothing over 2 MB); scratch in /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/<key>/.
- Git read-only; no commits. Label independent code versus reproduction. One verdict per claim
  (CONFIRMED / CORRECTED / REFUTED / NOT CHECKED / NEW) with your number beside Codex's. Claim labels [D] [C] [Q] [E]
  [A] [R]. Never claim validation, exclusion or new theory; data comparisons are retrospective.
- Compute budget: other jobs share the 4 cores. Use ONE thread (export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1), no
  3D lattice runs; ODE, BVP, algebra and small scans only. Run in the foreground.
- Return a digest of at most about 700 words.
