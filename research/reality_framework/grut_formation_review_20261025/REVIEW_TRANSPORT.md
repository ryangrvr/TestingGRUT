# Review of Codex 5f5513b: averaged classical source, infrared controls, coherent matching

1 November 2026. Claude, issue #13. Brief: [BRIEF_TRANSPORT.md](BRIEF_TRANSPORT.md). Four tracks, each with an
adversarial skeptic: `tracks/pt_rerun`, `pt_matching`, `pt_relax`, `pt_response`; plus `tracks/pt_m1tilt` (the same
transition under Codex's later M1 potential; its skeptic is still running, so that part is **preliminary**). The
synthesis below replaces the critic, whose run the usage limit stopped. Each track's NOTES.md and SKEPTIC.md give
methods, commands and full numbers.

This is same-session AI work on one supplied, reduced classical model. It is not outside replication. Nothing here is
observational evidence, a validation, an exclusion or a GRUT prediction.

## Result

**Everything Codex computed reproduces, and the bookkeeping is right.** The one IR12 history regenerates to 3×10⁻¹³;
the infrared preparation is a pure canonical Gaussian with the right normalisation; the averaged export keeps the
coherent outer product correctly; the Hann statistic returns the free-field antiparticle occupation; the coherent
matching numbers (J_r = 1.763725, lab statistic 1.737788) and the slow-ramp results are exact.

Three of Codex's interpretations need correcting, and the review adds the calculation Codex named as deciding:

1. **"Slow matching returns J_r"** holds for an externally imposed smooth mass ramp. It does not hold for the actual
   φ-condensation transition, which releases φ out of a Floquet-unstable state and does not conserve the radial
   action, whether the ramp is slow or not.
2. **The matrix-response range (0.44–1.70)** is a single-phase sample, and the canonical sum is an identity of the
   sampled data, so neither diagnoses the medium.
3. **The ≈1.75 budget is a prethermal plateau.** Chemical relaxation through the φ channels has 10–10⁴ e-folds to act
   before binding, so the population at binding is the equilibrium one at the actual energy per charge: about
   0.07–0.10 X̄ per net charge at the pilot's own energy, 2.3–2.9 if the e = 10⁻⁷ vector heat thermalises. Neither
   1.75 nor "about 0.1" is universal.

## Verdicts on Codex's claims

| Codex claim | Ours | Verdict |
|---|---|---|
| Eight-history budget 1.7505–1.7517 (s 14–24) and 1.7413–1.7454 (s 24–34) | IR12_s1_w1 regenerated: 1.7505527 / 1.7414468 (3×10⁻¹³); every window, spectrum and average array to ≤ 3.5×10⁻¹³. Adding the separately reported k > 8 tail gives 1.7530–1.7546 / 1.7501–1.7518, the window-to-window decline shrinks to −0.14…−0.17%, and the L12 refinement difference becomes −0.048% / −0.056%. The tail itself is grid dependent (0.0034 at N32 vs 0.0016 at N48) | CONFIRMED |
| Infrared preparation: pure Gaussian, Ω = √(k² + σ²) in the unstable directions; handoff RMS 4.3% / 2.2% | Symplectic eigenvalues ½ to 3×10⁻¹⁵; whitened covariance of the generated fields 0.999–1.003; handoff RMS 4.28–4.65% / 2.17–2.29%, continuum 4.44% / 2.24%, insensitive to σ from 0.1 to 30. New: the per-mode variance jumps ×2.4 at k = 1.414 and ×5.8 at 1.683 (an integrable 1/Ω squeeze), and the interval (1.683, 1.8), which holds about 11% of the continuum X variance, is undersampled by both boxes | CONFIRMED, with a note |
| Averaged export keeps the coherent outer product; charge and energy reconstruct | ⟨v₀v₀ᵀ⟩ from the trace to 1.5×10⁻¹⁵; charge 5.7×10⁻¹⁵, energy 7.8×10⁻¹⁶. The coherent part carries 113–129% of the net charge; squaring the mean would keep 0.4% | CONFIRMED |
| Hann statistic = free antiparticle occupation | Synthetic tones: ≤ 1.5×10⁻³ for m ≥ 1; leakage bias 5% at m = 0.4, 53% at m = 0.2. Exact orbit at 109 phases 1.737785–1.737794 | CONFIRMED |
| Coherent matching: lab 1.737788 vs J_r 1.763725 per coherent charge; slow ramp returns J_r; sudden 0.077–18.18 | Both numbers to 10⁻¹³. Derived: J_r − A₋ = Σ_{n≥1} n ω_n\|c_n\|² + Σ_{n≤−2}(\|n\|−1)\|ω_n\|\|c_n\|² ≥ 0, equal to the higher-harmonic weight (0.025937), zero only for the harmonic ellipse. C^∞ ramps converge super-polynomially (3×10⁻¹⁰ at T = 40). Sudden range at 64 phases is 0.013–18.30 and scales with the arbitrary m_out | CONFIRMED; sudden range CORRECTED |
| "No extra counting postulate is needed" | Agreed: in the free massive basis n₋ is fixed by the state. What remains is the dynamics that connects the orbit to that basis | CONFIRMED |
| "Slow matching to massive particles gives about 1.76" | True for an imposed m(s) ramp. For the actual transition (φ condensing out of φ = 0 while the orbit runs, in the untilted potential) the outcome is an ε-independent random split: n₋/q = 0.70 ± 0.02 X̄ plus 0.47 φ quanta per charge for the τ-20 orbit (J_r = 1.764); 0.197 ± 0.004 plus 0.124 for the late orbit (J_r = 0.48). Independent integrators agree run by run at ε = 10⁻² and diverge chaotically below 10⁻³; no trend from ε = 3×10⁻³ down to 3×10⁻⁶ | CORRECTED in scope |
| Matrix response: conjugate/normal 0.44–1.70 in the IR; canonical sum within 0.4%; free control 2.2×10⁻¹⁰ | Free control 2.1619489×10⁻¹⁰ (50-digit Yoshida map). The ratio spans 0.28–1.31 over one radial period, is orbit-locked (s₀ = 20.92 reproduces s₀ = 20) and depends on the delay window; it is not a medium constant. The canonical sum is an identity of the sampled data (≤ 3.4×10⁻⁵ on the Nyquist grid); the 0.1–0.4% is weight beyond \|ω\| = 20. The initial slopes equal 1 − D²(k² + ½ tr H̄_X)/6 + O(D³), so "consistent with ω² = k² + Hartree m²" is an identity, not a dispersion test | CORRECTED |
| What the response means | A circular-condensate Bogoliubov calculation gives the right order of the IR ratio (0.71–0.80 at k = 0.52, 0.36–0.45 at 1.57). A Floquet-driven HFB closure reproduces the UV poles (≤ 2.4% for k ≥ 3.1, ≤ 0.6% for k ≥ 4.7) and typical ratios within 5–15%. In the IR (k ≤ 2.1, 83% of the wave X̄ action) every homogeneous closure fails: the bare closure is tachyonic (Goldstone defect −6.0) and HFB gapped (+5.0); repaired closures fix the IR poles but overpredict the UV ratio ×1.9–2.7; momentum scattering ≥ 1.07× the forward amplitude within one period is model-independent evidence of inhomogeneity | NEW |
| Counting in the quasiparticle basis | A₊ − A₋ = Q holds to 10⁻³; S/Q = 4.48–4.50, so the stable budget is the total frequency action. J_r of the late zero mode exceeds the lab statistic by +0.5–1.5% like-for-like; quasiparticle-basis total 1.7625–1.7668. X ↔ X̄ conversion within a wave mode is exchange with the condensate (total conserved), so X̄ is not a quasiparticle quantum number of the wave sector but the sum is conserved | NEW |
| Chemical relaxation before binding is the deciding calculation | Agreed, and computed in part (below) | CONFIRMED |

## The deciding calculation: relaxation before binding (`tracks/pt_relax` and its skeptic)

- **Time available [D].** The reduced Hamiltonian has no explicit conformal-time dependence, so every rate is constant
  in conformal units. φ condenses when λ_φv²a² = g²⟨r²⟩, at a_onset = 5.0–9.2×10¹², i.e. s_onset = 2.5–4.6×10¹²,
  about 10¹¹ times Codex's s = 34. The budget survives only if Γ_chem < 2.2–4.0×10⁻¹³ per unit s.
- **Rates [E, corrected by the skeptic's explicit collision integral].** The track's parametric 2↔4 brackets were
  80–10⁵ times too high. With a tree-level Monte Carlo (Bose factors, thermal masses) the e-folds available before
  onset are: pure-X 2↔4 0.005 (negligible); g²-only 0.2–3 (marginal); mixed g²λ_φ 18–250; φ self 1.5×10³–2×10⁴;
  XX̄ ↔ φφ 5×10⁶. So relaxation runs through φ with a margin of 10–10⁴, and in a φ-free sub-model the budget would
  survive. Complex-X combinatorics (O(1–10)) and Bose-enhanced rates in the over-occupied phase are not computed.
- **Direct classical test [C].** Continuing IR12 to s = 1004: the total falls 1.751 → 1.651 (s ≈ 200) → 1.202 (≈ 500)
  → 0.832 (≈ 1000) at N32. The rate is cutoff dependent (N32's early rate is 3.2× the continuum estimate), and the
  classical UV catastrophe accelerates the late loss on finer grids, so no lattice rate bounds the physical one in
  either direction. The decay itself is real; the coherent part does not decay monotonically through s = 504.
- **Equilibrium end point [C].** At fixed n_Q and energy, with a condensate: 0.069–0.096 X̄ per net charge (0.110
  with μ = 0 on the X̄ branch); T_hidden ≈ 0.33–0.45 GeV at onset (T/m_X = 3.6–4.9). With the heavy-window release
  heat: R = 1.14 gives 0.19–0.26; R = 40.6–48.5 (e = 10⁻⁷) gives 2.3–2.9.

**Reading.** The ≈1.75 budget is what the eccentric orbit prepares; it is not protected. The population handed to the
annihilation era is set by the actual energy per charge at onset, which the vector heat dominates at e = 10⁻⁷.

## The transition itself (`tracks/pt_matching`; `tracks/pt_m1tilt`, preliminary)

- The untilted potential: φ = 0 is a Floquet-unstable state on the eccentric orbit (onset at 0.733 M_c,peri; with
  2g²/λ_X = 20 = n(n+1) the φ_k equation is exactly Lamé n = 4, with four gaps; gap 2 is the early pilot's k ≈ 7
  neutral band). The loss of radial action comes from φ growing out of that state, not from φ's inertia: a linear tilt
  β = 0.3 conserves J_r, β ≤ 0.1 gives the chaotic loss, β = 1 gives deterministic resonance capture.
- **Under Codex's M1 potential** (V_φ = λ_φ/4(φ² − v²)² − εv(φ − v)³, ε = 0.002; `tracks/pt_m1tilt`, 767 runs, own
  integrator; skeptic pending): the chaotic split disappears and the transition is deterministic. It is adiabatic for
  ε_rate ≥ 10⁻⁴ (τ-20 orbit) / 10⁻³ (late orbit), returning J_r to 2×10⁻⁴. Below a sharp threshold it **deterministically
  loses radial action, and the loss grows as the rate slows**: at ε_rate = 10⁻⁵–10⁻⁶ the τ-20 orbit gives 1.49–1.52 X̄
  (0.85 J_r) plus 0.16–0.18 φ quanta per charge; the late orbit 0.11–0.15 (0.23–0.30 J_r) plus 0.22–0.25. The
  mechanism is a tilt-driven φ oscillation slaved to the radial frequency up to ≈5.6 M_peri (τ-20) / 2 M_peri (late),
  with a free φ component released where that regime ends. A linear tilt of the same size (β = 0.2) releases nothing,
  so the loss is specific to the M1 shape. The physical rate is ε ~ 10⁻¹²–10⁻¹³, ten decades beyond the scan [E].
- Waves: with the actual φ(s) the wave modes are parametrically amplified (\|β_k\|² ≫ 1 for k from 0.25 to 8 by
  4–8 M_peri), so the wave switch-on is not adiabatic either and needs a lattice treatment of the transition.

## Corrections to our own statements

- REVIEW_EARLY's "conserved total ≈ 1.74": it is conserved over the pilot's 22 windows, but it is a prethermal plateau,
  not a protected number.
- REVIEW_EARLY's "count flat in the reference mass within 3%": about 4% on the infrared histories.
- `tracks/pt_relax`'s 2↔4 brackets (ours): corrected by its skeptic as above.
- The late coherent component's matching ambiguity is −0.055…+0.023 per net charge (−3.2%…+1.3% of the total),
  10–40× the eight-history control spread; the zero mode is not an adiabatic central orbit late (its instantaneous J_r
  swings ×3–5 within a window).

## What Codex's next calculation should compute

1. Γ_chem in conformal units against 1/s_onset, with the φ channels (XX̄ ↔ φφ, g²λ_φ mixed 2↔4) as the leading ones,
   including complex-X combinatorics and Bose enhancement.
2. Quantum kinetics (Boltzmann with a condensate, or 2PI) of X, X̄, φ and the eccentric condensate from the s = 34
   source, using classical runs only while occupations are ≫ 1; the IR wave sector needs a beyond-homogeneous
   closure or the field itself.
3. The timed Y/vector injections and their thermalisation, which set the energy per charge at onset.
4. The φ-condensation transition on the lattice (coherent orbit and waves together), under the actual potential.
5. Hand the resulting X̄ at T_h/m_X ≈ 3–8 to the case-L annihilation network, with the φ-transition heat.

For M1's frozen test (p ≤ 0.0047 at η = 1), the relevant number is the leftover after annihilation of an initial
0.07–0.1 (or 2–3 with the vector heat), not the 1.75 plateau.

## Status

Codex's packet is correct as computed. Its leading statistic is confirmed and now has a physical reading: a prethermal
budget that the φ channels relax before binding, and that the mass-generating transition itself redistributes.
Physical n̄, f_sym, the threshold, the full model and binding remain open. Candidate 1 is unchanged.
