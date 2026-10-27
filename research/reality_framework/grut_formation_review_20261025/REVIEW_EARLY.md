# Review of Codex 7d9c93b: matched early nonlinear X/φ pilot

27 October 2026. Claude, issue #13. Brief: [BRIEF_EARLY.md](BRIEF_EARLY.md). Four review tracks, each with an
adversarial skeptic, then a completeness critic: `tracks/en_rerun`, `en_indep`, `en_neutral`, `en_particles`,
`en_critic`. Each track's NOTES.md and SKEPTIC.md give methods, commands and full numbers.

This is same-session AI work on one supplied, reduced model. It is not outside replication. Nothing here is
observational evidence, a validation, an exclusion or a GRUT prediction.

## Result

**Every Codex number reproduces, and the code audit found no error.** All seven histories regenerate with Codex's
engine to at most 3.5×10⁻¹² relative. Three independently written integrators follow the published N48 history to at
most 7×10⁻⁶ through τ = 800:
- spectral FD4 with Blanes–Moan SRKN;
- spectral FD4 with our own Yoshida;
- PEFRL at N32.

The critic integrated the reduced-versus-full matching errors with its own code; they agree to 1.4×10⁻⁹. Codex's
qualitative conclusion stands: charged fragmentation at τ ≈ 60 cuts off the neutral resonance, and our fixed-background
"φ catches up at τ 665–697" does not describe the nonlinear continuation.

Three readings need correcting:
- the τ = 800 values are single orbit-phase samples;
- the neutral fraction depends on the unselected X seed amplitude;
- the "massless versus Hartree" particle question is the wrong one.

The particle analysis gives a new, conserved antiparticle budget, set out below.

## Verdicts

| Item | Codex | Ours | Verdict |
|---|---|---|---|
| Dense corrections: 10% warning; maximum through τ 80 | 56.45–56.475; 20.174% | bitwise; our earlier 56.46 and 20.2% | CONFIRMED |
| Reduced-versus-full matching | 0.547%; 0.98–1.03%; 1.48–2.72%; force 0.206%; frozen-Y ≤ 0.45% | 0.5471%; 0.977/1.031%; 2.717/1.482%; 0.2057%; exact/bound 0.38–1.00. Own integration agrees to ≤ 1.4×10⁻⁹. Two momenta, τ ≤ 50, e = 0 on both sides | CONFIRMED |
| Gauge omission | e = 0 | Coulomb energy of the measured charge fluctuations ≤ 2×10⁻¹⁰ of the total [E] | CONFIRMED as justified |
| initial.py / match.py | stated covariance and sampling | state pure to 3.6×10⁻¹¹. The generated field's whitened covariance has diagonal 0.983–1.009 (σ 0.012), and a χ² test per shell passes. No factor-2 error | CONFIRMED |
| early3d.cpp | adapted from 8159556 | forces, FD4, Yoshida and exact current conservation all correct. a(s) only sets the clock and outputs. τ(s) = 20 + a_m s + ½Cs² follows from dτ = a ds. After τ ≈ 40 the step is set by sampling (Δs = 0.05/a) | CONFIRMED |
| Fixed-charge energy bound (3/4)λ^{1/3}\|n_Q\|^{4/3} | homogeneous lower bound | holds for **every** lattice state at fixed charge (Cauchy–Schwarz, FD4 symbol ≥ 0, Jensen). It is 11.594% of the energy in all seven histories | CONFIRMED, stronger |
| Charge-contrast RMS = mean | τ 60.5–60.9 | 60.50–60.90 (three codes); its period mean is 1.14 over τ 57–70, so the clock is robust | CONFIRMED |
| X fluctuation variance = ½ | τ 69.2 | 69.21–69.28. This is an instantaneous crossing at pericentre, where \|χ̄\|² is 2% of its maximum; the fraction drops back to 0.06 by τ 75 | CONFIRMED, as a pericentre crossing |
| Energy outside the mean; Q₀/Q at τ 800 | 61–65%; 1.05–1.16 | 60.94–64.93%; 1.054–1.159. Over one late orbit: 46.7–72.7% and 1.05–1.21 (N48); all runs 45.8–74.5% | CONFIRMED as same-phase samples (see 1) |
| φ + interaction at τ 800 | 0.28–0.36% | 0.277–0.365%. N32 alone sets the low end (k ≈ 7 band under-resolved) | CONFIRMED, conditional on the X seed (see 2) |
| Fixed-orbit linear control | 2.09× | 2.0917128 (own code). Grid-dependent: 0.995 / 2.092 / 2.197 / 2.255 at N32 / N48 / N64 / continuum. Ensemble at N48: 2.18, unit crossing τ 698.5, the same physics as our 696.5–697 | CONFIRMED |
| References: massless 2.238, Hartree 1.180 | τ 800 | 2.237515 / 1.179962 (own code) | CONFIRMED at τ 800 |
| "Roughly 1.9-fold" difference | stated | phase sample: the ratio runs 0.96–2.54 over an orbit, mean about 1.6 | CORRECTED |
| Verification groups | energy drift 1.3×10⁻⁸; current drift 1.9×10⁻¹⁴ | 1.293×10⁻⁸; ≤ 4.3×10⁻¹⁴ (thread-dependent round-off) | CONFIRMED |
| Seeded zero point | not given | 5.3×10⁻⁴ X̄ per net charge, so negligible for the count. But only 21% of the late φ energy is above a half quantum (below) | NEW |

## Corrections

1. **The τ = 800 values sit on an eccentric coherent orbit.** The spatial mean still orbits with r_max/r_min ≈ 3–4
   and a period of Δη ≈ 0.89, which is 44–48 τ late. Every single-instant diagnostic swings with it. All seven
   histories end at the same orbit phase (0.25–0.27 of a cycle past the last maximum), so the control-to-control
   comparisons remain valid. Their spreads match the spreads of last-cycle means, with the extremes in the same
   places. But the headline ranges and the export are phase samples.
   - **Ask:** quote τ 700–800 cycle means and envelopes, not τ = 800 values. The N48 cycle mean of the energy outside
     the spatial mean is 60.7%.
   - **Ask:** publish a period-averaged spectral export (Δs ≳ 1.8, i.e. Δτ ≳ 100 near τ 800).
   - The isotropic Hartree mass includes the coherent \|χ̄\|², so it also swings, from 15 to 25.
2. **The neutral result depends on the X seed amplitude, which is not selected.** We scaled only the X fluctuations
   of the published τ = 20 field by a factor f (Codex's engine, N48):

   | f | Fragmentation τ | φ + interaction at τ 800 (τ 700–800 mean) |
   |---|---:|---:|
   | 1 | 60.6 | 0.31% (0.31%) |
   | 0.1 | 137.8 | 0.80% (0.81%) |
   | 0.01 | 230.8 | 2.64% (2.40%) |
   | 0 | 370.3 | 5.78% (6.67%), maximum 7.9% |

   - The qualitative claim holds: there is no catch-up toward 2.09× even with no X seed.
   - The value 0.28–0.36% belongs to the vacuum-like X seed and rises about 3× per decade of smaller seed.
   - Most of that value is not above-vacuum energy. At τ 792–800, φ energy above a half quantum per shell is only
     0.067–0.072% of ρ₀. The seed supplies 25–36% of the final value.
   - 64% of the late φ energy sits in unseeded k_lat > 8 shells with classical occupation 0.06–0.66. That is where
     the classical-statistical approximation fails (Epelbaum–Gelis–Wu, arXiv:1402.0115).
3. **Mechanism of the suppression.** Both mechanisms work, and either one alone suffices. We ran a 2×2 factorial of
   homogeneous same-seed proxies, giving φ + interaction at τ 800:
   - fixed orbit: 2.09;
   - X-fluctuation (Hartree) detuning only: 0.48%;
   - pump depletion only (a free orbit following the measured E₀, Q₀): 0.32–0.87%;
   - both: 0.21%.

   The curves separate from the linear control at fragmentation, τ ≈ 70–90. The forecast change is established by
   τ ≈ 200.
4. **The infrared is unseeded and box-limited.** The seed band starts at k = 1.8, but modes just below it are
   Floquet-unstable (μ 0.23–0.60, against a band maximum of 0.63). k = 1.737 grows 22.2× linearly, against 22.4× at
   the band maximum. Yet k < 1.8 carries most of both reference counts: 1.87 of 2.24 and 0.84 of 1.18. In the L = 16
   box, k < 1 holds 1.40. A single L = 16 run cannot separate volume from realisation variance.
   - **Ask:** run lower-cutoff and box-size controls.
5. **Small points.** The two codes never diverge chaotically. Long-wavelength perturbations grow only about 5× in X and
   20× in φ over τ 100–800 (≤ 0.15 per unit η), so a pointwise comparison works to τ 800. The 2.09× control exceeds
   the total active energy after τ 702.5, as Codex says.

## The antiparticle budget (new, [C] within the classical reduced pilot)

We measured the actual spectra (statistical and kicked-response unequal-time correlators) in late windows that are long
in conformal time:
- The response peaks follow ω² = k² + m² with m² ≈ 19, close to the Hartree 20.
- The count-weighted statistical branches give m²_eff ≈ 8–13.
- Most importantly, the time-averaged count is **flat in the chosen mass**: within 3% for any m² from 7 to 20, and
  within 6% from 5 to 26.

So the massless reference is wrong, because there is a mass gap. Any massive reference gives the same answer.

| Component | X̄ per net charge |
|---|---|
| Waves (nonzero modes), τ 340–800 | 1.12–1.24 (seed 2: 1.17–1.19; L16: 1.23–1.24) |
| Counter-rotating part of the coherent mean | 0.48–0.61 late (1.05 at τ 100–255) |
| **Total** | **1.74–1.76, constant over 22 windows from τ 20 to 1460, in L12, seed 2 and L16** |
| τ = 20 homogeneous orbit alone (exact Fourier series) | 1.7387 |
| Across X seed amplitudes | total 1.65–1.81 (the split shifts) |

- **The total X̄ number is conserved.** Fragmentation only moves antiparticles from the eccentric coherent orbit into
  waves, with an inverse cascade toward k ≈ 1. This is what number-conserving 2↔2 scattering does in classical
  turbulence (Micha–Tkachev).
- **The orbit's eccentricity sets the number.** As a check [E]: a harmonic orbit ψ = Ae^{−iωt} + Be^{iωt} with
  r_max/r_min = 8.06 holds B²/(A² − B²) = 1.55 X̄ per net charge. The quartic orbit gives 1.74.
- **Against our earlier estimates** (REVIEW_TIMING):
  - "pairs near k_peak" (1.4–2.9) contains the total but not the waves alone;
  - the equilibrated end point (0.10–0.12) is not approached by τ 1460;
  - the energy outside the mean is about 61%, not the 88% allowed;
  - every measured value lies inside the kin_source bracket of 0.13–2.6.

**The question is now relaxation, not reference frequency.** Within this pilot, about 1.74 X̄ per net charge are fixed
at preparation and conserved through τ 1460. Whether they relax toward the equilibrated value (about 0.1) before
binding depends on the number-changing rate with physical masses. That rate is bounded here only to ≳ 5×10³
conformal units, against roughly 10¹¹–10¹² until the binding transition. Two other things remain:
- a physical rule for whether the coherent mode's counter-rotating part counts as X̄;
- the infrared and box controls of item 4.

## Corrections to our own earlier statements

- REVIEW_TIMING's "φ catches up at τ 665–697" is a fixed-background result only. Nonlinearly, φ stays below 1% with
  the vacuum-like X seed, and below 8% with no X seed.
- REVIEW_TIMING's X̄ estimates: the waves sit at the low end of "gas at ω". The total is set by orbit eccentricity and
  conserved, so the equilibrated end point does not apply in this era.

## Status

Codex's packet is a correct, well-verified reduced pilot. It supplies the first computed, conserved antiparticle
budget (≈ 1.74 X̄ per net charge, classical, reduced, one history family). It replaces the five labelled alternatives
as the leading input to the kinetic source bracket, conditional on the items above.

The following remain open:
- the physical n̄;
- number-changing relaxation;
- f_sym and the threshold;
- the full X/Y/Σ/gauge and quantum evolution;
- the φ bias and binding transition;
- the vacuum density.

Candidate 1 is unchanged.
