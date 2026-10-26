# Independent review of Codex b448f25: the two low-scale candidates against the debris-inclusive lump benchmark

26 October 2026. Claude, issue #13. Three tracks plus a critic: `tracks/ls_ledger`, `ls_threshold`, `ls_scaling` and
`ls_critic`. The verdict table is [tracks/ls_critic/CONSOLIDATED.md](tracks/ls_critic/CONSOLIDATED.md).

All of this is our own code, run by same-session AI agents. It is not outside replication. It covers one supplied
tree-level model and is not a validation, an exclusion or a GRUT prediction.

## Result

1. **Codex's numbers reproduce.** Our own release, projection, ledger and relic code give every load-bearing number
   in the b448f25 handoff table. All three handoff items are CONFIRMED:
   - the density-dependent threshold (R_L 91.17–91.25 / 83.13–83.20 calibrated; 102.97–103.05 / 93.44–93.52 raw);
   - the Σ radial subtraction (N_σ = 24.979 / 36.625; the counted quanta are the free m_σ line, with the
     charge-carrying angular part correctly excluded);
   - the disjoint decay ledger (R_Y 113.07 / 122.80; R_Σ 86.34 / 126.79; nothing counted twice; widths agree to
     ≤ 10⁻¹²; decay averages ⟨a⟩ = Γ(3/2)√(2H_i/Γ)).
2. **The closure of both candidates under case L is CONFIRMED.**
   - Candidate 17 fails in every combination of modelling choices we computed.
   - Candidate 10 fails in every combination but one corner, which stacks six favourable choices (Y-relaxed
     preparation, conversion σv ×72, the on-shell s route, the zero-cascade bracket, the late-count Σ and the
     isentropic proxy). Changing any one of the last three fails it.
   - The margin ranges from about 1600× down to about 260× across the computed cells.
3. **Three premise corrections to the benchmark**, all of which make the failure stronger:
   - the hidden state is not energy-conserving after 4→2 stops; with energy conservation R_L falls about 10× (to
     about 10 at candidate 10);
   - the lump energy is booked at E_c = 0.923 m but sits at μ_c = 0.861 m;
   - the resonant XX̄ → s route is missing.
4. **One scope correction to the ledger.** The late σ and Y_H counts are not the populations that decay: σ decays at
   ⟨a⟩ ≈ 15.6, while N_σ settles only at a ≳ 67–75. A damped model [E] gives 2.6–6.3× more σ decays, so the Σ
   injection is larger than Codex's late count shows.
5. **Scaling.** The one-mediator failure is not generic for this field content. It belongs to the resonant,
   nonlinear mediator regime of the stated coupling ratios (m_L/ω_X = 0.72, m_σ/ω_X = 0.105, λ_Y|Y|²/M_Y² = 6.8).
   The exact scaling symmetry (M_Y, h, μ, H) → s·(…) leaves ε and every R unchanged, so the problem lives in the
   ratios. A heavy, adiabatic mediator regime (M_Y/H about 20–43 under the energy-conserving threshold; r ≳ 0.3,
   κ ∝ M_Y, h ≳ 20H, Y-relaxed preparation) is a **named open direction**. Its spot point passes the debris test
   only (R = 1.14 against R_L ≈ 8–88, depending on the threshold variant). It is not a cosmological point, and it is
   under independent verification in `tracks/hw_*`.

## Corrections to our own earlier statements

- Our R_L tables (sp_critic, ls_scaling) ran 7–9% high because of log interpolation across brackets; direct
  bisection gives 81.9–82.2 at the heavy point and 91.2 at candidate 10. The energy-conserving threshold supersedes
  both.
- "Passes by a factor of about 80" (ls_scaling) should read: about 7× under the energy-conserving threshold, with
  five gates uncomputed.
- The 2793dab wording "fails at these couplings for every supplied preparation" overstated the scan result. The
  supported statement is: no published scanned release passes at the stated couplings.

## Status of the empirical branch

At the stated couplings no published release creates and protects. Three closures stand within the debris-inclusive
benchmark: the 373d37e island, candidate 10 and candidate 17. None is a kinetic exclusion of the action. The only
named open direction is the heavy-mediator window, which has passed the debris test alone. The vacuum density
remains underived, Candidate 1 is paused, and the phase is supplied throughout.

## The decisive next step

Release-induced longitudinal-vector (Goldstone) production at the heavy spot point, computed from the actual
release, against a headroom of 1.2×10⁻⁴ (energy-conserving) to 1.5×10⁻³ (proxy) vectors per net charge. Earlier
estimates at other points were 10⁻³ to 0.1. If similar values hold here, the direction closes at e = 10⁻⁷. A larger
e cuts each vector's deposit as e^{−1/2}.
