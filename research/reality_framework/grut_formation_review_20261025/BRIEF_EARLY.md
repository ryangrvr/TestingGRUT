# Brief: review of Codex 7d9c93b (matched early nonlinear X/φ pilot), issue #13

27 October 2026. The Codex packet is `research/reality_framework/grut_early_nonlinear_20261027/` at commit 7d9c93b;
read it through git only. It accepts our review ec1eb7f (REVIEW_TIMING.md) and builds a reduced scalar-limit 3D pilot.

**What Codex did.**
- A Gaussian conformal state on the band 1.8 ≤ k/H_i ≤ 8 is propagated with the full X/Y/Σ six-field tangent system
  plus the φ system, from τ = 0 to 20, on the quiet (3, 0.5) history.
- The light marginal is then evolved nonlinearly in 3D with
  H = ½π² + ½|∇χ|² + (λ_X/4)|χ_X|⁴ + (g²/2)χ_φ²|χ_X|² + (λ_φ/4)χ_φ⁴, with λ_X = 10⁻³, g² = 0.01, λ_φ = 0.03.
- The lattice uses FD4 with Yoshida time stepping and a radiation clock a(s) = a_m + Cs, τ(s) = 20 + a_m s + ½Cs².
- Seven histories run to τ = 800: N = 32, 48, 64 at L = 12; a half time step; a second seed; a seed cutoff of 10;
  and L = 16 with N = 64.

**What Codex claims.**
1. Charge-contrast RMS reaches the mean at τ ≈ 60.5–60.9, and X fluctuation variance reaches half at τ ≈ 69.2.
2. At τ = 800, 61–65% of the active energy lies outside the spatial mean. The mean charge is 1.05–1.16× the net
   charge.
3. Neutral energy plus interaction is 0.28–0.36% at τ = 800. The same φ seed on a fixed homogeneous reduced X orbit
   instead extrapolates to 2.09× the initial active energy. So nonlinear continuation changes the neutral forecast,
   including our "φ catches up at τ 665–697".
4. Reference projections give 2.238 negative-charge quanta per net current (massless) and 1.180 (isotropic Hartree).
   Neither is a physical count.
5. Matching errors of the reduced against the full model through τ 50: 0.55% background, about 1% charged response,
   1.5–2.7% neutral transfer. The frozen-Σ / quadratic-Y estimate is about 0.45% or less.
6. Seven verification groups pass: energy drift 1.3×10⁻⁸, current drift 1.9×10⁻¹⁴, export reconstruction at
   round-off.

**Codex's request (CLAUDE_HANDOFF.md).**
- Check MATHEMATICS.md, match.py and initial.py, the early3d adaptation and the export conventions.
- Reproduce the short engine and seed checks first, then one main history if useful.
- Check the full-versus-reduced scope and the same-seed neutral comparison, especially the energy conventions.

The rules of [BRIEF.md](BRIEF.md) apply:
- write only in `tracks/<key>/`;
- git read-only, and no commits;
- label independent code versus reproduction;
- give a verdict per claim, with claim labels [D] [C] [Q] [E] [A] [R];
- never claim validation, exclusion or new theory;
- return a digest of at most about 700 words.

Credit the classical-statistical literature where it is used: Khlebnikov–Tkachev; Micha–Tkachev; Epelbaum–Gelis–Wu
arXiv:1402.0115.

## Prepared environment (already done)

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
# Codex tree (git archive of 7d9c93b, whole research/reality_framework): $SP/e7d/research/reality_framework/
export GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025
export GRUT_KINETIC_WORKDIR=$SP/kwork GRUT_ONSET_WORKDIR=$SP/owork
export GRUT_TIMING_WORKDIR=$SP/review/<key>/twork GRUT_MODE_WORKDIR=$SP/review/<key>/mwork
export GRUT_EARLY_WORKDIR=$SP/review/<key>/ework
cd $SP/e7d/research/reality_framework/grut_early_nonlinear_20261027/code && python3 -B -c "import en_common"
```

Codex scripts save JSON or NPZ into the packet directory.
- Only the track `en_rerun` may run them in place. It must first copy the published outputs to
  `$SP/review/en_rerun/published/`.
- Other tracks that need a Codex script, for example initial.py to build an initial field file, must copy the whole
  `grut_early_nonlinear_20261027` directory to `$SP/review/<key>/pkt/grut_early_nonlinear_20261027`. Sibling imports
  resolve through ROOT.parent, so symlink the sibling packets it imports (grut_source_timing_20261026,
  grut_reheating_modes_20261026, grut_reheating_preparation_20261026, grut_kinetic_continuation_20261026 and any
  others) next to the copy, or run from a copied parent tree.
- Never write inside `$SP/e7d` except from en_rerun.

Timing: the main N48 history to τ 800 takes about 93 s on 4 threads, N64 about 285 s. The machine has 4 cores shared
by four tracks, so use at most 2 threads each.

Useful own material:
- REVIEW_TIMING.md and `tracks/st_*`: notes, Floquet code, and the critic's neutral extension (φ/X = 1 at τ 665–697
  with Codex's 2cbc7f1 neutral control and its occupation and zero-point convention);
- `tracks/rm_review/code/quartic_floquet.py`;
- `tracks/kin_source/results/source_interface.json` (X̄ bracket 0.13–2.6).
