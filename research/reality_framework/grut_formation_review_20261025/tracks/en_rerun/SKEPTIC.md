# en_rerun_skeptic: adversarial check of the en_rerun findings on Codex 7d9c93b (matched early nonlinear pilot)

Track key `en_rerun_skeptic`, 27 October 2026. This note tries to refute the en_rerun verdicts in NOTES.md and
`results/en_rerun_verdicts.json`. The subject is Codex's packet `research/reality_framework/grut_early_nonlinear_20261027/`
at 7d9c93b. It was read from a fresh `git archive` in `$SP/review/en_rerun_skeptic/git7d`, which is identical to
`$SP/e7d` (`diff -rq` is clean, so en_rerun's restore held). SP = /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad.

**Bottom line.**
- **The headline numbers survive.** Every Codex number en_rerun confirmed comes out the same with our own code. That
  includes an independently written 3D integrator, which reproduces Codex's N48 trajectory to ≤ 8×10⁻¹³ through τ 85
  and to 6.5×10⁻¹¹ over τ 600–800.
- **We refute three of en_rerun's own statements, with recomputed numbers.**
  - Claim 3: "the curves agree within 4% to τ 100 and separate after about τ 150" is wrong. Separation starts at
    fragmentation, around τ 70–85.
  - Claim 4 and N2: the τ 600–800 ranges and means of the reference projections are biased.
    - They come from 20-τ snapshots, about two per orbit, so the orbit oscillation is aliased.
    - The mean of 1.98 also double-counts the duplicate τ-800 snapshot.
    - Dense sampling gives a massless orbit mean of 1.88 (range 1.18–2.77) and a Hartree mean of 1.166 (range
      1.07–1.24).
  - Claim 1b: "the cycle average never reaches ½ by τ 800" depends on the definition. The cycle mean of Codex's own
    instantaneous fraction crosses ½ at τ ≈ 455–535.
- **Two qualifications.**
  - The N3 "floor" holds only within the ω-only family. A common per-shell canonical shear lowers it from 1.060 to
    0.889.
  - The N1 range |χ̄|² = 314 ↔ 20 344 is the first orbit. The late orbit is about 650 ↔ 10 200.
- **We found five new issues with Codex's packet** (section 4).

Nothing here validates or excludes anything. It is all same-model numerical work on the reduced X/φ pilot.

## 1. What is independent of what

Labels: [C] computed, [D] derived, [A] argued, [E] estimate.

**Independent code.** All of it is in `skeptic/code/`. No Codex module is imported and no Codex binary is called.
- `own_engine.py` + `proj_own.py`: our own N48 integrator.
  - The FD4 Laplacian is applied spectrally through its exact symbol; Codex uses a real-space stencil.
  - The time stepper is our own Yoshida-4 implementation, on the same output schedule. τ is accumulated exactly as a
    running double sum from 20.05.
  - Our own diagnostics and projections use rfft with explicit half-space weights; en_rerun used a full complex fftn.
- `floquet_reduced.py`: Floquet exponents from the monodromy over one radial period, in the co-rotating frame. This is
  a different method from en_rerun's τ 20→55 singular values.
- `prep_checks.py`: a per-shell χ² test of the realised τ-20 field. This is a different test from en_rerun's whitening.
- `orbit_skeptic.py`: η-weighted cycle averages with both apocentre and pericentre conventions, and whole-orbit
  windows.
- `onset_compare.py`, `dense_projection_stats.py`, `controls_tau800.py`: comparisons.

**Inputs.**
- Codex's published files (FIELD_CURVES, FIELD_RESULTS, LINEAR_PHI_CONTROL, MATCH_L12).
- en_rerun's reproduced Codex-engine traces and snapshots. These were first checked against Codex's **published**
  FIELD_CURVES: the max relative difference over all seven runs is ≤ 1.5×10⁻¹¹ (`orbit_skeptic.json`).

**Own-engine runs [C].**
1. From the τ = 20 Codex-engine state to τ 85. Published N48 curves are reproduced to ≤ 7.7×10⁻¹³; the energy drift
   is 1.254×10⁻⁸, the same as Codex.
2. From the τ-600 snapshot to τ 800, with projections every Δτ 0.5. The trace is reproduced to ≤ 6.5×10⁻¹¹, and the
   τ-800 fields to 2.4×10⁻¹¹.
Each 100-τ chunk takes 131 s on 2 threads.

## 2. Verdicts on en_rerun's findings

| ID | en_rerun | Skeptic result | Refuted? |
|---|---|---|---|
| 1a | contrast = 1 at 60.50–60.90; first orbit-mean ≥ 1 cycle 56.65–69.6 (1.11–1.16) | Own engine: 60.60 (sample), 60.55 (interpolated). Pericentre cycle η-mean 1.10–1.15 (box16 1.103). With apocentre cycles the first cycle ≥ 1 is 62.95–76.55 (1.37–1.44) | No |
| 1b | X variance ½ at 69.20–69.25; cycle fraction 2.6–3.0%; cycle average never ≥ ½ (max 0.466–0.496) | 69.25 (sample) / 69.21 (interpolated), with a second crossing at 83.05. 2.65–3.01% holds only as a ratio of means over the pericentre cycle. Other conventions: mean of the fraction 9.0–9.7% (pericentre cycle), 13.5–14.4% (apocentre cycle). **The cycle mean of the instantaneous fraction exceeds ½ from τ 455–535 (max 0.509–0.534)** | Verdict stands; "never ½" is definition-dependent |
| 2 | τ 800 60.94–64.93%, 1.054–1.159; 700–800 means 60.1–63.0%, 1.07–1.18; last-orbit swing 47–75% | Reproduced: 700–800 means 60.14–62.98%, 1.068–1.179. Last two whole orbits 59.3–62.3%, 1.071–1.188. Since the last pericentre 46.9–74.5% | No |
| 3 | 0.277–0.365%; control 2.0917 (max 2.109, first ≥ 1 at 702.5); τ-20 8.032e-4 vs 8.038e-4; **agree within 4% to τ 100, separate after ~150** | Numbers confirmed (max 2.1094 at τ 796). **The pointwise ratio leaves 4% at τ 72.5, 10% at 86.5, 20% at 102 and 50% at 119.** Per-orbit mean ratios: 1.019 (56.65–69.6), 1.045 (69.6–84), 1.096 (84–100), 1.225 (100–118), 1.52 (118–138), 1.97 (138–159). en_rerun has no code behind the "4%/150" sentence | **Yes** (that sentence) |
| 4 | 2.237515 / 1.179962; **600–800 massless 1.25–2.64 (mean 1.98), Hartree 1.08–1.22 (mean 1.16)** | τ 800 confirmed with our own rfft code. Dense Δτ 0.5 history (403 samples): **massless 1.178–2.774, mean over the whole orbits 621–762.5 1.880** (per orbit 1.846 / 1.875 / 1.916); **Hartree 1.073–1.242, 1.166**. The en_rerun set has 12 rows with τ 800 twice (`state_s0000039/40`); without the duplicate the mean is 1.952 | **Yes** (range/mean) |
| 5 | 0.547%, 0.977/1.031%, 1.482/2.717%, 0.206%, 0.447%; scope 2 momenta, e = 0, τ ≤ 50 | Published = regenerated JSON. Read matching_check.py: k = 2.13, 7.1, e = 0 on both sides, τ 20–50 | No |
| 6 | 1.293e-8; ≤ 4.3e-14; halfstep ratio 16.04 | VERIFICATION diff: energy 1.293141e-8 in both; current 1.87e-14 vs 4.29e-14 (round-off) | No |
| D1 | bitwise | DENSE_CORRECTIONS / CORRECTED_CLOCKS identical except the wall-clock seconds | No |
| A-a | engine correct | An independent implementation reproduces the trajectory (see §1). a(τ) vs dense_b3 is 2.06e-6 at τ 240; no full-background comparison exists beyond 240 | No |
| A-b | preparation correct | Purity 3.56e-11. X-marginal symplectic eigenvalues 0.5000011–0.51999. Realised X variance 11.883 vs 11.440 (z = 0.75, N_eff = 375). Per-shell power over 186 shells: z mean 0.11, sd 0.92, KS p = 0.43 | No |
| A-c | variants insensitive (1.180 / 1.183 / 1.231) | Reproduced at τ 800 only. Over 600–800, dropping the zero mode changes the count by −6% to +34% (mean +6%); radial by −2% to +11% | No (qualified) |
| N1 | \|χ̄\|² 314 ↔ 20 344; common pericentres | Common pericentres confirmed, even at τ 762.2–762.5 (seed2 761.5). **314 ↔ 20 344 is the τ < 70 orbit.** The last whole orbit is 650 ↔ 10 195 (N48); across controls 649–846 ↔ 9 404–10 313 | No (numbers corrected) |
| N2 | IR split 1.874 / 0.836, k ≥ 1.8 "agree" 0.364 / 0.344; 600–800 swings | τ-800 split, k < 1 = 1.066 and box16 1.403 all confirmed. The time behaviour is corrected as in claim 4. For k ≥ 1.8 the massless/Hartree ratio over 600–800 is 0.96–1.36 (mean 1.135), not "agree". The total ratio is 0.96–2.54 (orbit mean 1.61) | **Yes** (time statements) |
| N3 | per-shell floor 1.060; 0.93–1.09 over 600–800; Hartree within ~11% | 1.060, and 1.056–1.142 across controls, confirmed. Dense 600–800: **0.917–1.116**. It is a floor only for one ω per shell. Add one U(1)-covariant canonical shear per shell (π → π − γχ, which leaves the charge term unchanged) and the floor is **0.889** at τ 800: 0.880–0.981 across controls, 0.844–0.912 over 600–800, and 0.666 per mode. Hartree is 1.20–1.43× this floor | No (qualified) |
| N4 | k = 1.737 grows 22.2× vs a band maximum of 22.4× at 1.9 | Own Floquet: μ(1.737) = 0.601 vs a maximum of 0.633 at k = 2.15. The whole range 0.1–5.95 is unstable, and the L = 12 box shells have μ = 0.23–0.60. Asymptotic growth over τ 20–55 is 5.9× vs 6.5×. en_rerun's 22× is the singular value of a non-canonical (q, p) map, inflated by transients. The conclusion stands | No |
| N5 | trajectory thread-independent; restarts ≤ 1e-12 | Our own engine, restarted at τ 600, agrees to 6.5e-11 at τ 800 | No |

**Side note.** "The machine has 2 cores" is wrong: lscpu shows 4 cores with 1 thread per core, and the affinity mask is
0–3. GNU `nproc` honours `OMP_NUM_THREADS=2`, which produced the "2".

## 3. Details behind the refutations

- **Claim 3 [C].** We aligned LINEAR_PHI_CONTROL (Δτ 0.5) with the published N48 φ + interaction curve.
  - The two curves agree at τ 20 (ratio 0.9993) and track each other through the first breakup cycle.
  - They part as soon as the charged field fragments.
  - Codex's own control energy oscillates with its homogeneous orbit:
    - The τ-800 value of 2.09 sits next to the orbital maximum (2.109 at τ 796).
    - The orbit mean first exceeds 1 only in the cycle τ 729–773 (mean 1.12; instantaneous crossing at 702.5).
  - The nonlinear-versus-control separation is 675× at τ 800 and about 360× on orbit means near 750. Codex's
    qualitative claim is unaffected.
- **Claim 4 / N2 [C].** The late radial period of the coherent mean is 44–48 τ: pericentres at 621.0, 666.5, 713.5
  and 762.5. The 20-τ snapshots therefore alias the orbit oscillation.
  - The massless count follows the orbit: correlation +0.94 with |χ̄|², peaking near apocentre. The Hartree count is
    anti-correlated (−0.86), because its m² includes the oscillating condensate (14.8–25.8).
  - Orbit means: massless 1.880, Hartree 1.166, ω-only shell floor 1.022, ω + shear floor 0.875.
- **Claim 1b [C].**
  - The instantaneous fraction (X2 − |χ̄|²)/X2 is close to 1 near pericentre and small near apocentre.
  - Its cycle mean and the ratio of cycle means differ by a factor 3–4 early on, and they straddle ½ late.
  - en_rerun chose the ratio of means without saying the alternative reverses the conclusion.

## 4. New issues with the Codex packet (missed by en_rerun)

1. **The massless value and the "roughly 1.9-fold" gap are single orbit-phase numbers.** REPORT.md's "2.238" and its
   "roughly 1.9-fold difference … much larger than the grid/time changes" are taken at τ 800. On orbit average the
   massless value is 1.88 and the ratio is 1.61. The instantaneous ratio runs 0.96–2.54, so the orbit-phase swing is
   larger than the gap itself. [C]
2. **The control's 2.09× is a near-peak sample.** REPORT.md says the control "extrapolates to 2.09 times the initial
   active energy by 800". That is the orbital peak; the orbit mean is about 1.1 at τ ~750. The control and the
   nonlinear run separate from fragmentation (τ ≈ 70–85), which ties the changed neutral forecast directly to breakup.
   [C]
3. **Both reference counts sit mostly in modes that were never seeded.**
   - k < 1.8 holds 1.87 of 2.24 (massless) and 0.84 of 1.18 (Hartree).
   - These modes start at exactly zero in the finite-band state. On the reduced orbit they are Floquet-unstable
     (μ = 0.23–0.60, against a band maximum of 0.63).
   - Their content is box-limited: in the L = 16 box, k < 1 holds 1.40 (80 modes) against 1.07 (26 modes) at L = 12,
     and the totals are 2.41 / 1.27.
   - So the exported "particle-like" content is mostly nonlinear transfer into box-scale modes whose initial state
     Codex did not specify. Any later matching (PARTICLE_MATCHING.md step 3) needs an IR seed and a box-size control.
     en_rerun's N2 and N4 each touch half of this but do not connect them. [C/A]
4. **One L = 16 run cannot separate volume from realisation variance.**
   - The L = 16 shifts relative to N48 are Q₀/Q +0.092, massless +0.17 and Hartree +0.086.
   - The L = 12 seed-to-seed shifts are +0.104, −0.18 and +0.032, which is comparable.
   - The L = 16 band also starts at k = 1.842, because n² = 21 gives 1.7996 < 1.8; at L = 12 it starts at 1.814.
   - Codex states that seed and volume dependence "remain visible". With one seed per box, the size of the volume
     effect is not measured. [C]
5. **The Hartree reference depends on the condensate phase.** Its m² contains the oscillating condensate, so
   variants that look interchangeable at τ 800 (with or without the zero mode, radial, tangential) differ by up to 34%
   at other orbit phases. [C]

## 5. Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
S=research/reality_framework/grut_formation_review_20261025/tracks/en_rerun/skeptic
G=$SP/review/en_rerun_skeptic/git7d/research/reality_framework/grut_early_nonlinear_20261027   # git archive 7d9c93b
W=$SP/review/en_rerun/ework; export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2; cd $SP/review/en_rerun_skeptic/work
python3 -I $S/code/orbit_skeptic.py $G/FIELD_CURVES.npz $W $S/results/orbit_skeptic.json N48 N64 N32 halfstep48 seed2_48 cut10_48 box16_64
python3 -I $S/code/own_engine.py $W/N48/state_s0000000.bin 85 1 20 own_20_85.json own_85.bin
python3 -I $S/code/onset_compare.py own_20_85.json $G/FIELD_CURVES.npz $S/results/onset_own_engine_N48.json
python3 -I $S/code/own_engine.py $W/N48/state_s0000029.bin 700 2 10 dense_600_700.json own_700.bin
python3 -I $S/code/own_engine.py own_700.bin 800 2 10 dense_700_800.json own_800.bin
python3 -I $S/code/dense_projection_stats.py $W/N48/state.tsv $W/N48/state_s0000039.bin own_800.bin $S/results/dense_projection_N48_600_800.json dense_600_700.json dense_700_800.json
python3 -I $S/code/controls_tau800.py $G/FIELD_RESULTS.json $W $S/results/controls_tau800.json N48 N64 N32 halfstep48 seed2_48 cut10_48 box16_64
python3 -I $S/code/floquet_reduced.py $G/MATCH_L12.npz $S/results/floquet_reduced.json
python3 -I $S/code/prep_checks.py $G/MATCH_L12.npz $W/N48/initial.bin $S/results/prep_checks.json
```

The results are in `skeptic/results/`; `skeptic_summary.json` holds the machine-readable verdicts. Everything ran on 2
threads, and total compute was about 7 minutes.

## 6. Limits

- The dense projection history covers N48 only, over τ 600–800 (three whole orbits). Other controls were checked only
  at their τ-800 snapshots.
- The own engine shares Codex's discretisation and time composition. It checks the implementation, not the
  discretisation. en_indep uses a different integrator.
- Matching beyond two momenta, gauge and heavy-field dynamics after τ 50, and any physical n̄ or f_sym are not checked.
- The shear-extended floor is still an equal-time projection family. It is not a physical count.

Credits: classical-statistical lattice method (Khlebnikov–Tkachev; Micha–Tkachev); UV limits of the
classical-statistical approximation (Epelbaum–Gelis–Wu, arXiv:1402.0115); Yoshida (1990) fourth-order composition.
