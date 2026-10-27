# en_particles_skeptic: adversarial re-check of the en_particles review of Codex 7d9c93b

27 October 2026. Claude (track en_particles_skeptic), issue #13. This file checks the findings in [NOTES.md](NOTES.md)
and `results/*.json`. For each claim I recomputed the key number with a different estimator or formulation, or reran
the check myself. Rules: [BRIEF.md](../../BRIEF.md) and [BRIEF_EARLY.md](../../BRIEF_EARLY.md).

Scope:
- Everything here concerns the reduced X/φ pilot only. Nothing validates or excludes anything. No physical n̄, f_sym or
  threshold is claimed.
- Code and results are in `skeptic/`.
- The Codex tree in `$SP/e7d` was only read. Git was used read-only, and nothing was committed.

## How each number was produced

| Check | What | Independence |
|---|---|---|
| s1 `sk_snapshots.py` | Per-mode counts n₋ = \|ωψ_k + iπ_k\|²/(4ω), built from separate real FFTs, with the FD4 symbol written from the stencil. Also Q₀/Q and the energy outside the mean, on the regenerated τ = 20, 360 and 600–800 snapshots | Own formulation. The fields are en_particles' regenerated history (Codex engine, reproduction-type) |
| s2 `sk_export.py` | Charge and energy rebuilt from SPECTRAL_EXPORT alone, bin-level references, τ = 20 occupation per mode by k | Own code on Codex's published export |
| s3 `sk_windows.py` | Hann, Slepian multitaper (NW 2.5), rectangular and per-mode virial (n + ½-type, √(⟨\|π\|²⟩⟨\|ψ\|²⟩)) counts. Per-shell branch centroids. Box-scale fractions. Full window and 100-sample sub-windows | Own estimators, run on en_particles' reproduction-type windows |
| s4 `sk_mass_scan.py` | Time-averaged equal-time n₋ for fixed m² from 0 to 60, the instantaneous Hartree mass and the Hann²-weighted virial count | Own code |
| s5 `sk_orbit.py` | τ = 20 homogeneous orbit split by the **exact quasi-periodic Fourier series**: radial period by quadrature, with no time window | Own method |
| s6 `sk_kick.py` | Kicked twin runs from **τ 700** instead of τ 500, with my own driver, kick seed 11 and ε = ±0.05. Analysed with the like-for-like transform and a **matrix-pencil pole fit** | Own driver and analysis. Codex's unmodified engine, compiled into scratch, run for reproduction |
| s7 `sk_synth_response.py` | The response transform applied to free fields of known m², and to an oscillator with Hartree-like mass modulation | Own code |
| s8 `sk_box.py` | Virial counts per shell over τ 547–800 for L12 (seed 1), seed 2 and L16 N64, and for L12 over τ ≈ 1100–1460. My own engine runs from en_particles' s = 14 states | Own driver and analysis. Engine runs are reproduction-type |
| s9 | L16 control regenerated to exactly τ = 800 and compared with Codex's published box16 values. A relabelled-clock restart τ 700 → 800 at h = 0.004 and h = 0.001 | Reproduction |
| s10 | Codex's published FIELD_CURVES.npz over τ 700–800 for all seven runs | Codex data, read only |

## Verdicts on the reviewer's claims

| ID | Reviewer | My recomputation | Refuted? |
|---|---|---|---|
| 1 massless 2.2375 | CONFIRMED | 2.237515 (amplitude form); 2.24277 bin-level | no |
| 2 Hartree 1.1800 | CONFIRMED | 1.179962; 1.17965 bin-level; m² 19.2636 | no |
| 3 export rebuild | CONFIRMED | 3.55e-15 / 2.44e-15. My full-space binning of the regenerated τ 800 state matches the published export to 1.76e-11 | no |
| 4 n± formula | CONFIRMED | Identical numbers from the amplitude form. At τ 20 the occupation per mode is 0.58–0.62 at k 3.5–7.5 and 0.49 at k 7.5–8, but 39–55 at k 1.8–2.6, which was already parametrically amplified. Unseeded modes are 0 | no |
| 5 1.9× split is phase-dependent | CORRECTED | Exactly reproduced: 2.174, 1.306, 2.643, 1.397, 2.161, 2.238; Hartree 1.153, 1.216, 1.081, 1.197, 1.173, 1.180. **The split itself is 1.07–2.44×** (time average 1.778 / 1.152 = 1.54). These are regenerated snapshots at Codex's output times, not published Codex snapshots | no (verdict upheld) |
| 6 Q₀/Q | CONFIRMED | 1.054625. Codex's own curve gives 1.0467–1.1694 over τ 340–800 | no |
| 7 count in the waves 1.145 etc. | NEW | **Number holds:** Hann 1.1451 (exact reproduction), multitaper 1.1233, virial 1.1186 (uniform weights) and 1.1515 (Hann²). Seed 2 gives 1.185 and L16 1.240 (virial, τ 547–800) | no |
| 7b "the measured dispersion picks the Hartree reference, to within 3–6%" | NEW | **Not supported.** A fixed-mass reference, time-averaged, reproduces the measured 1.145 to within 3% for **any m² from 7 to 20**, and to within 6% from 5 to 26. The count is smallest (1.152) at m² ≈ 12, not 20. The statistical spectrum that the count weighs has branch centroids giving m²_eff ≈ 8–13 at every k (Hartree 20.1). Only the response peaks follow Hartree, and they do so where the waves hold little X̄: k ≥ 4 holds 4%, k ≥ 3 about 10%. The agreement with Hartree reflects how flat the count is near the true frequency. Massless is the outlier because ω = k lies far below it | **yes** (the interpretation, not the number) |
| 8 zero mode + conserved total | NEW | τ 20 orbit by exact Fourier series: **1.738732**, which is the reviewer's number. The lines are +3.678 (63%) and −4.037 (37%). The total per sub-window is 1.741–1.753 (Hann) and 1.753–1.766 (virial). At τ 547–800 it is **1.761** for L12, seed 2 and L16 alike (virial). Conservation holds within each estimator, while the level depends on the estimator at about 1% | no |
| 9 zero point | NEW | 5.27e-4 per net charge; τ 20 waves 0.00185 / 0.00241 | no |
| 10 REVIEW_TIMING comparison | CORRECTED | Four minor fixes. (a) The (3,.5) equilibrated end point is 0.096–0.107; "0.10–0.12" includes (5,1). (b) The 1.75 total is already in the τ 20 orbit (1.739), before any pairs exist, so its lying inside "pairs 1.36–2.62" does not support the pairs mechanism. (c) "63% vs 88%" compares energy outside the mean, a phase sample spanning 47–73% over τ 700–800, with energy above the circular orbit, which is fixed at 88%. These are different quantities. (d) The count is flat to ±0.4% over Δs ≈ 34. That only bounds the relaxation time to ≳ 5×10³ conformal units, against ~10¹¹ available, so it does not support a statement about the binding era | no |
| 11 Δs = 1.83 | NEW | Δs = 1.8335 is confirmed. But "less than one period of the low-k modes" is wrong: 2π/ω₀ = 1.80, so τ 700–800 spans about 1.0 period, and about 2 periods of the 2ω₀ modulation that drives the reference swings. The restart at h = 0.004 reproduces τ 800 to 3–9e-7 in the charged fields but only **2e-5 in φ**. At h = 0.001 the differences are 4e-10 and 2.5e-8. The window step is about 4× Codex's effective late step (Δs ≈ 0.0009). Harmless for counts | no (wording) |
| 12 lab frame and Z± | NEW | At τ 700, with my own kick, Z₊ = 0.490–0.509 and Z₋ = 0.487–0.507. Low-k ρ peaks are symmetric about 0 to within 0.13. Note that Z₊ + Z₋ = 1 is the exact f-sum identity (G(0) = 0, G′(0⁺) = 1 for H = π²/2 + V(ψ)), so it checks only the numerics. Z± ≈ ½ also does not by itself show narrow poles: at k ≤ 1.6 a 4–6-pole pencil fit leaves 11–21% residual | no |
| 13 φ quanta | NEW | 0.00303 (virial) | no |

The response dispersion at k ≥ 3 also holds up:
- At τ 700, ρ peaks give m² = 17.7–18.8 (FD4 k), and pencil poles give 19.0–19.9 (one shell at 21.6).
- Synthetic test: the transform recovers free m² to ±0.16.
- A Hartree-like mass, 20.1 ± 8.3 modulated at 2ω₀, gives peak m² of 18.3 (k 1.05) to 19.5 (k 6.3).
- The reviewer's 18.9–19.4 against 20.1 is therefore consistent with an orbit-modulated Hartree mass.

## The box and the extension

At τ 547–800, 70% of the wave X̄ sits at k ≤ 1.6 in every run: 0.696 for L12, 0.707 for seed 2 and 0.703 for L16. The
peak is at k ≈ 0.8–1.2 in both boxes. So up to τ 800, "an inverse cascade toward k ≈ 1" is not a box artefact.

The split between waves and the zero mode does move with the box:

| Run | Waves | Zero mode | Total |
|---|---:|---:|---:|
| L12 | 1.170 | 0.591 | 1.761 |
| L16 | 1.240 | 0.521 | 1.761 |

In the τ ≈ 1100–1460 extension, the L12 box fundamental (18 modes, k = 0.52) holds **34%** of the wave X̄. The
sub-windows show it rising from 19% to 34%. The extension values (1.27–1.31) are therefore box-limited.

The L16 control is a faithful reproduction. Regenerated to τ 800 it gives 2.412006 / 1.265917 / Q₀/Q 1.146853, against
Codex's 2.412004 / 1.265917 / 1.146853.

## New issues with the Codex packet (missed by the reviewer)

1. **The τ = 800 headline fractions are phase samples (from Codex's own FIELD_CURVES).**
   - "61–65% of active energy outside the spatial mean": N48 oscillates **46.9–72.7%** over τ 700–800, with mean
     60.8%. Across all seven runs the range is 45.8–74.5%.
   - Neutral plus interaction "0.28–0.36%": N48 spans 0.277–0.344% (mean 0.311%), and all runs 0.246–0.405%.
   - Q₀/Q spans 1.046–1.211 across all runs.
   - Comparisons between controls at τ 800 are valid only because they share the phase, as in 62.722% versus
     62.734%.
   - Quote τ 700–800 means and envelopes instead.
2. **SPECTRAL_EXPORT is a single-phase snapshot.**
   - Every reference projection, and any kinetic input built from it, inherits a phase swing: massless 1.31–2.64 and
     Hartree 1.08–1.22 over τ 700–800.
   - Export a period average instead: Δs ≥ 1.8, which is Δτ ≳ 100 near τ 800.
3. **The isotropic Hartree mass swings with the orbit.**
   - It includes the coherent |ψ̄|², so m²_H is 15.0–25.4 across the τ 700–800 snapshots and 11.8–28.5 over τ 340–798.
   - MATHEMATICS.md presents m²_H as a single number without saying that it is evaluated at an arbitrary orbital
     phase.

## Minor points on the reviewer's write-up

- "ρ agrees to 7e-16" is weak evidence, because ρ is conserved. The 1.8e-11 export rebuild is the real state check,
  and I confirm it.
- The L16 control window used Δs = 0.1, not 0.05. This is unreported but harmless, since Nyquist 31 is well above
  ω_max ≈ 16.
- Wide-taper estimators fail on 5-unit sub-windows. A multitaper with NW 2.5 smooths by ±3 in ω and inflates the
  total to 1.86. A rectangular window leaks: the zero-mode n₋ is 0.64–0.90 in the late sub-windows, against 0.59–0.69 with Hann, and the total reaches 2.08. The reviewer's Hann choice and the
  virial estimator agree with each other. The others should not be used here.

## Limits

- These are classical single realisations.
- The response was checked at two start times: τ 500 (reviewer) and τ 700 (mine).
- Any particle and antiparticle split, whether of the zero mode or of the k ≲ 1.6 waves, is a basis choice.
- None of these numbers is a physical abundance.

## Reproduce

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; R=$SP/review/en_particles
S=$SP/review/en_particles_skeptic; K=tracks/en_particles/skeptic; export OMP_NUM_THREADS=2
g++ -std=c++17 -O3 -ffp-contract=off -fopenmp $SP/e7d/research/reality_framework/grut_early_nonlinear_20261027/code/early3d.cpp -o $S/bin/early3d
python3 -I $K/code/sk_snapshots.py $K/results/sk_snapshots.json $R/ework/N48/state_s00000{00,17,30,31,32,33,34,35,36,37,38,39}.bin
python3 -I $K/code/sk_windows.py $R/W_late/window.npz $K/results/sk_W_late_full.json
for w in W_first W_early W_late W_ext; do python3 -I $K/code/sk_windows.py $R/$w/window.npz $K/results/sk_sub100_$w.json 100 20; done
python3 -I $K/code/sk_mass_scan.py $R/W_late/window.npz $K/results/sk_mass_scan_late.json
python3 -I $K/code/sk_orbit.py $R/ework/N48/state_s0000000.bin $K/results/sk_orbit_tau20.json
python3 -I $K/code/sk_kick.py run $S/bin/early3d $R/ework/N48/state_s0000034.bin $S/kick_p 0.05     # and -0.05 -> kick_m
python3 -I $K/code/sk_kick.py analyse $S/kick_p $S/kick_m $K/results/sk_response_tau700.json
python3 -I $K/code/sk_box.py window $R/W_late/window.npz 99 $K/results/sk_box_L12_late.json            # W_ext 100 -> _ext
python3 -I $K/code/sk_box.py run $S/bin/early3d $R/ctrl/box16_s14.bin $S/box16 18.98331846877269 23.88331846877268 0.05 $K/results/sk_box_L16_late.json
python3 -I $K/code/sk_synth_response.py $K/results/sk_synth_response.json
python3 -I $K/code/sk_summarize.py $K/results $SP/e7d/research/reality_framework/grut_early_nonlinear_20261027/FIELD_CURVES.npz
```

Credits: Khlebnikov–Tkachev and Micha–Tkachev (classical lattice method and turbulent particle numbers); Epelbaum–
Gelis–Wu arXiv:1402.0115 (classical UV limits); Aarts–Berges (classical spectral functions) [R]; Hua–Sarkar (matrix
pencil) [R]; Thomson (multitaper) [R].
