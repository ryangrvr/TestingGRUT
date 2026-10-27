# en_neutral_skeptic: adversarial check of track en_neutral (Codex 7d9c93b, issue #13)

27 October 2026. Claude, track en_neutral_skeptic. Target: the findings in [NOTES.md](NOTES.md) about Codex's same-seed
neutral comparison in `grut_early_nonlinear_20261027/` (git 7d9c93b).

Rules: [BRIEF.md](../../BRIEF.md), [BRIEF_EARLY.md](../../BRIEF_EARLY.md).
- Everything here is conditional on Codex's reduced, conformal, classical X/φ pilot.
- Nothing validates or excludes anything. No commits were made.
- Code is in `skeptic/code/` and results are in `skeptic/results/`. Scratch files are in `$SP/review/en_neutral_skeptic/`.

## Method and labels

| Code | What it does | Label |
|---|---|---|
| `sk_lib.py` | Own library: FD4 symbol from the (−1, 16, −30, 16, −1)/12 stencil; seeds by full complex FFT of the τ-20 field; DOP853 orbits; cubic-spline measured pumps; RK4 transfer on a grid refined ×2 | independent of both Codex's and en_neutral's code |
| `factorial.py` | Same-seed linear δφ on eight pumps (2×2 depletion/detuning factorial, plus re-runs of en_neutral's L2/L3/L4) | independent code on reproduction-run data |
| `controls_extra.py`, `ensemble_grid.py`, `drift_check.py` | Floquet monodromy, exact ensemble mean and variance, grid dependence, concatenation check, fixed depleted orbits | independent |
| `nl_indep.py` | Own nonlinear lattice: spectral FD4 symbol and PEFRL 4th-order integrator. Codex and en_neutral both use the Yoshida triple jump. The N48 τ-20 field is resampled exactly onto N32 | independent |
| `snap_check.py`, `uv_occupation.py`, `xspec_check.py` | Own reading of en_neutral's reproduction snapshots and trace (Codex's unmodified engine, N48) | independent analysis of reproduction data |

The shared input is en_neutral's reproduction run of Codex's engine (`$SP/review/en_neutral/runs/N48`). Its final value
equals the published N48 curve to 3e-15, and our own real-space ρ₀ matches it to 4e-15. Published files were
re-extracted with `git show 7d9c93b:` and compared by md5.

## Verdicts on en_neutral's findings

| id | en_neutral's claim | Skeptic verdict | Ours | en_neutral |
|---|---|---|---|---|
| C3a [C] | φ + interaction 0.28–0.36% at τ 800 | stands | Own nonlinear N32 (PEFRL) **0.27733468%**; Codex 0.27733222% (8.9e-6). Published point values 0.2773–0.3648%; last-period means 0.2785–0.3669% | 0.2773–0.3648%; N32 0.27733026% |
| C3b [C] | φ Hamiltonian alone 0.19–0.27% | stands | Published 0.1882–0.2672%; own N32 0.18821681% | 0.188–0.267% |
| C3c [C] | Fixed orbit gives 2.09 | stands, plus NEW-1 | Own RK4 and DOP853: **2.0917128**; maximum 2.1094 | 2.0917128 |
| C3d [A] | Same conventions on both curves | stands; minor correction | Code agrees. U_φ max is **2.07e-6 ρ₀**, not ≤ 1e-6, which is still negligible | ≤ 1e-6 |
| N1 [C] | Seed is 25–35% of the final value; ≈ 0.042% zero point | stands; qualified | Period-averaged share 0.25–0.36 and rise 0.199–0.243%. The zero point depends on the convention: 0.024% (massless), 0.042% (instantaneous ω at τ 20), 0.044–0.046% (period mean), 0.063% (apocentre) | 25–35%; 0.042% |
| N2 [A] | Reconciles with the critic's 696.5–697 | stands, and closer | Ensemble at N48: **2.1837 ± 0.0775**, unit at 698.5. With continuum k: **2.264, unit at 696.5**. Critic: 2.22/2.17, 696.5/697.0 | 2.18, 698.5 |
| — | The single realization is 4% below the ensemble because of realization scatter | stands | z = **−1.19σ** (relative σ 3.5%, n_eff ≈ 794 pairs) | "scatter" |
| F1 [C] | μ_φ ≈ 0.206 at k ≈ 7 | stands | DOP853 monodromy: **0.20624 at k 7.008**; half-maximum band 6.74–7.27 | 0.2062 at 7.01 |
| **M_a** [C] | Pump depletion is REFUTED as the cause | **REFUTED** | See §1. Depletion alone gives **0.32% / 0.87%** against 2.09. The 2.71 figure sits **99.8% at k_lat < 2.5** | REFUTED; zero-mode pump → 2.71 |
| **M_b** [C] | Hartree detuning is CONFIRMED as dominant | **"dominant" REFUTED**; the numbers reproduce | Detuning alone gives **0.48%**, depletion alone 0.32–0.87%, both together 0.21%. Measured Hartree 0.19384%. The g²S series reproduces (0.17 … 6.81 … 14.97 … 47.7) | dominant; 0.194% |
| M_c [C] | Band leaves the lattice or the seed band: REFUTED | stands | S in k 1.8–4: 87–99% (τ 60–200). S in k < 1.8: 56% / 73% (τ 400 / 800). ≤ 0.08% above k 8 | same |
| N3 [C] | Band growth to τ ≈ 66, then decline | stands | Band energy 2.15e-4 → 6.28e-4 (×2.92) at τ 66; peak 9.77e-4 at τ 90; 1.62e-4 at τ 796; 64.3% of the energy above k 8 | ×2.9; 1.0e-3; 1.6e-4; 64% |
| N4 [C] | Without X noise, φ saturates at 4.3–5.9% | numbers stand; one support number corrected | Own N32 (PEFRL): 0.427% at τ 200, 3.73% at τ 400, maximum **5.905%** at τ 778, final 4.34%. "S/X2 = 0.59 at τ 500" is an instantaneous phase sample; the orbit-averaged value is **0.23** at τ ≈ 475 and 0.26–0.33 up to τ 730 | 0.43/3.7/5.9/4.3; S/X2 0.59 |
| C3e [A] | "Changes this forecast in the pilot"; not an exclusion | stands | Scope is fine. The mechanism wording should change (§1) | — |
| corr. 1 | The `$SP/e7d` tree is regenerated and differs from git | stale | The four files are now md5-identical to git. en_rerun restored the tree (its NOTES, line 56). No number is affected | — |
| corr. 4 | "Depletion is not the cause; the zero-mode pump alone → 2.71" | **REFUTED** | See §1 | — |

## 1. The depletion/detuning attribution does not survive a factorial test [C]

The test uses en_neutral's own framework: same-seed linear δφ, the same 14768 seeded modes (422 unique lattice k²),
energy as in Codex, and N48. Here ln σ is the largest singular value of the δφ transfer at τ 800; "band" means
6.74 ≤ k_lat ≤ 7.27.

| Pump (own code) | φ + int at τ 800 | Last-period mean | ln σ in band | Largest ln σ (k) |
|---|---:|---:|---:|---:|
| A: fixed orbit (Codex's control) | 2.0917 | 1.052 | 6.94 | 6.94 (7.04) |
| B: fixed orbit + measured g²S̄ (detuning only) | **0.477%** | 0.377% | 3.16 | 4.24 (1.81) |
| C: free orbit following the measured zero mode's per-period mean (E₀, Q₀) (depletion only) | **0.321%** | 0.311% | 3.15 | 3.68 (6.20) |
| Cp: same, with (E₀, Q₀) at pericentre (en_neutral's V1 convention) | **0.869%** | 0.827% | 3.17 | 4.41 (6.25) |
| D: C + g²S̄ (both) | 0.213% | 0.215% | 2.98 | 3.04 (3.80) |
| E: measured Hartree (en_neutral's L3) | 0.194% | 0.196% | 3.04 | 3.04 (7.14) |
| F: measured zero mode \|χ̄\|² only (en_neutral's L2) | 2.717 | 0.900 | 3.64 | **9.96 (1.81)** |
| G: measured zero mode + g²S̄ (≈ L4) | 0.217% | 0.219% | 3.01 | 3.01 (7.24) |
| Nonlinear N48 (reproduction) | 0.310% | 0.310% | — | — |

- **Depletion alone suppresses the resonance as strongly as detuning alone.**
  - C and Cp hold (E₀, Q₀) at the measured, decaying zero-mode values, from 0.998 ρ₀ down to 0.40 (period mean) or 0.53 (pericentre). Each period is an exact free quartic period.
  - They cut the forecast from 2.09 to 0.32% or 0.87%, a factor of 240–650.
  - Band amplification drops from 6.94 to 3.15–3.17, the same as detuning alone (3.16) and the measured Hartree pump (3.04).
  - Check of the construction: concatenating free periods at constant (E, Q) reproduces A to 3e-9 (2.0917128468).
- **Why the decay matters even though μ "only" falls from 0.206 to 0.15.**
  - The band moves with the orbit's energy (k ≈ 7.0 → 6.2 → 5.9), and the rate is lower on the weaker orbit. Over Δs = 23.9 the lower rate alone is a large factor.
  - A free orbit held at the late pericentre (E, Q) for the whole run gives 2.4%, 89× below 2.09. Held at the late period-mean (E, Q), it gives 0.29% (`drift_check.json`).
  - So both the lower rate and the band drift contribute. In the pericentre convention, drift adds a factor of about 2.7 (Cp 0.87% against 2.4%).
- **en_neutral's 2.71 × ρ₀ (L2) is not the matched resonance.**
  - We reproduce 2.717. Of its τ-800 energy, 93.3% sits at k_lat < 2.0 and 6.5% at 2.0–2.5; band amplification is only 3.64.
  - The cause is a low-k instability opened when the *Hartree-dressed* measured orbit is used without its Hartree mass in the φ equation. en_neutral's own `floquet_periods_N48.json` shows this: V2 has μ_max = 0.43–0.76 at k → 0 from τ 84, and the V2 matched band is 0 from τ 287.
  - L2 therefore says nothing about depletion.
- **Statement to use instead [C].**
  - In homogeneous same-seed proxies, X fragmentation ends the δφ resonance through two channels:
    - the depletion of the coherent pump (lower Floquet rate plus band drift);
    - the fluctuation mass g²S.
  - Each is sufficient alone (0.32–0.87% and 0.48%). Together they give 0.21%, close to the measured Hartree 0.19% and the nonlinear 0.31%.
  - Neither is "dominant" and neither is "refuted".
  - The homogeneous proxies cannot settle this further. S sits at k 1.8–4, only 2–4 times longer in wavelength than the k ≈ 7 φ modes, so it acts partly as an inhomogeneous background rather than a uniform mass [A].

## 2. NEW issues in the Codex packet that en_neutral missed

1. **The same-seed control depends on the grid [C].** Codex publishes only the N48 control.
   - Realization: N32 0.995, N48 2.092, N64 2.197, continuum k 2.255.
   - Ensemble: 0.996 (unit at τ 795.0), 2.184 (698.5), 2.220 (697.0), 2.264 (696.5).
   - The cause is FD4 dispersion. An on-axis continuum k = 7 maps to k_lat 5.91 / 6.72 / 6.90. The resonant lattice band is populated by continuum k 6.97–7.99 at N32, 6.83–7.57 at N48 and 6.75–7.26 in the continuum.
   - So N32/48/64 is not a pure resolution sequence for the φ resonance, and N32's own control is half of 2.09. The ≥ 300× gap is unchanged.
2. **The above-vacuum content of 0.28–0.36% is small [C].** This is a Hartree-ω reference projection, like Codex's own, and not a count.
   - At τ 792–800 (N48), the shell-averaged φ occupation above one half quantum carries **0.067–0.072% of ρ₀**, about 21% of the 0.32–0.34% quadratic φ energy.
   - The 64% of φ energy in unseeded k_lat > 8 shells has classical occupation 0.66 (k = 8) falling to 0.06 (k = 15). 28 of 32 shells are below ½. The vacuum half-quanta of those modes would be 0.39–0.41% of ρ₀.
   - Most of the raw number is therefore seeded half-quanta plus sub-half-quantum UV population. That is the regime where classical-statistical evolution is not a quantum population (Epelbaum–Gelis–Wu, arXiv:1402.0115).
   - en_neutral noted the UV share and the seed share, but not that the UV occupation is below the half quantum, nor the above-vacuum remainder.
3. **"X fluctuation variance reaches half at τ 69.2" is a pericentre artefact [C].**
   - The crossing happens when |χ̄|² is at 2.3% of its local maximum. The instantaneous fraction falls back to 0.06 at τ 75 and 0.008 at τ 80.
   - The orbit-averaged (|χ|²-weighted) share is 0.07 (τ 70–84), 0.26 (τ 100–118), and at most 0.47 up to τ 762.
   - Codex says these diagnostics are pericentre-dependent, but the headline reads as a breakup time. The charge-contrast clock is robust: its period mean is 1.14 over τ 57–70.
4. **τ-800 point values are single-phase samples [C].** N48 φ + int swings between 0.277% and 0.344% within τ 700–800. Differences between histories of ≲ 5% are within this swing.
5. **The outcome depends on the unselected X seed amplitude [A, from N4, reproduced independently].** With zero X noise, φ reaches 4.3–5.9%, against 0.28–0.36% with it. Codex's seven histories vary the seed realization, the cutoff and the box, but not the amplitude of the charged seed. Yet Codex states that the cosmic state is not selected.

## Limits

- All mechanism tests are homogeneous same-seed proxies on one N48 history, seed 1. The nonlinear checks are N32 only.
- The depleted orbit takes (E₀, Q₀) of the bare-potential zero mode in two conventions. A Hartree-dressed definition was not run.
- The occupations are reference projections with ω² = k_lat² + g²⟨|χ|²⟩ + 3λ_φ⟨φ²⟩. They are not populations.
- Classical, reduced, conformal model only. Nothing physical follows.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; K=$SP/review/en_neutral_skeptic; RUN=$SP/review/en_neutral/runs/N48
C=research/reality_framework/grut_formation_review_20261025/tracks/en_neutral/skeptic/code; export OMP_NUM_THREADS=2
mkdir -p $K/pub; git show 7d9c93b:research/reality_framework/grut_early_nonlinear_20261027/MATCH_L12.npz > $K/pub/MATCH_L12.npz
python3 -I $C/factorial.py $RUN $K/factorial_N48.json                    # 40 s
python3 -I $C/controls_extra.py $RUN $K/pub/MATCH_L12.npz $K/controls_extra.json   # 25 s
python3 -I $C/ensemble_grid.py $RUN $K/pub/MATCH_L12.npz $K/ensemble_grid.json
python3 -I $C/drift_check.py $RUN $K/drift_check.json
python3 -I $C/snap_check.py $RUN $K/snap_check.json; python3 -I $C/uv_occupation.py $RUN $K/uv_occupation.json
python3 -I $C/xspec_check.py $RUN $K/xspec_check.json
python3 -I $C/nl_indep.py $RUN $K/nl_N32_full.json 32 0.004 full 800      # 2 min
python3 -I $C/nl_indep.py $RUN $K/nl_N32_noXnoise.json 32 0.004 noXnoise 800   # 2 min
```

`$RUN` is en_neutral's reproduction of Codex's engine (see NOTES.md for how to regenerate it, about 3.5 min).
