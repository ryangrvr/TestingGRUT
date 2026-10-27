# pt_relax: does the ~1.75 classical X̄ budget survive to binding? (Codex 5f5513b, issue #13)

Claude review track, 27 October 2026. Brief: [BRIEF_TRANSPORT.md](../../BRIEF_TRANSPORT.md), key `pt_relax`. This is
same-session AI work on one supplied, reduced scalar model. It is not outside replication, and nothing here is
observational evidence, a validation, an exclusion or a GRUT prediction. Nothing was committed, and Codex's tree was
not edited.

## What is independent and what is reproduction

- **Reproduction tool.** The long evolutions use Codex's reviewed FD4/Yoshida engine (`fieldstep.cpp`, via
  `pt_common.library`). It runs from a private `git archive 5f5513b` copy in `$SP/review/pt_relax/rf/`. The shared
  `$SP/f55` tree had been rerun in place by pt_rerun, so it was not used.
  - The initial fields come from Codex's `prepare.initialize` in that copy.
  - nQ = 8112.665383318339 and ρ = 105481.24746711778 equal the published RUN_IR12_s1_w1.json.
  - The `initial.npy` SHA-256 (`69f758e2…`) differs from the published `170220d2…` but equals pt_rerun's
    regeneration. This is a last-bit platform difference (REPRODUCE.md anticipates it) and was not pursued.
- **Own code** (`code/`):
  - `long_run.py`: driver, restarts, Hann estimator, shell spectra and occupations;
  - `xcheck.py`: integrator check, using our independent BM6 integrator from `tracks/en_indep/code/indep3d.py`;
  - `analyze_long.py`, `phi_number.py`: analysis;
  - `scales.py`: onset, equilibrium, injection and rates. It imports nothing from Codex.

## 1. Conformal argument and the time available [D]

**Autonomy.**
- The reduced Hamiltonian in conformal variables has no explicit s dependence (`fieldstep.cpp`): there are no masses,
  the coupling is conformal and a″ = 0 in radiation.
- So any rate in conformal units is constant in s. A process with conformal rate Γ_c gives Γ_c·s_onset e-folds before
  the onset, and Γ_phys/H = Γ_c·a/C.
- Autonomy is broken only by:
  - the tachyonic φ term −κa²χ_φ²/2, with κ = λ_φv²/H_i² = 1.4407×10⁻²⁴, which is negligible until a ~ 10¹²;
  - logarithmic coupling running;
  - g* steps in the visible sector;
  - the timed injections (section 4b).

**Onset.** The onset is where κa² = g²⟨r²⟩, with s = (a − a_m)/C:

| ⟨r²⟩ source | a_onset | s_onset | Γ_crit = 1/s_onset |
|---|---:|---:|---:|
| brief, \|χ_X\| = 60 / 110 | 5.00×10¹² / 9.17×10¹² | 2.52×10¹² / 4.62×10¹² | 3.97×10⁻¹³ / 2.17×10⁻¹³ |
| measured, classical medium s 24–34 (⟨r²⟩ = 10104) | 8.38×10¹² | 4.22×10¹² | 2.37×10⁻¹³ |
| measured, N32 s 994–1004 (5844) | 6.37×10¹² | 3.21×10¹² | 3.12×10⁻¹³ |
| chemical equilibrium (r_c² + ⟨\|δψ\|²⟩ = 3927) | 5.22×10¹² | 2.63×10¹² | 3.80×10⁻¹³ |
| kin_source a_ref = 7.047×10¹² (T_vis 45.7 GeV) | 7.05×10¹² | 3.55×10¹² | 2.82×10⁻¹³ |

- **The ~1.75 budget survives only if Γ_chem < (2.2–4.0)×10⁻¹³ per unit conformal time**, i.e. (1.2–2.3)×10⁻¹⁴ T_eq.
- s_onset is 0.7–1.4×10¹¹ times Codex's s = 34.
- The conformal m_X at the onset is 3.5–6.4 [D].

## 2. Number-changing rates against Γ_crit [E]

T_eq = 17.56, r_c = 62.3 and Γ_crit = 3.8×10⁻¹³ (results/scales.json `item2_rates`).

The 2↔4 rates use Γ = c·V⁴·T with c = (10–10³)/(4π)⁵, where V is the Feynman vertex (4λ_X, 2g², 6λ_φ). These are
order-of-magnitude estimates, not computed collision integrals.

| Channel | Γ (conformal) | e-folds to onset |
|---|---:|---:|
| 2↔4, X only, thermal f ~ 1 | 1.4×10⁻¹³ – 1.4×10⁻¹¹ | 0.4 – 38 (marginal) |
| 2↔4 through two g² vertices | 9×10⁻¹¹ – 9×10⁻⁹ | 2×10² – 2×10⁴ |
| 2↔4 mixed g²·λ_φ (XX̄ → 4φ) | 7×10⁻⁹ – 7×10⁻⁷ | 2×10⁴ – 2×10⁶ |
| 2↔4 φ self (λ_φ = 0.03), once φ is populated | 6×10⁻⁷ – 6×10⁻⁵ | 10⁶ – 10⁸ |
| 2↔2 XX̄ ↔ φφ (number-conserving, but moves X̄ into φ) | 1.6×10⁻⁷ – 1.4×10⁻⁶ | 4×10⁵ – 4×10⁶ |
| condensate-assisted 1↔2 (vertex 2λ_X r_c, 2g² r_c) | 6×10⁻⁷ – 6×10⁻⁴ | 10⁶ – 10⁹ |
| **measured classical lattice drift (section 3)** | **2.4×10⁻⁴ – 2.2×10⁻³** | **10⁹ – 10¹⁰ × Γ_crit** |

- With a single-sign (circular) condensate, the condensate-assisted 1↔2 processes conserve total number including
  condensate quanta [A]. They act on the present eccentric condensate, not on an equilibrated one.
- Only the pure-X 2↔4 channel is marginal. Every channel through φ is fast by more than 10² e-folds.
- **φ does get populated** as the X cascade crosses the threshold set by m_φ ≈ g√⟨r²⟩ ≈ 10. The N32 continuation
  shows this (section 3).

**Literature** [Q, arXiv abstracts read this session]:
- Berges–Sexty 1201.0687: in relativistic theories, number-changing processes counteract condensation and drive the
  later decay of the condensate; the condensate forms through an inverse particle cascade in a dual cascade.
- Micha–Tkachev hep-ph/0403101: parametric resonance, then driven and then free turbulence, with self-similar
  distributions.
- Semikoz–Tkachev hep-ph/9507306: kinetic Bose condensation in λφ⁴ with 2↔2 Boltzmann.
- Kurkela–Moore 1207.1663: classical UV cascade, p_max ∝ t^{1/7}. Kurkela–Moore 1107.5050: parametric
  thermalization, t ~ α⁻²/T, for gauge theory.
- Epelbaum–Gelis–Wu 1402.0115: the classical-statistical approximation is not renormalizable.
- Epelbaum–Gelis–Jeon–Moore–Wu 1506.05580: quantum Boltzmann evolution with a condensate departs from the classical
  attractor.
- Jeon hep-ph/9409250: effective Boltzmann transport in scalar theory with thermal masses. The leading role of 2↔4
  there in bulk viscosity is [R].

## 3. Direct classical test [C]

**Runs.** Codex's engine was run with my estimator (Hann, windows of 10 in s, Δs = 0.04):
- IR12_s1_w1 (N32, L = 12) to s = 1004 (τ = 1.01×10⁶);
- IR12_fine_s1_w1 (N48, L = 12) as a cutoff control to s = 254 (τ = 6.6×10⁴).
- Compute time was about 61 min on 2 threads.

**Reproduction of Codex's windows** (coherent / waves k ≤ 8 / total, ours versus Codex):

| Window | Codex | Ours |
|---|---|---|
| IR12 14–24 | 0.670554 / 1.079999 / 1.750553 | 0.6705539 / 1.0799988 / 1.7505527 |
| IR12 24–34 | 1.741447 (total) | 1.7414468 |
| IR12_fine 14–24 | 1.751464 (total) | 1.7514644 |
| IR12_fine 24–34 | 1.744233 (total) | 1.7442332 |

**CONFIRMED.**

**The budget is a plateau, not a conserved number.** All-mode total A₋ per net charge (k ≤ 8 band plus all
k > 8 lattice modes):

| s window | N32 coherent | N32 waves k ≤ 8 | N32 k > 8 | **N32 total** | N48 total |
|---|---:|---:|---:|---:|---:|
| 24–34 | 0.411 | 1.330 | 0.010 | **1.751** | 1.750 |
| 94–104 | 0.606 | 1.064 | 0.047 | 1.717 | 1.726 |
| 194–204 | 0.559 | 1.025 | 0.067 | 1.651 | 1.703 |
| 244–254 | 0.511 | 0.982 | 0.073 | 1.566 | **1.690** |
| 494–504 | 0.357 | 0.734 | 0.112 | 1.202 | |
| 994–1004 | 0.152 | 0.518 | 0.162 | **0.832** | |

**Drift rates and fit.**
- N32 slopes per unit s: −4.9×10⁻⁴ (s 29–99), −6.5×10⁻⁴, −2.2×10⁻³ (199–299), −1.2×10⁻³, −9.3×10⁻⁴, −5.9×10⁻⁴
  (699–1000).
- An exponential fit for s ≥ 299 gives A₋ → **0.510** with **Γ = (1.540 ± 0.020)×10⁻³** (rms residual 0.0026).
- That asymptote lies inside the harmonic Rayleigh–Jeans lattice estimate of 0.41–0.61 [E].
- The fitted rate is a classical, cutoff-dependent thermalization rate, not a physical one.

**Cutoff dependence.**
- N48 drifts more slowly: −3.4×10⁻⁴ (s 29–99), then −2.4×10⁻⁴ (99–254); −3.4% by s = 254.
- N32 is at 1.566 at the same time.

**Where the action goes.**
- On N32 the energy above k = 8 rises from 1.2% to 31%.
- φ rises from 0.28% to 26% of the energy. Its frequency number per net charge rises from 0.0006 to 0.21.
- At fixed charge X + X̄ lose 1.84 per net charge, so the loss is mostly genuine number change, not conversion into φ.

**Checks.**
- Energy drift ≤ 1.4×10⁻⁷, charge drift ≤ 4×10⁻¹⁵, spectral charge 1 − 5×10⁻⁹.
- Integrator check over s 894–904 (`results/xcheck_integrators.json`): Codex's engine at h = 0.004 and h = 0.002, and
  our own BM6 (h = 0.004), give the same windowed A₋ to 1.9×10⁻⁸. The pointwise fields differ by at most 4.3×10⁻⁴ of
  the RMS.

**Where classical evolution fails** [C/A].
- N32 occupations stay ≥ 14 (s 200) and ≥ 97 (s 1000) in every shell. But the FD4 spectrum ends at |k| ≈ 10.7, below
  the quantum thermal scale (T_eq ≈ 18; the Bose spectrum peaks near 3T ≈ 50). The lattice therefore cannot represent
  the quantum end point, and its own end point (≈ 0.5) is a lattice artifact.
- On N48 the cascade front reaches shells with occupation 0.10–0.18 (X) and 0.12–0.19 (φ) by s ≈ 200–254, where the
  classical-statistical approximation is unreliable (1402.0115). Those shells hold only about 1.4% of the energy.

**Classical cascade time to the thermal scale** [E, using the literature scalings, not our fit]. With
p_max ∝ s^{1/7}–s^{1/5}, going from k ≈ 8 at s ≈ 200 to 3T ≈ 53 takes s ≈ 2.5×10⁶–1×10⁸, far below
s_onset ≈ 3×10¹². Our own k_rms exponent (0.3–0.4) is contaminated by the k ≤ 8 vacuum seed and is not used.

## 4. Equilibrium end point [E] (results/scales.json)

**Pilot state.** Fixed n_Q = 8112.7 and ρ = 1.0548×10⁵ (conformal), with all charge in a circular condensate.

| Model | T_conf | X̄ per net charge |
|---|---:|---:|
| free massless BE, X only (2 dof) | 19.40 | 0.110 |
| free massless BE, X + φ (3 dof) | 17.53 | 0.081 |
| Hartree, φ in equilibrium (μ = 2.04, r_c² = 3881, M_X = 2.84, M_φ = 6.40, φ share 29%, condensate charge 97.6%) | 17.56 | **0.069** |
| Hartree, φ absent | 19.40 | 0.095 |
| Hartree, φ decoupled holding 0.3% / 5% / 10% of the energy | | 0.095 / 0.091 / 0.086 |

- Range: **0.069–0.110**. This corrects our earlier 0.096–0.116 downward once φ shares the energy.
- Interactions shift the result by about 15%. The condensate is circular, so it contributes no X̄.

**At the onset** (m_X switched on at fixed conformal T and n_Q).
- With equilibrium ⟨r²⟩: a = 5.22×10¹², T_hidden = 0.448 GeV, T/m_X = 4.85.
- At a_ref: T_hidden = 0.332–0.367 GeV, T/m_X = 3.59–3.97.
- Thermal X̄ per net charge at μ/m = 0, 0.861 or 0.99: 0.056–0.107.
- The φ vacuum energy λ_φv⁴/4 is 12% (a = 5.2×10¹²) to about 40% (a_ref) of the hidden energy. The transition
  reheats the hidden sector, and that is not included here.
- After the onset, X̄ follows the massive-phase equilibrium until XX̄ annihilation freezes out. That is the case-L
  kinetic problem.

**4b. The end point depends on injected energy.** Heavy spot release (N_X = 2788.18, E_X = 38386.25) with injected
conformal energy R·E_X, R taken from kin_source; columns give φ in equilibrium / φ absent:

| R | T_conf | T/m_X at a_ref | X̄ per net charge |
|---:|---:|---:|---:|
| 0 | 13.8 / 15.2 | 2.82 / 3.12 | 0.099 / 0.136 |
| 1.144 (mediators) | 16.9 / 18.6 | 3.45 / 3.82 | 0.19 / 0.26 |
| 4.05 (vectors, e = 10⁻⁶) | 21.0 / 23.2 | 4.30 / 4.75 | 0.37 / 0.50 |
| 40.56 (vectors, e = 10⁻⁷ band) | 35.7 / 39.4 | **7.30 / 8.08** | **1.86 / 2.55** |
| 48.50 (broad) | 37.2 / 41.2 | 7.63 / 8.44 | 2.12 / 2.92 |

The R = 40.56 temperatures reproduce kin_thermo's T/m = 7.34 and 8.12 to 0.5% [C, independent route].

## Verdicts

| Item | Codex | Ours | Verdict |
|---|---|---|---|
| Classical budget s 14–24 / 24–34 (IR12, IR12_fine) [C] | 1.750553 / 1.741447; 1.751464 / 1.744233 | 1.7505527 / 1.7414468; 1.7514644 / 1.7442332 (Codex engine, own estimator) | CONFIRMED |
| Same budget beyond s = 34 [C] | not claimed conserved; "not a lifetime until binding" | all modes 1.751 → 0.832 by s = 1004 (N32, τ 1.0×10⁶); 1.750 → 1.690 by s = 254 (N48); classical Γ 2.4×10⁻⁴–2.2×10⁻³ | NEW (Codex's caution CONFIRMED) |
| "Only net charge is protected" [A] | stated | no second conserved number; it decays even classically | CONFIRMED |
| Chemical relaxation is the deciding calculation [A] | stated | Γ_crit = (2.2–4.0)×10⁻¹³ per unit s; φ channels ≥ 10² e-folds; only pure-X 2↔4 is marginal | CONFIRMED, quantified |
| Equilibrium "near 0.1" [E] | not adopted | 0.069–0.110 (pilot energy); 0.10–0.14 at the heavy release with R = 0; 1.9–2.9 with R = 40.6–48.5 | CORRECTED: depends on energy per charge |
| Slow matching ≈ 1.76; IR preparation; matrix response; export | | | NOT CHECKED (other tracks) |
| Our REVIEW_EARLY: "total X̄ number is conserved through τ 1460" | | true only for s ≲ 34–100; halves by s ≈ 10³ on N32 | CORRECTED (own) |

## Verdict on item 5 [A/E]

**The ~1.75 budget is very unlikely to survive to binding in the reduced model.**
- It is a prethermal plateau of the over-occupied classical medium. Even the classical lattice dynamics erodes it,
  at rates 10⁹–10¹⁰ above Γ_crit.
- Several quantum channels through φ give 10²–10⁸ e-folds before the onset.
- The relevant population at binding is the chemical-equilibrium value at the actual energy per charge:
  - 0.07–0.11 for the pilot's energy alone, which becomes 0.06–0.11 at T/m_X ≈ 3.6–4.9 once m_X switches on;
  - about 1.9–2.9 if the e = 10⁻⁷ vector heat (R ≈ 40–48) has thermalized.

The early classical number is not a lower bound on that population. Nor is the "~0.1" value a universal answer.

**What Codex's chemical-relaxation calculation should compute:**

1. **The right clock.** Compute Γ_chem in conformal units for the actual medium and compare it with 1/s_onset ≈
   (2–4)×10⁻¹³. The comparison is not with s = 34 or with H at τ ~ 10³.
2. **Quantum kinetics from the measured source.** Run quantum kinetic evolution (a Boltzmann equation with condensate,
   Hartree/Bogoliubov dispersion and 1 + f factors, or 2PI) of X, X̄, φ and the eccentric condensate. Start from the
   s = 34 source.
   - Include: 2↔2 for all three couplings; XX̄ ↔ φφ; condensate-assisted 1↔2; 2↔4 for λ_X, g² and λ_φ.
   - Output X̄/Q(s) and s_chem.
   - Use classical-statistical runs only while the front has f ≫ 1. On a lattice they cannot reach the
     T ≈ 18 (3T ≈ 50) end point, and their drift is cutoff-dependent (N32 versus N48).
3. **The end point at fixed (n_Q, ρ),** including φ's share: 0.069–0.110 here.
4. **The timed injections.** Include the Y and vector injections (R = 0.8–48) and their thermalization. These
   dominate the end point: 0.1 at R = 0, 2–3 at R = 40–48.
5. **The onset and after.** At the onset (a ≈ 5–9×10¹², T_h/m_X ≈ 3–8), include the φ-transition heat (12–40% of the
   hidden energy). Then hand the thermal X̄ to the case-L annihilation network.

Since the conformal problem is scale-free, a single quantum-kinetic run per energy-per-charge settles X̄/Q before the
onset [D].

## Limits

- The rates in section 2 are parametric brackets, not collision integrals. The pure-X 2↔4 channel is within the
  bracket of Γ_crit.
- The classical runs are one seed, one box (L = 12), and two cutoffs. N48 stops at s = 254. The classical drift rate
  is not a physical rate.
- The equilibrium uses averaged Hartree masses, not a full Bogoliubov spectrum. Thermal-mass effects on X̄ are about
  15%.
- No heavy, gauge, Y or vector dynamics are evolved. The injected-energy end points assume full thermalization.
- Nothing here sets n̄, f_sym, the threshold or the abundance.

## Reproduce

```bash
source tracks/pt_relax/code/env.sh   # creates no files; needs the private copy $RF (git archive 5f5513b research/reality_framework)
cd tracks/pt_relax/code
python3 -B long_run.py IR12_s1_w1 1004 550        # repeat until s = 1004 (restart files in $RUNDIR)
python3 -B long_run.py IR12_fine_s1_w1 254 560    # repeat until s = 254
cp $RUNDIR/IR12_s1_w1/restart.npz $RUNDIR/IR12_s1_w1_restart_s894.npz   # (taken at s = 894 in this run)
INDEP_WORKERS=2 python3 -B xcheck.py
python3 -B analyze_long.py; python3 -B phi_number.py
python3 -B scales.py $SP/review/pt_relax/r2.json  # r2.json: measured <r^2> {10104.0, 5844.2}
```

Results:
- `results/pt_relax_results.json`: summary;
- `long_run.json`: windows, equal-time series, occupations and spectra;
- `scales.json`, `xcheck_integrators.json`, `phi_number.json`.

Scratch, under `$SP/review/pt_relax/`: run parts, restarts and abstracts.
