# pt_response: the matrix response and what it means for the kinetic description (Codex 5f5513b, issue #13)

27 October 2026. Claude, track `pt_response`. Briefs: [BRIEF_TRANSPORT.md](../../BRIEF_TRANSPORT.md) and
[BRIEF.md](../../BRIEF.md).

This is same-session AI work on one supplied, reduced, classical X/φ scalar model. It is not outside replication.
Nothing here is a validation, an exclusion, a GRUT prediction or new theory.

Claim labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued.

## Which numbers are which

**Reproduction (Codex code, copied tree at `$SP/review/pt_response/rf`).**
- `code/regen_base48_states.py` regenerates the base48 field (N48, L12, seed 1) at s = 20, 20.24, 20.44 and 20.68. It
  uses Codex's initialiser and engine. The 27-column trace equals the published `FIELD_TRACES.npz` bitwise at every
  checked sample ([results/regen_base48_check.json](results/regen_base48_check.json)).
- Codex's `response.run` was then run for all three published cases. All arrays of `RESPONSE_s20_seed7`,
  `RESPONSE_s20.44_seed7` and `RESPONSE_s20_seed11` agree with the published files to 0.0 relative difference
  ([results/reproduction_response.json](results/reproduction_response.json)).

**Mixed (Codex engine, our diagnostics or probes).**
- `code/record_mean_hessian.py` advances the background with Codex's step and records our spatial means every
  Δs = 0.01. It agrees with the published trace to about 1e-9.
- `code/probe_study.py` runs Codex's exact tangent of the discrete map. Our probes, projections and estimators are
  used.
- Our own numpy integrator (`code/own_tangent_check.py`, drift-first Yoshida) agrees with that tangent to
  4.0e-8, 5.1e-8 and 1.0e-7 at delays 0.4, 1 and 2.

**Independent (our code only).**
- `reanalyse_published.py`, `hfb_driven.py`, `bogoliubov_circular.py`, `floquet.py`, `fit_poles.py`,
  `analyse_probes.py` and `coherent_counting.py`.
- Two cross-checks of our own code: the vectorised RK4 used in `hfb_driven.py` agrees with DOP853 to ≤ 5e-7, and the
  Floquet machinery reproduces the analytic circular-condensate value of A₋/W (0.9977 against 0.998).

## 1. Codex's response numbers [C]

These come from our own analysis of the published npz files, and match the bitwise reproduction.
- **Ratio.** The conjugate/normal RMS over delays 0–6, in shells 1–3 (k = 0.52, 1.05, 1.57), for the three probes,
  spans 0.4421–1.7043. This is identical to Codex.
- **Canonical sum.** The maximum of |Z₊ + Z₋ − 1| is 0.003966, against Codex's 0.4%.
  - Our independent transform reproduces Codex's ρ(ω) to 4e-16.
  - The sum is an identity [D]: Z₊ + Z₋ = ∫ωρ dω/2π = dG_N/dt at 0, which equals 1 for any background (equal-time
    commutator, f-sum).
  - It therefore checks the finite transform, not the quasiparticle content.
- **Free control.**
  - Our own Yoshida integration of the single FD4 mode (ω = 1.2541202921170822) gives a maximum error of 2.161928e-10,
    against Codex's 2.161931e-10.
  - There is no cross or anomalous response, by symmetry.

## 2. Definitions [D]

**The two channels.** Each tangent kicks one real momentum, π₁ or π₂, with the same real spatial pattern. The
response R_ab is the response of χ_a to a kick in π_b. For a complex kick p = δπ₁ + iδπ₂ of dψ/ds, with
ψ = χ₁ + iχ₂:

  δψ = G_N p + G_A p*,  G_N = ½(R₀₀ + R₁₁ + i(R₁₀ − R₀₁)),  G_A = ½(R₀₀ − R₁₁ + i(R₁₀ + R₀₁)).

So the normal channel is δψ responding to a δψ′ kick, and the conjugate channel is δψ responding to a δψ′* kick. The
identity holds to 1.1e-15, including for complex shell traces.

**The probe estimator.** It is a stochastic shell trace weighted by |src|²: Σ_k q_k src_k*/Σ|src_k|².
- Its expectation is the shell-averaged diagonal response.
- Off-diagonal (scattering) terms enter as noise. The initial slopes of 0.986–0.997 are consistent with
  ω² = k² + Hartree m².
- **Noise.** For the same field and a different seed (7 against 11, s0 = 20):
  - the normal channel differs by 61%, 31% and 24% in shells 1–3;
  - the conjugate channel differs by 67%, 39% and 28%;
  - in shells 9–15 the normal channel differs by 1.5–4%.

## 3. Low-noise retarded traces (new measurement) [C]

We used two kinds of probe:
- exact per-mode probes for shell 1 at s0 = 20: cos and sin of each of the 9 ± pairs, two components, 36 tangents;
- Hutchinson probes with unit modulus and random phase over shells 1–4: M = 8, giving 16 tangents per phase.

Results are in [results/probes_lownoise.json](results/probes_lownoise.json).

| s0 | Shell 1 | Shell 2 | Shell 3 | Shell 4 | Codex (one probe each) |
|---|---|---|---|---|---|
| 20 | **1.309** exact; 1.272 ± 0.136 Hutchinson | 0.884 ± 0.028 | 0.747 ± 0.018 | 0.619 ± 0.007 | 1.028 / 1.704; 0.932 / 1.028; 0.772 / 0.655; 0.690 / 0.650 |
| 20.44 | 0.501 ± 0.023 | 0.456 ± 0.013 | 0.588 ± 0.009 | 0.484 ± 0.010 | 0.633; 0.442; 0.607; 0.554 |

- **Corrected IR range.** It is 0.46–1.31, against Codex's 0.44–1.70.
  - Codex's 1.704 and 0.633 are single-probe excursions.
  - A white-noise estimator using only shell-1 sources, applied to the exact shell-1 matrix, scatters between 1.08
    and 1.61 (5–95%).
  - The orbital-phase dependence is real: ×2.6 in shell 1 and ×1.9 in shell 2. s0 = 20 is near apocentre; s0 = 20.44
    is just before pericentre.
- **Scattering.** The exact shell-1 data give the full G(k′, k). The amplitude scattered out of the kicked mode, into
  shells 1–4, divided by the diagonal amplitude, is:
  - 0.58 for t < 1;
  - 1.07 for t = 1–2;
  - 1.30–1.41 for t = 2–6.

  The Hutchinson incoherent remainder in shells 1–4 also exceeds 1 by t ≈ 1.5–2.5. For comparison, the IR oscillation
  period is about 1.4–1.7.

## 4. Bogoliubov comparison [D/C]

### (a) Circular condensate, analytic [D]

Set ψ_c = R e^{iω₀s}, η = e^{−iω₀s}δψ = a + ib (a radial, b tangential), and take a homogeneous HFB Hessian with
rotating-frame components H_rr and H_tt. Then:

  a″ − 2ω₀b′ + (k² + H_rr − ω₀²)a = 0,  b″ + 2ω₀a′ + (k² + H_tt − ω₀²)b = 0.

**Consistency.** A consistent circular condensate has H_tt = ω₀², which is the Goldstone condition. It also has
H_rr − H_tt = 2g, where g = λR² plus the pair term. The dispersion is

  Ω⁴ − 2Ω²(k² + g + 2ω₀²) + k²(k² + 2g) = 0.

- With a pure condensate (ω₀² = g): Ω±² = k² + 3g ± √(9g² + 4gk²).
- The phonon speed is c² = g/(g + 2ω₀²), which is 1/3 for a pure quartic.
- The radial gap is Ω₊(0)² = 2g + 4ω₀².

**Modes.** Write (a, b) ∝ (1, iβ)e^{−iΩs}, with β = 2ω₀Ω/(k² + H_tt − ω₀² − Ω²). The lab lines are
- u: ω₀ − Ω, with weight (1 − β)²;
- v: ω₀ + Ω, with weight (1 + β)²;

so |v/u| = |1 + β|/|1 − β|. The quasiparticle action is J = E_rot/Ω.

**Free limit.** Here Ω₊ = E + ω₀ is the antiparticle (β = −1) and Ω₋ = E − ω₀ is the particle (β = +1).

**Parameters at s ≈ 20.** These come from the recorded history ([results/bogoliubov_circular.json](results/bogoliubov_circular.json)):
- the measured mean angular velocity is 3.434;
- λ⟨|ψ_c|²⟩ = 5.27;
- the orbit-averaged rotating-frame mean Hessian is diag(24.21, 16.28).

Three models:
- C1 (HFB-consistent): g = 3.97, ω₀ = 4.03.
- C2: g = 5.27, ω₀ = 3.434.
- C3: g = 5.27, isotropic Hartree, ω₀ = 3.87.

**Results.**
- **Radial branch.** At k = 0.52, Ω₊ = 7.7–8.6 and the lab line is at −4.2 to −4.6. |v/u| is 0.019–0.050 at every k,
  so radial quasiparticles are almost pure lab antiparticles.
- **Phonon branch.** The minor-to-major amplitude ratio is 0.64 at k = 0.52, 0.47 at 1.05, 0.35 at 1.57, and falls to
  0.05 at 7.85. Both lab lines are positive while Ω₋ < ω₀.
- **Predicted RMS conjugate/normal ratio over 0–6** (phase-independent for a circular orbit): 0.71–0.80 at k = 0.52,
  0.54–0.60 at k = 1.05, 0.36–0.45 at k = 1.57, and 0.055–0.080 at k = 7.85.
- **Against the measurement.**
  - The IR order of magnitude is right, but the circular model has no phase dependence.
  - In the UV, the measured ratio is 0.25–0.32, against the circular 0.055–0.08. The circular model misses the static
    (non-rotating) anomalous term 2λab of an eccentric orbit ψ_c ≈ ae^{iωs} + be^{−iωs}. That term couples +ω to −ω
    resonantly at every k.

### (b) Floquet–Bogoliubov driven by the measured orbit [C]

We ran `hfb_driven.py` and `floquet.py` with three closures for H(s):
- **bare:** the condensate only, i.e. BdG around the eccentric orbit;
- **Hartree:** isotropic;
- **meanH:** the exact spatial mean of the Hessian, ⟨∂²V⟩(s). This is the Gaussian/HFB mean field, including the
  anisotropic fluctuation (pair) covariance, which reduces the orbit-averaged static anisotropy from 9.2 to 3.5.

**Bare BdG.** It is parametrically unstable for k² ≤ 7.1, with Floquet rates of 0.19–0.82 per unit s. It does not
describe the medium: its L2 errors are 2.4–18 times the signal.

**Isotropic Hartree.** The ratio is 1.10–1.24 for k ≤ 5.2 and 0.84–1.07 above, which is wrong in the UV. Its normal-channel error is
0.71–0.86.

**meanH, UV (k ≥ 2.6).**
- The pole positions match the measurement to ≤ 2.3%, and to ≤ 0.6% for k ≥ 4.7.
- The ratio matches to typically 5–15%, with a range of −27% to +21%. For example, shell 15 gives 0.343 / 0.261
  against measured 0.30–0.32 / 0.247.
- The early-time anomalous response follows the measured one closely.
- Two-pole fits of the measured response leave residuals of 6–20%. The fitted widths are Γ = 0.15–0.47, against a
  model baseline of 0.04–0.11, so Γ/ω = 0.017–0.088.
- These are narrow quasiparticles with a matrix (normal plus anomalous) structure.
- meanH is Floquet-stable, apart from a narrow band at k² ≈ 27.6 with rate 0.035.

**meanH, IR (k ≤ 2.1).** It fails after t ≈ 1:
- the normal-channel L2 error is 1.1–1.3;
- the measured dominant poles sit 8–20% below meanH (±3.5–4.5 against ±4.35–4.95);
- two-pole fits leave residuals of 24–38%, against 11–17% for the model's own trace; four poles are needed;
- scattering out of the kicked mode exceeds the forward amplitude within about one oscillation period (section 3).

**Answer.** The UV matrix response is quantitatively the Floquet–HFB Bogoliubov response of the eccentric condensate
plus its Gaussian fluctuations, but not of the condensate alone. The count-dominant IR is broad and does not conserve
momentum. It is not a set of momentum-diagonal quasiparticles at the scale of this measurement.

**Closure.**
- In the UV, a relativistic two-fluid closure of Zaremba–Nikuni–Griffin type is the natural one: a time-periodic
  condensate plus Floquet-HFB quasiparticles with anomalous terms.
- In the IR, and for the coherent sector, it is not justified. A matrix (Nambu) 2PI/Kadanoff–Baym treatment at NLO
  with an inhomogeneous condensate (Aarts–Berges; Berges), or the classical-statistical field itself, is required [A].

## 5. Consequence for counting [D/C/E]

**An identity [D/C].**
- A₊ − A₋ = Q to 2e-4 (the spectral charge identity). So A₋ = (S − Q)/2, where S = Σ|ω|P.
- Across all 18 windows, S per net charge is 4.483–4.504. The constancy of the X̄ budget is exactly the constancy of
  the total frequency action S.

**The coherent mode [C].** We computed J_r = (1/2π)∮p_r dr per radial period. It is the k = 0 radial-branch
quasiparticle action about the circular condensate with the same charge, and it is what slow matching preserves.
- **Same windows ([results/coherent_counting.json](results/coherent_counting.json)):**

  | Window (s) | J_r per net charge | Lab statistic |
  |---|---:|---:|
  | 14–24 | 0.6237 | 0.6075 |
  | 24–34 | 0.4798 | 0.4760 |
  | 0–6 | 1.6815 | 1.6657 |

  J_r is 0.8–2.7% larger than the lab statistic.
- Our own Fourier analysis of the published trace reproduces Codex's coherent window values exactly.
- This supports Codex's distinction between radial action and lab statistic (claim 3).

**Circular-reference (adiabatic) basis [E].** Invert the measured per-shell (A₋, A₊) of window s 14–24 onto the two
circular branches.
- The waves' radial-branch action is 1.1465–1.1508 (C1–C3), against a lab value of 1.1442.
- The total is 1.770–1.775 per net charge, against Codex's 1.7517: a 1.1–1.3% increase.

**Dynamical (Floquet–HFB) basis [C/E].** This basis reproduces the UV response.
- Its quasiparticles are particle–antiparticle mixtures:
  - lab A₋/W ranges over 0.30–0.70 for the two branches;
  - (A₋ + A₊)/W ranges over 1.00–1.16;
  - a pure harmonic-ellipse control gives 0.51 / 0.49 at k² = 57.8.
- The orbit-averaged static anisotropy is ΔH = 2.9–3.5, and its axis precesses at −0.045 rad per unit s.
  - The splitting is δω ≈ ΔH/2ω_k.
  - A single wave mode then converts X ↔ X̄ with a half-period π/δω of about 8 (IR) to 16 (k ≈ 7.9) s units.
  - That is shorter than the measured span, s 14–34.
- The branch inversion of the measured (A₋, A₊) is ill-conditioned: condition numbers reach 50 and 494, and some
  branch actions come out negative. **X̄ is not a quasiparticle quantum number for the waves in this basis.**
- What survives is the exact Q and the approximate S. A branch-resolved count would need frequency resolution better
  than δω ≈ 0.2–0.4, which means windows of Δs ≳ 15–30, against Codex's 10.

## Verdicts

| Item | Codex | Ours | Verdict |
|---|---|---|---|
| Conjugate/normal ratio, count-dominant IR | 0.44–1.70 (3 probes) | Published arrays reproduce 0.4421–1.7043 bitwise. Low-noise traces give 0.46–1.31 (s0 = 20: 1.309 / 0.884 / 0.747; s0 = 20.44: 0.501 / 0.456 / 0.588) | **CORRECTED** (probe noise; the phase dependence is real) |
| Canonical sum | within 0.4% | 0.397% | CONFIRMED. It is an identity and tests only the transform [D] |
| Free control | 2.2e-10, zero mixing | 2.161928e-10 (own Yoshida) against 2.161931e-10 | CONFIRMED |
| Channel definitions and normalisation | stated | δψ = G_N p + G_A p*; identity holds to 1e-15. The estimator is a \|src\|²-weighted shell trace | CONFIRMED |
| "Matrix response; does not exclude quasiparticles" | stated | UV: Floquet–HFB quasiparticles (poles ≤ 2.3%, ratio ±15%, Γ/ω ≤ 0.09). IR: broad, momentum-mixing, not mean-field. Condensate-only BdG is unstable | CONFIRMED + NEW |
| Coherent counting: radial action against lab statistic | 1.763725 / 1.737788 | J_r is the k = 0 radial-branch Bogoliubov action. Same-window J_r/lab is 0.6237/0.6075, 0.4798/0.4760, 1.6815/1.6657 | CONFIRMED (Bogoliubov reading) |
| X̄ per net charge in a quasiparticle basis | not given | circular/adiabatic basis 1.770–1.775 (+1.1–1.3%). Floquet basis: undefined (mixed branches, conversion half-period 8–16, ill-conditioned) | NEW [E] |
| X̄ budget = (S − Q)/2 | not stated | S per net charge 4.483–4.504 over 18 windows | NEW [D/C] |

## Limits

- **Coverage.**
  - Only base48 at s0 = 20 and 20.44 was probed, with one box (L12, N48).
  - Shells 1–4 were probed with our low-noise estimators; shells 5–15 use Codex's probes.
  - The response window is 6 units, so frequency resolution is about 1. The pole fits are descriptive.
- **Models.**
  - The Floquet idealisation repeats one measured radial period. The real background is not exactly periodic, and
    the fluctuations evolve.
  - φ mixing is neglected in the Floquet analysis. It is ≤ 4% of the normal RMS, except one noisy shell-1 probe at 11%.
  - The circular counting inversion assumes incoherent branch occupations and uses shell means of the FD4 k².
- **Physics not covered.**
  - The quantum half and renormalisation are untouched.
  - The physical n̄, f_sym, threshold, binding, heavy and gauge sectors, and the vacuum density remain open.
  - Candidate 1 is unchanged.
- **Credit.**
  - Bogoliubov/HFB and ZNG two-fluid kinetics: Zaremba, Nikuni and Griffin.
  - Nonequilibrium 2PI: Aarts and Berges (hep-ph/0103049); Berges and Sexty (1201.0687).
  - Classical-statistical validity: Epelbaum, Gelis and Wu (1402.0115).

## Reproduction

```bash
source $SP/review/pt_response/env.sh        # SP, GRUT_* workdirs, PKC = copied packet
T=research/reality_framework/grut_formation_review_20261025/tracks/pt_response
P=$SP/f55/research/reality_framework/grut_particle_transport_20261027
python3 -I $T/code/reanalyse_published.py $P $T/results/reanalysis_published.json
(cd $PKC/code && python3 -B $OLDPWD/$T/code/regen_base48_states.py $P/FIELD_TRACES.npz)        # 3.3 min
(cd $PKC/code && python3 -B -c "import response; response.run(20.,7); response.run(20.44,7); response.run(20.,11)")
(cd $PKC/code && python3 -B $OLDPWD/$T/code/record_mean_hessian.py $SP/review/pt_response/hessian_s20_26.8.npz)
# probes, chunked (about 13 s per sample for 36 tangents): repeat until FINISHED
(cd $PKC/code && python3 -B $OLDPWD/$T/code/probe_study.py 20 exact1 0 $SP/review/pt_response/probes 40)
(cd $PKC/code && python3 -B $OLDPWD/$T/code/probe_study.py 20 hutch 8 $SP/review/pt_response/probes 80)
(cd $PKC/code && python3 -B $OLDPWD/$T/code/probe_study.py 20.44 hutch 8 $SP/review/pt_response/probes 80)
cd $SP/review/pt_response
python3 -I $T/code/analyse_probes.py probes probes_analysis.json
python3 -I $T/code/hfb_driven.py hessian_s20_26.8.npz $P hfb_driven_v2.json probes
python3 -I $T/code/bogoliubov_circular.py hessian_s20_26.8.npz $P circular.json
python3 -I $T/code/floquet.py hessian_s20_26.8.npz $P floquet.json ; python3 -I $T/code/floquet.py --test
python3 -I $T/code/fit_poles.py probes hfb_driven_v2_series.npz $P poles.json
python3 -I $T/code/coherent_counting.py hessian_s20_26.8.npz $P coherent_counting.json
python3 -I $T/code/own_tangent_check.py pwork/base48/state_s20.npy probes/PROJ_exact1_s20.npz own_tangent_check.json
```

Large intermediate files (the PROJ_*.npz projections, 17–37 MB, and the series npz) stay in `$SP/review/pt_response/`.
