# en_critic: completeness critic for the review of Codex 7d9c93b (issue #13)

27 October 2026. Packet `research/reality_framework/grut_early_nonlinear_20261027/` at 7d9c93b, read through git and a
fresh `git archive` copy in `$SP/review/en_critic/pkt`. Inputs from four tracks (en_rerun, en_indep, en_neutral,
en_particles) and their skeptics. This track does not redo their work. It fills gaps in coverage, settles
contradictions between tracks, flags verdicts that go beyond the evidence, and states the main correction and the
main physics point.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.
Every number below says whether it comes from **own** code, **own integration of Codex's equations**, or
**reproduction**. Nothing here validates or excludes anything.

## 1. Coverage of Codex's review request

| Codex asked for | Covered by | Gap found | Settled here |
|---|---|---|---|
| MATHEMATICS.md | en_rerun (forces, FD4, Yoshida, clock), en_particles (n± formula) | Fixed-charge energy bound and scaling, frozen-Σ/quadratic-Y formula, gauge (e = 0) scope: unchecked | yes, §2 |
| match.py / initial.py | en_rerun, en_rerun skeptic (purity, χ² test) | none | – |
| early3d adaptation | en_rerun audit, two independent engines | adaptation versus parent 8159556 not diffed | yes, §2 |
| export conventions | en_rerun, en_particles (+ skeptic) | none | – |
| short engine and seed checks | en_rerun (bitwise) | none | – |
| one main history | en_indep, en_rerun skeptic (independent engines) | none | – |
| full-versus-reduced scope | en_rerun (reproduction only) | claim 5 **never independently integrated**; e = 0 not assessed | yes, §3 |
| same-seed neutral comparison and conventions | en_neutral, en_indep skeptic, en_rerun | **X-seed-amplitude dependence** raised by en_neutral skeptic but unresolved; counterfactual only at N32 | yes, §4 |

## 2. MATHEMATICS.md items nobody checked [D]/[C], own code (`code/math_checks.py`)

**Fixed-charge energy bound.** Codex calls (3/4)λ_X^{1/3}|n_Q|^{4/3}V a "homogeneous circular lower bound". It is in
fact a bound on **every** lattice state at fixed total charge:
- Cauchy–Schwarz gives |π_X|² ≥ j²/|χ_X|² pointwise, where j = χ₁π₂ − χ₂π₁.
- The FD4 symbol 4sin²(x/2)(1 + sin²(x/2)/3)/dx² is ≥ 0, and the φ terms are ≥ 0.
- Minimizing j²/(2r²) + λr⁴/4 over r gives (3/4)λ^{1/3}|j|^{4/3}, and Jensen gives ∫|j|^{4/3} ≥ V^{−1/3}|Q|^{4/3}.

Numerically:
- E_min/E = 0.11594 in all seven histories, constant to 2e-9.
- The spatial-mean energy exceeds the bound at its own charge Q₀ by ≥ 0.11 ρ at all times.
- At the τ-800 snapshot, the pointwise step is ≥ 0, and the local bound is 5.9 × the global one.

The infinite-space scaling (χ_R = χ(x/R)/R, π_R = π(x/R)/R²: Q fixed, E → E/R) is correct [D]. The remaining 88% of the
energy above this bound is the "allowance" in REVIEW_TIMING; it is not energy outside the mean.

**Frozen-Σ / quadratic-Y estimate.** From V ⊃ M²|y|² + 2Re(κ y* x³ + μ s* y²) + λ_X|x|⁴, the static quadratic Y is
M²y + 2μ s y* = −κx³. Its smallest eigenmass is M² − √2 μ h = 495. Then |ΔF_X|/|F_λ| ≤ 3κ²|x|²/(2λ_X(M² − √2μh)), and
with |x|² = |χ|²/(2a²) this is Codex's formula [D]. Over 2000 random X and Σ phases, the exact static solution divided
by the bound is 0.38–1.00, so the bound is valid and tight [C]. It decreases as max|χ|²/a² [D].

**early3d versus parent field3d.cpp (8159556).**
- Only the clock, output schedule and diagnostics changed; the kick, drift and FD4 code are the same.
- Latent, harmless: the kick keeps λ_φ(φ² − ν²)φ while the energy diagnostic uses λ_φφ⁴/4. Every input file writes
  ν = 0, so the two agree.

## 3. Claim 5 (full versus reduced), own integration [C]

`code/background_match.py` and `code/tangent_match.py`:
- Own scipy DOP853 integration (rtol 1e-12) of Codex's full homogeneous equations, `reheating.rhs`, plus the mode
  generator at e = 0 and the neutral mass. These are imported read-only; the equations are Codex's.
- Own reduced orbit and tangents from τ 20, with my own apocentre normalization.

| Quantity | Own | Codex |
|---|---|---|
| τ-20 state vs MATCH_L12 | χ 7.7e-13, χ′ 2.3e-13 relative; a 1.7e-13, C 1.0e-12 | – |
| background, orbit-scaled | 0.5470872% (maximum at τ 47.2) | 0.5470872% |
| charged tangent k = 2.13 / 7.1 | 0.9774883% / 1.0307871% | 0.9774883% / 1.0307871% |
| neutral transfer k = 2.13 / 7.1 | 2.7171576% / 1.4820005% | 2.7171576% / 1.4820005% |
| initial force mismatch | 0.2056958% | 0.2056958% |
| radiation clock, τ 20–50 | a within 1.3e-6 | – |

Verdict: **CONFIRMED**, now by own integration and no longer by reproduction alone. The scope is unchanged: two momenta,
e = 0 on both sides, through τ 50.

**Gauge omission [E].**
- The model's e_X is 1e-7 (1e-6 as a control).
- The Coulomb energy of the measured charge fluctuations, ½e²Var(j)/k_min², is ≤ 2.1e-10 ρ (e = 1e-7) or 2.1e-8 ρ
  (e = 1e-6) over τ 20–800.
- The gauge field matters only for k ≲ e|χ|_max ≈ 1.4e-5 (or 1.4e-4) H_i, far outside the L = 12 box.
- So e = 0 is a justified approximation inside the box. The cosmological-scale Gauss constraint is outside the pilot.

## 4. X-seed-amplitude control: the open neutral issue

`code/xseed_scan.py` and `code/analyze_xseed.py`.
- **Engine:** Codex's unmodified early3d.cpp, compiled from git with Codex's flags (reproduction-type engine).
- **Initial data (own):** the published N48 / L12 / seed-1 τ-20 field. My rebuilt f = 1 file has SHA-256 abe46838…,
  equal to the published one. Only the nonzero Fourier modes of (χ₁, χ₂, π₁, π₂) are multiplied by f. The X zero mode
  and the φ seed are unchanged.
- **Run settings:** threads 2, h_max 0.004, sampling every 0.05 τ.
- **Check:** an f = 1 run to τ 100 matches the published N48 curves to ≤ 1.3e-13.

| X seed f | Charge contrast 1 | X variance ½ | φ+int leaves linear control (−10%) | φ+int at τ 800 | φ+int mean, τ 700–800 | φ+int max |
|---|---|---|---|---|---|---|
| 1 (published) | 60.6 | 69.25 | 86.6 | 0.310% | 0.311% | 0.344% |
| 0.1 | 137.75 | 152.95 | 158.35 | 0.796% | 0.811% | 0.899% |
| 0.01 | 230.75 | 244.3 | 250.2 | 2.64% | 2.40% | 2.73% |
| 0 | 370.3 | 388.0 | 304.3 | 5.78% | 6.67% | 7.95% (τ 474) |

Energy drift is 1.27e-8 in every run. Data: `results/xseed_scan.json`.

**Reading [C]/[A].**
- **Claim 3's qualitative statement is robust to the X seed.** Even with no X seed at all, φ peaks at 7.9% and never
  approaches the fixed-orbit 2.09×. With no X seed, X fragments from φ-sourced fluctuations near τ 370.
- **The number 0.28–0.36% is not robust.** It is conditional on the vacuum-like half-quantum X seed: it rises ×2.6–3.3
  per decade of smaller seed, up to 6–8% with none.
- **Physical direction.** f < 1 needs a sub-vacuum (squeezed) growing quadrature. A more excited state, or seeding of
  the unseeded IR k < 1.8 modes, acts like f > 1. On the observed monotonic trend, 0.3% is then on the high side
  among unsqueezed preparations [A].
- **Fragmentation delay.** It is Δη ≈ 3.5–4.0 per decade of seed, giving μ_X ≈ 0.6, the Floquet value.
- **Correction to en_neutral.** Its N32 no-X-seed counterfactual gave 4.3–5.9%. At N48 the result is 5.5–7.8% over
  τ 700–800. N32 under-resolves the k ≈ 7 band: its linear control reaches only 0.995.
- **Remaining grid caveat.** The N48 linear control is 7% below the continuum value, so the no-seed value could be
  slightly higher on finer grids [E].

**Mechanism (settles the en_neutral contradiction).**
- For f = 1, 0.1 and 0.01, the nonlinear φ curve leaves the same-seed linear control 19.5–26 τ after X charge contrast
  reaches 1. Fragmentation timing therefore controls the termination.
- With no X seed, φ's own growth to about 1.5% ends the resonance first: φ+int leaves the control at τ 304, X contrast
  reaches 1 at τ 370.
- The en_neutral skeptic's factorial shows that pump depletion and Hartree detuning each suffice in homogeneous proxies.
  The mechanism is fragmentation-driven but not attributable to one channel. en_neutral's "detuning dominant" stays
  REFUTED.

## 5. Particle budget across X seeds

`code/zero_mode_refs.py`, own code.
- **Zero mode:** lab-frame windowed (Hann) spectrum of ψ̄ = χ̄₁ + iχ̄₂, uniform in η. Antiparticle content is
  n₀₋ = ∫_{ν<0}|ν|S dν/2π. The charge reconstructs to 1e-4.
- **Waves:** Codex's two reference projections at τ 800, plus a fixed m² = 12.

| f | Zero mode n₀₋/Q, τ 340–798 / 547–800 / 700–800 | Hartree, τ 800 | Massless, τ 800 | m² = 12, τ 800 | Rough total |
|---|---|---|---|---|---|
| 1 (reproduced trace) | 0.608 / 0.595 / 0.596 | 1.1800 | 2.2375 | 1.169 | ≈ 1.78 |
| 0.1 | 0.725 / 0.616 / 0.575 | 1.235 | 1.321 | 1.158 | ≈ 1.81 |
| 0.01 | 0.922 / 0.793 / 0.701 | 0.995 | 2.325 | 1.105 | ≈ 1.70 |
| 0 | 0.980 / 0.955 / 0.798 | 0.851 | 0.797 | 0.844 | ≈ 1.65 |

Before fragmentation the zero mode alone holds 1.72–1.73 (Hann; f = 0 over τ 25–330, f = 1 over τ 25–55). en_particles
gives 1.739 exactly, and the ~1% difference is the known estimator offset.

Reading [E]:
- The f = 1 zero-mode content reproduces en_particles' 0.59–0.61 with a different estimator.
- Across seeds the split moves: the zero mode goes 0.60 → 0.80–0.98 and the waves 1.18 → 0.85. The sum, 1.65–1.81, is
  consistent with a conserved ≈ 1.74 budget within the ±7% single-phase uncertainty of the equal-time wave count.
- The massless-versus-Hartree gap is not even a stable feature of the reference pair: at f = 0, massless < Hartree.

Harmonic check [D]: an ellipse with semi-axes A, B carries n₋/Q = (A − B)²/(4AB). With r_max/r_min = 8.06 this is
1.55. The quartic orbit gives 1.74, so the budget is set by the eccentricity of the τ-20 condensate orbit.

## 6. Contradictions settled

- **Claim 4.** The values are CONFIRMED. The "roughly 1.9-fold" gap is CORRECTED: it is one orbit phase (instantaneous
  ratio 0.96–2.54, orbit mean about 1.6). This resolves en_indep (CORRECTED) against its skeptic (CONFIRMED), and
  en_particles against en_rerun.
- **Claim 1b orbit average.** Convention-dependent (ratio of means ≤ 0.47–0.50; mean of ratio up to 0.51–0.53). Codex's
  69.2 is an instantaneous pericentre crossing, which Codex already caveats.
- **en_rerun N3 "floor".** It depends on the allowed transformation class: 1.06 with one ω per shell, 0.889 with shear,
  0.67 per mode. It is not a floor and should be withdrawn as one.
- **Shared tree.** `$SP/e7d` now matches git exactly (`diff -rq` against a fresh archive is clean); the stale-tree
  warnings no longer apply.
- **Machine.** It has 4 CPUs; nproc reports 2 when OMP_NUM_THREADS=2.

## 7. Verdicts stronger than the evidence

Track statements:
- en_particles 7b, "Hartree selected": refuted, because the count is flat in m².
- en_neutral's mechanism attribution: refuted.
- en_indep's CORRECTED on claim 4: too strong. Its N1 "no exponential growth": slow infrared growth (×5–20) is resolved.
- en_rerun's claim-3 sentence ("agree within 4% to τ 100") is wrong, and its τ 600–800 statistics are aliased.
- en_neutral's no-seed magnitude: an N32 artefact (§4).
- The en_neutral skeptic's "above-half-quantum φ energy 0.067–0.072% of ρ₀": single-track and not cross-checked.

Codex statements:
- "0.28–0.36%" holds only for the vacuum-like X seed.
- "1.9-fold" is one orbit phase.
- SPECTRAL_EXPORT is a single-phase snapshot.

## 8. Main correction and physics point for Codex

**Correction.**
- Every τ-800 number is one sample of a strongly eccentric coherent orbit: period Δη ≈ 0.89, and |χ̄|² swings about
  650–10 200 in the late orbit.
- The seven controls happen to share the phase, so the control spreads are valid. The orbit bias is about +2 points
  for energy outside the mean.
- The export and both reference projections inherit the swing: massless 1.18–2.77, Hartree 1.07–1.24, ratio 0.96–2.54.
- Publish period-averaged exports over several radial periods (Δτ ≳ 100–250 near τ 800) and state the 0.28–0.36% as
  conditional on the X seed (§4).

**Physics: what the particle count now needs.**
- **The massless/Hartree choice is settled classically.** Measured spectra are gapped (statistical m²_eff 8–13,
  response m² about 19), and the wave count is flat at 1.15–1.21 for m² from 5 to 26. The waves hold 1.12–1.24 X̄ per
  net charge.
- **The open question is a classical X̄ budget of about 1.74 per net charge.** It is conserved over Δη ≈ 34 and fixed
  by the eccentricity of the τ-20 condensate orbit. About a third of it (0.52–0.61) is still the condensate's
  counter-rotating content, and 70% of the wave part sits at k ≤ 1.6, in modes that were never seeded and are limited
  by the box.

The kinetic input therefore needs three things:
1. A physical prescription for the condensate's counter-rotating content: follow the zero mode until it decays, or
   show it is not X̄.
2. Infrared-seed and box-size controls.
3. The number-changing (chemical) relaxation rate in the measured medium, continued with the physical non-conformal
   masses. The τ-800 window bounds relaxation only to ≳ 5×10³ conformal units, against about 10¹¹–10¹² before binding
   at a ≈ 10¹²–10¹³ (REVIEW_TIMING [Q]). That rate decides between the conserved ≈ 1.7 and the equilibrated ≈ 0.1.

The classical-statistical method follows Khlebnikov–Tkachev and Micha–Tkachev (number-conserving cascade) [R]. Its
ultraviolet limits are given by Epelbaum–Gelis–Wu, arXiv:1402.0115 [Q].

## Not certified by this track

- The X-seed scan is N48, one seed realization, and τ ≤ 800. f < 1 is a sensitivity control, not a physical state.
- The particle totals in §5 mix a window-averaged zero mode with an equal-time wave count (±7%).
- No infrared-seed or box control was run.
- The matching (§3) uses Codex's full-model equations, so the full model itself is not re-derived.
- No physical n̄, f_sym, threshold or abundance is claimed.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; K=$SP/review/en_critic
T=research/reality_framework/grut_formation_review_20261025/tracks/en_critic; RF=$SP/e7d/research/reality_framework
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 GRUT_REHEATING_WORKDIR=$K/rwork GRUT_MODE_WORKDIR=$K/mwork
git archive 7d9c93b research/reality_framework/grut_early_nonlinear_20261027 | tar -x -C $K/pkt
P=$K/pkt/research/reality_framework/grut_early_nonlinear_20261027
g++ -std=c++17 -O3 -ffp-contract=off -fopenmp $P/code/early3d.cpp -o $K/bin/early3d
for f in 0 0.01 0.1; do python3 -I $T/code/xseed_scan.py $P $K/bin/early3d $K/runs/f$f $f 800; done   # ~3.3 min each
python3 -I $T/code/analyze_xseed.py $P $K/runs $T/results/xseed_scan.json
python3 -I $T/code/zero_mode_refs.py $T/results/zero_mode_refs.json f1_enrerun=$SP/review/en_rerun/ework/N48 f0.1=$K/runs/f0.1 f0.01=$K/runs/f0.01 f0=$K/runs/f0
python3 -B -I $T/code/background_match.py $RF $P $T/results/background_match.json
python3 -B -I $T/code/tangent_match.py $RF $T/results/tangent_match.json                                 # 80 s
python3 -B -I $T/code/math_checks.py $RF $P $SP/review/en_rerun/ework/N48/state.tsv $SP/review/en_rerun/ework/N48/state_s0000040.bin $T/results/math_checks.json
```

The runs in this session used the directory names f1em1, f1em2 and f0 for f = 0.1, 0.01 and 0.
