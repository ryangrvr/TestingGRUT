# st_physics_skeptic: adversarial re-check of the st_physics review of Codex b73de64

27 October 2026. Claude (track st_physics_skeptic). This file checks the findings in [NOTES.md](NOTES.md) and
`results/*.json`. For each claim I recomputed the key number with a different method where I could, rather than
re-reading the reviewer's prose. Rules: [BRIEF.md](../../BRIEF.md) and [BRIEF_TIMING.md](../../BRIEF_TIMING.md).
Everything is conditional on the supplied working model. Nothing here validates or excludes anything. Code and
results are in `skeptic/`. The Codex tree in `$SP/b73` was only read, and nothing was committed.

## How each number was produced

| Check | What | Independence |
|---|---|---|
| k1 | Floquet spectrum written in the **co-rotating radial/tangential frame**, monodromy taken pericentre to pericentre | Own code. A different formulation from the reviewer's lab-frame monodromy, which starts at apocentre and undoes the rotation afterwards. |
| k2 | Conformal time η by direct quadrature of ∫dτ/a. Late growth from **stroboscopic pericentre samples** of the published arrays and saved vectors. C, the radiation tail, L(τ) | Own analysis of Codex's published arrays, used as data |
| k3 | **Continuous between-sample maxima** around each published unit crossing: X-only propagation over every 0.5-τ interval from Codex's saved 12-vectors. Size of the neglected nonlinear terms | Own code on published data |
| k4 | **Full six-field linear response in physical variables** (δz, δż). It uses a numerical Jacobian (exact 4-point stencil) of the model force, our own DOP853 integration to τ = 300 and our own seeds and covariances. No `response.c`, no Codex Hessian or Ward projection, gauge terms dropped (e²/k² ≈ 2e-15). | Independent perturbation equations. The background `rhs`, `potential_force` and `prepared` are imported from Codex's Python as the *model definition*. |
| k5 | δφ parametric resonance on the matched X orbits, using the declared coupling g_X²φ²\|X\|² with g_X = 0.1 [Q] | Own code |
| k6 | Arithmetic re-derivations (P6, P7, P11–P13, Y_L clocks) | Own code |

## Verdicts on the reviewer's claims

| ID | Reviewer verdict | My recomputation | Refuted? |
|---|---|---|---|
| P1 μ_eff | CONFIRMED 0.6335 / 0.5599 | 0.63355 (0.63270–0.63390) / 0.55994 (0.55921–0.56026), with quadrature η rather than a/C | no |
| P2 late growth = quartic self-resonance | NEW | Holds. Saved-vector strobe at k = 2.13 gives 0.63357 against my Floquet 0.63315 (3,.5), and 0.54704 against 0.54743 (5,1). Charge-array strobe fits are within ≤ 0.1% for k = 1.2–2.8. Amplitude-array strobe fits deviate by up to 1.5% (at k = 1.2), so the "≤ 0.5% at every in-band k" wording depends on the method; ~1% is safe. | no |
| P3 auxiliary orbit | CONFIRMED | Co-rotating method: Codex's auxiliary orbit is the (3,.5) τ10 orbit (E = 105360.33, L = 8304.02), peak 0.631778 at k = 2.12780, μ/k = 0.43980, edge 3.13949, second band 5.57118–5.97645 with peak 0.16027 at 5.775. (5,1) τ10: 0.563131 at 1.88748, 0.44211, 2.78589, 4.92716–5.28682. | no |
| P4 k = 0 wording | CORRECTED | M−I at k = 0 has singular values (265, 5.6e-3, 3e-12, 1e-14). That is rank 2 with all multipliers equal to 1, so there are two 2×2 Jordan blocks (confirmed). The seed's k → 0 content (θ_X 0, θ_Y 0.169, θ_Σ 1.000) is not the 1:3:6 global rotation. Note that Codex's sentence sits in the *auxiliary-orbit* paragraph, so this is a context fix, not an error. | no |
| P5 unit samples at pericentre | CONFIRMED | r/r_min = 1.010, 1.143, 1.041, 1.086. **Continuous** reconstruction confirms that the first true unit crossing falls in the same pericentre passage. (3,.5) 0.03: sampled 1.485, continuous 1.585; the previous passage peaks at 0.945. (5,1), k = 2.13: 1.318 / 1.373; previous passage 0.829. Same pattern for the 1e-6 seed (previous passages 0.946 / 0.656). \|δχ_X\| ≥ r_max at τ 197.5 / 215.5. | no |
| P6 1e-9 seed | CONFIRMED | a = 77.017, τ = 1493.7; a = 66.257, τ = 1487.1. This linear extrapolation is moot physically (see P8). | no |
| P7 seed threshold | NEW | log10 s = −723.17 / −741.19 (arithmetic correct). Superseded by the X vacuum floor and by the δφ channel. | no |
| **P8** "X stays quartic-dominated; only nonlinear backreaction ends it; nothing in the declared model gives X a mass before a ~ 10¹²" | NEW | **Omits an in-model channel.** The declared g_X²φ²\|X\|² gives a conformal δφ mass² g²r² ∈ [3.1, 204] on the eccentric orbit. The Hill equation has a resonance band at k ≈ 6.7–7.3 with **μ_φ = 0.206** (3,.5) and 5.95–6.45 with **μ_φ = 0.185** (5,1). These rates are 0.33 of the X rate and strongly non-monotonic in g (0.81 at g = 0.05, 0.61 at g = 0.14). From a vacuum-like start, g²⟨δφ²⟩ reaches λr_min² after Δη ≈ 8 (τ ~ 150) and plausibly reaches O(1) by Δη ≈ 16 (τ ~ 400) [E]. That is before the 1e-6 crossing at τ 772/782. Codex's six-field system omits δφ. Smaller points: the "0.25 at a = 1" Y-sextic ratio is 0.141 for the actual a = 1 state, because 0.25 rescales the late apocentre. Higgs portals are zero in the model, so the thermal-φ caveat is correctly flagged as outside the model. | **yes (partial)** |
| P9 physical vacuum seed | CORRECTED | **Headline independently reproduced** with my own equations, upgrading it from reproduction-type. Seed per rad: dln\|X\| −0.04304, θ_X 0, θ_Y 0.16922, θ_Σ 1.00000. X-rotation weight 87.65 / 193.97 at k = 2.13. Vacuum A (the reviewer's covariance): 0.306 / 0.736 supplied rad per ln k. Covariance sensitivity: massless BD 0.255 / 0.603, Codex's regulator 0.241 / 0.565, X-only 0.253 / 0.601, giving **8–25×** the 0.03 seed rather than 10–25×. **Refuted sub-number:** "band RMS amplitude = 1 at τ 83.5 / 80.5" is a sampling artefact of Codex's pre-τ80 grid, which has gaps up to 3.2 τ. With 0.05-τ sampling it is **τ 69.3 / 66.9** (66.9–69.5 across covariances), and the charge crossing is 61.2 / 57.5. The (5,1) weight 86.5 at k = 2.8 is out of band (edge 2.786); in-band (5,1) weights are 194–221. | no (headline); sub-number yes |
| P10 no Q-balls at φ = 0 | NEW | \|x\| = √(λM_L²)/κ = 203.4 against 76.4 at a = 1. Arithmetic and logic hold. | no |
| P11 circular orbit | NEW | E_circ = 12611.9 (E_circ/E 0.11970, release 0.8803); 7011.1 (0.10854, 0.8915) | no |
| P12 X̄ estimate | NEW [E] | Pairs 1.361–2.625 / 1.509–2.853; gas 1.046–2.481 / 1.193–2.701; equilibrium 0.1067 / 0.1159 | no |
| **P13** cube charge | CORRECTED | Codex's number is exact at the k it states: the REPORT says "at k/H_i = 2.13", giving 2.0937e5 / 1.3440e5. The "primary peak" is the auxiliary orbit's peak (2.128). The own-peak values (2.100e5 / 1.932e5) are a useful clarification, not a correction of an error, so the verdict should be **CONFIRMED + note**. | **yes (verdict overstated)** |
| P14 C and tail | CONFIRMED | C = 1.985805 / 1.475868; tail matches to 6.08e-10 / 1.99e-9 | no |
| **P15** L −0.5% "inside a 1.6–1.8% oscillation, not a trend" | CONFIRMED + qualifier | The −0.50% endpoint is confirmed, but **the qualifier is wrong**. The 1.6–1.8% peak to peak comes from τ = 20–64, when L was 1.3–1.6% lower. For τ ≥ 200, L oscillates only **0.31 / 0.35% p-p** within [−0.58, −0.28]% / [−0.68, −0.33]% of the 64–80 mean. η-window means after τ 80 are flat to ±0.1%. The −0.5% is a level offset, because the 64–80 reference window sits ~0.45% above the later plateau. It lies outside the late oscillation band, is not a late trend, and Codex's sentence is accurate as written. | **yes (qualifier)** |
| P16 ordering stated correctly | CONFIRMED | The ordering holds. **Wording item 3 is refuted.** At the 0.03 unit samples \|δχ_X\|/\|χ_X\| = 1.52 / 1.26 (not ≈ r_min). The cubic remainder alone is 16% / 8% of k²\|δχ\|. The quadratic remainder 2χ\|δ\|² + χ̄δ², which the reviewer left out, is 31% / 18%. The total nonlinear remainder is **39% / 23% of the linear force** at the sample and 49% / 28% at most over the preceding radial period, not "≤ 7%". Codex's "certainly lost small-perturbation control" is defensible. Item 8 ("unless an unmodelled X mass intervenes") should name the modelled-but-omitted δφ channel. | no (overall); item 3 yes |

## New issues with the Codex packet that the reviewer missed

1. **Pre-τ80 threshold "brackets" are not brackets [C].** Before τ = 80 the response is saved on the inherited
   MODE_PROFILES grid, with gaps of up to 3.2 τ. The radial period there is ~5–15 τ, and pericentre spikes are
   ~0.3 τ wide. With 0.05-τ sampling (k4, independent equations), the 0.03-seed 0.1-amplitude threshold on (3,.5)
   is first reached at **τ = 56.5 (a = 14.95)**, at k = 2.0/2.13/2.4. RESPONSE_RESULTS gives the bracket
   **[82.5, 83.0]**, and SOURCE_CLOCKS propagates τ = 83. On (5,1) the published 0.1 crossing is at k = 0.02,
   which I did not test. All unit crossings (after τ 80) are unaffected: dense sampling gives exactly 134.0 and 147.5.
2. **δφ parametric resonance is omitted [C/E].** See P8. The declared coupling gives μ_φ ≈ 0.21 / 0.18 at
   k ≈ 7.0 / 6.2. Codex says in general terms that "neutral production … may modify" the growth. The rate is now
   quantified, and it can act before the 1e-6 and 1e-9 crossings. It does not change the time ordering: it adds an
   early nonlinear channel.
3. **Positive finding [C]: independent-equation reproduction.** My physical-variable linearisation reproduces
   Codex's published amplitude and charge responses to ≤ 3.4e-9 (3,.5) and ≤ 6.6e-9 (5,1) relative through τ = 300
   at k = 2.0, 2.13, 2.4 and 2.8. The background matches to 7.5e-12 / 2.2e-14. My Jacobian's initial eigenvalues
   are −2.832, −2.000, 15.68, 529.06, 1320.08 and 3201.44. Neither Codex nor the reviewer had an independent check
   of the response arrays from the equations.
4. **Y_L clocks [C].** My own trapezoid plus analytic tail gives F_N = 1.3479e-5 / 1.4837e-5 and
   F_E = 3.7161e-8 / 4.3066e-8 at the 0.03 unit samples. These agree with Codex's numbers to 2e-6 relative.

## Limits

- k4 uses the model's background equations and potential from Codex's Python. It is independent of the
  perturbation code, not of the model definition.
- The δφ result is a Floquet rate on the matched late orbit [C]. Its backreaction time is an order-of-magnitude
  estimate [E]: vacuum normalisation, the pre-τ20 era and ξ_φ are not computed.
- The vacuum covariance at a = 1 remains a supplied choice. I tested four and found a factor of ~1.3 spread.
- Nothing nonlinear was evolved.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; R=$SP/b73/research/reality_framework
T=research/reality_framework/grut_formation_review_20261025/tracks/st_physics/skeptic; P=$R/grut_source_timing_20261026
export OMP_NUM_THREADS=1 GRUT_REHEATING_WORKDIR=$SP/review/st_physics_skeptic/rwork; D=$SP/review/st_physics_skeptic
python3 -I $T/code/k1_corot_floquet.py $P $R/grut_reheating_modes_20261026/QUARTIC_ORBIT.json $T/results/k1_corot_floquet.json  # 60 s
python3 -I $T/code/k2_strobe_growth.py $P $T/results/k1_corot_floquet.json $T/results/k2_strobe_growth.json
python3 -I $T/code/k3_continuous_crossing.py $P $T/results/k3_continuous_crossing.json
python3 -B $T/code/k4_phys_tangent.py $R quiet3 300 $D/k4_quiet3.npz   # ~2 min each
python3 -B $T/code/k4_phys_tangent.py $R quiet5 300 $D/k4_quiet5.npz
python3 -B $T/code/k4b_analysis.py $R $D $T/results/k4b_seed_vacuum.json
python3 -I $T/code/k5_phi_resonance.py $T/results/k1_corot_floquet.json $T/results/k5_phi_resonance.json  # 40 s
python3 -I $T/code/k6_arith.py $P $T/results/k2_strobe_growth.json $T/results/k1_corot_floquet.json $T/results/k6_arith.json
```
