# st_physics: physics and wording audit of Codex b73de64 (growth versus source clocks)

27 October 2026. Claude (track st_physics). Packet: `research/reality_framework/grut_source_timing_20261026/` at
b73de64, read from the git archive in `$SP/b73`, which was not modified. Rules: [BRIEF.md](../../BRIEF.md) and
[BRIEF_TIMING.md](../../BRIEF_TIMING.md). Everything is conditional on the supplied working model. Nothing here
validates or excludes anything.

**Independence.**

- **Own code:** `code/orbit_tools.py`, `s1`, `s2`, `s3` and `s5`. The Floquet routine is adapted from our own
  `tracks/rm_review/code/quartic_floquet.py`. `s2` and `s3` read Codex's *published* arrays (backgrounds, response
  arrays and saved tangent vectors) as data.
- **Codex engine with our own seeds (reproduction-type, labelled):** `s4` and `s4b`. These use Codex's
  `response.c`, compiled by `timing_common.build()` in our own work directory, with Codex's preparation and Hessian.
  The seeds, vacuum covariance and post-processing are ours.

## 1. Radiation mapping and the late growth exponent [D, C]

With a² = 2C(τ−τ\*), we have dη = dτ/a = da/C. After τ = 80 the published a(τ) follows this form to 6e-10 (3,.5) and
2e-9 (5,1). The published unit-amplitude samples give:

| History | a(.03) | a(1e-6) | C | μ_eff = ln(3e4)·C/Δa | spread over the sampling bracket | from the local-charge crossings |
|---|---:|---:|---:|---:|---:|---:|
| (3,0.5) | 23.0527 | 55.3653 | 1.98581 | **0.6335** | 0.6327–0.6339 | 0.6338 |
| (5,1) | 20.8777 | 48.0495 | 1.47587 | **0.5599** | 0.5592–0.5603 | 0.5601 |

**Matched quartic orbits (own code).** For the conformal X invariants on the published backgrounds after τ = 20,
E is constant to 0.21% and 0.25% peak to peak. L oscillates by 1.6% and 1.8% through exchange with Y and Σ.

| Orbit | E | L | r_max/r_min | T_η | peak μ at k | primary upper edge | second band | μ/k as k→0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| (3,.5), τ10 | 105360 | 8304.0 | 7.89 | 0.81404 | **0.6318 at 2.128** | 3.1395 | 5.571–5.976 | 0.4398 |
| (3,.5), late mean | 105298 | 8118.4 | 8.07 | 0.81434 | 0.6331 at 2.129 | 3.1418 | 5.569–5.975 | 0.4405 |
| (5,1), τ10 | 64596 | 5346.2 | 8.49 | 0.92061 | **0.5631 at 1.887** | 2.7859 | 4.927–5.287 | 0.4421 |
| (5,1), late mean | 64615 | 5207.8 | 8.72 | 0.92076 | 0.5645 at 1.889 | 2.7886 | 4.927–5.287 | 0.4429 |

For the τ10 (3,.5) orbit, the period, the precession and Codex's QUARTIC_ORBIT exponents are reproduced to 2e-10.
Symplectic error is ≤ 4e-13, and a Radau control agrees to 2e-13. Pure-quartic scaling (χ → sχ, with
s = (E5/E3)^{1/4} = 0.885) predicts 0.559 at k = 1.88 for (5,1). The small remaining difference comes from the
eccentricity, which is 0.204 for (3,.5) and 0.189 for (5,1).

**Is the late growth the quartic self-resonance? Yes, quantitatively [C].**

- **Per-momentum late growth.** We fitted the envelope of Codex's published per-radian amplitude and charge
  responses over τ = 200–790, using one radial period per window. For every in-band momentum the fitted rate matches
  the late-orbit Floquet exponent:

  | k | (3,.5) fit / Floquet | (5,1) fit / Floquet |
  |---|---:|---:|
  | 1.2 | 0.4748 / 0.4743 | 0.4644 / 0.4620 |
  | 2.0 | 0.6304 / 0.6294 | 0.5633 / 0.5611 |
  | 2.13 | 0.6339 / 0.6331 | 0.5494 / 0.5474 |
  | 2.8 | 0.4879 / 0.4881 | outside its band; fit ≈ 0 |
  | 5.8 | 0.1592 / 0.1590 | outside its second band; fit ≈ 0 |

  Codex's other k samples lie outside the bands and give fits of about 0. For k ≤ 0.1, the late response grows
  secularly, at about 1/η (fit 0.04–0.05 against Floquet 0.009–0.04).
- **Time-domain check.** We took Codex's saved 12-component tangent at τ = 200, kept only its X part, and propagated
  the X-only quartic linearisation on the published background. At τ = 800 it reproduces Codex's full six-field and
  longitudinal tangent at k = 2.13 to 0.17% (3,.5) and 0.15% (5,1). Over that interval the full response grows by
  1.3e4 and 4.6e3. On our own pure-quartic background the agreement is 0.6% and 1.4%.

The brief's hand estimates (0.634 and 0.560) are CONFIRMED. The (5,1) crossings use the sampled k = 2.0, whose growth
is 0.6% below the true peak at k = 1.89. This changes the timing negligibly: e^{0.06} over the interval.

## 2. Seed extrapolation and validity [D, E, A]

The relation is a_break(s) = a_ref + (C/μ) ln(s_ref/s), with a_ref taken from the 1e-6 samples.

- **The 10⁻⁹ seed.** The unit crossing would fall at a = 77.0 (τ ≈ 1494) and a = 66.3 (τ ≈ 1487). The 0.1 crossing
  would fall at τ ≈ 1197 and 1199. The largest 1e-9 response through τ = 800 is 1.6e-3 and 1.4e-3. This is
  CONFIRMED as consistent with "not reached by 800".
- **Seeds that would break down only after the Y_L median.** The Y_L medians fall at a = 5231 and 4510. A seed would
  need s ≤ 10^−723 (3,.5) or 10^−741 (5,1) to cross unity only after them.
- **Why nothing in the model shuts the resonance off before a ~ 10¹².**

  | Term | Size relative to the quartic λ\|χ\|² | Effect |
  |---|---|---|
  | Y-induced sextic −κ²\|x\|⁶/M_L² (strongest channel, M_L² = 495; Y_H gives 0.38× this) | 0.25 at a = 1; 6.3e-3 at τ10; 4.6e-4 at τ134; 8e-5 at a = 55; 9e-9 at a = 5.2e3 (apocentre, falls as 1/a²) | negligible after τ ≈ 10 |
  | Direct Y terms (B00 ∝ κy, and A01 at second order through Y_L), relative to k²+λ\|χ\|² at k = 2.13 | ≤ 7e-3 (τ 20–100), ≤ 1.6e-3 (τ 100–400), ≤ 7e-4 (τ 400–800) | negligible |
  | Gauge (longitudinal mass and Coulomb term), e²/λ | 1e-11 (e = 1e-7), 1e-9 (e = 1e-6) | consistent with Codex's e-insensitivity |
  | ξR and a''/a | R = 0 in radiation, up to an η²-suppressed X trace | negligible |

  The quartic λr² itself ranges from 0.31 (pericentre) to 20.4 (apocentre) for (3,.5), and from 0.21 to 16.0 for
  (5,1).
- **When could X stop being quartic-dominated?** The neutral field stays pinned while \|X\| > x_c = 8.5e-12 H_i.
  The pericentre falls below x_c at a ≈ 1.2–1.5e12 and the apocentre at a ≈ 1.1–1.2e13. The vacuum mass a·m_X equals
  √λ·r_c at a ≈ 2.5–2.9e12. **Within the declared model the conformal quartic dominates until a ~ 10¹², so only
  nonlinear backreaction ends the resonance [E].**
- **Caveat [A].** Suppose a φ population were thermalised with the visible bath (conformal T_vis = 3255 and 2807
  H_i). Then a²m_T² = g²(aT)²/12 ≈ 8.8e3 and 6.6e3, which is about 4e2–3e4 times λr² (apocentre to pericentre). That would remove the quartic
  resonance. Whether φ is populated at all is outside this model and is not computed.

## 3. Physical seed estimate [E]: the brief's premise needs correcting

- **The naive estimate [D].** The X-phase Bunch–Davies amplitude per ln k is (k/2π)/|χ_X(a=1)|, with
  |χ_X| = √2·76.3998 = 108.05. This gives 2.95e-3 at k = 2 and 3.14e-3 at k = 2.13.
- **What the supplied seed actually is (s4b).** It is *not* an X phase. Per radian it is:
  - a Σ rotation of 1.000;
  - an X phase change of exactly 0, because X is held real;
  - δ|X|/|X| = −0.0430;
  - a Y phase change of 0.169.
- **Its growing-mode weight.** In the late growing mode at k = 2.13, one supplied radian has the weight of only
  1/87.6 (3,.5) and 1/194 (5,1) of a genuine 1-rad local X rotation (with δp = δq). The ratio is 78–121 and
  87–221 across k = 1.8–2.8. The naive 3.1e-3 therefore corresponds to about 0.27 and 0.61 supplied radians.
- **Full six-field adiabatic vacuum at a = 1.** We used ⟨δq²⟩ = 1/2ω and ⟨δp²⟩ = ω/2 in the eigenbasis of M(0),
  where the eigenvalues are −2.83, −2.00, 15.7, 529, 1320 and 3201, and propagated with Codex's engine. The result is
  equivalent to **0.31 (3,.5) and 0.74 (5,1) supplied radians per ln k at k = 2.13**. X carries 97–98% of that power.
  This is 10–25 times the 0.03 seed, not 10 times below it.
- **Implied breakdown (supplied seeds = Codex's numbers).**

  | History | Vacuum band RMS (k = 1.8–2.8): charge = 1 | amplitude = 1 | Extrapolated unit crossing | 0.03 rad (Codex) |
  |---|---:|---:|---:|---:|
  | (3,.5) | τ 62.5 | τ 83.5 | τ ≈ 63 | 134 |
  | (5,1) | τ 57.6 | τ 80.5 | τ ≈ 52 | 147.5 |

  The 1e-6 and 1e-9 seeds lie 5–8 decades below this zero-point floor.
- **Caveats.**
  - The preparation at a = 1 is not an inflationary history.
  - For k ≲ 2 the modes are non-adiabatic (ω = 0.64–1.6 in the light directions), so the vacuum normalisation is
    uncertain by a factor of a few. That shifts the times by roughly ±15 in τ.
  - The k band below 1.8 is omitted, so the band RMS is a lower bound.
  - We use the standard classical-stochastic treatment of amplified zero-point modes (Khlebnikov–Tkachev). It is
    not a renormalised quantum calculation.

## 4. Early-era interactions, lumps and the X̄ estimate [D, E]

- **No bound lumps at φ = 0 [D/E].**
  - X is massless here, with a repulsive quartic plus a Y-induced attractive sextic. Q-balls need min 2V/ρ² < m² = 0,
    which means V < 0 somewhere.
  - With Y integrated out at tree level, that requires |x| > √(λM_L²)/κ = 203 H_i (Y_L) or 330 H_i (Y_H). At that
    point Y ≈ 2.6 M_L, so the sextic truncation has already failed there.
  - The histories have |x| ≤ 76.4 (a = 1) and ≤ 16 after τ = 10.
  - The Coulomb term is repulsive (e²), and φ is pinned at 0.
  - So this epoch has no localized bound states. Self-resonance and the k → 0 instability (μ ≈ 0.44k) produce waves
    and inhomogeneity, not Q-balls. No object charge is set before the binding transition at a ~ 10¹².
- **Circular orbit at fixed charge [D].** E_circ = ¾λ^{1/3}L^{4/3}:

  | History, orbit | E_circ | E_circ/E | Charge fraction ε | Energy above circular |
  |---|---:|---:|---:|---:|
  | (3,.5), τ10 | 1.2612e4 | 0.1197 | 0.2035 | **88.0%** |
  | (3,.5), late | — | — | — | 88.4% |
  | (5,1), τ10 | 7011 | 0.1085 | 0.1891 | **89.1%** |
  | (5,1), late | — | — | — | 89.5% |

  These percentages are the upper bound on what the coherent orbit can hand to waves at fixed charge.
- **X̄ per net charge [E].** We put that energy into X X̄ pairs near each history's own k_peak. ω ranges from
  massless up to √(k²+3λr_c²), with r_c = 64.0 and 55.3.

  | History | Excess energy into pairs | Everything a gas at ω | Equilibrium end-point |
  |---|---:|---:|---:|
  | (3,.5) | 1.36–2.62 | 1.05–2.48 | 0.096–0.107 |
  | (5,1) | 1.51–2.85 | 1.19–2.70 | 0.105–0.116 |

  The equilibrium end-point assumes the waves equilibrate before the binding transition: a circular condensate
  holding the charge plus a massless Bose gas at conformal T ≈ 19.4 and 17.2. About 10¹¹ conformal times are
  available, against a period of 0.8.
- **Comparison with kin_source [A].** The kin_source bracket is 0.13 (R5 thermal) to 2.59 (R2 adiabatic). The early
  wave count sits at the bracket's upper end. The equilibrated end-point sits just below its lower end. Which one
  applies depends on the unsolved wave evolution over the roughly 10¹¹ conformal times before binding.
- **Wavelength-cube charge [D].** At k = 2.13 it is 2.094e5 and 1.344e5, which is CONFIRMED. At (5,1)'s own peak
  (k = 1.887) it is 1.93e5, so the range at each orbit's primary peak is 1.93–2.10e5.

## 5. Wording audit (with corrected wording)

1. **REPORT title**, "Conditional breakup response precedes the prescribed heating". This overstates: nothing broke
   up. → "Conditional linear growth reaches order unity before the prescribed heating."
2. **REPORT, HANDOFF and RESULT_TABLE**, "A supplied .03rad phase reaches unit X-amplitude response by the sampled
   tau134/147.5". The numbers are confirmed. Two qualifiers are missing:
   - The seed is the Σ preparation phase, with the X phase unchanged and δ|X|/|X| = −0.043 per rad. Readers,
     including our brief, take it as an X phase.
   - Both unit samples fall at pericentre passages, where |X| = r_min ≈ r_max/8 (r/r_min = 1.007 and 1.038).

   → "The linearized response to a supplied 0.03-rad shift of the preparation (Σ) phase reaches unit fractional |X|
   modulation at a pericentre passage by sampled τ 134/147.5. The perturbation equals the orbit radius r_max by
   τ ≈ 198/216."
3. **MATHEMATICS**, "When a seeded linear radius modulation reaches1, that linear solution has certainly lost
   small-perturbation control". This slightly overstates. At that sample |δχ_X| ≈ r_min ≈ 0.12 r_max, and the
   neglected cubic terms are ≲ 7% of the k² term near pericentre. → "has reached O(1) local |X| modulation at
   pericentre; the perturbation becomes comparable to the orbit within about one further e-fold of a." The 1e-6
   equivalent is τ ≈ 916/930 by extrapolation, so ordering conclusions are unchanged.
4. **REPORT and MATHEMATICS**, "These statements concern the isolated quartic diagnostic; the actual coupled
   reheating background is not periodic and is not assigned a constant Floquet exponent" and "not constant Floquet
   exponents for the nonperiodic coupled history". These *understate* the result. → "After τ ≈ 20 the coupled history
   is radiation dominated. Its conformal X orbit conserves E to 2e-3, and the late full-system growth equals the
   matched orbit's Floquet exponent to ≤ 0.5% at every in-band k (0.633 at k = 2.13 and 0.561 at k = 2.0). The X-only
   quartic tangent reproduces the full tangent to 0.2%."
5. **MATHEMATICS**, "Exactly k=0 remains a symmetry/zero-growth limit". This is imprecise. → "At k = 0 the auxiliary
   multipliers are all unity, with Jordan blocks: there is no exponential growth, but E and L changes dephase
   linearly. In the full system the supplied seed's k → 0 limit is a physical relative X/Σ phase change, not a
   symmetry; its late response at k ≤ 0.1 is this secular growth." The value μ/k ≈ 0.440 is CONFIRMED: 0.4398 and
   0.4405, flat down to k = 0.002. By τ = 800 its growth is μη ≤ 0.6 for k ≤ 0.05 (at most about 2×).
6. **REPORT and MATHEMATICS band edges** (3.135–3.140 and 5.57–5.98). These are CONFIRMED for the (3,.5) τ10 orbit
   only. The text should say so, because (5,1) has 2.786 and 4.93–5.29.
7. **REPORT**, "That kinematic number has the order of the reference objects". This invites a link that the physics
   rules out in this era. → "…is not a fragment charge; with φ = 0 the X sector has no localized bound states, so
   this era sets no object charge."
8. **REPORT**, "A 1e-9rad seed does not reach either selected amplitude threshold before the chosen stop800. These
   are supplied seed controls…". This *understates* how strong the ordering is. Add: "An adiabatic zero-point
   estimate corresponds to ≈0.3–0.7 rad of this seed per ln k, so the 1e-6 and 1e-9 controls lie far below the
   physical floor. Any seed above ~10^−720 crosses before the Y_L median unless an unmodelled X mass intervenes."
9. **REPORT**, "the collisionless homogeneous description loses linear control" and "must therefore carry an
   inhomogeneous release". This is acceptable as a conditional. → Add "within this model (no X thermal or vacuum mass
   before a ~ 10¹²)".
10. **REPORT**, "The unperturbed light charge changes by about−0.5%…". The −0.50% endpoint value is confirmed. However,
    L oscillates by 1.6–1.8% peak to peak (exchange with Y and Σ), so the endpoint value is not a trend.
11. **The user's summary**, "A 0.03-radian seed reaches order-one modulation by τ≈134–148". → "The linearized
    response to a 0.03-rad preparation-phase seed reaches order-one…". The time ordering stated there is correct and,
    if anything, conservative.

## Limits

- Everything is linear or estimated. No nonlinear evolution, wave spectrum, φ backreaction or thermalisation is
  computed.
- The vacuum estimate depends on an ill-defined non-adiabatic initial vacuum for k ≲ 2 and on a preparation that is
  not an inflationary history.
- The X̄ numbers are energy-budget estimates, not populations.
- The s4 and s4b numbers use Codex's engine and are not independent of its equations.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; R=$SP/b73/research/reality_framework
T=research/reality_framework/grut_formation_review_20261025/tracks/st_physics; P=$R/grut_source_timing_20261026
export OMP_NUM_THREADS=2 GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025 \
  GRUT_KINETIC_WORKDIR=$SP/kwork GRUT_ONSET_WORKDIR=$SP/owork GRUT_TIMING_WORKDIR=$SP/review/st_physics/twork \
  GRUT_MODE_WORKDIR=$SP/review/st_physics/mwork
python3 -I $T/code/s1_orbits_floquet.py $P $R/grut_reheating_modes_20261026/QUARTIC_ORBIT.json $T/results/orbits_floquet.json  # 35 s
python3 -I $T/code/s2_late_growth.py $P $T/results/orbits_floquet.json $T/results/late_growth.json
python3 -I $T/code/s3_xonly_tangent.py $P $T/results/xonly_tangent.json
python3 -B $T/code/s4_vacuum_seed.py $P/code $T/results/vacuum_seed.json        # Codex engine, own seeds, 45 s
python3 -B $T/code/s4b_seed_projection.py $P/code $T/results/seed_projection.json
python3 -I $T/code/s5_estimates.py $P $T/results/late_growth.json $T/results/orbits_floquet.json $T/results/vacuum_seed.json $T/results/estimates.json
```
