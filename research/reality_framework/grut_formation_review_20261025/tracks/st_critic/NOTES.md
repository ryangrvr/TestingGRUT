# st_critic: completeness critic for the review of Codex b73de64 (conditional growth versus source clocks)

27 October 2026. Claude (track st_critic). Inputs: the structured findings and skeptic checks of st_rerun,
st_response, st_clocks and st_physics (their NOTES.md, SKEPTIC.md and results/), Codex's packet
`grut_source_timing_20261026/` at b73de64 (read from `$SP/b73`, not modified) and, for one check, Codex's earlier
`grut_reheating_modes_20261026/` (2cbc7f1) code in the same tree. Rules: [BRIEF.md](../../BRIEF.md),
[BRIEF_TIMING.md](../../BRIEF_TIMING.md). Everything is conditional on the supplied working model. Nothing here
validates or excludes anything, and nothing was committed.

**My own computations (two quick checks).**

| File | What | Independence |
|---|---|---|
| `code/settle.py` → `results/settle.json` | Late growth rate and orbit-scale crossings from Codex's published npz arrays, under the different window/normalisation choices the tracks used | Own code; Codex arrays read as data |
| `code/neutral_extend.py` → `results/neutral_extend.json` | Codex's own 2cbc7f1 neutral (δφ) control, extended from τ 80 to τ 800 on the two quiet histories | **Reproduction-type**: Codex's `modes.c -DNEUTRAL_SECTOR`, seeds and energy-ratio definition; only the time range and bookkeeping are ours. It reproduces Codex's published τ ≤ 80 maxima (0.0017169 / 0.0016724) to 3e-9 relative |

## 1. Coverage of the items Codex asked to be checked

| Codex request | Who checked it | Status |
|---|---|---|
| MATHEMATICS.md | Hamiltonian, A0 elimination, seed consistency (p0 = q0) and δ\|X\|/\|X\| re-implemented in four different formulations (st_response, three skeptics); charge (st_rerun A3); clock formulas and tail (st_rerun, st_clocks, skeptic); Q_cube units (st_clocks); wording (st_physics) | Covered. One sentence was not checked by anyone: "The neutral binding field is not backreacted in this extension". Codex's 2cbc7f1 neutral control stops at τ 80; the timing claims run to τ 800. Filled here (§3) |
| verify_response.py | st_rerun only, by reproduction (all 7 groups pass; numbers identical). The skeptic hash-checked but did not rerun | Covered as reproduction only. These are same-author checks; the independent evidence for the response arrays is the agreement of four other formulations (≤ 3e-8 of the envelope at in-band k, ≤ 5e-7 at low k where the e-term differs) |
| Direct vector versus full transfer | st_rerun A1 (code diff, Ward identity 3.9e-16, generator 6e-16, 2.6613e-9 reproduced); independent codes agree | CONFIRMED |
| Phase-current sign | st_rerun A3, st_response R7 (independent), skeptic (−1.0000000017) | CONFIRMED |
| Matched radiation clocks | st_rerun R9, st_clocks C1 (independent engine), skeptic (third background) | CONFIRMED. Codex's C has a run-to-run spread of about 1e-9 (quiet vs gauge runs), larger than its quoted late ranges; harmless |
| Statement of the conditional time ordering | st_physics P16, st_clocks T1, st_rerun | CONFIRMED with qualifiers (§4). st_clocks' T1 covers only the clock side |

The 42 histories: all reproduced bitwise (st_rerun); all 18 k on both quiet histories checked by an independent code
(st_response skeptic); gauge5 brackets by st_rerun's skeptic (e = 0, which the e-comparison shows is adequate).

## 2. Contradictions between tracks, and how they resolve

1. **Late growth rate.** st_rerun N2 gives 0.45–0.58 per unit conformal time (its skeptic "confirmed" it with the same
   construction). st_physics gives 0.6335 / 0.5599, matching the matched-orbit Floquet exponent. **Settled [C]:**
   st_rerun takes the running maximum in 20-τ windows, but the late radial period is 24 τ (τ ≈ 200) to 45 τ
   (τ ≈ 780), so many windows contain no pericentre spike. With one-radial-period windows on the same published
   arrays, the max-over-k amplitude slope is **0.633 (3,.5) and 0.563 (5,1)** in every sub-range (150–400, 400–800,
   200–790), against 0.48–0.57 with 20-τ windows. The charge response, which is not spike-dominated, gives 0.63 / 0.56
   either way. So st_rerun N2 is CORRECTED, and its [E] extrapolation for the 1e-9 seed (τ 1704 / 1429) is superseded
   by st_physics' τ ≈ 1494 / 1487 [E]. Its skeptic's confirmation was stronger than its evidence (same method).
2. **Orbit-scale crossing.** st_rerun's skeptic: 196.3 / 234.0, and the 1e-6 seed reaches 0.41 / 0.34 by τ 800.
   st_clocks' skeptic and st_physics: 197.5 / 215.5, with 0.34 / 0.19. **Settled [C]: a matter of definition.** The
   first uses the radial part δ\|X\| over max\|X\| within ±8 τ. Late in the run that window covers only about a third
   of a radial period, and it under-reads r_max by a factor of 0.49–0.51 at τ > 700, which inflates the late ratios
   about 2×. On the published arrays: full \|δχ_X\| / r_max (radial period) gives **197.5 / 215.5** (k = 2.13), and
   radial only / r_max gives 219 / 234.5. The 1e-6 seed reaches 0.34 (3,.5) and 0.11 (5,1, k = 2.13) by τ 800. The
   extrapolated orbit-scale crossing for 1e-6 is τ ≈ 900–965 [E].
3. **"Certainly lost small-perturbation control" (MATHEMATICS.md).** st_physics: slightly overstated (cubic terms
   ≲ 7%). Its skeptic: defensible (including the quadratic remainder, the nonlinear terms are 39% / 23% of the linear
   force at the sample). st_rerun's and st_clocks' skeptics: overstated relative to the orbit scale. **Resolution:**
   these measure different things, and all are right on their own terms. Locally, at the pericentre passage, linear
   theory has failed (\|δχ_X\|/\|χ_X\| = 1.26–1.52; nonlinear remainder 23–39%). The orbit-scale perturbation is
   0.12–0.20 of r_max there and reaches 1 about 60–70 τ later. The sentence stands if it names the pericentre
   criterion. st_physics' "≤ 7%" omitted the leading quadratic term and is withdrawn.
4. **Q_cube.** st_physics said CORRECTED; the others said CONFIRMED. → **CONFIRMED at the stated k = 2.13**, with a
   note: on (5,1) the unit crossings are carried by k = 2.0 (1.62e5), and the (5,1) orbit's own peak is k = 1.887
   (1.93e5).
5. **Vacuum floor (st_physics P9).** "10–25× the 0.03 seed" becomes **8–25×** over four covariances. "Band-RMS
   amplitude = 1 at τ 83.5 / 80.5" is the same pre-τ80 sampling artefact as the 10% warning; with dense sampling it is
   **τ 69.3 / 66.9** (charge 61.2 / 57.5). Skeptic values adopted.
6. **Light charge −0.5% (P15).** st_physics' "inside a 1.6–1.8% oscillation" is wrong for the late epoch (late
   oscillation 0.31–0.35% p-p; the −0.5% is a level offset of the τ 64–80 reference). Codex's sentence is accurate.
7. **Condensate-assisted Y_L → XX [E].** 1.1e-3 (st_clocks) versus 0.77–1.76e-3 (skeptic, lower-limit dependence).
   Use the range; ≪ 1 either way.
8. **Minor, no conflict.** 10% warning 56.46 / 56.458 / 56.457 / 56.5 across four codes; fine-k shift of the
   local-charge crossings 0.28 versus 0.31 τ (use 0.31); st_rerun's τ 40–80 spacing is 1.49–3.22, not 1.89–3.22.

## 3. A skeptic's partial refutation that was left unresolved: the neutral (δφ) channel

st_physics P8 said nothing in the declared model ends the quartic resonance before a ~ 10¹². Its skeptic found an
omitted in-model channel: the declared g_X²φ²\|X\|² (g_X = 0.1) gives δφ a parametric band with μ_φ ≈ 0.21 / 0.18 at
k ≈ 7 / 6.2, and estimated O(1) backreaction "by τ ~ 400" [E]. Two corrections to that framing, then the result:

- Codex did not ignore δφ. 2cbc7f1 integrated 816 neutral histories (136 k) through τ 80 with conditional vacuum seeds
  and found φ energy ≤ 0.17–0.19% of the X conformal energy. b73 removed that branch from response.c (exact at linear
  order, because φ = 0 decouples it) and did not extend the control.
- **Extension to τ 800 with Codex's own engine and definitions [C, reproduction-type]:**

| History | φ/X energy ratio 1% | 10% | 1 (Codex's "linear approximation invalid") | Late μ_φ (per η), peak k |
|---|---:|---:|---:|---:|
| (3,.5) | τ 202.5 | τ 427.5 | **τ 696.5–697.0** | 0.202, k 7.07 |
| (5,1) | τ 194.5–195.0 | τ 409.5 | **τ 664.5–665.0** | 0.185, k 6.19 |

The μ_φ values independently match the skeptic's Floquet result (0.206 / 0.185). The skeptic's "O(1) by τ ~ 400" was
early by about 270 τ; at τ 400 the ratio is 5–9%. **Consequence:** in Codex's own conditional vacuum treatment, the
prescribed-background approximation fails through the neutral sector at τ ≈ 665–697. That is after the 0.03-rad
crossings but **before the 1e-6 crossings (772 / 782)**, and long before any 1e-9 crossing. The 1e-6 and 1e-9 numbers
are therefore linear extrapolations on a background already invalid in this model. The ordering against heating is
unaffected and, if anything, strengthened.

## 4. Verdicts stronger than their evidence

- st_rerun N2 / N3 and its skeptic's confirmation: window artefact (§2.1).
- st_rerun's skeptic, orbit-scale 196.3 / 234.0 and 0.41 / 0.34: window and normalisation artefact (§2.2).
- st_clocks T1 "CONFIRMED": clock side only. The crossing side rests on the other tracks.
- st_physics P8 (§3), P13 (§2.4), P15 (§2.6), wording item 3 (§2.3) and the P9 sub-number (§2.5).
- st_response R8: "the residual is Codex's integration noise" is not established. A floor of 5e-9 to 1e-8 is shared
  by three codes. R11: its "~10%" caveat is too pessimistic; the gauge term is tested to about 1e-4 of its size in
  the charge channel.
- st_rerun's skeptic, "upper bound narrower than stated": Codex already writes that in-medium widths are not
  computed. The 7e4× / 1.3e4× enhancement needed is useful, but it is not a missed caveat.
- R7, the 10% CORRECTED: correct as to value. But Codex's PROTOCOL asks for sampled brackets, and MATHEMATICS warns
  that excursions can fall between samples. The defect is that [82.5, 83.0] was then propagated as the warning time
  into SOURCE_CLOCKS.json (8.35e-6 → 5.68e-6) and TIMING.png (pre-80 peak drawn at 0.089 versus 0.202), and the
  inherited 2cbc7f1 table (8.88% → 20.2%).

## 5. The single most important point for Codex

**The supplied-seed linear X timing is not the controlling clock.** At vacuum level, the same model leaves the
homogeneous linear regime earlier through two channels:

- **X itself.** The six-field adiabatic vacuum at a = 1 is equivalent to 0.24–0.74 supplied rad per ln k, which is
  8–25× the 0.03 seed. Its band RMS reaches unit charge and amplitude modulation at τ ≈ 57–69 (dense sampling).
- **The binding field φ.** In Codex's own neutral control extended to τ 800, the φ energy reaches the X energy at
  τ ≈ 665–697, before the 1e-6 crossings.

The supplied 1e-6 and 1e-9 seeds lie 5–8 decades below the zero-point floor. Codex's next nonlinear release
calculation should therefore:

- start from vacuum-seeded X/Y/Σ **and δφ** fluctuations by τ ≲ 50;
- sample at ≤ 0.05 τ, because pericentre spikes last about 0.3–1 τ.

Both points already appear in Codex's own handoff ("neutral field/bias and quantum matching"). This makes them
quantitative. The standard treatment of amplified zero-point modes as classical fields is Khlebnikov–Tkachev;
condensate fragmentation follows Kusenko–Shaposhnikov, Enqvist–McDonald and Kasuya–Kawasaki. Caveats: the vacuum
normalisation at a = 1 is a supplied choice, uncertain by about 1.3× across covariances and more for non-adiabatic
k ≲ 2; the φ seeds are Codex's conditional Gaussian seeds; and nothing nonlinear is computed.

## 6. Consolidated verdict table

| Item | Verdict | Ours | Codex | Method |
|---|---|---|---|---|
| 7 verification groups | CONFIRMED | 2.6613e-9; Friedmann 3.73e-9; Gauss 6.94e-10; tol 1.789e-8 | 2.6613e-9; 3.8e-9; 7e-10; 1.8e-8 | reproduction |
| Direct vector = full transfer | CONFIRMED | Ward 3.9e-16; 4 formulations ≤ 3e-8 of envelope in band (≤ 5e-7 low k) | stated | diff + independent codes |
| Phase-current sign | CONFIRMED | old-archive ratio −0.9999999986 / −1.0000000017 | stated conversion | independent |
| 0.03 rad, unit \|X\| | CONFIRMED (pericentre criterion) | 133.60 / 147.32; orbit-scale 197.5 / 215.5 | by 134 / 147.5 | 3 independent codes + dense |
| 1e-6 rad, unit \|X\| | CONFIRMED (linear; after φ failure) | 771.94 / 781.81 | 772 / 782 | independent |
| Local charge = 1 | CONFIRMED | 119.58 / 134.41; 737.42 / 751.07 (fine-k ≤ 0.31 earlier) | 120 / 134.5 | independent |
| 1e-9 none by 800 | CONFIRMED | max 1.6e-3 / 1.4e-3; [E] ~1494 / 1487 | none | independent + [E] |
| 10% warning, (3,.5), 0.03 rad | **CORRECTED** | 56.46 (k 2.13); gauge3 same | [82.5, 83.0] (k 2.4) | 4 codes, dense |
| Other 10% warnings | CONFIRMED | 55.17; 563.47 / 576.81 | in brackets | independent |
| C, τ* | CONFIRMED | 1.98580528 / 1.47586779; 0.193674 / −0.169062 | 1.9858052758 / 1.4758677865 | independent |
| Widths, medians | CONFIRMED | 6.890794e6; 4.355172e17 / e20; a 5231 / 4510 | same | independent |
| Y_L fractions at unit samples | CONFIRMED | 1.3479 / 1.4837e-5; 3.7161 / 4.3067e-8 | same (PCHIP 2e-7) | independent |
| Auxiliary quartic spectrum | CONFIRMED ((3,.5) orbit only) | 0.4398; 3.1395; 5.571–5.976; 0.63178 @ 2.128 | 0.440; 3.135–3.140; 5.57–5.98 | independent |
| k = 0 "symmetry" wording | CORRECTED (wording) | Jordan blocks; the seed's k → 0 limit is a relative X/Σ phase | symmetry limit | independent |
| Q_cube at k = 2.13 | CONFIRMED + note | 2.0937e5 / 1.3440e5 ((5,1) k 2.0: 1.62e5) | 1.34–2.09e5 | independent |
| Late growth = quartic resonance | NEW | 0.633 / 0.563 per η; X-only tangent to 0.17% | not claimed | own analysis |
| Neutral φ beyond τ 80 | NEW | φ/X energy 1 at τ 697 / 665; μ_φ 0.202 / 0.185 | ≤ 0.17% through τ 80 (2cbc7f1) | reproduction-type extension |
| Vacuum floor | NEW [C/E] | 8–25× the 0.03 seed; band RMS = 1 at τ 57–69 | supplied controls | independent + Codex engine |
| Conditional ordering | CONFIRMED with qualifiers | first 1% of X-sector heat at τ 9.99e4 versus every breakdown ≤ 782 | stated | derivation |

Qualifiers on the ordering: it concerns X-sector heat (the Y_H → Y_L + A_L cascade is 99% done by τ 134 but heats no
X); the unit criterion fires at pericentre; the 1e-9 seed is extrapolated only; the neutral channel fails before the
1e-6 crossing; and in-medium or stimulated widths are not computed (they would need ≥ 1.3e4× enhancement to matter).

## Limits

- Only two quick checks were run here. All other numbers are the tracks' and skeptics'.
- The φ extension uses Codex's engine and conditional seeds (reproduction-type), with no backreaction. Its unit ratio
  marks failure of the prescribed background, not a computed φ population.
- The time grid after τ 80 is 0.5 τ. The φ energy ratio is smooth there, so its first-sampled crossings are accurate
  to about 0.5 τ.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; R=$SP/b73/research/reality_framework
T=research/reality_framework/grut_formation_review_20261025/tracks/st_critic
export OMP_NUM_THREADS=2 GRUT_MODE_WORKDIR=$SP/review/st_critic/mwork \
  GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025
python3 -I $T/code/settle.py $R/grut_source_timing_20261026 $T/results/settle.json            # < 10 s
python3 -B $T/code/neutral_extend.py $R $T/results/neutral_extend.json                        # ~ 6 s
```
