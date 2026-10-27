# Review of Codex b73de64: conditional growth versus source clocks

27 October 2026. Claude, issue #13. Brief: [BRIEF_TIMING.md](BRIEF_TIMING.md). Four review tracks, each followed by an
adversarial skeptic, then a completeness critic: `tracks/st_rerun`, `st_response`, `st_clocks`, `st_physics`,
`st_critic`. Each track's NOTES.md and SKEPTIC.md give methods, commands and full numbers.

This is same-session AI work on one supplied model. It is not outside replication. Nothing here is observational
evidence, a validation, an exclusion or a GRUT prediction.

## Result

**Every Codex number we checked reproduces, and no skeptic refuted any of them.** That covers the crossing times,
the clocks and the auxiliary quartic spectrum. The conditional time ordering is correct as a linear statement, and the
review makes it stronger. Most numbers were checked with our own code:
- The linear response was re-derived in four formulations that share no code with response.c:
  - st_response: our reheating background with our own Wirtinger Hessian;
  - a cosmic-time second-order form with dual-number Hessian-vector products;
  - a sympy-Hessian physical-variable form;
  - a numerical-Jacobian form.
- They match Codex's published arrays to 3×10⁻⁹ – 4×10⁻⁸ of the response envelope. A 5×10⁻⁹ – 10⁻⁸ floor is shared by
  all codes and is immaterial.
- The clocks and widths come from our own formulas, our own phase-space quadrature and our own reheating engine
  extended to τ = 800.

Several statements are corrected or sharpened:
- one sampling artefact;
- the scope of the crossing criterion;
- what the "phase" seed is, and where the vacuum sits relative to it;
- the neutral channel past τ = 80.

## Verdicts

| Item | Codex | Ours | Verdict |
|---|---|---|---|
| verify_response.py, seven groups | direct vector vs transfer 2.66×10⁻⁹; Friedmann 3.8×10⁻⁹; Gauss 7×10⁻¹⁰ | identical; 3.73×10⁻⁹; 6.94×10⁻¹⁰ (reproduction; all four npz bitwise identical) | CONFIRMED |
| response.c is the full generator applied to one column, with the Ward projection kept | stated | diff with modes.c: only the stride, seed input and removed neutral branch change. M g = ward to 3.9×10⁻¹⁶. DOP853 error control covers the tangent (relative once grown); tolerance and seed-rescaling tests move the response ≤ 10⁻⁷ | CONFIRMED |
| Phase-current sign | old archive = −new | least-squares ratio −1.0000000017 | CONFIRMED |
| Local charge with canonical p is gauge covariant | stated | equals the X part of the Gauss-law current to 4×10⁻¹⁵ | CONFIRMED |
| 0.03 rad, unit δ\|X\|/\|X\| | by τ 134 / 147.5 | first continuous crossing 133.60 / 147.32 (three independent codes) | CONFIRMED, see criterion note |
| 10⁻⁶ rad, unit | by 772 / 782 | 771.94 / 781.81 | CONFIRMED (linear; see neutral channel) |
| Local charge = 1, 0.03 rad | by 120 / 134.5 | 119.58 / 134.41 | CONFIRMED |
| 10⁻⁹ rad | none by τ 800 | none (largest 1.6×10⁻³); extrapolated τ ≈ 1490 [E] | CONFIRMED |
| 10% warning, (3, 0.5), 0.03 rad | [82.5, 83.0] | **56.46** (k = 2.13), in four codes | **CORRECTED** (sampling; see below) |
| C = a²H/H_i; τ\* | 1.986 / 1.476 | 1.98580528 / 1.47586779; 0.193674 / −0.169062; a² = 2C(τ − τ\*) fits to 10⁻¹¹ | CONFIRMED |
| Γ_YL = 3κ²m_L/256π³; Γ_A = e²m_Aβ³/48π with m_A = 6ehH_i | as used | re-derived (Y_L real eigenmode, m_L = 22.25 H_i; one charge-1 pair channel). The 1000× median ratio is Γ_A ∝ e³ at fixed h | CONFIRMED [D] |
| Medians τ and a | 6.89×10⁶ (a 5231 / 4510); 4.36×10¹⁷; 4.36×10²⁰ | same to ≤ 5×10⁻¹⁵ and 3.5×10⁻¹⁰ | CONFIRMED |
| Y_L fractions at the unit samples | number 1.35–1.48×10⁻⁵; energy 3.72–4.31×10⁻⁸ | 1.3479 / 1.4837×10⁻⁵; 3.7161 / 4.3067×10⁻⁸. The 2×10⁻⁷ offset is from PCHIP interpolation, not quadrature | CONFIRMED |
| Auxiliary quartic spectrum | μ/k 0.440; edge 3.135–3.140; band 5.57–5.98 | 0.4398; 3.1395; 5.571–5.976; peak 0.63178 at 2.128 | CONFIRMED, (3, 0.5) orbit only. The (5, 1) orbit: peak 0.5631 at 1.887, edge 2.786, band 4.93–5.29 |
| "k = 0 is a symmetry / zero-growth limit" | stated | there is no exponential growth, but the monodromy has two Jordan blocks (secular dephasing). For the full system the supplied k → 0 seed is a physical relative phase, not a symmetry | CORRECTED (wording) |
| Q_cube = \|N_X\|(2π/k)³ at k = 2.13 | 1.34–2.09×10⁵ | 1.3440 / 2.0937×10⁵ exactly (N_X is a count per comoving H_i⁻³) | CONFIRMED. On (5, 1) the crossings come from k = 2.0, where Q_cube = 1.62×10⁵ |
| Light-charge change −0.5% to τ 800 | stated | −0.497% / −0.504%, a level offset from the τ 64–80 window (late oscillation 0.3% peak to peak) | CONFIRMED |

## Corrections

1. **The pre-τ 80 grid is too sparse.** The inherited sample spacing is 1.5–3.2 τ between τ 12 and 80. The response
   peaks at pericentre passages that last about 0.3 τ.
   - So the (3, 0.5) 10% bracket [82.5, 83.0] is not a bracket. The first crossing is τ 56.46, and the same holds for
     gauge3.
   - The error carries into SOURCE_CLOCKS.json (evaluated at τ 83). The Y_L number fraction there is 8.35×10⁻⁶,
     against 5.68×10⁻⁶ at τ 56.46.
   - It also carries into TIMING.png, whose pre-τ 80 curve under-reads the response 2.3×.
   - The same artefact affects 2cbc7f1's "amplitude through τ 80" for (3, 0.5): 8.88% published, 20.2% on a dense
     grid. We did not check that number in bca4927.
   - All unit crossings after τ 80 are unaffected, because the spacing there is 0.5 τ.
   - The correction strengthens the ordering.
2. **The unit criterion fires at pericentre.** The X orbit is eccentric (\|X\|max/\|X\|min ≈ 8). At the unit samples,
   linear theory has failed locally: \|δχ_X\|/\|χ_X\| = 1.26–1.52, and the neglected terms are 23–39% of the linear
   force. Yet the perturbation is only 0.12–0.20 of the orbit amplitude. It reaches the orbit scale at τ 197.5 / 215.5
   for 0.03 rad.
   - So "certainly lost small-perturbation control" should name the pericentre criterion.
   - The quoted times depend on the criterion at the tens-of-percent level.
   - The ordering does not depend on the criterion.
3. **The supplied "phase" seed is a Σ preparation-phase rotation with X held real** (δθ_X = 0; δln\|X\| = −0.043 per
   rad). In the growing mode, one supplied radian weighs as much as 1/88 (3, 0.5) or 1/194 (5, 1) of a 1-rad X rotation.
   It should be called a preparation-phase seed, not an X phase seed.
4. **The late growth is quantitatively the quartic self-resonance.** Codex's wording, "not assigned a constant Floquet
   exponent", understates this.
   - In radiation, a² = 2C(τ − τ\*) makes conformal time η = a/C. The crossing pairs then imply μ_eff = 0.6335 / 0.5599
     per unit η.
   - Isolated quartic orbits matched to each history give 0.6318–0.6331 / 0.5631–0.5645.
   - Strobed late growth of Codex's saved tangents matches to about 0.1% (charge) and about 1% (amplitude).
5. **Clock-side scope.** The sentence "Before substantial prescribed decay heating" is right for X-sector heat. One
   prescribed decay, Y_H → Y_L + A_L, is 99% complete by τ 134–148, but it puts nothing into X. The first 1% of the
   X-sector heat (Y_L → 3X) arrives at τ ≈ 1.0×10⁵. A condensate-assisted Y_L → XX channel (not in the dictionary)
   would decay 0.8–1.8×10⁻³ of Y_L by τ 134 [E], still ≪ 1. To reverse the ordering, in-medium enhancement would have
   to be 10⁴–7×10⁴.
6. **Minor points.**
   - The e = 10⁻⁶ versus 10⁻⁷ comparison is nearly a null test, because the gauge term is about 10⁻⁶ of the response.
   - The τ 64–80 mean uses 6 samples, a 10⁻⁵ effect.
   - C differs run to run by about 10⁻⁹ (quiet versus gauge runs of the same background).
   - Drop the comparison of Q_cube with the reference objects (next section).

## Two new physics points (they strengthen the ordering)

**The vacuum floor sits well above the small supplied seeds.** The k ≈ 2 modes are inside the horizon at a = 1, so
their initial state is close to the adiabatic vacuum. From the six-field adiabatic vacuum we estimate an equivalent
seed of 0.24–0.74 supplied radians per ln k [E]. That is 8–25× the 0.03 seed, under four reasonable covariances; the
normalisation is uncertain by a factor of order one. Its band RMS reaches unit charge and amplitude modulation at
τ ≈ 57–69 (a ≈ 13–17), with dense sampling. The 10⁻⁶ and 10⁻⁹ controls lie 5–8 decades below this floor.
In radiation the growth is exponential in a. For the breakdown to come after the Y_L median (a ≈ 5×10³), the seed
would have to be below about 10⁻⁷²⁰. So the ordering is not a matter of seed choice.

**The neutral field φ is driven parametrically and catches up with X.** The declared coupling g_X²φ²\|X\|² gives δφ a
Floquet band on the eccentric orbit: μ_φ = 0.202 at k ≈ 7.1 on (3, 0.5) and 0.185 at k ≈ 6.2 on (5, 1). We extended
2cbc7f1's own neutral control from τ 80 to τ 800, with Codex's engine, seeds and energy definition. This is
reproduction-type, and it matches the published τ ≤ 80 values to 3×10⁻⁹. The φ/X energy ratio reaches:
- 1% at τ 202.5 / 195;
- 10% at τ 427.5 / 409.5;
- 1 at τ 697 / 665.

So in this model the fixed-background approximation fails through φ before the 10⁻⁶ crossings (772 / 782). Those
crossings are linear extrapolations on a background that is no longer valid. The ordering against heating is
unaffected.

**There are no bound lumps in this era [D/E].** With φ = 0, X has a repulsive quartic self-interaction. The attractive
sextic from Y exchange is 0.6% of it at τ 10 and 10⁻⁸ by a ≈ 5×10³. The potential would turn negative only for
\|x\| > 203 H_i, against at most 76 at a = 1. The wavelength-cube charge is therefore kinematic only; this epoch
selects no object charge.

**The expected end state of the early release [E].** The charge goes into a rotating condensate. At fixed charge, the
circular orbit holds only 12.0% / 10.9% of the X energy, so up to 88–89% can pass to waves. Turned into X/X̄ pairs near
the peak k, that gives about 1.4–2.9 X̄ per net charge. The fully equilibrated massless end point gives 0.10–0.12.
Our kin_source bracket of 0.13–2.6 sits between these limits.

## What this changes

Codex's conclusion stands and is strengthened. A homogeneous release cannot be carried into the heating era, at any
physically sensible seed. The successor nonlinear calculation should:
- start from vacuum-seeded X, Y, Σ **and δφ** fluctuations by τ ≲ 50;
- sample at ≤ 0.05 τ, because pericentre spikes last 0.3–1 τ;
- keep the charge and energy ledger;
- export the wave spectrum (expected: a charged rotating condensate plus a symmetric X/X̄ wave population carrying up
  to about 88% of the X energy) into the later heating and binding stages.

Its X̄ per net charge would replace the five labelled alternatives in our kinetic source bracket.

The following remain open:
- the nonlinear history, thermalisation and the φ population;
- the binding transition at a ≈ 10¹²–10¹³;
- f_sym and the kinetic threshold;
- the selected state, the vacuum density and any distinctive signal.

Candidate 1 is unchanged.
