# Skeptic pass on track st_clocks (Codex b73de64, source timing), issue #13

27 October 2026 · Claude (st_clocks_skeptic) · adversarial AI review in the same session, not outside replication.

**Task.** Try to refute every st_clocks finding. For each one I recompute the key number by a method different
from st_clocks' own, and I look for unit, sign, convention, sampling or indexing errors and for overstated verdicts. I
also list new issues with Codex's packet that st_clocks missed.

**Provenance.**
- All code is in `skeptic/code/` and is own code. No Codex or st_clocks code was imported or run.
- Codex files (`SOURCE_CLOCKS.json`, `RESPONSE_RESULTS.json`, `quiet*.npz`) and st_clocks' JSON and npy files were read
  as data only.
- Two analyses below are re-readings of **Codex's saved linear solutions**: the crossing criterion and the 1e-9
  extrapolation. They are labelled that way and are not new dynamics.
- Labels: [D] derived, [C] computed, [E] estimate, [A] argued, [Q] quoted.

## Result in one line

No st_clocks number is refuted. Every value was recomputed by a different method and agrees.

The verdict on T1 ("linear control is lost well before substantial prescribed decay heating": CONFIRMED) is
**overstated in scope**. st_clocks confirmed only the clock side. The crossing side depends on Codex's choice of
criterion, which st_clocks did not examine (§3). The ordering itself survives every variant I tried.

## 1. Methods (different from st_clocks)

| Item | st_clocks method | Skeptic method (`skeptic/code/`) |
|---|---|---|
| Widths | hand-written complex-basis amplitude, (κ/√2)·3! per channel; dblquad of the Dalitz area | `widths_skeptic.py`: sympy expansion of the potential in real canonical fields. Numeric Y mass-matrix eigenvectors; every real-field vertex from 3rd and 4th derivatives along the eigen-directions; sum over all real final states with their own symmetry factors. Phase space by recursive factorisation (scipy quad). Two Σ phases (−π/2 as in the engine preparation, and 0). |
| Vector width | Σ\|ε·(p₁−p₂)\|²/3 in closed form | explicit polarisation vectors, Minkowski contraction, 40×64 angular quadrature in the A rest frame; m_A from \|DΣ\|² with \|Σ\|² = h²/2 |
| Background C, τ*, a(τ) | rpreheat.c: full 19-variable system, own DP5 to τ = 800 | `background_skeptic.py`: inflaton + radiation only (the hidden sector enters Friedmann at η² = 3e-15, about 1e-11 [A]), scipy DOP853 at rtol 1e-13 (control 1e-11). After τ₁ = 160 (inflaton/radiation < 2e-34) the radiation era is exact and closed form. |
| F_E | cubic spline of the integrand on a 0.05 grid | the energy integrals are ODE quadrature states (no interpolation); tail by mpmath gammainc at 50 digits, checked by mpmath quadrature in u = γ(s − t) |
| PCHIP attribution | asserted | tested directly: Codex's own a(τ) data with PCHIP, and the same data with a cubic spline |
| Condensate-assisted [E] | rp data, integrated from τ = 1 | coefficient re-derived by real-field vertex sum (`widths_skeptic.py`); integral on **Codex's** \|X\| data and on rp data; lower limit 0, 1, 5 (`ordering_skeptic.py`) |
| N1 populations | rp engine diagnostic | own free-quanta decomposition applied to **Codex's** saved fields (`population_skeptic.py`) |

## 2. Per-claim results

| ID | st_clocks verdict | Codex | st_clocks | Skeptic (own) | Skeptic call |
|---|---|---|---|---|---|
| W1 Γ_YL | CONFIRMED | 13413.681870567 GeV | same | **13413.6818705671** (real-basis vertices 3κ, 3κ for Y_L → X1X1X2 and X2X2X2; phase-space ratio 1 − 1e-16; rel −1.1e-16) | not refuted |
| W2 Γ_A | CONFIRMED | 2.12233e-10 / 2.12233e-7 GeV | same | **2.1223251086668e-10 / 2.1223251086668e-7** (rel −1.0e-14 / −2e-16) | not refuted |
| W3 1000× | CONFIRMED | 1000 | 1000.0000000000 | 1000.0000000000153 | not refuted |
| M1 medians | CONFIRMED | 6.890794e6, 4.355172e17, 4.355172e20 | same | same (≤ 1e-14) | not refuted |
| C1 C | CONFIRMED | 1.9858052758 / 1.4758677865 | 1.9858052754 / 1.4758677875 | **1.9858052753537 / 1.4758677874894** (vs st_clocks −1.35e-11 / −1.19e-11, the expected hidden-sector size; vs Codex −2.3e-10 / +6.9e-10); τ* 0.1936739888 / −0.1690615009; C within 1e-6 after τ 22.45 / 11.4 | not refuted |
| A1 a at medians | CONFIRMED | 5231.40 / 4509.97 … | same | 5231.400162 / 4509.966852 (rel −1.13e-10 / +3.44e-10); vectors same | not refuted |
| F1 number fraction | CONFIRMED | 1.34790e-5 / 1.48370e-5 | same | 1.3479013e-5 / 1.4836963e-5 | not refuted |
| F2 energy fraction | CONFIRMED, offset from PCHIP | 3.716111e-8 / 4.306670e-8 | 3.716112e-8 / 4.306671e-8 | **3.71611179e-8 / 4.30667095e-8** (vs st_clocks 9e-11 / 2e-10; vs Codex +2.28e-7 / +1.94e-7) | not refuted; **attribution verified** (see below) |
| F3 tail, late limit | CONFIRMED | formulas | 1.9e-8 / 1.7e-8 | closed form = quadrature to 1e-50; limit rel +1.95e-8 / −1.70e-8 [D: u = γ(s − τ*)] | not refuted |
| Q1 Q_cube, units | CONFIRMED (geometric) | 2.0937e5 / 1.3440e5 | same | same. N_X = a³n/H_i³ is confirmed from the engine source (n = −2Im(X̄Ẋ)a³, physical complex-canonical X in H_i units, τ = H_i t) | not refuted (note: see §4.4) |
| T1 ordering | CONFIRMED, qualified to X-sector heat | statement | clocks only | clock side reproduced; crossing side is criterion-dependent (§3); the 1e-9 seed is unresolved; the [E] number depends on its lower limit | **not refuted, verdict overstated in scope** |
| N1 populations | NEW (diagnostic) | — | N_YL 0.00682 / 0.0178, R_Y 0.0080 / 0.029 | on Codex's fields: N_YL 0.0068156 / 0.0177683, N_YH 9.69e-6 / 4.33e-6, E_X 1.05202e5 / 6.4619e4, R_Y 0.008038 / 0.029377 | not refuted |

**F2 attribution test [C].** Codex's own a(τ) samples with PCHIP, integrated by my 12-point Gauss rule on each
interval, give 3.71611094e-8 and 4.30667012e-8, equal to Codex's numbers to 1e-10. The same samples with a not-a-knot
cubic spline give 3.71611179e-8 (own +1.7e-9). The 2e-7 offset is therefore PCHIP interpolation on Codex's grid, as
st_clocks said; the grid's largest step is 3.217 at τ ≈ 77. Codex's "orders 4 and 8 agree to 3.4e-16" does not test
this.

**Cascade and σ [D, C].**
- The Y_H Y_L G vertex is exactly √2μ: 10.125 against 10.125. This gives γ_casc = 0.0350420, Γ = H at τ = 14.5 / 14.1
  (a = 7.538 / 6.490), and 99.09% / 99.43% complete by τ = 134 / 147.5.
- σ → AA: 1.12540. σ → Y_LY_L: 0.0111306. These agree.
- Not in the dictionary: Y_H → 3X has γ = 1.633e-7. Its branching against the cascade is 4.7e-6, so the X heat from
  the cascade epoch is ≲ 1e-9 of E_X [E]. "No X heat" holds to that level.

**Condensate-assisted Y_L → XX [E].**
- The coefficient c₂ = 9κ²/(8πm_L) is exact in the real-field sum, independent of the condensate phase to 7e-16.
- Integral to the unit sample:

  | Lower limit | (3, 0.5): Codex data | (3, 0.5): rp data | (5, 1) at τ = 147.5 |
  |---|---|---|---|
  | τ = 1 | 1.121e-3 | 1.121e-3 | 1.173e-3 |
  | τ = 0 | 1.764e-3 | — | 1.679e-3 |
  | τ = 5 | 0.767e-3 | 0.779e-3 | 0.803e-3 |

- 52–56% of the integral accrues at τ < 5, before the homogeneous Y_L count has settled (N_YL settles by τ ≈ 3–5).
- So "1.1e-3, about 80× the three-body" is really **0.8–1.8e-3, 57–131×**. It is still ≪ 1.
- Γ₂/γ₃ at the unit sample is 18.3 / 17.1.

## 3. New issue: the crossing times depend on the criterion

This bears on T1 and on the quoted "134–148" and "772–782".

**What Codex's criterion measures.** The amplitude criterion is the radial projection divided by the **instantaneous**
radius, (χ_X·δχ_X)/|χ_X|² (`timing_common.responses`). The homogeneous X orbit is very eccentric:
|X|max/|X|min = 7.7–8.1 on (3, 0.5) and 8.2–8.7 on (5, 1) within one radial period. The criterion therefore fires at
pericentre passages. The sampled per-radian response at k = 2.13 is 49.5 at τ = 134 but 34.1 at τ = 200 on (3, 0.5).

**What the field perturbation actually is [C].** Re-reading Codex's saved tangent vectors (k = 2.13), I normalised the
full |δχ_X| by the local orbit maximum, using windows of ±0.9a, ±1.5a and ±3a in τ; the results are identical.

- **At Codex's own unit samples** the field-level perturbation is only:
  - 0.03 rad: |δχ_X|/max|χ_X| = **0.19 / 0.15** at τ = 134 / 147.5;
  - 1e-6 rad: **0.18 / 0.12** at τ = 772 / 782.5.
- **With the criterion |δχ_X| = orbit amplitude**:
  - the 0.03 rad crossing moves to **τ = 197.5 / 215.5** (a ≈ 28 / 25);
  - the 1e-6 seed reaches only **0.34 / 0.19 by τ = 800**. Extrapolating the stable late growth (0.32 / 0.37 per unit
    a) gives about τ ≈ 900 / 965 [E].
- **Why this is the more relevant measure [A].** For the cubic force, second-order terms relative to the linear ones
  scale as |δχ|/A, with A the orbit amplitude. The local-charge criterion δn/n_X is amplified by about 1/ε relative to
  that (ε ≈ 0.15, the orbit charge fraction), which is why it fires earlier still.
- **What follows.**
  - MATHEMATICS.md's "When a seeded linear radius modulation reaches 1, that linear solution has certainly lost
    small-perturbation control" is overstated. What is certain is that the linearised *radius* at pericentre fails; the
    field perturbation is then 12–20% of the orbit amplitude.
  - The quoted crossing times are criterion-dependent at the tens-of-percent level: +64–68 τ for 0.03 rad, and
    ≳ +120–180 τ [E] for 1e-6 rad.
- **Ordering unchanged.** The first 1% of prescribed X heat (Y_L → 3X) comes at τ ≈ 9.99e4, which is still ≥ 100×
  later under every criterion.

## 4. Other new issues st_clocks missed

1. **The 1e-9 seed is unresolved.** It is a supplied seed, but it crosses no threshold by τ = 800: Codex's own envelope
   times 1e-9 is 1.58e-3 / 1.44e-3.
   - The headline "for the supplied phase seeds … loses linear control long before …" is computed for two of the three
     seeds.
   - Extrapolating the stable late growth of Codex's envelope gives τ ≈ 1.46e3 (Codex criterion) or about 1.67–1.75e3
     (orbit criterion) [E]. That is still ≪ 1e5, but it is not computed.
2. **Codex's C carries a run-to-run spread of about 1e-9.** The quiet and gauge runs integrate the same homogeneous
   background but give C = 1.9858052758 against 1.9858052791 (1.7e-9) and 1.4758677865 against 1.4758677874 (6e-10).
   That is larger than Codex's quoted "late ranges" of 4e-11 and 1.3e-10. It is immaterial at the four quoted digits;
   st_clocks' "tighter" comparison is right, but Codex's C accuracy is about 1e-9, not about 1e-10.
3. **The k that carries the crossings.** On (5, 1) both unit crossings are carried by k = 2.0, not 2.13. At k = 2.0,
   Q_cube = 1.62e5 instead of 1.34e5. Also on (5, 1), the 0.03 rad 10% warning comes from the infrared k = 0.02 mode
   (τ 53–55, k/aH ≈ 0.17). That mode's running-maximum response per radian grows only from 1.8 to 14 between τ = 20 and 800, not like a band mode. "At the
   primary peak" therefore names the auxiliary-orbit peak, not the mode that crossed.
4. **Scope of st_clocks' T1 verdict.** The verdict line says CONFIRMED, but st_clocks lists the crossing times as not
   certified. It should read "CONFIRMED for the clock side; the crossing side is not checked." See §3 for why the
   crossing side matters.

## 5. Limits

- The background here omits the hidden sector (an effect of about 1e-11 [A]) and checks C only to that level.
- The criterion analysis and the 1e-9 extrapolation reuse Codex's saved linear solutions. They do not recompute the
  mode equations, the direct-vector adaptation or the phase-current sign. Those belong to other tracks.
- Widths are tree-level vacuum widths. The condensate-assisted channel is an order-of-magnitude estimate.
- Nothing here validates or excludes the model. This is a supplied working model with supplied fields.

## Reproduction

```bash
cd research/reality_framework/grut_formation_review_20261025/tracks/st_clocks/skeptic/code
python3 -B widths_skeptic.py       # ~2 s  -> ../results/widths_skeptic.json
python3 -B background_skeptic.py   # ~5 s  -> ../results/background_skeptic.json (reads Codex JSON/npz as data)
python3 -B ordering_skeptic.py     # <1 s  -> ../results/ordering_skeptic.json (needs st_clocks' rp npy in $SP/review/st_clocks/rp/runs)
python3 -B population_skeptic.py   # <1 s  -> ../results/population_skeptic.json
python3 -B criterion_skeptic.py    # ~1 s  -> ../results/criterion_skeptic.json
python3 -B verdicts_skeptic.py     #       -> ../results/skeptic_verdicts.json
```
