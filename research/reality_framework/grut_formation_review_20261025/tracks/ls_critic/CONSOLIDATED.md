# Consolidated review: Codex b448f25, the two low-scale candidates against the debris-inclusive lump benchmark (track ls_critic)

26 October 2026 · Claude (issue #13) · consolidates tracks ls_ledger, ls_threshold and ls_scaling, and adds
cross-track recomputations (§3). Packet: b448f25 (`origin/codex/grut-low-scale-debris-20261025`,
`research/reality_framework/grut_low_scale_debris_20261025/`). It was read only through `git archive b448f25` into
the scratchpad. No Codex code was run by any track. Nothing was committed.

**Scope.** This is one supplied tree-level working model. Nothing here is evidence that X, Y, Σ, φ or the vector
exist. It is not a validation, an exclusion or a GRUT prediction. All four tracks are same-session AI work, so this is
not outside replication. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued.

Couplings [Q]: M_Y = 1e12, h = 1e11, μ = 5e11 GeV, κ = 1e-3, e = 1e-7, λ_Y = λ_Σ = 1; light target g_X = 0.1,
λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV, Q = 2e5; Y_Q = 4.8294e-9; allowance f = 0.01. Candidates: c_H = 10,
H_i = 434867798206.6179 GeV, δ = π/6 (10) and 5π/24 (17), σ = 0.001 ensembles (48 members).

## Bottom line

- **Codex's numbers reproduce.** Our own release, projection, ledger and relic code give every load-bearing
  number in the handoff table. All three handoff items are CONFIRMED.
- **The closure of both candidates under case L is CONFIRMED.**
  - Candidate 17 fails in every combination of modelling choices we computed.
  - Candidate 10 fails in every combination but one corner (§2). That corner stacks four favourable choices: the
    Y-relaxed preparation, conversion σv ×72, the on-shell s route and the zero-cascade bracket. Its lower ledger then
    sits 6% below the benchmark threshold. The cascade vectors, the damped σ count or energy conservation each remove
    that corner.
- **Three premise corrections to the benchmark.** The benchmark's arithmetic is right. Its premises have three gaps
  (ls_threshold):
  1. The hidden state is not energy-conserving after 4→2 stops.
  2. The lump energy is booked inconsistently.
  3. The resonant XX̄ → s route is missing.

  With energy conservation the case-L threshold falls about 10× (R_L ≈ 10 at candidate 10). The benchmark therefore
  favours the candidates, so their failure is stronger, not weaker.
- **One scope correction to the ledger.** The late σ and Y_H counts are not the populations that decay (ls_ledger)
  [C]. A damped model makes the Σ injection 2.6–6.3× larger [E].
- **Scaling (ls_scaling).** The one-mediator link is not generic. It belongs to the resonant, nonlinear mediator
  regime of the stated coupling ratios. The heavy-mediator (EFT) regime is a named open direction, not a viable
  point. Under the energy-conserving threshold it shrinks (§3c, §4):
  - the M_Y/H window narrows from about 12–43 to about 20–43;
  - the release-induced vector headroom drops from 1.5e-3 to about 1.2e-4 vectors per net charge.

## 1. Verdicts on Codex's handoff items

| Item | Codex (10 / 17) | Consolidated finding (10 / 17) | Verdict |
|---|---|---|---|
| **Density-dependent threshold** | R_L calibrated 91.174–91.246 / 83.115–83.181; raw 102.956–103.037 / 93.393–93.467 | Own Bessel-series and Radau code: 91.17–91.25 / 83.13–83.20 and 102.97–103.05 / 93.44–93.52. Our earlier run_inj is 1–3% lower because its endpoint is 10⁻⁶ GeV [C] | CONFIRMED |
| Final f, conservative / central (raw) | 10.8493, 12.0800 / 13.9423, 15.2658 | 10.849, 12.080 / 13.942, 15.265, unchanged at N = 800 and 1600 [C] | CONFIRMED |
| Post-tail bound | f ≥ 9.90 / 12.80 | 9.90 / 12.80 (tail depth 0.0916 / 0.0854) [C] | CONFIRMED |
| Threshold premises | the supplied equilibrium lump proxy | (a) Energy-conserving hidden state: f = 0.013 at R = 10 (candidate 10, raw), so R_L ≈ 10, not 103. (b) Lumps sit at μ_c = 0.861 m but are booked at E_c = 0.923 m, a 16% energy-work violation; with μ_c = E_c = 0.923 the proxy R_L is 99.9. (c) On-shell s route omitted; it halves f [C, in proxies] | CORRECTED (premises; the arithmetic stands) |
| **Σ radial subtraction** | N_σ = 24.9767 / 36.6526 | 24.979 / 36.625 (16 members; subset 24.979 / 36.627). The fit gives a free m_σ line of 24.97950 against the formula's 24.97946; static offset 3.7e-5; Y-forced lines 1e-8; the charge-carrying angular part is correctly excluded [C] | CONFIRMED |
| Use of the late N_σ as the decaying population | "late projected numbers conserved to decay" | σ decays at ⟨a⟩ = 15.6, but N_σ settles only at a ≳ 67–75 [C]. The damped model has 2.6–6.3× (10) and 3.0–6.3× (17) more σ decays, giving R_Σ = 227–543 / 379–797 [E] | CORRECTED (scope; magnitude [E]) |
| **Disjoint decay ledger** | R_Y 113.069 / 122.806; R_Σ 86.335 / 126.884; R_lower 199.404 / 249.690; R_central 228.417 / 280.531 | 113.07 / 122.80; 86.34 / 126.79 (sampling); 199.41 / 249.59; 228.42 / 280.43 [C]. Widths agree to ≤ 1e-12. ⟨a⟩_L = 43299.7 and ⟨a⟩_A = 4.1434e11. Chain and time-dilation corrections are ≤ 4e-5. Nothing is counted twice [D, C]. The closed form reproduces 113.07 / 86.34 / 29.00 (ls_scaling) | CONFIRMED |
| Cascade bracket [0, N_Y] | "uncertain, bracketed" | Not a strict upper bound, because Y_H is refilled by H–L mixing before the counts settle. This affects the central and upper ledgers only [C] | CONFIRMED lower end; upper end CORRECTED (scope) |
| Ensemble ε, N_Y | 0.1574128 / 0.1470462 | 0.1574113 ± 7e-6 / 0.1470609 ± 5e-5 (16 members; same-member deviation ≤ 1.2e-7) [C] | CONFIRMED |
| Primordial vectors | R ≈ 0.011 | Quoted, not rechecked (sp_critic confirmed the evolving-mass method) [Q] | NOT CHECKED |

## 2. Is the closure under case L robust to the benchmark's modelling choices?

All entries are for the conservative (zero-cascade) ledger unless marked. The first three threshold columns are
the benchmark's own isentropic proxy; the last two are ls_threshold's energy-conserving (ke) state. Except where
marked [E], all entries are [C].

| Ledger (R) | Proxy | Proxy, σv ×72 | Proxy, ×72 + s | ke | ke + s |
|---|---|---|---|---|---|
| 10, Codex preparation (199.4) | R_L 103.0; f 10.85 | R_L 148.2 | — | f 6.45–7.84 | f 3.73 |
| 17, Codex preparation (249.7) | R_L 93.4; f 13.94 | R_L 134.7 | — | f 8.41–10.15 | f 4.97 |
| 10, Y relaxed (148.86) | R_L 110.0–110.3; f 6.58 | R_L 158.0–158.4 | **R_L 158.9–159.4; f 8e-7 (passes)** | — | ke + ×72 + s: f 2.64–3.11 |
| 10, Y relaxed, central (178.95) | — | — | f 2.26 | — | — |
| 17, Y relaxed (192.69) | — | — | R_L 149.5–150.0; f 5.23 | — | — |
| Damped Σ, Codex preparation (333–649 / 491–909) [E] | above every proxy threshold here | | | | |

**Reading.**

- **Candidate 17** fails in every cell.
- **Candidate 10** fails in every cell except one. That corner needs all six of the following:
  - the Y-relaxed preparation (not Codex's);
  - conversion σv ×72, ls_threshold's hot-daughter sensitivity value;
  - the on-shell s route;
  - the zero-cascade bracket;
  - the late-count Σ;
  - the isentropic proxy.

  Changing any one of the last three fails it: the central ledger gives f = 2.26, the damped Σ count [E] fails, and
  energy conservation gives f = 2.6–3.1.
- **Single choices.** No single choice rescues either candidate. The verdict is robust. The size of the margin is
  not: it runs from about 1600× down to about 260× in the computed cells, and to zero only in the corner.
- **Y-only control.** The candidate-10 Y-only failure is thin: 113.07 against 103.0, which is 1.10×, or 1.03× with
  backreaction [E]. It is not load-bearing.

## 3. Cross-track recomputations (this track)

- **(a) Y-relaxed preparation, re-checked [C].**
  - Method: our own analytic stationarity solve, with the sp_yledger yrelease.c release, ls_project projection and
    our own closed-form ledger.
  - Minimum: x = 75.33750, Y = (−5.97966, −0.22171), |Σ| = 2.65000. The gradient is 1e-13 and the Hessian is
    positive. This agrees with ls_scaling to about 1e-8.
  - Release results:

    | Candidate | ε | R_lower | ls_scaling R_lower | R_central |
    |---|---:|---:|---:|---:|
    | 10 | 0.16487 | 148.86 | 148.92 | 178.95 |
    | 17 | 0.15783 | 192.69 | 192.76 | 224.61 |

  - Raw yields are 5.970 / 5.701.
  - Note: arg X carries a torque of −284 / −341 (H units). The phase remains supplied.
- **(b) Our own R_L tables were biased high [C].**
  - ls_scaling's R_L values came from the sp_critic grid, log-interpolated across brackets where f changes by up to 8
    decades. A direct proxy bisection gives:

    | Point | Direct bisection | Table |
    |---|---:|---:|
    | Heavy point (ε = 0.1457, normalised yield) | 81.9–82.2 | 87.9 |
    | Candidate 10 (calibrated) | 91.2 | 99.5 |

  - The tables therefore run 7–9% high. This also corrects our sp_critic "R_L ≈ 20–175", which ke supersedes in any
    case.
- **(c) Energy-conserving threshold at ls_scaling's heavy point [C, ke proxy].** Same light sector, ε = 0.1457,
  normalised yield.
  - Crossings:
    - μ_c = E_c = 0.861: R_L,ke ≈ 8.2. f = 0.0057 at R = 8 and 0.045 at 8.9; the N = 1000 control gives 0.041.
    - μ_c = 0.923: ≈ 7.0 (f = 3.7e-4 at R = 5, 0.25 at 8.9).
    - With the s route: between 15 (f = 5e-14) and 20 (f = 0.093).
  - Base point, R = 1.14: erased in every variant (f ≤ 1e-196; the root-solver fallbacks sit at the p floor).
  - m = 20 (R = 8.9): fails ke at 4× the allowance and passes ke + s.
  - m = 15 (R = 34.2): fails both (f = 1.36 / 0.51).
- **(d) Vector headroom at the base point [C, E].** Headroom = (R_L − 1.1445)/(19.642 × 2788.2) vectors per net
  charge:

  | Threshold | Headroom |
  |---|---:|
  | Table (ls_scaling) | 1.58e-3 |
  | Direct proxy | 1.48e-3 |
  | ke | 1.1–1.3e-4 |
  | ke + s | 2.5–3.4e-4 |

## 4. What ls_scaling implies for the one-mediator link

- **Within the benchmark assumptions the link is not generic [D, C].** Those assumptions are: homogeneous classical
  instant release, a supplied phase, the late-decay ledger and the case-L proxy.
  - At the stated ratios (m_L/ω_X = 0.72, m_σ/ω_X = 0.105, λ_Y|Y|²/M_Y² = 6.8) the charge forms resonantly, and the
    same resonance makes real Y and Σ quanta. That is why ε and debris rise together, and why no single-coupling move
    (×0.1–×4) passes.
  - The scaling symmetry (M_Y, h, μ, H) → s(·) leaves ε and every R unchanged. So the problem lives in the ratios.
- **Named open direction: heavy, adiabatic mediators.** The charge then comes from virtual Y exchange (α6).
  Conditions:
  - r = √2μh/M_Y² ≳ 0.3;
  - κ ∝ M_Y;
  - h ≳ 20H;
  - a Y-relaxed preparation, which decides the outcome there: Codex's preparation gives R ≈ 1414.
- **The window.** It is M_Y/H from about 12 (proxy) or about 20–21 (ke) up to 43 λ_Y^{1/4} (tree stability). The ke
  edge is interpolated in m and is [E]. Spot point: H = 1.334e11, M_Y = 4.00e12, h = 5.33e12, μ = 9.55e11 GeV,
  κ = 3.46e-3. It gives ε = −0.146, Y = 1.00 Y_Q, R = 1.14 and c6 = 1.25e6.
- **What "passes" means here.** The window passes the debris-against-threshold test only. It is not a cosmological
  point.

## 5. Overreach flags

1. **ls_scaling, "passes by a factor of about 80" and "a passing window".**
   - The 80 is measured against an interpolated proxy threshold that runs 7% high. Under ke the base-point margin is
     about 7× and the window narrows.
   - Five gates are uncomputed: the vector, phase holding, strings, κ running and the case-L premise.
   - Read the result as an open direction.
2. **ls_scaling, "Conservative R_lower CORRECTED".** Codex's "conservative" refers to the cascade bracket within its
   stated preparation. A Y-relaxed release is a changed history, which Codex scoped out. Relabel: CONFIRMED for the
   stated preparation. It is not a lower bound across preparations (§2 uses it as a robustness probe).
3. **ls_ledger, "fail by a wider margin than Codex reports".** This rests on a phenomenological friction model [E].
   What is demonstrated [C] is only that the late counts are not the decaying populations.
4. **ls_threshold, "R_L ≈ 10".**
   - Its basis is a five-point scan at one density, in an equilibrium proxy with instant condensation, and runs near
     zero hit about 200 root fallbacks.
   - (c) reproduces the same edge (about 7–8) at another point and grid.
   - It is a proxy threshold, not a kinetic one. "×370–1600" holds for Codex's preparation only (see §2).
5. **Codex, "Even Y alone fails".** True, with a 1.10× margin for candidate 10. The end point matters there: our
   earlier code gives f = 2.18 and Codex 0.57.

## 6. Status of the empirical branch

**At the stated couplings no published release creates and protects.** Three closures stand:

- the 373d37e island (f ≈ 15, sp_critic);
- candidate 10;
- candidate 17.

All three closures hold only within the debris-inclusive benchmark (homogeneous classical releases, late vacuum
decays, equilibrium lump proxy). None is a kinetic exclusion of the action.

**The only named open direction is ls_scaling's heavy-mediator EFT window.** It has passed only the debris test.

**Not touched by any of this:**

- the vacuum density remains underived;
- candidate 1 is paused;
- the phase is supplied throughout.

**Single most decisive next step: release-induced longitudinal-vector (Goldstone) production at the heavy spot
point.** Compute it from the actual release, and compare it with the headroom of 1.2e-4 (ke) to 1.5e-3 (proxy)
vectors per net charge. Earlier estimates at other points were 1e-3 to 0.1 [Q, REVIEW_SCANPOINT]. If similar values
hold here, the direction closes at e = 1e-7. A larger e cuts each vector's deposit as e^{−1/2}.

## 7. Open, not certified

- **Early decays.** A kinetic treatment of σ → A_L A_L and Y_H → Y_L A_L at a ≈ 5–70, including backreaction
  (ε −4%, N_Y −6 to −9% [E]).
- **Case-L premise.** A kinetic or Boltzmann replacement for both the proxy and ke; the two differ by about 10× in
  R_L. Also missing:
  - the D relic, which can only add;
  - the s route beyond Maxwell–Boltzmann statistics.
- **Primordial vectors.** R ≈ 0.011, quoted only.
- **Heavy direction.** Phase holding (relative-phase mass about 1.84 H), strings (Gμ ≈ 9e-12 if they form), κ running
  (0.70 of the tree bound), and survival of instant release.

## Reproduction (this track, about 45 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/ls_critic/code
OMP_NUM_THREADS=1 python3 -B $T/heavy_threshold.py $T/../results/heavy_threshold.json
OMP_NUM_THREADS=1 python3 -B $T/ke_refine.py       $T/../results/ke_refine.json
python3 -B $T/relaxed_check.py                     $T/../results/relaxed_check.json
OMP_NUM_THREADS=1 python3 -B $T/stack_check.py     $T/../results/stack_check.json
OMP_NUM_THREADS=1 python3 -B $T/stack_ke.py        $T/../results/stack_ke.json
```

These import our ls_threshold `lsrelic.py`, sp_yledger `ymodel.py`/`yrelease.c` and ls_ledger `ls_project.py`
unchanged. The binaries are built in the scratchpad.

Credits:
- Q-ball formation and thermodynamics: Kusenko–Shaposhnikov; Enqvist–McDonald; Laine–Shaposhnikov.
- Cannibal freeze-out: Carlson–Machacek–Hall.
- Bose integrals: Kolb & Turner.
- Asymmetric-relic annihilation: Graesser–Shoemaker–Vecchi.
- Mixed and decaying dark-matter structure constraints are the context for f = 0.01.
