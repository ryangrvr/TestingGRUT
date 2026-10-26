# Track ls_scaling: how the debris failure scales with the couplings (scope question after Codex b448f25)

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication.

Packet: b448f25 (`origin/codex/grut-low-scale-debris-20261025`), read only through `git archive` into the scratchpad.
No Codex code was run. Codex JSON was used only as quoted input [Q] for one ledger check. Nothing was committed.

**Scope.** One supplied tree-level working model. Benchmark assumptions are unchanged throughout:
- a homogeneous classical release, instant at a = 1 with instant radiation;
- a supplied CP phase;
- the late-particle decay ledger;
- the case-L proxy with the 1% allowance.

Nothing here is evidence that X, Y, Σ or the vector exist. It is not a validation, an exclusion or a viable point.

Labels: [D] derived, [C] computed, [E] estimate, [Q] quoted, [A] argued.

## Question

At the stated couplings every published release fails, because its own Y and Σ debris R exceeds R_L. Is that failure
structural for this field content, or is there a direction in coupling space where:
- R < R_L;
- ε and the yield stay adequate;
- c6 stays below 1.486e13?

## Answer

**The failure is not structural, under the stated benchmark assumptions.** It belongs to the stated coupling *ratios*.

At the stated couplings, charge forms in a resonant, nonlinear mediator regime that also produces real Y and Σ
quanta. In the opposite regime the mediators are adiabatic, and the charge comes from virtual Y exchange (the α6
operator). Real quanta are then suppressed by powers of ω_X/m. Our full three-field releases find a window there that
passes by a factor of about 80. It needs a coordinated move, a different (physical) preparation, and gates that are
still open (§5).

## 1. Exact relations [D] (numbers checked [C])

**Scaling symmetry.** Scale (M_Y, h, μ, H_inf) → s·(M_Y, h, μ, H_inf), with κ, λ's, e, c_H and δ fixed. Then:
- the release in H units is unchanged, so ε, N_Y, N_σ, E_X and every R are unchanged;
- yield ∝ s^{3/2};
- c6 ∝ s^{−2}.

The check at s = 0.3066 reproduces ε and R to machine precision. The yield ratio is 0.16981 = s^{1.5} and the c6 ratio
is 10.635 = s^{−2}. **The debris problem therefore lives entirely in dimensionless ratios. H_inf only sets the yield
and c6.**

**Closed-form ledger.** These reproduce Codex's 113.07, 86.34 and 29.00 at candidate 10:
- R_Y = 64.47 κ⁻¹ √(m_L/H) · N_Y/E_X
- R_Σ = 24π √(h/(eH)) · N_σ/E_X
- R_casc = 12π √(h/(eH)) · N_YH/E_X

R_Y is linear in κ at fixed forcing, because N_Y ∝ κ². The vector deposits scale as e^{−1/2} h^{1/2}.

**Charge, EFT regime.**
- ε ≈ (7–10) η, where η = α6 X0²/λ_X and α6 = rκ²/(2M_Y²(1−r²)), with r = √2 μh/M_Y² [C, using our gp_critic
  companion code].
- The companion caps η ≤ r/12.
- At the stated r = 0.0707 this gives ε ≲ 0.07–0.10. **ε ≳ 0.15 needs r ≳ 0.3**, i.e. μh ≳ 0.2 M_Y².

**Yield.** Y ∝ ε c_H^{3/2} H^{3/2}/λ_X.

**Survival.** c6 = 2ηλ_X² M_red²/(c_H H²).
- It is fixed by ε and H alone, and is independent of the mediator direction.
- Values: 1.18e5 at H_i, 1.25e6 at the yield-normalised H. It binds only for H ≲ 4e7 GeV.
- Stated point: 2.11e5.

**Vacuum stability (tree level, large field).** κ < (λ_X/1.19)^{3/4} λ_Y^{1/4} = 4.93e-3 λ_Y^{1/4}.

**Heavy-limit debris from the instant release** [D, magnitude E]:
- N_Y ≈ 9 c_H² κ² X0⁶/m_L⁷. At fixed ε (κ ∝ M_Y) this gives N_Y ∝ M_Y^{−5} and R_Y ∝ M_Y^{−5.5}.
- Measured slopes: N_Y 4.2→4.8 and R_Y 4.9→5.3. The magnitude is high by a factor of 2.1–2.6 for m ≥ 20.
- N_σ falls faster than the h^{−9} power law once m_σ ≫ ω_X (resonant-to-adiabatic crossover).

## 2. Why the stated point fails [C]

| Ratio at the candidates | Value |
|---|---|
| m_L/ω_X | 0.72 (resonant) |
| m_σ/ω_X | 0.105 (light radial Σ, displaced by μY² during inflation) |
| λ_Y \|Y\|²/M_Y² | 6.8 (nonlinear Y, no EFT) |
| κX0²/M_Y² | 1.07 |

**Preparation.** Codex leaves Y at its δ = 0 value while rotating Σ. Letting Y relax to its conditional minimum at
the same phases (which a heavy Y does during inflation) gives:

| Candidate | Preparation | ε | R_lower | R_central | R_L (ours) |
|---|---|---:|---:|---:|---:|
| 10 | Codex | 0.1574 | 199.5 | 228.4 | 99.5 |
| 10 | Y relaxed | 0.1648 | 148.9 | 179.0 | 107.4 |
| 17 | Codex | 0.1471 | 249.5 | 280.3 | 89.2 |
| 17 | Y relaxed | 0.1579 | 192.8 | 224.6 | 99.9 |

Both still fail. For δ > π/2 the relaxed |Σ| collapses to zero, because the μY² forcing (s ≈ 2.75) dominates the vev
(0.16).

**One-at-a-time moves (×0.1–×4, Codex preparation, 6 phases each).** None passes.

| Move | Effect |
|---|---|
| M_Y ×2, ×4 | ε ≤ 0.006 |
| κ ×0.5 | ε ≤ 0.013 |
| μ ×0.5 | ε ≤ 0.0045 (the wrong direction) |
| λ_Y ×10 | ε ≤ 0.005 |
| h ×3 with μ ×1/3 | ε ≤ 0.006 |
| κ ×3 (best case) | ε 0.237, R 379 vs R_L 170 |
| h ×3 | ε 0.135, R 1042 vs 79 |
| μ ×2 | R 2954 vs 106 |
| λ_Y ×0.1 | R ≥ 8700 |

In this regime ε and debris rise together. That is the "structural" link, but it holds only while the mediators are
resonant.

## 3. Heavy-mediator direction [C]

Setup: c_H = 10, Y-relaxed preparation, δ = π/2, a_end = 40 (a control at 100 agrees to 4 digits), and R_L from our
normalised table.

| Point | ε | R_central | R_L | Yield/Y_Q |
|---|---:|---:|---:|---:|
| m = M_Y/H = 10 | −0.113 | 145 | 60 | 4.2 (fails) |
| m = 15 | −0.134 | 34 | 77 | 5.2 |
| m = 20 | −0.142 | 8.9 | 85 | 5.7 |
| **m = 30** | **−0.146** | **1.14** | **88** | 5.9 |
| m = 40 | −0.145 | 0.25 | 88 | 5.9 |

The m scan holds r = 0.45, h = 40H and η = 0.5 η_edge, so κ = 1.153e-4 m.

**Base point at the yield-normalised scale (from the symmetry).**
- H = 1.334e11 GeV, M_Y = 4.00e12, h = 5.33e12 and μ = 9.55e11 GeV.
- κ = 3.46e-3, which is 0.70 of the stability bound; λ_Y = λ_Σ = 1, e = 1e-7.
- Results: ε = −0.146, yield 1.00 Y_Q, R = 1.14 against R_L = 88, c6 = 1.25e6.

**Robustness around the base point.**
- Phases 15°–165°: R ≤ 1.5. The best ε is −0.154 (R_L 96).
- Edge fraction 0.25–0.9: ε 0.074–0.29, R ≤ 1.4, always below R_L.
- r = 0.2 gives ε = 0.075. r = 0.0707 gives ε = 0.030, which is too small.
- h: 20H gives R_Σ = 0.36; 10H gives R_Σ = 38 (passes with a 2.8× margin); ≤ 3H is an uncontrolled release (μ ≥ 95H
  makes Y_L tachyonic).
- λ_Y = 10 passes. λ_Y = 0.1 passes the release but violates boundedness.

**Window.** It runs from m ≈ 13 (debris) to m ≲ 43 λ_Y^{1/4} (stability, at ε ≈ 0.15, r = 0.45). That is a factor of
about 3 in M_Y/H.

**Preparation control.** Codex's preparation at the same base point gives N_YH = 77 and R_central ≈ 1414. In the
heavy regime the preparation is decisive.

## 4. Verdicts

| Claim | Verdict |
|---|---|
| Codex deposits at candidate 10: R_Y 113.07, R_Σ 86.34, R_casc 29.00 | CONFIRMED (closed form; our release gives 113.07 / 86.40 / 28.92) |
| No create-and-protect point among the published releases at these couplings | CONFIRMED, including with Y relaxed and under ×0.1–×4 moves |
| R_lower 199.4 / 249.7 "conservative" | CORRECTED: it depends on the preparation (148.9 / 192.8 with Y relaxed); still fails |
| The failure is generic for this field content | REFUTED within the benchmark assumptions |
| A smaller μ helps | REFUTED: r must increase |
| Survival constrains the direction | REFUTED: c6 is fixed by ε and H; it is 1e7 below the ceiling |

## 5. Not certified (needed before the window means anything)

1. **Release-induced longitudinal-vector/Goldstone production.** Not computed by Codex or here.
   - One vector per comoving H³ costs 19.6 E_X at h = 40H.
   - The headroom allows only about 1.6e-3 vectors per net charge.
   - Earlier estimates elsewhere were 1e-3–0.1 per net charge [Q, REVIEW_SCANPOINT]. **This is the decisive next
     gate.** A larger e cuts the deposit as e^{−1/2}.
2. **Phase holding.** A supplied phase is still needed: the relative-phase mass is about 1.84 H during inflation [D].
3. **Strings.** Gμ ∝ h²: about 9e-12 at h = 5.3e12 GeV, if strings form [E, scaled from 3.1e-15].
4. **Model assumptions.** Instant release and reheating, the homogeneous classical treatment, early-decay backreaction,
   and the case-L premise itself.
5. **Running.** κ is at 0.7 of its tree-level bound, and running was not checked.

## Reproduction (~40 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/ls_scaling
python3 -B $T/code/scan.py $T/results/scan_B.json B
python3 -B $T/code/scan.py $T/results/scan_A.json A
python3 -B $T/code/prep_compare.py $T/results/prep_compare.json
python3 -B $T/code/analytic.py $T/results
```

**Code provenance.**
- `scalrel.c` is our sp_yledger `yrelease.c` with its couplings made runtime arguments. It is validated at candidate
  10: ε = 0.15741241, N_YH = 16.779, N_YL = 22.184, N_σ = 24.978, the same as ls_ledger.
- R_L comes from our sp_critic tables. They are about 9% above Codex's calibrated values at the same ε; this is
  immaterial at the margins found here.
