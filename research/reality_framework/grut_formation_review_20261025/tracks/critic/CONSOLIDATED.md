# Consolidated review: Codex packets ff49e5d, 9636950, bb211d1 (capture and retention) and 59a1a71

Claude review, critic track, 25 October 2026. This consolidates five track records: `reproduce`, `independent_solver`,
`spectrum`, `operators` and `conversion` (their NOTES.md files are the record). `xcharge` covers a different packet
(8986d34) and `affleck_dine` has no NOTES.md yet; neither is consolidated here. Read-only on Codex's branches;
nothing committed. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued.

Every number below is one model with supplied fields and one prepared spherical condensate. Nothing here is
evidence that the fields exist, a validation, an exclusion or a cosmological source.

**Critic recomputations** (`code/critic_checks.py`, `results/critic_checks.json`; own analysis code, reading the
independent_solver raw histories, the spectrum and conversion result files and Codex's JSON via `git archive`):

- **K1.** The old-damper deficit, Q60(old damper) − Q60(box), from the independent_solver raw runs at h = 0.5 and
  dt = 1/64, in percentage points of Q0:

  | τ | 1400 | 2000 | 2500 | 3000 | 4000 | 6000 | 9600 | 10800 | 12000 |
  |---|---|---|---|---|---|---|---|---|---|
  | Deficit (pp) | +0.0005 | +0.037 | −0.588 | −0.415 | −0.309 | −0.225 | −0.280 | −0.194 | −0.059 |

  The 9600–12000 mean is −0.188 pp. Under a 1e-6 amplitude perturbation the same window means move by at most
  0.027 pp.
- **K2.** Trajectory sensitivity of the escape through r = 240 by τ = 12000, under the same 1e-6 amplitude
  perturbation: ΔQ = +13.6 and ΔE = +11.7.
- **K3.** Five converged runs from three independent codes. Retention spans 35.3565–35.3631%, so 35.360 ± 0.003
  pp. Escaping charge through r = 240 by τ = 12000 is 23519 ± 6, and escaping energy 25482 ± 4.
- **K4.** The late-window offset, worked by three routes, gives 34.735, 34.740 and 34.741% [E].
- **K5.** Codex's 32-node rate table, interpolated with PCHIP in log space, overestimates the dense-scan integral
  by 1.431× over p/m_X = 0.45–0.56 and by 1.097× over the full band. The full-band ratio matches the optical-depth
  ratio 54.61/49.74 = 1.098.
- **K6.** In the free_soft packet, 0.300 of the energy lies below k = 0.0163 and 0.079 below k = 0.01.

## 1. Verdict table (all handoff items)

### ff49e5d: causal receiver

| Item | Codex | Consolidated value | Verdict |
|---|---|---|---|
| 1 CFS realization and zero memory; the C++ code implements it | stated | Re-derived; resolvent agrees to 8.6e-14 (operators [D,C]) | CONFIRMED |
| 1 Maximum reflection, k ≥ 0.01 | 1.82e-8 charged / 6.516e-4 neutral | 1.8222e-8 / 6.5161e-4 from three methods (operators, reproduce, independent_solver analytic) | CONFIRMED, for k ≥ 0.01 only |
| 1 Linear stability | none growing | max Re λ from −2.3e-8 to −1.9e-7 | CONFIRMED (linear only) |
| 1 Nonlinear matching, empirical | untested theorem | Absorber-free runs agree with the CFS runs to ≤0.007 pp in retention and ≤0.05% in escape (K3) | CONFIRMED empirically for these observables |
| 2 Summation-by-parts pair, 4th order, lead roots | stated | Exact; observed order 3.998; |z| ≥ 25.96 | CONFIRMED (O(1) error where U_rr ≠ 0 at the wall; no reported value changes) |
| 3 Packet errors 1.1e-12, 4.32e-15, 3.78e-13, 6.88e-6 | stated | Bitwise, by reproduction only | CONFIRMED numbers; scope CORRECTED: free_soft tests k ≳ 0.016 only (K6) |
| 4 Retention 35.36%; two finest grids within 0.019% | 35.3631% | 35.360 ± 0.003 pp numerically (K3); ±0.03 pp under 1e-6 perturbations | CONFIRMED |
| 4 Time-step control, layer move, coarse-step failures | +0.0305%; 4.144e-7 | identical (reproduction) | CONFIRMED |
| 4 Late coarse run, τ 19200–24000 | 34.574% | 34.74% (K4) [E] | CORRECTED |
| 5 Exact discrete currents | stated | Exact for the semi-discrete scheme [D]; closure is RK4 error of order ≈4.9 | CONFIRMED |
| 5 I025 closure 3.3e-7 / 6.28e-8; bitwise identical to P025_time | stated | identical (reproduction only) | CONFIRMED |
| 5 I025 escaping E and Q through r ≈ 240 | 25618.6 / 23680.9 | 25482 / 23519 (K3). Codex's values hold for its h = 0.25 trajectory | CORRECTED as physical values (+0.54% / +0.69%); mQ/E 0.9230, not 0.9244 |
| 6 Snapshot bound-mode numbers | 11182.655, ω 0.864632 | reproduced by a re-implementation | CONFIRMED |
| 6 "Time/resolution dependence" of the split | REPORT table | Split oscillates over 11063–13954 with period ≈35, which is the core breathing 2π/0.181 measured independently by spectrum. So 11183 vs 13324 is a phase difference, not a resolution effect [E]. On excited cores the split errs by 20–36% | CORRECTED (interpretation) |
| REPORT: mass shift 0.00174 leaves outgoing energy "uncontrolled" | stated | The maximum is non-resonant neutral radiation; the slow part is 3.9e-5 at r = 240; 97.1% of the X number passes | CORRECTED (interpretation) |
| REPORT: the 35.15 → 35.36 shift came from a boundary "that failed its soft-wave test" | stated | The deficit is set at τ ≈ 1800–2500 (K1). Only p ≳ 0.4 waves can return by then; the old layer reflects 4.9–14.8% of them. Soft returns (p ≈ 0.05) shrink the deficit inside the window | CORRECTED |
| Gate: certified particle spectrum | open | N_X̄ = 255.5 ± 0.6 (0.656% of Q0); N_s 109–110, E_s ≈ 354; X spectrum at p/m_X ≥ 0.03 to ±1.6% | PARTLY CLOSED (leading-order classical) |
| 9/9 checks; 48 arrays | stated | 9/9; 48/48 | CONFIRMED (reproduction of author checks) |

### 9636950: soft receiver

| Item | Codex | Consolidated | Verdict |
|---|---|---|---|
| Old damper reflects 57.1 / 31.7 / 16.7 / 3.41 / 14.8% | stated | 31.73% at p = 0.02, from two independent codes | CONFIRMED |
| Antiparticle estimate is absorber-sensitive (29%) | 28.95% | Bitwise; but +127% with a Hann window and ×10.6 with a late window | CONFIRMED as computed; not a stable measure |
| Conditional antiparticle fraction | 1.13–1.16% | 1.06–1.07%, absorber-free, two codes | CORRECTED |
| 25.1% passes the local cut | stated | Codex's cut was not recomputed; with the slow shift 97.1% passes | Reinterpreted [A+C] |
| N0125_box and N0125_far bitwise identical | stated | trace SHA values match | CONFIRMED |

### bb211d1: capture channel and retention

| Item | Codex | Consolidated | Verdict |
|---|---|---|---|
| 1 Three-channel matrix, √2, charge convention, π(2ℓ+1)/p², closed s, unitarity, reciprocity | stated | Re-derived; unitarity 3e-10 | CONFIRMED |
| 1 σv, cold and at e = √1.5 | 1.13039e-16 / 8.12250e-15 | 1.13038e-16 / 8.12250e-15 | CONFIRMED |
| Table between its nodes | PCHIP | ℓ = 5 resonance at p = 0.508 (7.69e-14, FWHM 0.012); integral overestimated 1.43× (K5) | CORRECTED |
| 2 Portal spectral kernel | stated | — | NOT CHECKED |
| 3 Transport in ν | stated | Same table: +0.18%. Corrected table: forward events −3% (z = 10⁶), −8% (1 keV), −13% (20 MeV) | CONFIRMED; diagnostics CORRECTED |
| 3 Late anti-X tail | 1.08e-6 | ≲2.6e-7, order 1e-7 (grid diffusion) | CORRECTED |
| 3 Frozen profile | assumed | Relaxation raises the release by 0.17% and u_esc by up to 0.25% | CONFIRMED, bounded |
| 4 Retention 35.15%, spread 0.069% | stated | 35.36% (old-boundary bias −0.21 pp) | CORRECTED |
| 4 Nodal Hamiltonian and endpoint derivation | stated | reproduced only | NOT CHECKED |
| Antiparticle bound ≤ 2.69% | bound | 1.06–1.07% | CONFIRMED |

### 59a1a71: conversion streaming

| Item | Codex | Consolidated | Verdict |
|---|---|---|---|
| 1 Γ/H ≈ 27; half first-converted by z ≈ 974500; survival kept separate | stated | 27.0066; 974383–974487; survival 1e-22 vs net ≲2.6e-7 | CONFIRMED |
| 1 Net anti-X today | 7.9e-7 | ≲2.6e-7 | CORRECTED |
| ω/m_X 0.861195 vs E/(Q m_X) 0.923316 | stated | Two independent profile solvers (conversion, xcharge) agree to 2e-13 | CONFIRMED |
| 2 Absorbing boundary 1+2Ω; u_esc; momentum at crossing | 2.532076 | identical; MC q matches to 0.02%. The δ-grid code itself was not reproduced | CONFIRMED |
| 2 Median exit | 890452 | 889170–889455 | CORRECTED (minor) |
| 3 Horizons | 194.95 pc / 2.36919 / 8.37832 Mpc | ≤0.02%; about 199 pc at 20 MeV with g*(T); elastic scattering roughly halves the effective length [E] | CONFIRMED (collisionless upper bounds) |
| 3 Path-moment estimator | excluded | — | NOT CHECKED |
| 4 Mass fractions | 3.175 / 3.359 / 6.718 / 5.469% | 3.17513 / 3.35882 / 6.71764 / 5.46881% | CONFIRMED |

## 2. Disagreements between tracks and how they resolve

1. **Retention uncertainty.** Two tracks bound the 0.019% agreement differently:
   - reproduce puts it at about 0.1% relative;
   - independent_solver says agreement below 0.03 pp is "not convergence evidence".

   Five converged runs from three independent codes agree to ±0.003 pp (K3). So the numerical solution is
   converged at Codex's level. The 0.02–0.03 pp scale is the physical robustness to 1e-6 perturbations. Reproduce's
   0.1% comes from old-executable runs whose remote damper still reflects fast waves. It is a bracket, not an error
   bar.
2. **bb211d1 corrected value.** reproduce gives 35.38% (old executable, damper at 720); independent_solver gives
   35.36%. The old executable gives 35.335–35.383% depending on where the damper sits, so 35.36% is the
   consolidated value.
3. **Bias mechanism.** independent_solver says "soft-wave return accounts for essentially all"; reproduce says
   "fast waves". K1, computed from independent_solver's own runs, settles it:
   - the deficit appears at τ ≈ 1800–2500, which is a fast-wave return;
   - the in-window recovery from −0.28 to −0.06 pp fits p ≈ 0.05 returns, whose return time is about 9600.

   Reproduce is right about the offset. Independent_solver is right about the slope flattening and wrong about the
   attribution.
4. **Late window.** spectrum's "CONFIRMED to grid error" for 34.574% should read CORRECTED. All three tracks imply
   34.74%. reproduce's notes still quote 34.57%.
5. **Escape totals.** independent_solver and spectrum differ by 9.3 in Q and 7.3 in E, which is below the
   trajectory sensitivity (K2). Both values are consistent.
6. **Antiparticle number.** spectrum's "256.1 ± 0.4" leaves out the other converged values (255.0–255.4). The
   consolidated value is 255.5 ± 0.6.
7. **Escape by τ = 48000.** These numbers come only from coarse-dt runs:
   - independent_solver's 25503 (65.4%) carries its dt bias of +201 at τ = 12000, giving about 25300 [E];
   - spectrum's net charge is 25447, or about 25385 after its offset [E].

   They still disagree by 0.3% after correction. This is unresolved and minor.

## 3. Overreach in the tracks

- **operators 6c.** "REFUTED" targets a separation that Codex explicitly disclaimed. What the track actually found
  is a correction of the interpretation, plus a new negative result. Its "own solver" agrees with Codex to 2e-11,
  as do reproduce #1 and operators 1b: these are re-implementations of the same scheme, which test coding fidelity
  rather than the method.
- **reproduce.** Its two "independent discretizations" are both Codex executables. That makes them a cross-family
  control, not independent evidence.
- **independent_solver.** The mechanism claim is wrong (§2.3). The 65.4% at τ = 48000 is not corrected for dt
  bias.
- **spectrum.** Its certification criteria were stated before use, which is sound. Three problems remain:
  - "antiparticle" means classical negative-frequency content at leading order [A];
  - the ± on N_X̄ is too narrow;
  - on the 25.1% claim, it swapped in its own criterion for Codex's rather than recomputing Codex's.
- **conversion.**
  - The labels "excluded (approx.)" break the brief's no-exclusion rule; they should read "in tension (screen)".
  - The z ≈ 10⁶ call rests on DESI DR1 n_* = −2.309 ± 0.019, read from an abstract or search summary, with the
    correlation ignored, for a feature that does not have a ΛCDM shape. It moves to "marginal" under elastic
    coupling. The structure verdicts go beyond what the approximations support.
  - The net anti-X value is an upper envelope, so it should be quoted as ≲2.6e-7.

## 4. Corrections Codex should apply

1. Quote the I025 totals as values for the h = 0.25 trajectory. The physical values are Q = 23519, E = 25482 and
   mQ/E = 0.9230.
2. The late window is 34.74%, not 34.57%. Q60 is still falling: about 34.45% at τ = 38400–48000 [E].
3. Attribute the old-boundary bias to fast-wave returns at τ ≈ 1800–2500. The bb211d1 value becomes 35.36%.
4. The reflection figure "< 2e-8" holds for k ≥ 0.01 only. At k = 0.005 it is 1.35e-4 (charged) and 2.6e-2
   (neutral). The free_soft packet covers k ≳ 0.016 only.
5. Label the bound-mode table as dependent on the breathing phase, not on resolution. Do not use the frozen split
   for excited cores.
6. Base the linearity statements on the slow mass shift. The antiparticle fraction is 1.06–1.07%. Drop "29%" as a
   measure of sensitivity.
7. Add rate-table nodes across p = 0.49–0.53. The optical depth becomes 49.7. The anti-X tail is ≲2.6e-7, the median
   exit z ≈ 889400, and the 20 MeV horizon about 199 pc with g*(T).
8. State two omissions:
   - **The s → XX̄ channel.** Its width is Γ = 0.0188 g_X² m_X, and the classical runs leave it out. The escaping X̄
     count is 256 without it and up to 365 with it [A].
   - **Post-escape elastic scattering.** Because of it, the horizons are upper bounds.
9. Code hygiene:
   - create `data/` in refinement.py;
   - check the ofstreams;
   - give the late window as 9610–12000;
   - pass the constants as parameters.

## 5. Certified now, and still open

**Certified, by the stated criteria and for this prepared condensate only:**
- the causal layer's linear reflection (k ≥ 0.01) and its stability;
- the summation-by-parts operators;
- retention of 35.36% over τ = 9600–12000;
- escape totals by τ = 12000;
- N_X̄, N_s and the X spectrum at p/m_X ≥ 0.03, at classical leading order;
- the profile ω, the node rates, Γ/H, u_esc, the collisionless horizons and the mass fractions.

**Open:**
- the X spectrum below p/m_X = 0.02 (1530 X), the halo, future core emission and asymptotic retention;
- a bound/coherent split for breathing cores;
- reflection below k = 0.01;
- the s → XX̄ decay;
- the portal kernel and bb211d1's nodal derivation;
- post-escape transport;
- a structure likelihood;
- 3D and cosmological formation, general hot-core leakage and the vacuum density.

## 6. The most decisive next calculation

**Formation-sourced conversion.** Run the absorber-free, certified free-mode extraction (spectrum criteria C1–C6;
r ≥ 480; T ≥ 48000; two independent codes) on a small family of prepared condensates spanning density and
charge-to-energy ratio. Then replace 59a1a71's imposed s-decay source with the per-object X̄ yield and momentum
spectrum this produces: the direct prompt X̄, plus the conditional s → XX̄ channel. Run that through the corrected
conversion network on breathing, charge-losing cores.

The prototype's X̄ yield is 1.9% of retained core charge, or about 2.7% if the s quanta decay to XX̄. 59a1a71
imposes 6261 anti-X per Q = 200000 object, which is 3.1%. These are different objects, so the comparison sets only
the order of magnitude. Even so, the 3.4–6.7% converted-X fraction that the structure screen tests scales directly
with this yield. The calculation decides whether that fraction exists, and it removes the imposed source and the
frozen Q = 200000 assumption in one step.

## Reproduce

    python3 -I code/critic_checks.py results/critic_checks.json

It reads the scratchpad runs of `independent_solver` (review/independent_solver/runs/), the spectrum and conversion
result JSONs, and the Codex JSON extracted with:

    git archive ff49e5d research/reality_framework/grut_causal_receiver_20261024 | tar -x -C $SCR/critic/ff49e5d
    git archive bb211d1 research/reality_framework/grut_core_conversion_20261023 | tar -x -C $SCR/critic/bb211d1

Credits: CFS/PML absorbers (Berenger; Kuzuoglu–Mittra; Roden–Gedney); Q-ball formation (Kusenko–Shaposhnikov;
Enqvist–McDonald; Kasuya–Kawasaki); mixed and decaying dark-matter structure constraints, as cited by the conversion
track.
