# Skeptic of track pt_m1tilt: the φ-condensation transition under Codex's GRUT-M1 potential

1 November 2026. Claude review of Codex 5f5513b / d4747cb (issue #13). Code: `skeptic/code/`, results:
`skeptic/results/`, per-run files (with traces): `$SP/review/pt_m1tilt_skeptic/runs/`. Same-session AI work on one
supplied, reduced, classical model: **everything here is [C]; nothing is a cosmological result, a validation or an
exclusion.**

## Method (independent of the track's)

- **Integrator** (`sk_sym.c`): a symplectic Kahan–Li 17-stage order-8 composition of drift–kick–drift steps in the
  *extended* phase space (s is a coordinate, p_s its momentum, H_ext = p²/2 + V(q, s) + p_s), so the drift M(s) is
  exact and E(s) − E(0) − W(s) = ΔH_ext is the integrator check (W = −Δp_s is the work integral). The step depends on s
  only (h = 2π/(32 Ω_max(M(s)))), so every step is a symplectic map. Not the track's adaptive DOP853; no code shared.
  Force on χ and ∂V/∂s re-derived from V = λ_X/4 r⁴ + g²/2 χ² r² + λ_P/4 (χ² − M²)² − ε_t M (χ − M)³ and checked against
  central differences (χ-force 3×10⁻⁸, ∂V/∂s 2×10⁻¹⁰ relative; `sk_verify.json`). Convergence order measured on the
  coupled system: slope 8.0 (KL8) and 6.0 (Yoshida-6) against scipy DOP853 at rtol 10⁻¹³; H_ext conserved to 10⁻¹² in a
  test run.
- **Orbits**: own time-domain evaluation (DOP853 rtol 10⁻¹³, J_r = (1/2π)∫p_r² ds over one radial period):
  J_r/q = 1.76372528 (τ-20), 0.48049185 (late); T_r = 0.81446568, 1.06160355. Same as the track to 8 digits.
- **Measures**: own implementation of n₋/q = [(p_X²/2 + m_X²r²/2)/m_X − q]/(2q) (m_X = gM, and g|χ|) and
  N_φ/q = ⟨(δχ′²/2 + m_φ²δχ²/2)/m_φ⟩/q, 60 X periods at 5, 10, 20, 40 M_c,apo, plus a 400-point log grid in M with
  4-radial-period short windows (n₋, N_φ, E_X, E_φ, E_int, quartic term, χ/M, H_ext, charge).
- **Same start as the track**: χ₀ at the instantaneous tilted minimum near 0, χ′₀ = 0, M₀ = 0.5 M_c,peri,
  M = M₀e^{εω_r s} (growth 0) or M₀(1 + κs) (growth 1).
- Scripts: `sk_common.py`, `sk_run.py` (restartable runner, ≤ 2 workers), `sk_freq.py` (spectra of χ and r² in
  32-radial-period windows along M), `sk_compare.py`, `sk_aggregate.py`, `sk_verify.py`.

## Verdict table

| # | Track claim | Track | Skeptic | Verdict |
|---|---|---|---|---|
| 1 | Runs reproduce (ε = 10⁻³, 10⁻⁴, both orbits, 10 phases) | per-run files | run-by-run \|Δn₋/q\| ≤ 1.2×10⁻⁷ (10⁻³), ≤ 1.3×10⁻⁶ (10⁻⁴); H_ext ≤ 4×10⁻¹¹, charge ≤ 2×10⁻¹² | **CONFIRMED** |
| 2 | Lossy runs at ε = 10⁻⁵ reproduce; independent of integrator | τ-20: 1.5155, 1.5420, 1.5147, 1.5194; late 0.1425, 0.1401 (ph 0, 0.5) | τ-20: 1.5148, 1.5422, 1.5148, 1.5194 (\|Δ\| ≤ 7×10⁻⁴); late 0.1409, 0.1314 (Δ 1.5×10⁻³, 8.7×10⁻³); step 32 → 48 per period: τ-20 unchanged to 10⁻⁵, late 0.1409 → 0.1446; Yoshida-6: 1.5148 / 0.1414 | **CONFIRMED** (late orbit weakly sensitive, as the track noted) |
| 3 | ε = 10⁻⁶ | τ-20 1.4896 (ph 0.358), late 0.1054 (ph 0) | 1.4880 / 0.1411; (J₀ − n₋)/N_φ = 1.5001 | CONFIRMED for τ-20; late differs by 0.036 (sensitivity) |
| 4 | Adiabatic for ε ≥ 10⁻⁴ (τ-20) / 10⁻³ (late) to ≤ 2×10⁻⁴ | | loss ≤ 2×10⁻⁴ (τ-20, 6 phases) and ≤ 10⁻³ (late, 6 phases) | **CONFIRMED** |
| 5 | Late window is the free massive regime; measure is the free-basis antiparticle number | | at 40 M_c,apo: λ_X r⁴/4 ÷ E_X ≤ 4×10⁻⁷, E_int/(m_X q) ≤ 2×10⁻⁵, rms δχ/M ≤ 10⁻³, n₋(gM) = n₋(g\|χ\|) to 10⁻⁵; n₋ = (E_X/m_X − q)/2q is the n₋ of ψ = Ae^{−imt} + Be^{imt} (E/m = n₊ + n₋, q = n₊ − n₋) | **CONFIRMED** |
| 6 | Work and charge ledgers | W ≤ 3×10⁻⁷ (10⁻⁵), 10⁻⁵ (10⁻⁶) | H_ext ≤ 10⁻⁸ (10⁻⁵), 7×10⁻⁹ (10⁻⁶), ≤ 4×10⁻¹⁰ (3×10⁻⁷–3×10⁻⁸ at reduced M_end); charge ≤ 5×10⁻¹¹ in every run | CONFIRMED, improved |
| 7 | "Sharp" ε threshold (3×10⁻⁵–10⁻⁵ τ-20; 3×10⁻⁴–10⁻⁴ late) | | fine scan, 3 phases: τ-20 captured phases 0/3 at ε ≥ 3×10⁻⁵, 2/3 at 2.5×10⁻⁵ and 2×10⁻⁵, 3/3 at ≤ 1.5×10⁻⁵; late 0/3 at ≥ 3×10⁻⁴, 1/3 at 2.5 and 2×10⁻⁴, 3/3 at ≤ 1.5×10⁻⁴ but with losses rising 19–49 % (1.5×10⁻⁴) → 49–66 % (10⁻⁴) → 69–73 % (10⁻⁵) | **CORRECTED**: a phase-dependent capture window of a factor 2 (τ-20) to 4 (late) in ε, not a sharp threshold |
| 8 | Growth law M ∝ e^{εω_r s} (radiation era is linear) | not tested | linear M = M₀(1 + κs) with κ matched to εω_r at the release M (6 M_c,peri τ-20; 2.5 late): τ-20 loss 0.05 % (10⁻⁴), 0.2–0.4 % (3×10⁻⁵), 13–15 % (10⁻⁵); matched at M_c,peri instead (6× slower at release): 12–13 % already at ε = 10⁻⁴ | **CONFIRMED**: the control parameter is the local d ln M/ds at the release point; the law itself does not matter |
| 9 | Mechanism: slaved φ swing, free φ released at 5.6–7 / 2.2–3 M_c,peri; "resonance-capture-like, characterised, not derived" | | the released χ component is born at **1.501 f_r** (τ-20, 5.66 M_c,peri) and **1.497 f_r** (late, 2.03 M_c,peri) and then drifts to m_φ/ω_r = 1.2253; in the adiabatic run no half-integer component ever appears; the exchange conserves **2J_r + 3N_φ**: (J₀ − n₋)/N_φ = 1.4999 ± 0.0001 in the track's own 48 τ-20 runs at 10⁻⁵, 1.499 ± 0.006 (late, 64 runs at 10⁻⁴), 1.500 in every lossy run of mine | **IDENTIFIED**: the 2ω_φ = 3ω_r subharmonic (parametric) resonance of φ driven by the X radial motion; see §2 |
| 10 | Slow limit: loss grows as ε → 0 and "saturates slowly"; τ-20 n₋/q = 1.49–1.52 (0.85 J_r), late 0.11–0.15, quoted as the slow-limit values | 10⁻⁵: 1.520 ± 0.012; 10⁻⁶: 1.489 (8 runs) | 3×10⁻⁷: 1.4877, 1.4856 (saturation-like); **10⁻⁷: 0.593, 0.403, 1.294 (losses 66, 77, 27 %)**; **3×10⁻⁸: 1.301 with (J₀ − n₋)/N_φ = 1.25, not 1.5** (a second, 5:4 resonance active). Late: 3×10⁻⁷: 0.142, 0.101; 10⁻⁷: 0.117, **0.024** | **REFUTED as a limit**: loss(ε) is non-monotonic and widens below 10⁻⁶; the quoted numbers are the values of one ε decade, not a slow limit; nothing can be said about ε ~ 10⁻¹² |
| 11 | Sensitivity to the supplied tilt ε_t (0.002) | not tested | τ-20, ε = 10⁻⁵: loss 0 (ε_t = 0.0015, 0.0017, 0.0018), 0.2–0.4 % (0.0019), 14 % (0.002), 48–52 % (0.0022), 52–54 % (0.0025), 63–70 % (0.003), 0.7–0.9 % (0.004); ε_t = 0.001 (two vacua): 47–79 % and one of three runs ends in the false vacuum χ = −0.77 M. Late, ε = 10⁻⁴: 0 (0.001, 0.0015), 0–6 % (0.0017), 0–28 % (0.0018), 49–66 % (0.002), 55–69 % (0.0022), 73–95 % (0.0025), 94 % (0.003), 92–93 % (0.004) | **NEW**: the headline numbers are specific to ε_t = 0.002 within ±5 %; ±10 % in ε_t moves the τ-20 loss between 0 and 50 % |
| 12 | Deterministic, phase-smooth, no chaotic split | | true for τ-20 at ε ≥ 3×10⁻⁷ (phase spread ≤ 0.03); at 10⁻⁷ the three phases give 0.40, 0.59, 1.29 — a broad phase dependence returns; late at 10⁻⁷: 0.117 vs 0.024 | REFINED |
| 13 | ε_t = 0 reproduces the untilted split; β = 0.2 alone keeps J_r | | not re-run | not tested (ε_t = 0.001 result above is consistent with the Z₂-like regime returning below the one-vacuum threshold) |

## 1. Reproduction and robustness (A)

`sk_reproduction.json`. Ten phases at ε = 10⁻³ and 10⁻⁴ agree with the track's DOP853 runs to ≤ 1.3×10⁻⁶ in n₋/q
and N_φ/q, so the two codes integrate the same equations. The lossy ε = 10⁻⁵ runs reproduce run by run (τ-20 to
7×10⁻⁴; late to 10⁻²) and are unchanged under 32 → 48 steps per fastest period and under a different composition
(Yoshida order 6, 64 steps). The step-count limit of the track (work identity 10⁻⁵ at ε = 10⁻⁶) is not an issue here:
my extended-energy check is ≤ 7×10⁻⁹ at 10⁻⁶ and ≤ 4×10⁻¹⁰ in the longest runs, with the same n₋. The ε = 10⁻⁵ loss
is real.

## 2. Mechanism: a 2ω_φ = 3ω_r resonance, not a generic "slaved-to-free" release (D)

`sk_track_ratio.json`, `results/freq/*.json`, `sk_freq.py`.

- **Action bookkeeping.** In every lossy run, of the track and of mine, the radial action lost and the φ quanta gained
  satisfy ΔJ_r = −(3/2) ΔN_φ, i.e. **2n₋ + 3N_φ = 2J_r** (track: 1.00001 ± 0.00001 of 2J_r over 48 τ-20 runs at
  10⁻⁵; late 1.0004 ± 0.002 at 10⁻⁴, drifting to 0.990 at 10⁻⁶). That is the conserved combination of a resonance with
  3ω_r = 2ω_φ: three radial quanta become two φ quanta.
- **Spectra.** In 32-radial-period windows the χ spectrum is pure harmonics of the measured f_r (harmonic power fraction
  1.000, half-integer fraction 0.000) until the release point. There a new χ line appears at exactly 1.50 f_r — τ-20,
  ε = 10⁻⁵: 1.501 f_r at 5.66 M_c,peri (power ratio 0.04, then 0.09, 0.07, 0.12, 0.17, 0.27…); late, ε = 10⁻⁴: 1.497 f_r
  at 2.03 M_c,peri (0.013 → 0.28 → 0.39 → 0.55 → 0.91) — and then detaches, its ratio to f_r drifting to the free value
  m_φ/ω_r → 1.2253. In the adiabatic τ-20 run (ε = 10⁻⁴) no half-integer line appears at any M (half-integer fraction
  0.000 throughout 4.5–14 M_c,peri) although the same frequencies are crossed. The coupling that does this is
  g²δχ²r²/2 (parametric in δχ, driven by the third harmonic of r²); the linear tilt alone has the wrong structure, which
  is why β = 0.2 conserves J_r.
- **Lock.** At ε_t = 0.004 (τ-20, ε = 10⁻⁵) the line sits at 1.5002–1.5004 f_r from 4.4 to 7.1 M_c,peri while its power
  grows to 0.77 of the main line and n₋ falls to 0.70 — and then the exchange **reverses**, the line leaves (1.47, 1.41,
  1.37…) and n₋ recovers to 1.749. So the 3:2 resonance can capture, hold and release; the net loss is the
  difference between the forward and backward exchange and depends on the details of the crossing (ε, ε_t, phase).
  That is why the loss is so irregular (§3, §4). At ε = 10⁻⁷ the release at 5.96–6.6 M_c,peri is violent: f_r doubles,
  n₋ drops from 1.18 to 0.57 within one window.
- **Energy budget (C).** Actions: J_r → n₋ + (3/2)N_φ. Energies at the end: E_X = m_X(q + 2n₋), E_φ = m_φN_φ =
  2.449 m_X N_φ against the X deficit 2m_XΔJ_r = 3m_XN_φ; the exchange itself conserves energy (2ω_φ = 3ω_r), and the
  0.55 m_XN_φ difference at the end is less work done by the drift afterwards (φ quanta cost m_φ ∝ M, radial quanta
  ω_r → 2m_X ∝ M, but the ratio moves from 1.5 to 1.22). E_int and the quartic term are negligible in the window (claim 5).
  Nothing is in the circular part: the charge is conserved to 5×10⁻¹¹.

## 3. The ε dependence is not a threshold plus a slow saturation (D, F)

`sk_scan.json`, `sk_slow.json`. Loss 1 − n₋/J_r per phase (τ-20: phases 0, 0.3, 0.7; late the same):

| ε | τ-20 | ε | late |
|---|---|---|---|
| 10⁻⁴ … 4×10⁻⁵ | ≤ 0.2 % | 10⁻³ … 4×10⁻⁴ | ≤ 1.3 % |
| 3×10⁻⁵ | 1.1, 0.1, 0.1 % | 3×10⁻⁴ | 0.5, 0.2, 3.0 % |
| 2.5×10⁻⁵ | 0.2, **10.1, 10.1** % | 2.5×10⁻⁴ | 0.1, 4.3, 5.3 % |
| 2×10⁻⁵ | **14.0, 14.5**, 0.4 % | 2×10⁻⁴ | 2.7, 2.1, **20.0** % |
| 1.5×10⁻⁵ | 12.5, 15.0, 13.1 % | 1.5×10⁻⁴ | 19.0, 19.8, 48.9 % |
| 10⁻⁵ (4 phases) | 12.6–14.1 % | 10⁻⁴ (6 phases) | 48.5–65.5 % |
| 10⁻⁶ | 15.6 % | 10⁻⁵ (4) | 68.6–72.7 % |
| 3×10⁻⁷ (M_end 10) | 15.7, 15.8 % | 10⁻⁶ | 70.6 % |
| **10⁻⁷ (M_end 5)** | **66.4, 77.2, 26.6 %** | 3×10⁻⁷ (M_end 10) | 70.4, 78.9 % |
| **3×10⁻⁸ (M_end 5)** | **26.2 %, with ΔJ_r/ΔN_φ = 1.25** | 10⁻⁷ (M_end 5) | 75.6, **95.1** % |

(M_end 10 and 5 M_c,apo bias n₋ by −4×10⁻⁴ and −3.5×10⁻³ for τ-20, +4×10⁻⁴ for late, from the 40-M_c,apo runs'
intermediate windows.) Reading:

- The onset is a **capture window**, not a threshold: at a given ε some phases are captured and lose ~10–14 % (τ-20)
  while others lose < 0.5 %. The captured fraction rises from 0 to 1 over a factor 2 (τ-20) or 4 (late) in ε.
- Between 10⁻⁵ and 3×10⁻⁷ the τ-20 loss is flat at 14–16 %, which is what the track read as a slow saturation. At
  **10⁻⁷ it jumps to 27–77 %** with a broad phase dependence, and at **3×10⁻⁸ the exchange ratio changes to 1.25**:
  a second resonance (4ω_φ = 5ω_r, conserved 4J_r + 5N_φ; m_φ/ω_r crosses 1.25 at ≈ 17 M_c,peri in the massive regime)
  takes part. The late orbit at 10⁻⁷ loses 76–95 %.
- So the quantity the track extrapolates, loss(ε), is **non-monotonic and discontinuous in ε**, as expected for
  successive resonance captures with partial reversal. Its ε = 10⁻⁵–10⁻⁶ values are not "slow-limit values"; they are
  the values on one plateau. **About ε ~ 10⁻¹² nothing can be concluded from this model's scan**, in either direction
  (the τ-20 loss at 10⁻⁷ ranges from 0.27 to 0.77 J_r over three phases; a run at ε = 10⁻¹² would take ~10⁵ CPU-hours
  with either code). The only robust statements are: (i) the slow transition under the M1 potential is *not*
  adiabatic for any ε below the capture window, and (ii) the loss is quantised by the resonance order (2:3, then 5:4).

## 4. The result is a property of ε_t = 0.002 (E)

`sk_epst.json`. ε_t enters only through the tilt (the one-vacuum threshold is ε_t = 0.001716). With everything else
fixed (τ-20, ε = 10⁻⁵; late, ε = 10⁻⁴), the loss as a function of ε_t is: zero up to 0.0018, 14 % at 0.002, 50 % at
0.0022–0.0025, 65–70 % at 0.003 and back to ~1 % at 0.004 for τ-20; for the late orbit it rises from zero (≤ 0.0015)
through 50–65 % (0.002) to 93–95 % (≥ 0.0025). Below the one-vacuum threshold (ε_t = 0.001) the τ-20 outcome is
chaotic again (47 % and 79 % at phases 0.5 and 0; the neighbouring phase 0.001 ends in the false vacuum χ = −0.77 M,
where the +M measures are meaningless). Codex's ε_t = 0.002 sits 17 % above that threshold. Since ε_t is a supplied
parameter with no derivation in the packet, the headline "0.85 J_r + 0.16 φ quanta" has no robustness: ±10 % in ε_t
moves it between J_r and 0.5 J_r. All exchanges at ε_t ≥ 0.0019 are 2:3 to 10⁻³.

## 5. Issues with NOTES.md (F)

1. **§3.1 / §4.2 "sharp threshold … the slower the transition, the more radial action it loses … still creeping up at
   10⁻⁶"**: the threshold is a capture window (§3), and the monotone reading rests on two ε decades; the next decade
   breaks it (10⁻⁷: 0.27–0.77 J_r; 3×10⁻⁸: a different resonance). The sentence "the ε → 0 limit is therefore not J_r"
   is right; "the loss grows as ε falls" as a trend is not established.
2. **§4.3 "numbers per orbit charge, taken at ε = 10⁻⁵–10⁻⁶ … slow-limit values of this scan"**: these are plateau
   values at ε_t = 0.002; §4 above shows they are not stable against ε_t and §3 that they are not stable against ε.
   Any downstream use (e.g. "≈ 1.67–1.70 X̄ per charge if φ → XX̄") inherits this.
3. **§3.4 mechanism**: "characterised, not derived" undersells the data. The track's own per-run files contain the
   identification: 2n₋ + 3N_φ = 2J_r to 10⁻⁵, i.e. the 2ω_φ = 3ω_r resonance, confirmed by the 1.5 f_r line at the
   release point. The "grazing 1:2 commensurability" candidate is not it (the linear m_φ/ω_r ≈ 2 at the release is the
   frozen-orbit number; the resonant frequency is the nonlinear one of the large swing, which is 1.5 ω_r there).
4. **§3.2 "deterministic at every ε"**: holds for the ε scanned by the track; at 10⁻⁷ the τ-20 outcome spreads over
   0.40–1.29 across three phases (still deterministic, but as broad as the untilted split).
5. **§2 "the same drift law and start as pt_matching"**: fine, and the radiation-era linear law gives the same result
   at matched local rate (claim 8), so this is not a weakness; but the *location* of the matching matters by a factor 6
   in the effective ε, which should be stated when comparing with the cosmological ε.
6. **Labels**: the track's "[C]" is right; its §4.3 "[E]" extrapolation should say that the model's own ε dependence is
   irregular below 10⁻⁶, not only that "the loss is still rising".

## Limits of this skeptic pass

- ε = 10⁻⁷ rests on 3 (τ-20) + 2 (late) runs at M_end = 5 M_c,apo, ε = 3×10⁻⁸ on one run; these show the irregularity
  but do not map it. ε_t was scanned at one ε per orbit and two phases.
- The 5:4 identification at 3×10⁻⁸ rests on the exchange ratio alone (no spectrum taken).
- Modes 1 (β = 0.2) and ε_t = 0 were not re-run.
- All limits of pt_m1tilt apply (homogeneous zero mode, classical, no gauge/Y/vector fields, pure-quartic late orbit).

## Reproduce

```bash
cd research/reality_framework/grut_formation_review_20261025/tracks/pt_m1tilt/skeptic/code
python3 -B sk_verify.py
python3 -B sk_run.py tau20:1e-3:0.0 tau20:1e-3:0.25 tau20:1e-3:0.5 tau20:1e-3:0.75 late:1e-3:0.0 late:1e-3:0.25 late:1e-3:0.5 late:1e-3:0.75   # likewise 1e-4, 1e-5
python3 -B sk_run.py tau20:1e-5:0.0:nper=48 tau20:1e-5:0.0:method=6,nper=64 late:1e-5:0.0:nper=48 late:1e-5:0.0:method=6,nper=64
for e in 1e-4 8e-5 6e-5 5e-5 4e-5 3e-5 2.5e-5 2e-5 1.5e-5; do python3 -B sk_run.py tau20:$e:0.0 tau20:$e:0.3 tau20:$e:0.7; done
for e in 1e-3 8e-4 6e-4 5e-4 4e-4 3e-4 2.5e-4 2e-4 1.5e-4; do python3 -B sk_run.py late:$e:0.0 late:$e:0.3 late:$e:0.7; done
for e in 1e-4 3e-5 1e-5; do for rm in 6 1; do python3 -B sk_run.py tau20:$e:0.0:growth=1,rm=$rm,mend=10 tau20:$e:0.5:growth=1,rm=$rm,mend=10; done; done
for e in 1e-3 3e-4 1e-4; do for rm in 2.5 1; do python3 -B sk_run.py late:$e:0.0:growth=1,rm=$rm,mend=10 late:$e:0.5:growth=1,rm=$rm,mend=10; done; done
for t in 0.001 0.0015 0.0017 0.0018 0.0019 0.0022 0.0025 0.003 0.004; do python3 -B sk_run.py tau20:1e-5:0.0:epst=$t tau20:1e-5:0.5:epst=$t late:1e-4:0.0:epst=$t late:1e-4:0.5:epst=$t; done
python3 -B sk_run.py tau20:1e-6:0.358016 late:1e-6:0.0; python3 -B sk_run.py tau20:3e-7:0.0:mend=10 tau20:3e-7:0.5:mend=10 late:3e-7:0.0:mend=10 late:3e-7:0.5:mend=10
python3 -B sk_run.py tau20:1e-7:0.0:mend=5 tau20:1e-7:0.25:mend=5 tau20:1e-7:0.5:mend=5 late:1e-7:0.0:mend=5 late:1e-7:0.5:mend=5
python3 -B sk_run.py --budget 540 tau20:3e-8:0.0:mend=5    # repeat the call until 'done' (restart file); ~25 min
python3 -B sk_freq.py tau20 1e-5 0.0; python3 -B sk_freq.py late 1e-4 0.0 0.002 1.2 30 44; python3 -B sk_freq.py tau20 1e-4 0.0 0.002 4.5 14 16
python3 -B sk_freq.py tau20 1e-5 0.0 0.004 4 14 14; python3 -B sk_freq.py tau20 1e-7 0.0 0.002 4 12 12
python3 -B sk_compare.py; python3 -B sk_aggregate.py
```

Runtime about 2.5 h on two workers (3×10⁻⁸ alone 25 min; 10⁻⁷ runs 6 min each).
