# Consolidated review: the heavy-mediator window after independent verification, and Codex 2908e91's question (track hw_critic)

26 October 2026 · Claude (issue #13) · consolidates tracks hw_verify, hw_envelope and hw_gates with ls_scaling,
ls_threshold and ls_critic/CONSOLIDATED.md, and adds cross-track recomputations (§2; own code in `code/`, numbers in
[NOTES.md](NOTES.md)). Packet: Codex 2908e91 (`origin/codex/grut-concept-checkpoint-20261025`,
`research/reality_framework/grut_concept_checkpoint_20261025/`), read only through `git archive` into the scratchpad.
No Codex code was run by any track. Nothing was committed.

**Scope.** One supplied tree-level working model. All tracks are same-session AI work, so this is not outside
replication. Nothing here is evidence that X, Y, Σ, φ or the vector exist. It is not a validation, an exclusion, a
viable point or a GRUT prediction. Domain throughout: homogeneous classical release, supplied CP phase, late-decay
ledger, case-L proxy at f = 0.01; light target g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV, Q = 200000,
Y_Q = 4.8294e-9. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued.

## Bottom line

**Codex's question** — within the declared domain, does producing sufficient lasting charge necessarily generate
excessive processed late debris? — **No.** The heavy-mediator window **exists after independent verification**: five
integrators (ls_scaling, hw_envelope, hw_verify, hw_gates, this track) give ε = −0.146, 1.00 Y_Q and mediator debris
R = 1.14 at the spot point (H = 1.334e11, M_Y = 4.00e12, h = 5.33e12, μ = 9.55e11 GeV, κ = 3.46e-3, r = 0.45, c6 = 1.25e6),
and R = 1.14 passes every case-L threshold variant on record: ×71 against the direct isentropic proxy (81.5–82), ×6–7
against the energy-conserving proxy (7.0–8.3), ×17 with the on-shell s route (19.6). Codex's point (i) is realised
dynamically: the forced Y carries the X⁶ torque and leaves 0.35 free quanta per comoving H_i³ (1.3e-4 per unit charge),
0.007 for a slow release.

**What does threaten the window is not produced by the charge source.** The dominant processed late debris at the spot
point is the inflationary longitudinal vector (hw_gates piece A: R_A = 48.5 from own modes, 64 from the GMR formula at
e = 1e-7). It passes the isentropic proxy (0.61–0.82 of the direct R_L) and **fails the energy-conserving proxy by
6–9** (2.5–3.3 with the s route) at the stated e = 1e-7; since R_A ∝ 1/e, it passes every variant for e ≳ 3e-6
(margin ≥ 2.3) and by ≥ 4.7 for e ≥ 1e-5. Which proxy is right is the decisive open calculation (§7).

**Where the implication does hold** (all [C]): the resonant regime m = M_Y/H ≲ 12 (proxy) or ≲ 20 (ke); the stated
couplings at every release duration tried (ε falls to 0.06 by Γ_φ = H_i and reverses sign below 0.3 H_i while R rises
2–6×); and any preparation or holder that leaves |Y − Y_*| > 0.37 H at release (Codex's preparation: 77 quanta,
R = 1415; a Y-coupled holder quenched in ≪ 0.3 H⁻¹: R = 160–1415).

## 1. Cross-check of shared numbers

| Quantity | ls_scaling | hw_verify | hw_envelope | hw_gates | ls_critic | this track | Status |
|---|---|---|---|---|---|---|---|
| Spot ε | −0.14572 | −0.14567 | −0.14576 | −0.14568 | — | −0.14568 | agree ≤ 0.03 % |
| Spot R_central / yield | 1.144 / 1.00 | 1.144 / 1.0000 | 1.144 / 1.0007 | 1.144 / 1.0000 | — | 1.144 / 1.0000 | agree |
| N_YH, N_YL, N_σ | 0.01725, 0.3345, 3e-9 | same | same | same | — | same | agree |
| Codex-preparation R | 1414 | 1415 | 1415 | 1415 | — | 1415 | agree |
| R_L at the spot (proxy) | 88 [table] | **81.5** direct | 88 [Q] | 87.85 [table] | **81.9–82.2** direct | — | table 7–8 % high; direct codes agree to 0.6 % |
| R_L energy-conserving (ke) | — | — | — | — | ≈ 8.2 / 7.0; 15–20 with s | 8.26 / 6.97 / 19.6 (interpolated) | one track; proxy 10× above it |
| Relative-phase mass | 1.84 H [D] | — | — | 2.53 H | — | 2.531 H, f_eff 17.10 | CORRECTED to 2.53 |
| Y relaxation / tolerance | "physical prep" [A] | 3H/2; pass needs \|Y−Y_*\| ≤ 0.37 H; R < R_L after 1.02 held e-folds | — | holder-dependent; D1/D2 need switch-off ≳ 0.3 H⁻¹ | — | held N = 0.5/1/2 → R 385/88/7.9; D1 series 1415/161/11.1/2.4/4.6 | both CONFIRMED |
| Holder-off drift | — | ε to −0.24 at 0.5 e-fold, +0.04 at 2 | — | ε flat to 0.5 e-fold, −0.10 at 1 | — | both reproduced; the difference is the charge accumulated during the unheld stage (§2.3) | hw_gates CORRECTED |
| Codex envelope bound at the spot | — | ladder at 0.20–0.24 of bound | margin 2.65 (R_lead) / 4.07 (R_actual) | — | — | — | consistent (1/0.24 = 4.2) |
| Codex yield table / extrapolation 11.87 | — | reproduced; at M_Y/H_i = 5.96 | reproduced to 1e-12 | — | — | — | CONFIRMED |
| Vector headroom per net charge | 1.6e-3 | 1.3e-3 heavy-quanta equiv. | — | 1.58e-3 | 1.48e-3 direct; 1.1–1.3e-4 ke; 2.5–3.4e-4 ke+s | 1.47e-3; 1.06–1.30e-4; 3.4e-4 | agree |
| Inflationary vector N_A | — | — | — | 2.468 (modes); 3.26 (GMR) | — | 3.2598 (GMR, own conversion) | GMR conversion CONFIRMED |
| Vector e-scaling | e^{−1/2} | — | — | per vector e^{−1/2}; piece ∝ 1/e | — | ∝ 1/e from the formula | CORRECTED to 1/e |
| Finite release (heavy) | open | Γ_φ/H = 1: ε −0.279, R 0.032; 0.1: −0.380, 0.029 | — | — | — | identical to printed digits | CONFIRMED |
| Window m (mediator only) | 13–43 | 12.1–12.3 to 42.8 | 10 outside matching | — | 20–43 under ke | 11.5–11.8 proxy; 20–21 ke; 17 ke+s | consistent |
| Scaling symmetry | machine precision | s = 0.1–10, ≤ 4e-13 | — | — | — | run at the normalised H directly | CONFIRMED |
| Strings, running, light sector | estimates | — | — | computed | — | not redone | NOT CHECKED here |

## 2. Recomputed disputed numbers (this track, own code)

1. **Phase mass** [C]: gauge-projected 6×6 Hessian at the aligned minimum, light angular mode m² = 6.4057 H²,
   m = 2.531 H, f_eff = 17.10 H. hw_gates CONFIRMED; ls_scaling's 1.84 H (leading |X| = 70.7, no Σ kinetic weight)
   CORRECTED. Along δ at the Y-relaxed state the projected curvature is −0.83 H² (no mass; the state is not stationary).
2. **R_L** [Q, two direct codes]: 81.5 (sp_yledger code, hw_verify) and 81.9–82.2 (lsrelic code, ls_critic) replace
   the table's 88. The ke crossing (ls_critic's runs) is 8.26 (μ_c = E_c = 0.861 m), 6.97 (0.923 m), 19.6 with the s
   route (log-interpolated here).
3. **Holder-off drift** [C]: free de Sitter evolution of the Y-relaxed state for N e-folds, then release. With the
   acquired velocities carried (hw_verify's method): ε = −0.163, −0.178, −0.204, −0.239, −0.185, −0.067, +0.039 at
   N = 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2 (hw_verify identical). With the velocities zeroed: −0.146, −0.147, −0.148,
   −0.151, −0.109, −0.016, +0.040 (hw_gates' scan reading identical to 4 %). The phase trajectory δ(N) agrees with
   hw_gates' 1D model to 2 %. The difference is physical in the homogeneous model: the unheld phase rolls under the X⁶
   torque during inflation and the charge density it accumulates (−2179 H³ after 0.5 e-fold, 78 % of the release
   charge) is carried into the release. R stays 1.2–1.5 in every case. Resolution: the holder sets ε to 10 % only if it
   stays on to within ≈ 0.05 e-fold of the end of inflation; |ε| ≥ 0.10 survives a switch-off up to ≈ 1.2 e-folds before
   the end; the debris verdict does not depend on it.
4. **Finite release** [C]: Γ_φ/H_i = 10 / 1 / 0.1 → ε −0.172 / −0.279 / −0.380, N_Y 0.087 / 0.011 / 0.007, R 0.20 / 0.032 /
   0.029 (0.21 / 0.037 / 0.045 with the reheating ⟨a⟩), diluted yield 1.24 / 1.16 / 0.55 Y_Q. hw_verify CONFIRMED.
5. **Adiabatic switch-off of a Y-coupled holder (D1)** [C]: instant R 1415; Δa = 0.1 / 0.3 / 1: R 161 / 11.1 / 2.4 at
   ε −0.117 / −0.092 / −0.065; ∝ H: R 4.6 at ε −0.075. hw_gates CONFIRMED to three digits.
6. **Full debris against every threshold** [C, arithmetic on hw_gates' and ls_critic's numbers]: table in NOTES.md. At
   e = 1e-7: 0.61–0.82 of the direct proxy, 6.0–9.3× the ke proxy, 2.5–3.3× ke + s. Passes all variants for e ≥ 3e-6.
7. **GMR conversion** [C]: N_A = 3.2598 H_i³, identical to hw_gates' formula value; their own mode integration gives
   0.757 of it, which is the 24 % bracket they quote.

## 3. Verdicts on Codex 2908e91

| Codex statement | Verdict | Basis |
|---|---|---|
| Exact sources (−6 Im K; 2 Im K − 4 Im L; 2 Im L), gauge sum 0, n_D sourced by −12 Im L; 48-sample checks | CONFIRMED [D, C] | hw_envelope (sympy + dynamical release, 6e-7) |
| Forced response Y_* = (−AJ + BJ*)/D, ΔV, zero displacement action with nonzero torque | CONFIRMED [D, C]; added: the n_Y source vanishes on the forced branch | hw_envelope |
| ρ_phase = n_X²/2f² + n_X²/72h² | CONFIRMED [D] | hw_envelope |
| (i) A nonzero source does not by itself imply residual heavy quanta | CONFIRMED, now dynamical: 1.3e-4 quanta per unit charge (relaxed), 0.007 per H_i³ at 5416 units (Γ_φ = 0.1 H_i); the condition is adiabatic preparation (holder change ≫ 1/M_Y) | hw_envelope, hw_verify, hw_gates, this track |
| (ii) Envelope bound \|a³n_X\| ≤ 12\|α6\|R⁶/H_i; table 1.22e-8 … 1.74e-4 Y_Q; ∝ r_i^{7/2}; extrapolations 18.24 / 11.87 | CONFIRMED to the printed digits [C]; the 11.87 extrapolation sits at M_Y/H_i = 5.96 | hw_envelope, hw_verify |
| Quiet regime unproductive at the fixed couplings | CONFIRMED and sharpened: ε ≲ 0.005, ≲ 3.5e-3 Y_Q for M_Y/H_i ≥ 30 on Codex's line | hw_envelope |
| The bound "excludes neither all preparations nor changed couplings" | CONFIRMED: the spot point satisfies the same bound with margin 2.65–4.07 on the line η = α6X0²/λ_X fixed | hw_envelope |
| "No published scanned release passes the stated leftover benchmarks under its stated assumptions" | CONFIRMED (ls_scaling one-at-a-time scan; hw_verify finite-Γ at the stated point: worse) | ls_scaling, hw_verify |
| "Neither all possible preparations nor the complete field content has been excluded" | CONFIRMED: the window is the counterexample for the field content | all hw tracks |
| Question: enough lasting charge ⇒ excessive processed late debris? | REFUTED within the domain for the mediator debris at every threshold variant; the window's late debris is the inflationary vector, whose verdict depends on the threshold variant (pass proxy / fail ke at e = 1e-7; pass all at e ≥ 3e-6) | this track |

No numerical error was found in Codex's packet.

## 4. The window after verification

**Exists** (debris test, every threshold variant, mediator debris). **Requires**, all supplied:
1. *Coordinated coupling change* from the stated ratios: r = √2μh/M_Y² 0.0707 → ≳ 0.3 (0.45 used), κ ∝ M_Y at fixed
   η = α6X0²/λ_X ≈ r/24 (κ ×3.46; 0.70 of the stability bound, 0.68–0.79 under one-loop running), M_Y/H_i 70.7 → 30
   (window 12–43 proxy, 20–43 ke, mediator only), h/H_i ≥ 20 (40 used), λ_Y = λ_Σ = 1; H_i free by the exact scaling
   symmetry (1.334e11 GeV for 1.00 Y_Q; c6 ∝ H⁻² binds only below 4e7 GeV); and, if the energy-conserving threshold is
   adopted, e ≳ 3e-6 instead of 1e-7 (no Floquet obstruction up to 0.1 found; longitudinal modes and backreaction not
   computed).
2. *Preparation*: Y within 0.37 H of its forced response at release — the Y-relaxed state (the held-phase attractor at
   rate 3H/2, reached from the Codex state in 1–3 e-folds, for a holder outside the Y equation) or a Y-coupled holder
   switched off over ≳ 0.3 H⁻¹ ≈ 10/M_Y (then ε = 0.08–0.10 and 0.6–0.8 Y_Q, restored by H ×1.2–1.4). A sudden quench of
   a Y-coupled holder fails by 6.5–14.
3. *Phase*: a CP phase supplied through a holder term (D1 |μ_P| = 10 H, a_H ≈ 3–10; D2 κ_P = 0.77 κ; D3 c7 with
   Λ₇ = 2.3e14 GeV); the stated action holds δ = 0 where ε = 0 (e^{−43,000}). A held phase is quiet (spread ≤ 1e-4 rad,
   S/S_max ≈ 1e-30, P(band) = 1). The value of ε is fixed to 10 % only if the holder stays on to within ≈ 0.05 e-fold of
   the end of inflation; |ε| ≥ 0.1 survives a switch-off up to ≈ 1.2 e-folds before the end; the debris does not care.

| Gate | Status at the spot point |
|---|---|
| Mediator debris vs R_L | PASS ×71 (direct proxy), ×6.1–7.2 (ke), ×17 (ke + s) [C, five codes + ls_critic] |
| Codex envelope bound | SATISFIED, margin 2.65–4.07 [C] |
| 1 Inflationary longitudinal vector (piece A; piece B negligible) | PASS under the isentropic proxy (0.61–0.82 of R_L); **FAIL under ke at e = 1e-7** (×6–9; ×2.5–3.3 with s); conditional PASS for e ≥ 3e-6 (≥ ×2.3), e ≥ 1e-5 (≥ ×4.7) [C + arithmetic]; longitudinal Floquet modes and backreaction at larger e OPEN |
| 2 Phase holding | holder with a supplied CP phase REQUIRED; held phase quiet PASS; Y-blind holder PASS for any switch-off; Y-coupled holder PASS only for switch-off ≳ 0.3 H⁻¹ [C, two codes]; the ε value is set by the switch-off epoch — OPEN (supplied) |
| 3 Strings | not formed (U(1) broken through inflation and release, T_h,max = 0.25 h, 0.04 windings) → not a constraint [C, A]; Gμ = 8.9e-12 if a network formed (LVK O3 model B/C-1 would exclude it retrospectively) |
| 4 Running | PASS at one loop (κ/κ_max 0.70 → 0.68 at \|X\|_inf; λ_eff > 0) [C] |
| 5 Light sector | Candidate 1 unchanged; c6 = 1.25e6 unchanged by Planck operators to 0.2 %; m_A = 3.2e6 GeV, vector decay at 12 TeV [C] |
| Instant release | finite Γ_φ favourable (R 1.14 → 0.03, \|ε\| 0.146 → 0.38) [C, two codes, one release model] |
| Scaling symmetry | CONFIRMED to 1e-12 [C] |
| Case-L premise (kinetic), lattice/fragmentation, backreaction, two loops, preheating, holder origin | OPEN / NOT CHECKED |

## 5. Corrections

**For Codex 2908e91.**
1. Build the envelope R from the actual prepared r_i: the leading radial relation is 5–8 % low where the β6 companion
   matters, and the envelope |X| ≤ R/a then fails by 8 % in amplitude (hw_envelope).
2. The bound is a factor ≈ 2.6 estimator (triangle inequality 98 % saturated; slack ⟨|X|⁶⟩/R⁶ ≈ 0.4), not only an
   inequality (hw_envelope).
3. State the fixed-coupling conclusion with its domain: on the line η fixed, r ≳ 0.3 and h/H_i ≥ 20 the same bound and
   the full release reach Y_Q inside the matching regime (every diagnostic < 0.05) (hw_envelope).
4. Point (i) becomes dynamical only with an adiabaticity condition: the forced Y carries the torque with ≤ 0.35 free
   quanta when the preparation or holder changes on time scales ≫ 1/M_Y; a quench faster than ≈ 0.3 H⁻¹ leaves 50–78
   quanta, and Codex's prescribed preparation is that quench (hw_gates, hw_verify, this track).
5. The processed late quantity Codex asks for is, at the heavy point, dominated by a species the packet does not name,
   the inflationary longitudinal vector, not by Y or Σ (hw_gates, this track).

**For our own statements.**
1. R_L = 88 at the spot → 81.5–82.0 (direct proxy); "passes by a factor of about 80" (ls_scaling) → ×71 (direct proxy),
   ×6–7 (ke), ×17 (ke + s) for the mediator debris, and ×1.2–1.6 (proxy) / **fails** (ke) once the inflationary vector is
   added at e = 1e-7.
2. ls_scaling's relative-phase mass 1.84 H → 2.53 H.
3. ls_scaling's "Y-relaxed is what a heavy Y does during inflation" → true for held phases and a holder outside the Y
   equation, or for a Y-coupled holder switched off adiabatically; the tolerance is |Y − Y_*| ≤ 0.37 H.
4. ls_scaling's "a larger e cuts the deposit as e^{−1/2}" → the vector piece falls as 1/e.
5. hw_gates' gate-1 "PASS by 1.2–1.8" → under the isentropic proxy only (0.61–0.82 of the direct R_L); under the
   energy-conserving proxy it fails by 6–9 at e = 1e-7 and needs e ≳ 3e-6.
6. hw_gates' "the holder must stay on until ≲ 0.7 e-folds before the end" rested on a zero-velocity reading of the
   phase; carrying the charge generated during the unheld stage, ε changes by 10 % within 0.05 e-fold and by 60 % at
   0.5 e-fold, while |ε| ≥ 0.1 holds to ≈ 1.2 e-folds and R ≤ 1.5 throughout.
7. ls_critic's "the window narrows to 20–43 under ke" → for the mediator debris; with the vector at e = 1e-7 the ke
   window is closed, and re-opens at 24–43 (e = 1e-6) or 21–43 (e = 3e-6).
8. hw_verify's "phase-sector mass² 1.0 and 5.0 H²" are diagonal entries of the 2×2 phase block, not masses; the
   projected curvature along δ at the relaxed state is −0.83 H².
9. hw_verify's "Codex preparation fails at every phase" → at every phase sampled (5° steps, 20°–145° where R_L is
   defined) at the heavy point.
10. The earlier scan-point review wording (REVIEW_SCANPOINT.md / sp_critic CONSOLIDATED, cited by the owner as the
    2793dab wording) "the one-mediator link fails at these couplings for every supplied preparation" overstated: only
    two preparations (Codex's and Y-relaxed) and six phases were tried there. Supported, as the owner corrected: **no
    published scanned release passes at the stated couplings**. Not supported and not stated: "every possible initial
    phase fails" (at the changed ratios the same field content passes the mediator-debris test for every sampled phase
    with |ε| ≥ 0.068, and the stated-coupling phase scan covered six values) and "the field content is excluded" (the
    window is the counterexample within the declared domain).

## 6. Codex's envelope-bound statements that are consistent with the window

All of them. (a) The bound and its table are reproduced and the spot point satisfies the bound (2.65–4.07). (b) The
r_i^{7/2} scaling and the extrapolation's exit from matching hold on Codex's line (M_Y, h, μ, κ fixed, H_i varied).
(c) "Excludes neither all preparations nor changed couplings" is exactly the window: the same bound on the line η fixed
with r 0.0707 → 0.45, κ ×3.46, M_Y/H_i 70.7 → 30, h/H_i = 40 allows 2.65–4.07 Y_Q, and the full release gives 1.00.
(d) "A nonzero source does not imply residual heavy quanta" has a dynamical instance (1.3e-4 quanta per unit charge).
(e) "No published scanned release passes" and "neither all preparations nor the field content has been excluded" are
both confirmed. (f) "A proof restricted to sudden prepared radiation releases should retain that restriction": the
window's verdict survives finite release, favourably. (g) "An integral of bare heavy-field energy is not the late
quantity": agreed, and at the heavy point the late quantity is the vector ledger. The one phrase that needs a domain
label is "the quiet limit is quantitatively weak" → "at the fixed couplings".

## 7. The single most decisive next calculation

**A kinetic case-L threshold at the heavy point.** The window's only failing gate (the inflationary vector, R = 49–65
at the stated e) sits between the two equilibrium proxies (energy-conserving 7–20; isentropic 80–88), which differ by
10× and decide whether the stated e = 1e-7 passes or e ≳ 3e-6 is required. A Boltzmann treatment of the hidden X–X̄
plasma from onset (T_h ≈ 12 m) with energy conservation, the actual X̄ number, the threshold 4→2 rates, the on-shell s
route and lump formation kinetics (not instant condensation), evaluated at ε = 0.146 and Y_Q for injected R from 1 to
90, gives R_L directly. Runner-up, if the energy-conserving value stands: release-induced longitudinal-vector
production at e = 3e-6–1e-5 with longitudinal Floquet modes and backreaction.

## 8. Open, not certified

- Both case-L thresholds are equilibrium proxies; the D relic, Boltzmann statistics for the s route and lump kinetics
  are missing; the ke runs have root fallbacks at R ≤ 6 and a 1–2 % energy-work residual.
- Homogeneous, classical, tree level everywhere: no lattice, fragmentation, finite-k mediator dynamics or backreaction;
  one loop without thresholds; one finite-release model; instant reheating in the ledger.
- The phase value and the holder (its CP phase, size and switch-off) are supplied by hand; no mechanism selects them.
- Piece A spans 0.76–1.0 of the GMR formula; longitudinal Floquet modes and backreaction at e ≥ 3e-6 not computed;
  string formation argued, not simulated; GW limits loop-model dependent and retrospective.
- The vacuum density remains underived; Candidate 1 is paused; no cosmological point is claimed.

## Reproduction (this track, ~8 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/hw_critic/code
cd $T; python3 -B run_checks.py ../results/checks.json; python3 -B run_holder.py ../results/holder.json
python3 -B thresholds.py ../results/thresholds.json; python3 -B window_edges.py ../results/window_edges.json
python3 -B summary.py ../results
```
Inputs reused unchanged: the hw_verify ladder, hw_gates' vector numbers, ls_critic's ke runs (all [Q, own tracks]).
Credits: Dine–Randall–Thomas; Affleck–Dine; Graham–Mardon–Rajendran; Starobinsky–Yokoyama; Kusenko–Shaposhnikov;
Enqvist–McDonald; Kasuya–Kawasaki; Carlson–Machacek–Hall; Kolb & Turner; Codex 2908e91 and b448f25.
