# Track hw_critic: notes (consolidation of hw_verify, hw_envelope, hw_gates with ls_scaling, ls_critic, ls_threshold; Codex 2908e91)

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication. The verdicts, corrections and the
answer to Codex's question are in [CONSOLIDATED.md](CONSOLIDATED.md); this file records method, numbers, limits and
reproduction. Codex 2908e91 (`origin/codex/grut-concept-checkpoint-20261025`) was read only through `git archive` into
the scratchpad; no Codex code was run. Nothing was committed; only this directory was written.

**Scope.** One supplied tree-level working model (homogeneous classical release, supplied CP phase, late-decay ledger,
case-L proxy at f = 0.01). Nothing here is evidence that X, Y, Σ or the vector exist; it is not a validation, an
exclusion or a viable point. Labels: [D] derived, [C] computed, [E] estimate, [Q] quoted, [A] argued, [R] recalled.

## Method

1. Read the NOTES.md of hw_verify, hw_envelope, hw_gates, ls_scaling, ls_threshold and ls_critic/CONSOLIDATED.md,
   their results JSONs, Codex 2908e91 (CLAUDE_HANDOFF, REPORT, SOURCE_CHECK, CHECKS.json), the owner's scope
   corrections in REVIEW_LOWSCALE.md, and the earlier wording in REVIEW_SCANPOINT.md / sp_critic/CONSOLIDATED.md.
2. Cross-checked every number that two or more tracks share (CONSOLIDATED §1).
3. Recomputed the disputed or single-track load-bearing numbers with own code (`code/`):
   - `hwc.py` [C]: own compact model of the spot point (H_i = 1 units): potential, complex gradient, own minimisers for
     the aligned and Y-relaxed preparations, finite-difference canonical 6×6 Hessian with the exact gauge mode
     projected out, de Sitter stage with the two phases held or free, release in ln a as the independent variable
     (scipy DOP853, rtol 1e-10, atol 1e-12) in a radiation background or with an inflaton of width Γ_φ decaying into
     radiation (X Hubble mass ∝ ρ_φ), the free-Y split about the forced response, the Σ radial count, the
     closed-form decay ledger [D, re-derived], the prompt-radiation yield and the sp_critic case-L table [Q].
     Four other integrators exist for this point (ls_scaling cosmic-time DP5 C; hw_envelope Python DOP853 cosmic time;
     hw_verify conformal Yoshida-8 C; hw_gates DP5 C); this is a fifth, in a different independent variable.
   - `run_checks.py` [C]: spot point (both preparations, window 30–60 and 30–100); held-phase relaxation from the
     Codex state (N = 0.5, 1, 2 e-folds); the disputed holder-off drift: free de Sitter evolution from the Y-relaxed
     state for N = 0.05–2 e-folds, then release carrying the acquired velocities (hw_verify's method) and with the
     velocities zeroed (hw_gates' zero-velocity reading); finite release Γ_φ/H_i = 10, 1, 0.1.
   - `run_holder.py` [C]: hw_gates' D1 holder (μ → μe^{iδ} during inflation) switched off instantly, exponentially over
     Δa = 0.1, 0.3, 1 and ∝ H (μ_P ∝ a⁻²).
   - `thresholds.py` [C, arithmetic on quoted track numbers]: the full late debris (mediator ledger + the inflationary
     longitudinal vector of hw_gates) against every threshold variant on record, as a function of e with hw_gates'
     scalings (R_A ∝ 1/e, R_Σ, R_casc ∝ e^{−1/2}); an independent conversion of the Graham–Mardon–Rajendran formula to
     a comoving vector number per H_i³; vector headroom per net charge under each threshold.
   - `window_edges.py` [E]: lower edge m of the window for each threshold, from the hw_verify ladder R(m) (power-law
     interpolation; constant R_L at the spot ε, so 0.5 lower in m than hw_verify's ε-dependent 12.1–12.3).
   - `summary.py`: prints every quoted number from `results/*.json`.

## Numbers (own code; comparison values [Q, own tracks])

**Spot point** (H = 1.3335e11 GeV, M_Y = 30 H, h = 40 H, μ = 7.1595 H, κ = 3.4587e-3, c_H = 10, δ = π/2, Y-relaxed):
aligned minimum |X|, |Y|, |Σ| = 80.2021, 3.44245, 28.3107 H (hw_gates identical); relaxed 76.39978,
(−2.11151, −0.93933), 28.29314 (hw_envelope identical); ε = −0.14568, a³n_X = −2788.17, E_X = 38399.6,
N_YH = 0.01725, N_YL = 0.3345, N_σ = 3.2e-9, R_Y = 0.805, R_Σ = 1.2e-7, R_casc = 0.339, **R = 1.144**, yield 1.0000 Y_Q,
max a|X| = 79.96 (hw_envelope 79.97), Gauss residual 4.6e-12; α6 = 3.750e-9 H⁻² = 2.109e-31 GeV⁻², c6 = 1.2504e6.
Window 30–100 control: ε −0.14567, R 1.1440. Codex preparation: ε = −0.15883, N_YH = 76.81, R = 1415.
All five integrators agree to ≤ 0.03 % in ε and ≤ 0.1 % in R.

**Phase mass.** Gauge-projected 6×6 Hessian at the aligned minimum: m² = 6.4057, 14.77, 572.9, 1334.3, 3209.9 H²
(gauge eigenvalue 9e-15); light angular mode m = **2.531 H**, f_eff = 17.10 H (hw_gates 2.53, 17.1: CONFIRMED;
ls_scaling 1.84 H [D]: CORRECTED). Along δ at the Y-relaxed state the projected curvature is −0.827 H² (hw_gates
−0.83); hw_verify's "phase-sector mass² 1.0, 5.0 H²" are the diagonal entries of the 2×2 phase block, a different
projection, not a mass.

**Held relaxation from the Codex state** (phases held, release with velocities): N = 0.5 / 1 / 2 e-folds → R = 384.8 /
88.3 / 7.87, ε = −0.138 / −0.142 / −0.1455 (hw_verify 385 / 88.3 / 7.87: CONFIRMED).

**Holder-off drift (the disputed item).** Free de Sitter evolution from the Y-relaxed state for N e-folds, then release:

| N (e-folds) | 0.05 | 0.1 | 0.2 | 0.5 | 1 | 1.5 | 2 |
|---|---|---|---|---|---|---|---|
| δ at the end (rad; starts at π/2 = 1.571) | 1.565 | 1.547 | 1.483 | 1.147 | 0.467 | 0.023 | −0.090 |
| hw_gates 1D drift [Q] | — | 1.548 | 1.489 | 1.169 | 0.476 | — | −0.095 |
| n_X at the end of the stage (H³; release charge is 2788) | −403 | −746 | −1293 | −2179 | −2096 | −925 | −3 |
| ε, velocities carried (own) | −0.163 | −0.178 | −0.204 | −0.239 | −0.185 | −0.067 | +0.039 |
| ε, hw_verify [Q] | −0.163 | −0.178 | −0.204 | −0.239 | −0.185 | — | +0.039 |
| ε, velocities zeroed (own) | −0.146 | −0.147 | −0.148 | −0.151 | −0.109 | −0.016 | +0.040 |
| ε, hw_gates scan reading [Q] | — | −0.146 | −0.148 | −0.151 | −0.104 | — | −0.028 |
| R_central (own, both variants) | 1.18–1.20 | 1.30–1.32 | 1.30 | 1.33–1.35 | 1.43–1.48 | 1.47 | 1.48–1.49 |

Both tracks are reproduced: the phase trajectory of hw_gates' 1D model is right (δ agrees to 2 %), and so is
hw_verify's ε(N). The difference is that the unheld phase rolls under the X⁶ torque during inflation and the charge
density it accumulates (−2179 H³ after 0.5 e-fold, 78 % of the release charge) is carried into the release; the
zero-velocity reading drops it. So: ε is set to 10 % only if the holder stays on to within ≈ 0.05 e-fold of the end of
inflation; |ε| ≥ 0.10 survives a switch-off up to ≈ 1.2 e-folds before the end; the debris is unaffected (R ≤ 1.5,
N_YH ≤ 0.016 throughout). hw_gates' "≲ 0.7 e-folds" is CORRECTED in that sense.

**Finite release** (Y-relaxed, inflaton width Γ_φ): Γ_φ/H_i = 10 / 1 / 0.1 → ε = −0.1720 / −0.2791 / −0.3804,
N_Y = 0.087 / 0.0107 / 0.0074, R = 0.204 (0.209 with the reheating ⟨a⟩) / 0.0322 (0.0370) / 0.0292 (0.0446),
diluted yield 1.244 / 1.159 / 0.545 Y_Q (hw_verify: identical to the printed digits; their 7.32 / 6.83 / 3.21 at the
ladder scale divide by 5.889). CONFIRMED: slower release gives more charge and less debris at this point.

**D1 holder switch-off** (own): instant R = 1415 (ε −0.159); exp Δa = 0.1: ε −0.1173, R 161; 0.3: −0.0923, R 11.1;
1: −0.0653, R 2.38; ∝ H: −0.0750, R 4.57 (hw_gates 1415 / 161 / 11.1 / 2.4 / 4.6: CONFIRMED).

**Thresholds at ε = 0.1457, normalised yield** (all [Q] from the named tracks; crossings log-interpolated here):
sp_critic table 87.85; direct proxy 81.5 (hw_verify, sp_yledger code) and 82.0 (ls_critic, lsrelic code; the 0.6 %
is the known endpoint difference); Codex calibration ≈ 80; energy-conserving (ke) proxy 8.26 (μ_c = E_c = 0.861 m),
6.97 (0.923 m); ke + on-shell s route 19.6. Vector headroom per net charge: 1.47–1.48e-3 (direct proxy),
1.06–1.30e-4 (ke), 3.4e-4 (ke + s); hw_gates' piece A is 8.85e-4 (own modes) to 1.17e-3 (GMR formula).

**GMR formula, own conversion:** Ω_A/Ω_DM = √(m_A/6e-6 eV)(H/1e14 GeV)² = 4.11e4 with ρ_DM/s = 4.4e-10 GeV gives
n_A/s = 5.65e-12 and N_A = **3.2598 H_i³** (hw_gates "3.26": CONFIRMED; their own Proca modes give 2.468 = 0.757 of it).

**Full late debris against the thresholds** (mediator ledger 1.144 at e = 1e-7 plus R_A = 48.5 [own modes] / 64.0 [GMR];
R_A ∝ 1/e, R_casc ∝ e^{−1/2} [D, hw_gates]):

| e | R_total | / direct proxy 81.5 | / ke 8.26 | / ke 6.97 | / ke + s 19.6 |
|---|---|---|---|---|---|
| 1e-7 | 49.6 / 65.1 | 0.61 / 0.80 | **6.0 / 7.9** | **7.1 / 9.3** | **2.5 / 3.3** |
| 3e-7 | 17.2 / 22.3 | 0.21 / 0.27 | 2.1 / 2.7 | 2.5 / 3.2 | 0.88 / 1.14 |
| 1e-6 | 5.8 / 7.3 | 0.07 / 0.09 | 0.70 / 0.89 | 0.83 / 1.05 | 0.29 / 0.37 |
| 3e-6 | 2.5 / 3.0 | 0.03 / 0.04 | 0.30 / 0.36 | 0.36 / 0.43 | 0.13 / 0.15 |
| 1e-5 | 1.3 / 1.5 | 0.02 | 0.16 / 0.18 | 0.19 / 0.21 | 0.07 / 0.08 |

Mediator debris alone (1.144) passes every variant: ×71 (direct proxy), ×6.1–7.2 (ke), ×17 (ke + s).

**Window lower edge in m = M_Y/H** [E]: mediator only 11.5–11.8 (proxy variants; hw_verify's ε-dependent value
12.1–12.3), 20.3–21.0 (ke), 16.9 (ke + s); with the vector at e = 1e-7: 14.4–15.3 (proxy), **closed** (ke, ke + s,
R_A alone exceeds R_L); at e = 1e-6: 24–27 (ke), 18 (ke + s); at e = 3e-6: 21–22 (ke), 17 (ke + s). Upper edge 42.8
(κ = κ_max) [Q hw_verify].

## Limits

- Everything is homogeneous, classical and tree level; the de Sitter stage holds the phases by fiat or not at all (no
  holder term except the D1 test); one finite-release model (perturbative inflaton decay, Hubble mass ∝ ρ_φ, instant
  thermalisation, reheating factors applied as (Rr/3)^{∓3/4, 1/4}); the field energy is neglected in H.
- The ledger and both case-L proxies are reused from our earlier tracks; this track adds no new threshold physics. The
  ke values rest on ls_critic's 13 + 10 runs (N = 500–1000, energy-work residual 1–2 %, root fallbacks at R ≤ 6 that
  mean only "erased"); ke is itself an equilibrium proxy, not a kinetic calculation.
- The vector piece, its e-scaling, strings, running and the light sector are hw_gates' numbers [Q]; only the GMR
  conversion and the threshold arithmetic were redone here. R_L for ε ≠ 0.1457 uses the sp_critic table (7–9 % high).
- Window edges use the spot-point R_L for every rung (the ε spread along the ladder is 0.113–0.146).
- No lattice, fragmentation, finite-k mediator dynamics, backreaction, two loops or preheating.

## Credit

Dormand–Prince 8(5,3) (scipy DOP853); Dine–Randall–Thomas Hubble-induced terms and Affleck–Dine [R]; Graham, Mardon &
Rajendran PRD 93 (2016) 103520 [R]; Starobinsky–Yokoyama equilibrium [R]; Q-ball formation and fragmentation
(Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki) [R]; Codex 2908e91 and b448f25 [Q]; our tracks ls_scaling,
ls_threshold, ls_critic, sp_critic, sp_yledger, hw_verify, hw_envelope, hw_gates [Q where cited].

## Reproduction (~8 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/hw_critic
cd $T/code
python3 -B run_checks.py   ../results/checks.json        # ~5 min, 23 releases
python3 -B run_holder.py   ../results/holder.json        # ~1 min, 5 releases
python3 -B thresholds.py   ../results/thresholds.json    # 1 s
python3 -B window_edges.py ../results/window_edges.json  # 1 s (reads thresholds.json)
python3 -B summary.py      ../results                    # prints every quoted number
```
