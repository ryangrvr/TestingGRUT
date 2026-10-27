# en_neutral: the same-seed neutral comparison, its conventions, and why the nonlinear run suppresses φ

27 October 2026. Claude, track en_neutral, reviewing Codex 7d9c93b (`grut_early_nonlinear_20261027/`), issue #13.

Inputs were read only through git 7d9c93b: MATCH_L12.npz, FIELD_CURVES.npz, LINEAR_PHI_CONTROL.json and
early3d.cpp. They were extracted with `git show` / `git archive` to `$SP/review/en_neutral/`. Note: the shared tree
`$SP/e7d` is no longer the published state. en_rerun has regenerated LINEAR_PHI_CONTROL.json, FIELD_CURVES.npz,
RUNS.json and MATCH_L12.npz there, so every number below uses the git copies.

Rules: [BRIEF.md](../../BRIEF.md), [BRIEF_EARLY.md](../../BRIEF_EARLY.md).
- Everything is conditional on Codex's reduced, conformal, classical X/φ model.
- Nothing here validates or excludes anything.
- Nothing was committed.

## What is independent and what is reproduction

| File | What | Label |
|---|---|---|
| `code/en_lib.py` | Own pieces: the Fourier-realization generator written from MATHEMATICS.md, the real-space FD4 energy (np.roll stencils), the RK4 for the orbit and the δφ transfer matrices, the turning points and radial period by quadrature, and the Floquet monodromy | independent. The only borrowed ingredient is the documented random-key recipe, so the *same* realization is drawn |
| `code/fixed_orbit.py` → `results/fixed_orbit.json` | Codex's fixed-orbit control recomputed, in realization, ensemble and zero-point-subtracted versions, plus the Floquet exponent on the matched reduced orbit | independent |
| `code/indep_nl.py` → `results/indep_N32_*.json` | Own numpy nonlinear lattice: FD4 via np.roll and a fixed-step 4th-order composition. Runs the N32 history and the counterfactual without X seed noise | independent |
| `code/run_codex_engine.py` → `results/reproduction_N48_run_info.json` | Codex's unmodified early3d.cpp, compiled from a git-archive copy and run on an initial file from our own generator. Snapshots every 2 τ (32-bit); the trace is the engine's own, at Δτ = 0.05 | **reproduction** |
| `code/floquet_periods.py`, `code/hartree_control.py`, `code/mechanism.py`, `code/xspec.py` | Own analysis of that reproduction history: period-by-period Floquet exponents, same-seed linear φ on the *measured* masses, φ spectra against the linear control, and the X fluctuation spectrum | independent analysis of reproduction data |
| `code/conventions.py` → `results/conventions.json` | Arithmetic on Codex's published curves | own arithmetic |

Units: Codex's H_i code units. s is conformal time, with τ(s) = 20 + a_m s + C s²/2. k is the FD4 lattice momentum √(symbol). The matched band means 6.74 ≤ k ≤ 7.27, the half-maximum band of the matched-orbit Floquet exponent (1291 of the 7384 seeded mode pairs at N48).

## 1. Conventions: what exactly is compared

From the code (`analyze.py`, `linear_phi_control.py`, `early3d.cpp`) and MATHEMATICS "Active energy and exchange":

- **Nonlinear curve.**
  - Numerator: E_φ + E_mix = K_φ + G_φ(FD4) + λ_φ⟨φ⁴⟩/4 + (g²/2)⟨φ²|χ_X|²⟩, with the local |χ_X|².
  - Denominator: ρ(t), the total active energy. It is conserved to 1.3e-8, so it equals ρ₀.
- **Linear control.**
  - Numerator: Σ_modes ½(|π_k|² + (k_lat² + g² r_orbit(s))|φ_k|²). It uses the same N48 discrete modes and the same FD4 k_lat, with the homogeneous reduced orbit r = |χ̄|² in the interaction term.
  - Denominator: `trace['rho'][0]` of the N48 run, so ρ₀ is identical.
  - Both curves are raw, with no vacuum subtraction.
- **Residual differences.** The control has no λ_φ φ⁴ term (≤ 1e-6 of ρ₀), no φ zero mode, and uses the homogeneous rather than the local |χ|². At τ 20 the two curves are 0.08032% and 0.08038% of ρ₀. **Verdict: conventions identical [C].**
- **The seed is not negligible [C, Codex curves].**
  - At τ 20 the seeded φ + interaction energy is 0.079–0.083% of ρ₀ (0.128% for cut10). That is 25–35% of the τ 800 values; for N48 it is 0.0804%, or 26% of 0.3098%.
  - In an ensemble occupation convention (instantaneous-ω half quantum), about 0.042% of ρ₀ is zero-point.
  - The net rise over the seed is 0.198–0.243% in all seven histories. cut10's higher final value (0.365%) is its larger seed (rise 0.237%).
- **Reconciling with our critic's "φ/X = 1 at τ 665–697" [C].** These are the same fixed-background physics in two conventions.
  - The critic used Codex's 2cbc7f1 control: full-model background, continuum k from 1e-5 to 40, ensemble covariances, zero-point subtracted, X-energy denominator. It gave 2.22 / 2.17 at τ 800 and a unit crossing at 696.5–697.0 on (3,.5).
  - Codex's box control, as an ensemble and zero-point subtracted (own code), gives **2.18 at τ 800 and a unit crossing at 698.5**.
  - The single N48 realization gives 2.09 and 702.5. That is 4% below the ensemble, which is realization scatter of the narrow band.
  - Subtraction is irrelevant late on (2.1837 versus 2.1831), and the denominators differ by 0.2% (E₀/ρ₀ = 0.99778).
  - The (5,1) value of 665 belongs to a history not run nonlinearly here.

## 2. Independent fixed-orbit control and Floquet rate

Own generator: ρ₀ = 105467.43730088249, against Codex's ...255 (relative difference 5.5e-16); nQ = 8113.477086664831.
Own RK4 at h = 1e-3: orbit energy drift 1.2e-11, transfer determinant error 1.9e-9.

| Quantity | Ours | Codex | Verdict |
|---|---:|---:|---|
| φ + interaction / ρ₀ at τ 800, same realization | **2.0917128** | 2.0917129 | CONFIRMED (max relative difference over all 1561 samples 2.4e-8) |
| Maximum / first sampled unit crossing | 2.10944 / 702.5 | 2.10944 / 702.5 | CONFIRMED |
| Ensemble expectation of the same control | 2.184; unit at 698.5 | n/a | NEW |
| μ_φ on the matched reduced orbit | **0.2062 per s at k = 7.01**; half-maximum band 6.74–7.27; half-grid 0.2062 | n/a (earlier reviews 0.202 / 0.206 at k ≈ 7.07, full model) | CONFIRMED |
| Late growth fit of the control, τ 300–800 | 0.192 per s | n/a | consistent |

The matched zero mode is very eccentric: E₀/E_circ(Q₀) = 8.60, r_max/r_min = 8.06, radial period 0.8145 in s, and
E₀/ρ₀ = 0.9978. At τ 800 the control's energy sits at k_lat ≈ 6.98 (energy-weighted).

## 3. Why φ stays near 0.3% in the nonlinear run

**Nonlinear number [C].**
- Reproduction (Codex engine, our initial file): N48 final 0.0030979302403044, against the published ...3043476. Energy drift 1.254e-8; charge drift 2.2e-14.
- Independent code: N32 final **0.27733026%**, against Codex's 0.27733222% (7e-6 relative). E_φ is 0.18821366% against 0.18821514%. Energy drift 1.35e-8.

**(a) Pump depletion (the zero mode's eccentricity decays): REFUTED as the main cause [C].**
- The coherent mean stays an eccentric pump. At mid-period samples:
  - r_max/r_min goes 8.06 (τ 20) → 7.6 (τ 60–70) → 4–6 (τ 100–200) → 2.7–3.7 (τ 300–800);
  - E₀/E_circ goes 8.6 → 2.1–3.0;
  - |χ̄|² still has a modulation depth of 0.88–0.97 throughout.
- Period-by-period Floquet on the free orbit with the zero mode's own (E₀, Q₀) gives μ_max = 0.206 → 0.147–0.150. The band drifts 7.0 → 5.9 as the mean loses energy.
- Same-seed linear φ driven by the *measured* zero-mode pump g²|χ̄(s)|² alone reaches **2.71 × ρ₀ at τ 800**, more than the fixed-orbit 2.09.

**(b) Detuning by the X fluctuation mass g²⟨|δχ_X|²⟩: CONFIRMED as the dominant mechanism [C].**
- Period means of g²S, with S = ⟨|χ_X|²⟩ − |χ̄|², compared with the band's k² ≈ 49:

  | τ | 26–35 | 57–70 | 70–84 | 84–100 | 100–317 | 317–762 |
  |---|---:|---:|---:|---:|---:|---:|
  | g²S | 0.2 | 2.8 | 6.8 | 14.9 | 24–33 | 36–48 |

- Floquet exponent, per period, on the measured homogeneous mass g²⟨|χ_X|²⟩(s) + 3λ_φ⟨φ²⟩(s):
  - in the matched band: 0.205 (τ 26–57) → 0.185 (57–70) → 0.099 (70–84) → **0 from τ 84 onward**;
  - at any k: ≤ 0.10 (84–100), ≤ 0.06 (≥ 118), ≤ 0.035 (≥ 258), ≤ 0.013 (≥ 418). The residual bands sit at k ≈ 3–5.
- Accumulated amplitude e-folds after τ 84: 0.18 at any fixed k on the measured Hartree mass, against 1.59 for the free zero-mode orbit.
- **Same-seed linear φ on the measured homogeneous mass gives 0.194% at τ 800** (mean over τ 700–800: 0.196%). Compare the nonlinear 0.31%, the fixed orbit's 2.09, and the zero-mode-only 2.71.
- Replacing S by its one-period running mean gives 0.20%. So the *mean* shift matters, not S's oscillation.
- Timing:
  - X fluctuations reach half the variance at τ 69.2 and charge contrast reaches 1 at τ 60.6 (Codex; reproduced);
  - S/⟨|χ|²⟩ is 0.09 (τ 60), 0.20 (τ 80) and 0.48 (τ 100) in 10-τ means;
  - the shutoff falls at τ ≈ 70–84.

**(c) The band leaves the lattice or the seeded band: REFUTED [C].**
- The residual Hartree bands, at k ≈ 3–5, lie inside both the lattice and the seed band.
- The nonlinear φ field is populated at every lattice k.
- S is carried by high-occupation infrared modes, not by lattice UV. At τ 60–200, 87–99% of it is in k_lat 1.8–4. At τ 400 and 800, 56% and 73% are in k < 1.8 (box-scale modes). The share above k = 8 stays at or below 0.08%.
- S/⟨|χ|²⟩ agrees across N32/48/64, seed2, cut10 and box16 to within about 0.04 (10-τ means; N32/48/64 within 0.02).

**Inhomogeneous scattering adds φ energy rather than suppressing it [C].**
- The nonlinear run has 0.31%, against 0.19% for the homogeneous-Hartree control.
- φ energy migrates to the unseeded k_lat > 8 range: 0% at τ 20, 26% at τ 140, 64% at τ 790.

**φ's own spectrum: the band growth stops around τ 70 [C].**
- Nonlinear energy in the matched band (10-τ running mean, fraction of ρ₀) tracks the linear control within about 5% until τ 66 (ratio 0.95). It grows from 2.1e-4 to 6.3e-4.
- From there it falls behind: ratio 0.92 at τ 70, 0.87 at τ 76.
- It peaks at about 1.0e-3 near τ 90 (linear: 1.4e-3), then declines to 1.6e-4 at τ 800, below its starting value. The linear control reaches 2.06 by then.
- The mean band occupation (Hartree ω) peaks at about 5.3 (τ 88–100) and falls to 0.9–1.0 by τ 640–800.

**Counterfactual: no X seed noise (independent code, N32) [C].**
- φ follows the linear resonance: 0.43% at τ 200, 3.7% at τ 400.
- X then fragments from fluctuations sourced by φ (S/⟨|χ|²⟩ = 0.59 at τ 500), and φ saturates at **4.3–5.9%** of ρ₀ (maximum 5.9% at τ 778).
- So in this model the 0.3% result requires the charged field's own seeded fluctuations to fragment the condensate first. The X instability (μ_X ≈ 0.63 [Q], our earlier review) is three times faster than μ_φ ≈ 0.21.

**Literature.**
- Classical-statistical treatment of amplified fluctuations: Khlebnikov–Tkachev, and Micha–Tkachev for the turbulent cascade [R].
- Its UV limits: Epelbaum–Gelis–Wu, arXiv:1402.0115 [Q].
- Termination of parametric resonance by backreaction through ⟨χ²⟩ and by rescattering is standard in preheating, for example Kofman–Linde–Starobinsky 1997 [R]. What is computed here is the size of that effect in this pilot.

## 4. Scope of "nonlinear continuation changes this forecast in the pilot"

**CONFIRMED as scoped [A].** Codex says "in the pilot" and "not a full-EFT exclusion". Our numbers support it, and they supply the mechanism: X fragmentation detunes the δφ resonance after τ ≈ 70–84. The cause is not depletion of the coherent pump.

Suggested additions:
- The same-seed fixed-orbit curve and our earlier τ 665–697 crossing are one linear extrapolation in two conventions.
- About a quarter of the 0.3% is the classical seed itself.
- If X fragmentation is delayed (no X seed), φ saturates near 4–6% instead.

What remains open:
1. **Heavy and gauge fields.** Y/Σ/gauge dynamics after τ 20. Gauge charge can change the fragmentation and hence S, which sets the detuning.
2. **Quantum effects.**
   - The seed's half-quanta are 25–35% of the final value.
   - 64% of the late φ energy lies above k_lat = 8 in unseeded modes populated by the cascade. That part depends on lattice and seed; the N32 value is 0.277%.
   - The detuning itself is carried by high-occupation X modes, so it is less exposed.
3. **Non-conformal physics.** The reduced model has no masses, so S and |χ̄|² never redshift relative to each other. In the full model, m_X, v, the φ potential and bias, and the later φ transition and wall history can change the ratio between detuning and pump at larger a.
4. **Box.** Late S sits in the lowest box modes. box16 gives a similar φ value (0.323%), but the infrared is finite.
5. **Single history.** Only (3,.5) was run nonlinearly.
6. **Energy, not populations.** These are energy fractions, not particle numbers.

## Verdicts

| Codex claim / item | Verdict | Ours | Codex | Method |
|---|---|---|---|---|
| [C] φ + interaction 0.28–0.36% at τ 800 | CONFIRMED | 0.2773–0.3648% (published); N48 0.30979302403% reproduced to 3e-15; independent N32 0.27733026% | 0.28–0.36% | reproduction + independent |
| [C] φ Hamiltonian alone 0.19–0.27% | CONFIRMED | 0.188–0.267% | 0.19–0.27% | own arithmetic |
| [C] Same seed on the fixed orbit gives 2.09 × initial active energy | CONFIRMED | 2.0917128 (relative difference 2.4e-8); unit at 702.5; ensemble 2.18 / 698.5 | 2.0917; 702.5 | independent |
| [A] Same conventions in both curves (interaction included, no subtraction, same momenta) | CONFIRMED | same ρ₀ and modes; at τ 20, 0.08032% versus 0.08038%; φ⁴ ≤ 1e-6 | stated | code + independent |
| [C] The seed's share of 0.28–0.36% | NEW | 25–35% (N48 26%); rise over seed 0.20–0.24% | not stated | own arithmetic |
| [A] Earlier τ 665–697 crossing versus 2.09 | NEW (reconciled) | same physics; ensemble/subtracted box control 2.18, unit at 698.5, against the critic's 2.22/2.17 and 696.5–697 | "fixed-background diagnostic" | independent |
| [C] μ_φ ≈ 0.20 at k ≈ 7 on the matched orbit | CONFIRMED | 0.2062 at k 7.01; band 6.74–7.27 | n/a | independent |
| Mechanism (a): pump depletion | REFUTED as cause | free-orbit μ 0.147–0.206; zero-mode-only pump → 2.71 × ρ₀ | n/a | own analysis |
| Mechanism (b): fluctuation (Hartree) detuning | CONFIRMED (dominant) | g²S 7–48 after τ 70; matched-band μ → 0 from τ 84; measured-mass linear φ 0.194% | n/a | own analysis |
| Mechanism (c): band leaves lattice/seed | REFUTED | residual bands at k 3–5, inside; S in IR modes (≤ 0.08% above k 8) | n/a | own analysis |
| [C] φ band growth before τ ≈ 70, then stop | NEW | within about 5% of linear to τ 66; peak ≈ 1.0e-3 ρ₀ near τ 90; 1.6e-4 by τ 800 | n/a | own analysis |
| [C] Counterfactual without X seed | NEW | φ saturates at 4.3–5.9% | n/a | independent |
| [A] "Changes this forecast in the pilot"; not a full-EFT exclusion | CONFIRMED (scoped) | see §4 | stated | argument |

## Limits

- The mechanism tests use homogeneous (Hartree-type) proxies built from the measured trace. Inhomogeneous effects are measured only through the gap between the nonlinear result and the homogeneous-Hartree control (0.31% versus 0.19%).
- One seed and one lattice (N48) for the mechanism analysis. The independent nonlinear check and the counterfactual are N32 only.
- The fixed-step RK4 transfer on the trace grid is validated: the free-orbit control on that grid gives 2.0909, against 2.0917 with the fine step.
- A first-pass finite-time Lyapunov estimate (ln σ_max / W) was biased to 0.4–0.7 by the swing within each period. It was dropped, and no conclusion uses it.
- Everything is classical-statistical within the reduced model. No physical population is implied.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; K=$SP/review/en_neutral
T=research/reality_framework/grut_formation_review_20261025/tracks/en_neutral; D=research/reality_framework/grut_early_nonlinear_20261027
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2; mkdir -p $K/published $K/pkt
for f in LINEAR_PHI_CONTROL.json FIELD_CURVES.npz MATCH_L12.npz; do git show 7d9c93b:$D/$f > $K/published/$f; done
git archive 7d9c93b $D | tar -x -C $K/pkt
python3 -I $T/code/conventions.py $K/published $T/results/conventions.json                       # seconds
python3 -I $T/code/fixed_orbit.py $K/published $K/fixed_orbit.json 1e-3                         # 45 s
python3 -I $T/code/run_codex_engine.py $K/published $K/pkt/$D/code/early3d.cpp $K/runs/N48 48 12 2 800 2   # 3.5 min, 1 GB
python3 -I $T/code/floquet_periods.py $K/runs/N48 $K/floquet_periods_N48.json                  # 45 s
python3 -I $T/code/hartree_control.py $K/published $K/runs/N48 $K/hartree_control_N48.json     # 40 s
python3 -I $T/code/mechanism.py $K/published $K/runs/N48 $K/mechanism_N48.json                 # 80 s
python3 -I $T/code/xspec.py $K/runs/N48 $T/results/xspec_N48.json
python3 -I $T/code/indep_nl.py $K/published $K/indep_N32_full.json 32 0.004 800 full           # 2.5 min
python3 -I $T/code/indep_nl.py $K/published $K/indep_N32_noXnoise.json 32 0.004 800 noXnoise   # 2.5 min
```

The JSON files in `results/` are copies of these outputs (the npz side files stay in scratch).
