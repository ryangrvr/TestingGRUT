# Track gp_leftovers: symmetric X/X̄ leftovers at the 7737487 gauge point (handoff item 5 of b5174e5)

25 October 2026 · Claude (independent review track) · Codex packets read only through `git archive` of b5174e5
(`grut_gauge_point_test_20261025/`, `grut_z6_gauge_bridge_20261025/`, and the portal-debris ledger
`grut_portal_debris_window_20261023/` as it stands at b5174e5) into the scratchpad. Own code in `code/`, outputs in
`results/`. No Codex code was run. No commits.

**This is a conditional calculation inside one supplied model. It is not a GRUT prediction, and nothing here is
evidence that X, φ, D, Y, Σ or the inflationary inputs exist.** Same-session AI agents are not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 1. Inputs

- Light target [Q]: g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV, h_D = 0.01, m_D0 = 0.1 v, so
  m_X = 0.092409, m_s = 0.226355 and m_D = 0.092870 GeV.
- Reference object [Q]: Q = 200000, ω/m_X = 0.861195, E/(Q m_X) = 0.923316.
- Portal κ = 9.006e-10 [Q]. Gauge point [Q]: e = 1e-7, h = 1e11 GeV, so m_A' = 6eh = 6e4 GeV.
- Calibration [Q]: Y_Q = 4.8294e-9 and B = 4.4e-10 GeV per visible entropy. The ledger is the portal-debris witness:
  objects plus binding-energy scalars, with s → X X̄.
- Conversion [Q, bb211d1 REPORT §2]: σv_cold = 1.13039e-16 cm³/s and σv_daughter = 8.1225e-15 cm³/s on the cold
  Q = 200000 profile.
- Onset [Q, tracks/affleck_dine]: the φ transition starts when the condensate energy reaches ρ_c = 3λ_X X_c⁴, with
  X_c = 1.1318 GeV.
- Sibling results used for comparison [Q]:
  - gp_phase: net ⟨ε⟩ is consistent with 0 (|⟨ε⟩| ≤ 0.03 at 2σ), with local rms 0.14–0.23;
  - gp_matching: median local |ε| = 0.11.

## 2. Method (own code)

- `budget.py` [C, D]:
  - Reproduces the pair budget.
  - Derives the remainder law three ways: (a) the quoted 0.9(ε^-4/3 − 1); (b) the non-relativistic (NR) count
    f = (m_X Y_Q/B)(1/ε − 1); (c) an adiabatic map from the quartic-era fraction ε = n/n_max(ρ) to the NR count. The
    map uses the invariants a³n and a³J_r of the planar orbit, with J_r = (1/2π)∮p_R dR. A harmonic control is exact
    to 3e-14, and near ε = 1 it gives j = J_r/n = (1 − ε)/√6 [D].
  - Builds the transition-energy ledger: the sudden-onset energy (ρ_c + V0)/n_Q; the homogeneous fixed-charge φ
    potential E(φ²; n), which gives the order of the transition; and the zero-pressure coexistence point.
- `amp42.py` [C]: tree amplitudes of the light sector by Berends–Giele recursion, with vertices taken as exact sympy
  derivatives of U. Checks:
  - 2→2 X X̄, X X and D D → X X̄ agree with the analytic amplitudes to 1e-16;
  - charge-violating amplitudes vanish (3e-20);
  - the λ_X-only 4→2 amplitude equals the hand diagram sum 12λ_X²/m_X²;
  - the result is rotation-invariant.
- `channels.py` [C]: thresholds and Γ/H from T_vis = 86 GeV to 2.5 keV, with n_Q = Y_Q s.
- `catalysis.py` [D, C]: the charge budget of conversion, checked by a mass-action integration.
- `equilibrium.py` [C]: maximal-equilibrium benchmark.
  - The post-onset hidden sector is fully thermalized at fixed energy and charge, using Bose–Einstein thermodynamics
    with a condensed phase at chemical potential μ_c.
  - Case H: a homogeneous condensate, μ_c = m_X, with X̄ removed by 4→2 only.
  - Case L: lumps with μ_c = ω_Q and E/Q as in the reference object, with X̄ removed by 4→2 plus conversion at
    n_c/Q · σv_cold.
  - Evolution is isentropic, with X̄ relaxation integrated by exact exponential steps.
  - Controls (`eq_controls.py`): doubling the steps changes p_f by ≤ 0.8%; the 4→2 coefficients are scaled by ×10 and
    ×0.1.
- `vector_source.py` [E]: the gp_matching flag. The longitudinal A' is estimated with the GMR scaling, followed by
  A' → X X̄ and the case-H plasma model.

## 3. Results

### 3a. Codex's budget [C]

- Y_pair,max = δB/(2m_X) = 2.38072e-11 at δ = 1%, which is 4.92964e-3 pairs per net charge.
- At δ = 10% and 100% the values are 2.381e-10 and 2.381e-9.
- The quoted law gives ε_min = (1 + δ/0.9)^-3/4 = 0.991747. All of these match POINT_RESULTS.json.

The ledger already holds 0.031305 free pairs per net charge from binding debris. Calibrated matter per net charge
is 0.985928 m_X: objects 0.923316 plus pairs 0.062611. **The ledger contains no φ-transition (V0) energy.**

### 3b. Remainder law [D, C]

| | f_extra(ε = 0.9) | f(0.5) | f(0.1) | f(0.01) | ε_min at 1% | 10% | 100% |
|---|---:|---:|---:|---:|---:|---:|---:|
| quoted 0.9(ε^-4/3 − 1) | 0.136 | 1.37 | 18.5 | 417 | 0.99175 | 0.9240 | 0.5710 |
| NR count, 1.01427(1/ε − 1) | 0.113 | 1.01 | 9.13 | 100 | **0.99024** | 0.9103 | 0.5035 |
| adiabatic quartic → NR | 0.093 | 0.86 | 8.09 | 90 | **0.98808** | 0.8930 | 0.4640 |

The 4/3 law counts the excess quartic-era *energy* at a common time. It overstates the frozen remainder by 1.2–1.45×
near ε = 1, rising to 4.2× at ε = 0.01. The conclusion that a frozen remainder needs ε ≈ 0.99 is unchanged.

### 3c. Every channel for X X̄ at the light point [C unless marked]

| Channel | Status for cold pairs (√s = 2m_X = 0.1848 GeV) | Rate |
|---|---|---|
| X X̄ → s s | closed (needs 0.4527; KE 0.134 GeV per X) | — |
| X X̄ → s | closed (needs 0.2264; KE 0.0208 GeV per X); in a hot plasma s ↔ X X̄ only equilibrates (Γ_s→XX̄ = 1.73e-5 GeV) | — |
| X X̄ → D D | closed by 0.46 MeV per X; D is Z₂-stable and D D → X X̄ is open with σv = 4.78e-9 GeV⁻² | D stores dark energy and cannot drain it |
| X X̄ → A' | closed: m_A'/(2m_X) = 3.2e5; no light U(1)' charged state and χ = 0 | — |
| X X̄ → SM (s–h mixing θ = 1.31e-11) | open (e⁺e⁻, γγ; μμ closed) | σv = 2.2e-34 GeV⁻² (2.6e-51 cm³/s); Γ/H ≤ 3e-20 at all epochs; s → SM width 1.3e-32 GeV |
| 3→2 (X X X̄ → X + neutral) | forbidden: 2m_X < m_s and a single D is Z₂-odd | — |
| 4→1 (→ s) | forbidden kinematically | — |
| 4→2: XXX̄X̄ → XX̄, XXXX̄ → XX, XXX̄X̄ → DD | open | \|M\| = 0.0534, 0.0126, 1.0e-3 GeV⁻². A = \|M\|²Φ₂/(2m)⁴/S_f = 0.0842, 2.36e-3, 1.6e-5 GeV⁻⁸. s-exchange enhances \|M\| 38× over λ_X alone. Γ = H at T_vis = 6.5 GeV (X̄ in an asymmetric background), and 3.1 / 1.4 / 0.56 GeV for p = 1 / 10 / 100 X̄ per net charge |
| X X̄ bound states | none: α m_X/m_s = 3.2e-4 (α = 8e-4); a bound pair could not decay in the hidden sector anyway | — |
| Conversion X̄ + core(Q) → X + core(Q−2) (bb211d1) | open: E_out = 0.2516 GeV. The inverse needs E_X ≥ m_X + 2ω | n_core σv_cold/H > 1 down to T_vis = 0.62 keV if all charge sits in Q = 2e5 cores; costs 2 bound charges per X̄ |
| X̄ + core → core(Q−1) + s | closed: m_X + ω = 0.172 < m_s | — |
| X capture, X + core → core(Q+1) (+ s) | closed for 2 bodies; XX + core → core(Q+1) + X is open (3-body; rate not computed) | — |
| Object + anti-object | geometric σ = π(2r₉₀)², r₉₀ = 62.9/GeV | collides within a Hubble time for v_rel > 3e-13 (86 GeV) up to 5e-5 (2.5 keV). It turns bound pairs into free hidden quanta and erases nothing |

**No channel removes hidden-sector energy.** Annihilation into lighter states is closed, apart from a portal channel
that is 20 orders of magnitude too slow. Erasure can only be number-changing inside the hidden sector. Rest mass then
becomes kinetic energy that redshifts. The two channels that do this are 4→2 against a dense (Bose-condensed)
majority, and conversion on same-sign cores.

### 3d. Relic of the symmetric component versus ε

Here p_f is the frozen X̄ per net charge, and f = 2p_f m_X Y_Q/B.

1. **Frozen remainder, no thermalization** [C]: f is given by §3b. The 1% threshold is ε ≥ 0.990 (0.988 adiabatic).

2. **Conversion only, non-thermal bookkeeping** [D, C]. Take a fraction β of each sign bound in same-sign lumps right
   after fragmentation.
   - The minority is fully erased iff ε ≥ 1 − β. The final bound fraction is c_f = 1 − (1 − β)/ε.
   - Otherwise at least p_f = ((1 − β)/ε − 1)/2 remains.
   - Late matter is (c_f·0.9233 + 1 − c_f) m_X per charge. That is within 1% of calibrated for c_f ≥ 0.055, giving
     ε ≥ (1 − β)/0.945.
   - Mass-action control: exact. β is not computed for this model. bb211d1's prepared-collapse retention of 0.351 gives
     ε ≥ 0.686; β = 0.9 gives ε ≥ 0.106.

3. **Maximal equilibrium, homogeneous condensate, 4→2 only (case H)** [C]. The thermalized plasma starts at
   T_h0 = 2.5–3.8 m_X.
   - At onset, most of the charge is condensed for ε ≥ 0.05 (72% at ε = 0.1). While condensed, μ_X = m_X and
     n_X̄ ∝ e^(−2m_X/T_h).
   - At freeze-out the condensate survives only for ε ≥ 0.2.
   - Freeze-out is at T_vis = 5.6–6.2 GeV for ε ≥ 0.1, and 4.2 / 3.4 GeV at ε = 0.05 / 0.03.

   | ε | 0.7 | 0.5 | 0.3 | 0.2 | 0.15 | 0.1 | 0.05 | 0.03 | 0.01 | 0.003 | 0.001 |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | entropy per charge σ | 0.13 | 0.32 | 0.74 | 1.25 | 1.76 | 2.8 | 5.8 | 9.8 | 30 | 100 | 300 |
   | f | 1e-6 | 1.2e-4 | 2.2e-3 | 8.6e-3 | 0.018 | 0.050 | 0.44 | 1.23 | 5.4 | 20 | 60 |

   - The allowed range is **ε ≥ 0.19 at 1%** (0.12–0.30 for 4→2 coefficients scaled ×10 or ×0.1), ε ≥ 0.080 at 10%,
     and ε ≥ 0.033 at 100%.
   - The adiabatic-energy start changes f by ≤ 16% for ε ≤ 0.3. Larger relative changes occur only where
     f < 1e-3.
   - For ε ≥ 0.88, a sudden onset leaves no thermal component (energy deficit, §3e).

4. **Maximal equilibrium with lumps and conversion (case L)** [C]:
   - **erased (p_f < 4e-8) for ε ≥ 0.01;**
   - **not erased for ε ≤ 0.003** (p_f = 9.8 and 30, the same as case H, because lumps never form);
   - lumps persist throughout only for ε ≥ 0.1;
   - for 0.01–0.05 they evaporate and must re-form by solitosynthesis, which needs the uncomputed 3-body capture.

5. **Longitudinal A' source** [E, conditional on GMR scaling with the vacuum mass]. At H_inf = 4.7e13 GeV:
   - ρ_A/s = 0.31 GeV, and A' → X X̄ (Γ = 4.0e-12 GeV) at T_vis = 1.7 TeV;
   - that gives 1053 pairs per net charge;
   - the thermalized plasma has σ = 1.5e5, so p_f = 1.3e4 and f = 2.6e4, independent of ε.
   - f = 2.2 at 1e11 GeV. Reaching 1% needs H_inf ≲ 2e10 GeV [E]. That is below the minimum H_inf of the one-field
     c_H = 1 and 10 windows in tracks/affleck_dine.

**Against the sibling charge fractions**:
- Mixed microscopic domains make ε_eff = net ≤ 0.03. That fails options 1–3. Option 4 holds only above a threshold
  between 0.003 and 0.01, and only with solitosynthesis.
- Large unmixed domains (|ε| ≈ 0.17) are marginal in option 3 (f ≈ 1.3%) and leave matter/antimatter domains.

### 3e. Sudden-quench debris and the ledger [C, D]

- **Reproduction.** Codex's six values (2m_X/m_s)(λ_φv⁴/4)/s(T) = 0.3408–4.1914 are reproduced exactly.
- **Self-consistent onset.** With n_Q = ε n_c and n_c = 0.12967 GeV³ (T_on = 86.1 GeV at ε = 1), the debris is
  **0.378/ε** of calibrated matter.
- **The ledger does not contain V0, but V0 cannot be added to it.**
  - V0, the binding debris and the remainder all come from one budget. In a sudden onset that budget is
    (ρ_c + V0)/n_Q = **0.8672 m_X/ε per net charge**.
  - Ledger + quench (1 + 0.378/ε) exceeds that energy for **ε > 0.502**, so adding the two is double counting.
- **Energy deficit at high ε.**
  - For ε ≥ 0.880, the onset energy is below the calibrated 0.986 m_X per charge.
  - For ε ≥ 0.939 it is below the reference object's own E/Q (0.9233 m_X). Q = 2e5 objects cannot all form from a
    sudden ε ≈ 1 onset at the calibrated v.
  - E/Q(1e6) = 0.845 m_X [Q]. By log interpolation, Q ≳ 6e5 would fit [E].
- **Evolutions that avoid the debris:**
  1. Quasi-static phase separation at zero pressure. P = λ_X R⁴ − V0 vanishes at R_coex = 1.529 GeV, which is
     1.35 X_c, before the spinodal. Q-matter then coexists with the +v vacuum, and no V0 is released. E/Q = ω =
     0.740 m_X. This requires bubble nucleation (not computed).
  2. Homogeneous tracking. At fixed charge the transition is only weakly first order: d²E/d(φ²)² =
     λ_φ/2 − g⁴/(6λ_X) = −1.67e-3 at the spinodal, the jump is at R = 1.001 X_c, and it releases 1.5e-4 m_X per
     charge. φ then follows its minimum and V0 goes into X rest mass, which is the adiabatic law of §3b.
     tracks/affleck_dine finds finite-k growth faster than H near X_c, so fragmentation pre-empts this route.

## 4. Verdicts

| Claim | Codex | Mine | Verdict |
|---|---|---|---|
| Extra cold-pair budget at 1%: Y_pair ≤ 2.381e-11, 0.0049296 per net charge | stated | 2.38072e-11, 4.92964e-3 | CONFIRMED |
| ε ≥ 0.99175 "under the review's remainder model" | quoted law | arithmetic 0.991747 confirmed; the law itself overstates the frozen remainder; thresholds 0.9902 (NR count) and 0.9881 (adiabatic) | CONFIRMED (arithmetic); law CORRECTED |
| Our REVIEW_XCHARGE remainder 0.9(ε^-4/3 − 1) "unless hidden processes erase it" | — | erasure exists only through in-sector number changing, under the conditions of §3d | CORRECTED (self) |
| Quench debris 0.34–4.2× calibrated | 6 values | same to all digits; self-consistent 0.378/ε | CONFIRMED |
| Not a lower bound; adiabatic evolution may avoid it; the ledger may already include it; double counting is wrong | stated | the ledger does not include V0; the two are alternative uses of one onset energy, and adding them exceeds it for ε > 0.50; avoided by coexistence-point separation or homogeneous tracking (weak first order) | CONFIRMED in substance; "already included" CORRECTED |
| No symmetric population erased by assumption | open | channels computed (§3c–d) | CONFIRMED (not assumed) |
| A' channel closed if m_A' > 2m_X | task | closed by 3.2e5 | CONFIRMED |
| Number-changing / anti-object fate | NOT CHECKED | 4→2 coefficients, conversion budget, equilibrium relic (§3d) | NEW |

## 5. Limits

- Homogeneous thermodynamics. Case L treats lumps as a condensate (E/Q = ω_Q proxy).
- 4→2 uses threshold amplitudes without Bose enhancement; enhancement would favour erasure.
- s and D are left out of the plasma degrees of freedom. That raises σ by ≤ 19% [E] and shifts the case-H threshold
  to about 0.22 [E].
- Full thermalization at onset is an optimistic assumption. Coherent lumps and non-thermal cascades are not simulated.
- The conversion rate is the cold-ground-profile value; excited or larger lumps are not covered.
- Solitosynthesis (3-body capture) is not computed. β is not computed. Domain mixing is not computed.
- The GMR estimate ignores the changing Higgsing (X-dominated during inflation).
- The g* table is coarse (about 10%).

## 6. Credit

- Cannibal dark matter: Carlson, Machacek & Hall, ApJ 398 (1992) 43 [R].
- Solitosynthesis: Griest & Kolb, PRD 40 (1989) 3231; Frieman, Gelmini, Gleiser & Kolb, PRL 60 (1988) 2101;
  Postma, PRD 65 (2002) 085035 [R].
- Q-ball formation and Q/anti-Q at small ε: Kusenko & Shaposhnikov, PLB 418 (1998) 46; Enqvist & McDonald, PLB 425
  (1998) 309; Kasuya & Kawasaki, PRD 61 (2000) 041301 and PRD 62 (2000) 023512 [R].
- Bose condensation kinetics and turbulent thermalization: Semikoz & Tkachev, PRL 74 (1995) 3093; Micha & Tkachev,
  PRD 70 (2004) 043538 [R].
- Berends & Giele, NPB 306 (1988) 759 [R].
- Graham, Mardon & Rajendran, PRD 93 (2016) 103520 [R].
- Conversion channel: Codex bb211d1 [Q].

## 7. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/gp_leftovers/code
git archive b5174e5 research/reality_framework/grut_portal_debris_window_20261023/data | tar -x -C <scratch>
python3 budget.py              # ~40 s  -> results/budget.json
python3 amp42.py               # 1 s    -> results/amplitudes.json
python3 -I channels.py <scratch>/research/reality_framework/grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv
python3 catalysis.py           # -> results/catalysis.json
python3 equilibrium.py         # ~10 min, 4 cores -> results/equilibrium.json ('full' adds the adiabatic-energy start)
python3 eq_controls.py         # ~3 min -> results/equilibrium_controls.json
python3 vector_source.py       # ~8 min -> results/vector_source.json
python3 summary.py             # -> results/summary.json
```

`results/equilibrium_adiabatic_energy.json` comes from the first corrected pass. It used the same rates; case L
stopped at the first freeze-out.
