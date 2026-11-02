# Track wt_critic: consistency of the repair tracks, the corrected chain, and the repaired model M1r

1–2 November 2026 · Claude, consistency critic of the repair pass (BRIEF.md of this directory; issue #13 context).
Same-session AI work on one supplied model; not outside replication. Nothing here is observational evidence, a
validation, an exclusion of the real world or a new theory. **Every repair is a model choice, not a consequence of
closure** (closure supplies no law: grut_closure_forcing_20261027/SYNTHESIS.md). Git was read-only; nothing committed.
Labels: [D] derived, [C] computed with our own code, [E] estimate, [Q] quoted from a reviewed track, [A] argued,
[R] recalled.

Inputs read: the four tracks (fx_holder, fx_chain, fx_formation, fx_observables), their NOTES.md, SKEPTIC.md and
results; REVIEW_M1, REVIEW_TRANSPORT, REVIEW_KINETIC, REVIEW_REHEATING(_MODES), kin_source and hw_gates notes.
Own code: `code/chain_consistency.py` (uses `code/hartree_copy.py`, an unchanged copy of
`tracks/fx_chain/skeptic/code/hartree_indep.py`) and `code/obs_freeX.py`. Results: `results/chain_consistency.json`,
`results/obs_freeX.json`, `results/m1r_summary.json`.

## Result in one paragraph

The four tracks agree on every shared number:
- hidden couplings, the object family E(Q), the survival threshold S/Q < 3.4716, the onset density
  (151 m_X³ at T_h ≈ 8.3 m_X) and p at the spot point;
- with one exception, which is the provenance of the heat budget.

fx_chain and fx_formation, and through them fx_observables, ran the thermal history on the **old instant-release spot
point** (kin_source: λ = 1, no holder, ε = 0.146, Y_L debris R = 0.806). fx_holder's repaired release (inflaton
trilinear, λ_Y = λ_Σ = 0.5, (3, 0.5) history) is a different state: ε = 0.127, R_med = 0.40. I settled this with a
normalisation-free invariant:
- E/E_c(N) = |ε|^{-4/3}, checked to 4×10⁻⁵ across 20 releases;
- a transfer of the inflationary (GMR) vector heat, which does not depend on the release;
- our own Hartree solver.

The repaired chain gives S/Q = 4.68–5.13 at e = 10⁻⁶ against the 3.47 survival threshold [C/E]. So the corrected
physical picture of the supplied model is:
- **every object made at the φ crossover evaporates by T_h ≈ 0.8 m_X;**
- X̄ freeze out in a pure gas (p ≈ 6×10⁻⁴–1.3×10⁻³);
- **the dark matter today is a free, cold, asymmetric 92 MeV scalar gas**, not Q = 2×10⁵ Q-balls.

This holds for every release computed, Codex's own M1 (Option 1A, λ = 1) included. The repaired model M1r is
self-consistent: the relic mass per charge is computed to 0.15%, giving H_i = 1.293×10¹¹ GeV. But that number rests on
a supplied renormalisation condition: portals vanishing at M_P instead of at M_Y give 1.08×10¹¹. M1r has **no
observable signature** beyond CDM gravity. The single decisive calculation is S/Q of the actual M1r release against
3.47. It fails if S/Q < 3.47, i.e. vector heat below 0.52 of the transferred value. Then the relic is macroscopic
Q-matter and H_i moves +17%.

## 1. Are the tracks mutually consistent?

| Quantity | fx_holder | fx_chain | fx_formation | fx_observables | Settled |
|---|---|---|---|---|---|
| g, λ_X, λ_φ, v, ε (tilt) | M1 | M1 | M1 | M1 | consistent |
| λ_Y = λ_Σ | **0.5** (repaired) | heat budget from the λ = 1 spot point | same | same | **inconsistent** → checks C1–C4 |
| release state (N_X, E_X, ε) | (3, 0.5) trilinear: ε = 0.1274 | kin_source instant release: ε = 0.1457 | from fx_chain | – | **inconsistent** → C1–C4 |
| e | ledger at 10⁻⁷ and 10⁻⁶ | 10⁻⁶ | 10⁻⁶ | vapour branch uses S/Q = 1.235 and 2.07, i.e. **no vectors** (e ≫ 2×10⁻⁶) | stated per branch; S/Q = 1.235 is not reachable in M1r (floor 1.85 [C]) |
| H_i, Y_Q | 1.351×10¹¹ for the **Q = 2×10⁵ relic** (Y_Q 5.0794×10⁻⁹) | 1.343×10¹¹, 5.0794×10⁻⁹ | same; skeptic: free-X relic needs ×0.957 | same | M1r: **1.293×10¹¹, Y = 4.756×10⁻⁹** [D on C] |
| onset n_Q, T_h | – | 151.0 m_X³, 8.28 m_X | from fx_chain | – | M1r recomputed: 151.5, 8.15 [C] (2%) |
| E(Q), Q_min, binding Q | – | binding 71263, ω 0.94017 | Q_min 56709, binding 71263, E/Q(2×10⁵) 0.937400, Ω_b 0.775865 | same solver; skeptic's collocation agrees to 10⁻⁶ | consistent |
| survival threshold | – | skeptic: σ_min = 3.4716 (bulk), 2.9284 (2×10⁵) | Y_crit = 0.288054 → 3.4716 | uses 3.47 | identical |
| S/Q at e = 10⁻⁶ (spot) | – | skeptic 4.936 | 4.95 (5.4 is an artefact, skeptic) | – | 4.936 |
| p at e = 10⁻⁶ | – | 9.10×10⁻⁴ (spot) | used 2.26×10⁻³ (wrong case, skeptic) | used 2.3×10⁻³ and 5×10⁻³ | spot 9.10×10⁻⁴; **M1r 5.6×10⁻⁴ (band) / 1.3×10⁻³ (broad)** [E] |
| cascade-vector scaling | – | "1/e" | skeptic "1/e" | – | interface ratio 0.31623 = 10^{−1/2}: **e^{−1/2}** (fx_chain skeptic right) [C]; e_crit shifts < 2% |
| primordial-vector scaling | – | 1/e | 1/e | – | interface ratio 0.09994: 1/e [C] |
| object population | calibrates on Q = 2×10⁵ | η ≥ 0.85–0.92 "remaining" | η ≈ 0 at e = 10⁻⁶; Q-matter 10¹²–10²⁰ at low heat | §1, §3a, §4 assume Q = 2×10⁵ cores | **no objects in M1r** (C4); fx_observables §1, §3a, §4 and fx_holder's target do not apply |
| clock sensitivity | 1.71 | – | – | – | 3.34 at λ = 0.5 (skeptic) |
| decisive 3D run | – | – | 48³ bath at S/Q = 5 | – | infeasible (skeptic); needed only in the Q-matter branch |

**C1. Release-state invariant [C].** Before onset the X sector (massless quartic X plus φ) is conformal. Its relaxed
state at fixed comoving charge N and conformal energy E therefore depends only on E/E_c(N), where E_c is the energy of
the circular orbit of charge N. Both engines record ε = N_X/N_c(E_X) with N_c = c E^{3/4}. Values of c:
- kin_source spot point: 6.97713;
- the 19 fx_holder runs: 6.97709–6.97735.

So **E/E_c = |ε|^{−4/3}** compares releases across engines with no normalisation of a. Values:
- spot point: 13.04;
- unheld (3, 0.5): 8.60;
- Option 1A, λ = 1 (Codex's M1): 12.31;
- trilinear, λ = 1: 17.43;
- **trilinear, λ = 0.5 (M1r): 15.61**.

**C2. Heat per charge in E_c units [C/E].** W = |ε|^{−4/3}(1 + R_med) + ν W_vec,spot(e).
- R_med is fx_holder's own decay ledger for that release (Y_L + cascade + σ at the same e) [E ledger].
- W_vec,spot = 13.04 × R_vec. That is 52.8 (band) or 62.5 (broad) at e = 10⁻⁶, scaling as 1/e [C kin_source].

**C3. Vector transfer [D/E].** The hidden vectors are the inflationary Graham–Mardon–Rajendran longitudinal modes
(hw_gates gate 1), not release products. Their abundance per entropy, ρ_A/s ∝ m_A^{1/2}H_i² [R], is set at H = m_A in
the radiation era (a ≈ 65 here), long after reheating, so it is the same for any release history. The derivation:
- Decay heat per vector at fixed T goes as m_A/T_dec ∝ H_i^{1/2} (m_A = 240 e H_i, T_dec ∝ (Γ_A M_P)^{1/2}).
- The circular energy per charge scales as e_c ∝ n_X^{1/3} ∝ (Y s)^{1/3}.
- Together: ν = (H_i/H_spot)²(Y_spot/Y)^{4/3}.

For M1r (H_i = 1.293×10¹¹, Y = 4.756×10⁻⁹ against 1.3335×10¹¹ and 4.8294×10⁻⁹): **ν = 0.960**. With portals zero at M_P
(H_i = 1.076×10¹¹): ν = 0.664. The inflationary Higgs norm differs only slightly between λ = 1 and 0.5 (Σ 28.4 → 28.6
H_i). hw_gates found the light-vector abundance insensitive to the mass profile (P0 = P1). Bracket: band → broad
×1.18; our modes are 0.76 of the GMR formula.

**C4. S/Q and survival [C].** The copied Hartree solver runs at R_eff = W/13.04 − 1. Its controls reproduce the
skeptic's 1.22792 at R = 0 and 4.93631 at R = 4.9665. The thresholds in W units are:
- bulk (S/Q 3.4716): 49.06;
- Q = 2×10⁵ (2.9284): 39.31;
- binding (2.5531): 32.91.

| release (history, λ, holder) | E/E_c | R_med (10⁻⁶) | S/Q, no vectors | S/Q at e = 10⁻⁶, ν = 1, band / broad | ν needed for bulk survival at 10⁻⁶ |
|---|---:|---:|---:|---:|---:|
| spot instant release (used by fx_chain, fx_formation) | 13.04 | 0.913 | 2.06 | 4.94 / 5.39 | < 0.456 |
| (3, 0.5), λ = 1, unheld | 8.60 | 0.008 | 0.88 | 4.13 / 4.61 | < 0.76 |
| (3, 0.5), λ = 1, Option 1A (Codex's M1) | 12.31 | 2.76 | 3.32 | 5.94 / 6.37 | < 0.05 |
| (3, 0.5), λ = 1, trilinear | 17.43 | 0.214 | 1.81 | 4.75 / 5.21 | < 0.53 |
| **(3, 0.5), λ = 0.5, trilinear (M1r)** | **15.61** | **0.397** | **1.85** | **4.78 / 5.25** | **< 0.515** |
| (3, 0.5), λ = 0.5, Option 1A (fx_holder's fallback) | 11.33 | 5.48 | 4.72 | 7.13 / 7.53 | never |
| (5, 1), λ = 0.5, trilinear | 13.38 | 4.95 | 5.02 | 7.39 / 7.79 | never |
| opposite inflaton branch, λ = 1 | 9.21 | 9.67 | 5.89 | 8.16 / 8.55 | never |

**M1r with the transferred ν = 0.960** [C/E]:
- S/Q = **4.68 (band) / 5.13 (broad)**, margins ×1.35 / ×1.48 over 3.4716;
- T_h at onset 8.15 m_X, n_Q 151.5 m_X³, X̄ 0.40 per charge;
- objects survive only for e ≥ 1.86×10⁻⁶ (band) / 2.20×10⁻⁶ (broad), bulk; e ≥ 2.9×10⁻⁶ for Q = 2×10⁵;
- portals zero at M_P: S/Q ≈ 3.89 (release state of that run not recorded, so the same state is assumed [E]), margin ×1.12.

Readings:
1. The evaporation verdict at e = 10⁻⁶ does not depend on which release the tracks used. It holds for every release
   computed, including Codex's M1.
2. The object branch exists only for the trilinear or unheld (3, 0.5) histories with e ≳ 2×10⁻⁶. For Option 1A at
   λ = 0.5, and for the (5, 1) trilinear, the Y debris alone exceeds the threshold, so objects never survive. That
   removes fx_holder's fallback as a route to objects.
3. The verdict is one-sided robust. Every omitted process raises S/Q: entropy production at the crossover,
   fragmentation, the evaporation cycles, finite-k Y heating.

**C5. p(S/Q) [E].** fx_chain's gas-only scan, mapped through the same Hartree S/Q(R), gives d ln p/d ln(S/Q) =
7.5–9.9, matching the skeptic's 9.3–9.8. The mapping relies on the skeptic's finding that p depends only on S/Q.
M1r: **p = 5.6×10⁻⁴ (band), 1.3×10⁻³ (broad)**.

## 2. The corrected chain, inflation to today (M1r, e = 10⁻⁶)

1. **Inflation and preparation** (H = H_i).
   - Physics: quadratic-inflaton surrogate. Y, Σ, X sit in the conditional minimum held by
     μ_eff = μ + g_I e^{iα_I}φ_I = 10i H_i. All hidden mass² ≥ 10.19 H_i² > 9/4 H², so they stay in the vacuum.
     GMR longitudinal vectors are generated.
   - Numbers:
     - prep (90.79987, 6.88773, 28.57624) H_i, phase −π/2, unique minimum [C fx_holder; skeptic 2×10⁻¹⁶];
     - κ/κ_max 0.834, λ_Y,eff 0.2585 [C];
     - g_I = 8.00×10⁻⁷ at H_i = 1.293×10¹¹ (g_I ∝ H_i), α_I = 2.192144, φ_I(slice) = 0.8165 M_r [D];
     - vectors N_A ≈ 0.65–0.77 H_i³ (spot), ν = 0.96 [C hw_gates, E transfer].
   - Supplied: the inflaton model and the sign of φ_I; α_I; inflaton sequestering (other portals ≲ 5×10⁻¹⁵,
     Planck-suppressed ≲ 10⁻⁸, shift symmetry); an inflaton quartic ≥ 1.42×10⁻²³; the renormalisation condition of
     the induced portals (c(M_Y) = 0); ξ_H ≈ 0.1 for Higgs survival [Q].
   - Open: the earlier inflation is represented only by the static minimum (tracking smooth over s ∈ [−1.5, 4] [C]);
     the GMR abundance on this background.
2. **Reheating and charge release** (H_i t ≲ 2000, T_RH ≈ 2–5×10¹⁴ GeV).
   - Physics: the holder switches off with the inflaton amplitude (u → μ). X's first roll turns the CP phase into
     charge, and Σ exchanges charge with X slowly.
   - Numbers:
     - Y/Y_Q(H_M1) = 0.990733 (window 0.9889–0.9926) [C fx_holder, skeptic];
     - N_X = 8045.5, E_X = 1.887×10⁵, ε = 0.1274, N_YL = 0.567 [C].
   - Supplied: the (3, 0.5) history (H_i moves −16% / +21% with (5, 1)); the clock (H_i ∝ φ_I^{−2.23}).
   - Open: finite-k Y production (FX17); the long-wavelength condensate instability (REVIEW_REHEATING_MODES); X ↔ Σ
     freeze-out beyond H_i t = 2000.
3. **Debris decays.** Y_H → Y_L A at a ≈ 5; Y_L → 3X at ⟨a⟩ ≈ 4×10³ (T_vis ≈ 8×10¹⁰ GeV); vectors at a ≈ 10⁹
   (T_vis ≈ 3×10⁵ GeV).
   - Numbers: R_med = 0.397 E_X [E ledger]; vectors 52.8ν (band) to 62.5ν (broad) E_c per charge; total
     W ≈ 72.5–81.8 E_c per charge [C/E].
   - Supplied: e = 10⁻⁶; ε_BX = 0.
   - Open: ν on the actual background.
4. **Pre-onset relaxation** (a ≈ 10⁴ → 3.8×10¹²; conformal X sector).
   - Physics: 2→2 thermalisation (~1.5×10³ collisions [E]); φ channels relax X̄ with 10–10⁴ e-folds to spare [Q].
   - Numbers: circular condensate (92%) plus Hartree gas; **S/Q = 4.68–5.13** [C]; spot chain 4.936 [C].
   - Assumed: thermalisation and chemical relaxation [E/Q].
5. **φ crossover and spinodal decomposition** (T_vis ≈ 120 → 80 GeV, T_h ≈ 12 → 8 m_X).
   - Physics: the tilt makes it a crossover with a unique φ root [D]. The condensate crosses the bulk density 449 and
     the spinodal 233 m_X³ at a/a_on = 0.68 and 0.84, on a slow passage (Γ*/H ≈ 2–5×10⁷).
   - Numbers:
     - the circular condensate makes no X̄ (true-equilibrium start: n₋/q = 8.9×10⁻¹⁰) [C];
     - first lumps Q ≈ 4×10¹¹–1.2×10¹³ [E];
     - the M1 cubic raises p by 5.8% [C].
   - Open: entropy produced here, which only raises S/Q.
6. **Two-phase era and complete evaporation** (T_h ≈ 8 → 0.8 m_X, T_vis ≈ 80 → 7.5 GeV).
   - Physics: exchange runs 10⁵–10¹¹ times faster than H, so thermodynamics decides. No light hidden neutral exists, so
     the vapour must carry all the entropy, and a condensed phase can coexist only if S/Q < 3.47.
   - Numbers: η = 0 from T_h ≈ 0.80 m_X onward. Ripened lumps (≤ 7×10¹⁹) are far below what the net flux removes
     (10²⁴–10³⁴) [C threshold, E kinetics].
7. **X̄ freeze-out in a pure gas** (T_h ≈ 0.1–0.2 m_X, T_vis ≈ 1.2–1.6 GeV).
   - Physics: ss → XX̄ carries > 99.9% of the pair change.
   - Numbers: **p = 5.6×10⁻⁴–1.3×10⁻³** [E; network C] (spot 9.10×10⁻⁴).
   - Open: the network is a one-temperature moment recipe [Q].
8. **Supersaturated free gas to today** (T_h < T_s ≈ 0.10 m_X, T_vis ≲ 1 GeV).
   - Physics: objects are favoured again but cannot nucleate. The gas is non-degenerate (z = e^{5/2 − S/Q} = 0.11,
     fixed under adiabatic expansion) and the barrier is W ≈ 10⁵ T [E].
   - No channel is open at rest: XX̄ → φ needs 2.45 m_X, φφ needs 4.9 m_X, XX̄X → Xφ needs 3.45 m_X [D]. So the X̄ are
     inert.
   - Numbers:
     - λ_fs ≈ 0.3–0.5 pc comoving, M_fs ≈ 4×10⁻⁹–2×10⁻⁸ M_⊙ [E skeptic; D];
     - σ_XX/m = 6.0×10⁻⁸ cm²/g [D];
     - ΔN_eff ≈ 0 [C].
9. **Calibration (closing the loop).**
   - M = 1 + 2p = 1.0011–1.0026 m_X per charge, so Y = D/(M m_X) = 4.756×10⁻⁹.
   - **H_i = 1.293×10¹¹ GeV** (1.292 broad) [D on C]. With portals zero at M_P: 1.076×10¹¹.
   - Coordinated scaling leaves R(e) in H_i units unchanged; ν and the 2% T_vis shift are included [D/E].

Alternative model choice **M1r′** (e ≳ 1.9–2.2×10⁻⁶, S/Q 1.85–3.47):
- the relic is macroscopic Q-matter (⟨Q⟩ 10¹²–10²⁰) plus a 2.6–11% warm X vapour;
- M = 0.782–0.802 and H_i = 1.50–1.52×10¹¹ GeV;
- η alone moves the result by ±1.2%, so it cannot be calibrated to 1% until the sticking and ⟨Q⟩ are computed
  [Q fx_formation and skeptic].

## 3. The repaired model M1r

**Inputs** (all supplied):
- hidden light sector: g = 0.1, λ_X = 10⁻³, λ_φ = 0.03, v = 0.924 GeV, tilt ε = 0.002;
- heavy sector: M_Y = 30 H_i, h = 40 H_i, μ = 7.159 H_i, κ = 3.4587×10⁻³, **λ_Y = λ_Σ = 0.5**, ξ_X = −5/6;
- **holder**: g_I e^{iα_I}φ_IΣ*Y² with α_I = 2.192144 and g_I fixed by μ_eff = 10i H_i at the slice, plus:
  - the sign of φ_I during inflation;
  - a reflection-even inflaton decay sector;
  - inflaton sequestering and shift symmetry;
  - an inflaton quartic ≥ 1.42×10⁻²³;
  - **the induced-portal renormalisation condition** (zero at M_Y);
- inflaton (m_I, Γ_I) = (3, 0.5) H_i;
- **e = 10⁻⁶**, ε_BX = 0, zero SM portals;
- ξ_H ≈ 0.1, g* = 106.75, D = 4.4×10⁻¹⁰ GeV (0.6% above Planck), the 1% tolerance.

**Derived:** H_i = 1.293×10¹¹ GeV, g_I = 8.0×10⁻⁷.

**Retired:** Q = 2×10⁵ (never formed), ξ = 30.75 (Option 1A), the 10% pulse and the 20 MeV clock (no cores, so no
pulse gas), and η = 1, p = 0 (now computed).

**Fixed** relative to M1 (each a model choice):
- holder EFT (Λ_eff = M_Y/ξ ≈ H_i): renormalizable dimension-4 trilinear, technically natural, switch-off built in, 12–24×
  less light-Y debris [C fx_holder];
- Landau pole: λ = 0.5 < 0.566 [Q cf_scale];
- relaxation justification: the relaxed start is isentropic and p depends only on S/Q, so placement and the
  prethermal 1.75 plateau no longer matter. Still conditional on φ-channel relaxation [E] [C fx_chain skeptic];
- coherent-orbit X̄ inventory and ε-sensitive resonance capture: absent for the circular relaxed condensate [C];
- the supplied object size and η: replaced by a computed relic [C/E];
- the opposite-branch "trap": not a separate minimum [C skeptic].

**Not fixed:**
- the CP sign (α_I plus the sign of φ_I) is still supplied, and the reflection-even inflaton premise is fragile, because
  A φ_I|H|² makes the Higgs tachyonic on one branch [A];
- the portal renormalisation condition: ±1.13% in H_i per unit ln Q, −16.8% at M_P;
- the clock (H_i ∝ φ_I^{−2.23}) and the inflaton history (−16% / +21%);
- the vector abundance on the real background; finite-k Y; X ↔ Σ freeze-out;
- e and ε supplied; D rounded; ρ_vac underived; ξ_H supplied;
- **the evaporation verdict makes the Q-ball sector cosmologically idle.** Its only role is a transient between
  T_h ≈ 10 and 0.8 m_X.

**Frozen mass test (p, η).**
- As written, the test assumes Q = 2×10⁵ objects. M1r has η ≈ 0 [C absence, E nucleation] and p = 5.6×10⁻⁴–1.3×10⁻³
  [E]. The reading is 0.0667805 + 2.13356p = **0.068–0.070, a failure by ×6.8–7.0**. This is a calibration failure, not
  a model failure.
- In M1r the test is replaced by recalibration on the computed relic. The mass side then closes to **0.15%**
  (band against broad p).
- The 1% benchmark now depends on the yield side: windowing ±0.13%, D rounding 0.4%, the portal renormalisation
  condition (±1.1% per unit ln Q), the clock and the history. At fixed stated inputs it is about ±0.5%. Over the
  admissible choices of the renormalisation condition it is −17%.

## 4. Observable status (honest)

M1r is secluded dark matter with **no non-gravitational signature and no departure from CDM on any accessible scale**:
- **Clustering.** λ_fs ≈ 0.3–0.5 pc comoving (M_fs ~ 10⁻⁸ M_⊙), 10⁵ below the Lyman-α and satellite reach
  (k ≈ 25–50/Mpc [R]). There are no nuggets and no warm tail. Codex's parsec dip is absent because no cores exist.
- **Self-interaction.** σ/m = 6.0×10⁻⁸ cm²/g [D], 10⁶–10⁸ below SIDM interest.
- **Radiation.** There is no hidden state lighter than X, so ΔN_eff ≈ 0. The vector decays to XX̄ before BBN.
- **Annihilation.** Every channel is closed at rest. The 0.06–0.13% X̄ admixture is inert.
- **Direct detection.** Zero in M1r. With a minimal portal added, single 92 MeV X scatter with
  σ_Xn = 4.4×10⁻³⁶ θ² cm² (φ–h mixing) or 8.3×10⁻³⁹ λ_HX² cm². At the couplings fx_observables costed (A = 2.7×10⁻¹⁰ GeV,
  b = 2.2×10⁻¹⁰, λ_HX = 9×10⁻¹¹) that is 7×10⁻⁵⁹–2×10⁻⁵⁸ cm² [D]. Freeze-in at those couplings would add ΔS/Q ~ 10³,
  and p ∝ (S/Q)^{~9} would then wreck the relic [E].
- **The composite DD signature** (form-factor zeros at 27.6/34.1/41.7 keV, the light-WIMP mimic) needs Q ≲ 10⁶
  objects. Neither computed branch has them: M1r has none, M1r′ has ⟨Q⟩ ≥ 10¹² with every recoil past the first zero.
- **M1r′'s only signature** is the latent-heat vapour in P(k) [E, generic in shape]:
  - ζ ≈ 10⁻⁴: a 7–11% deficit at k = 50/Mpc;
  - ζ ≳ 5×10⁻⁴: needs the skeptic's two-fluid or kinetic treatment;
  - the size depends on the uncomputed sticking s and ⟨Q⟩.

So M1r makes no prediction that any current or planned probe can test. Its content is internal consistency.

## 5. The single decisive calculation

**The hidden entropy per charge S/Q of the actual M1r release, computed end-to-end with no input borrowed from the
instant-release spot point**, compared with 3.4716. It needs three pieces:
1. the GMR longitudinal-vector abundance (hw_gates gate-1 Proca solver) on the λ = 0.5 trilinear-held inflationary
   Higgs norm, followed through the (3, 0.5) matter era;
2. finite-k Y_L production by the trilinear-modulated light-Y mass on the same background (FX17);
3. the relaxed Hartree S/Q from (1 + R)E_X/E_c, under the adopted portal renormalisation condition (repeat for zero
   at M_P, margin ×1.12).

**Failure condition: S/Q ≤ 3.47**, equivalently a vector heat below 0.515 of the transferred estimate at e = 10⁻⁶.
- In that case the condensed phase survives and the relic is Q-matter (M ≈ 0.78–0.80 m_X per charge).
- H_i moves from 1.29 to 1.50–1.52×10¹¹ GeV (+17%).
- The calibration cannot then be stated to 1% until the sticking and ⟨Q⟩ are measured (a grand-canonical bath
  measurement, not fx_formation's infeasible 48³ box).

A computed S/Q > 3.47 settles the free-X relic, because everything this calculation still omits raises S/Q.
Second-ranked, and not a calculation: the H_i error budget is dominated by the supplied renormalisation condition of
the induced inflaton portals, which no computation fixes.

## 6. Limits and reproduction

Limits:
- The transfer factor ν is a scaling argument [E]. It is not a Proca solve on the M1r background.
- R_med is fx_holder's homogeneous decay ledger [E].
- p(S/Q) assumes the skeptic's result that p depends on S/Q alone.
- The Planck-renormalisation S/Q assumes the c = 0 release state.
- No new evolution runs were made. This track re-uses the four tracks' computed states with our own solver and
  arithmetic.

```
export OMP_NUM_THREADS=1
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/wt_critic/code
python3 -B -W ignore chain_consistency.py ../results/chain_consistency.json   # 8 s
python3 -B obs_freeX.py ../results/obs_freeX.json                             # < 1 s
```
