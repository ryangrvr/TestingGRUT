# Track fx_chain: the corrected thermal history at e = 10⁻⁶ and the antiparticle leftover p for GRUT-M1

1 November 2026 · Claude, repair pass (BRIEF.md of this directory; issue #13 context). Same-session AI work on one
supplied model. It is not outside replication. Nothing here is observational evidence, a validation, an exclusion
or a GRUT prediction. **Every repair below is a model choice, not a consequence of closure**: the relaxed-state
start, the formation rules, e = 10⁻⁶, ε = 0.002 and Q = 2×10⁵ are all supplied. Labels: [D] derived, [C] computed,
[Q] quoted, [E] estimate, [A] argued. Nothing was committed.

## Code and reuse

| File | What | Reuse |
|---|---|---|
| `code/heat_onset.py` | items 1–2: R budget and clocks, relaxed Hartree plasma, onset, M1 homogeneous condensate with thermal masses, φ_eq along the history | Hartree solver copied from `tracks/pt_relax/skeptic/code/equilibria.py`; R and clocks read from `tracks/kin_source/results/source_interface.json`; M1 conventions of `tracks/m1_vacuum` (re-implemented; bulk and spinodal reproduced) |
| `code/frag_m1.py` | item 2b: slow passage through the M1 spinodal, first-lump charge [E] | dispersion functions copied from `tracks/m1_vacuum/code/coldfluid_m1.py` |
| `code/transition/` | item 3: circular and nearly circular X + φ transition | copy of `tracks/pt_m1tilt/code` (own DOP853, mode 3 = M1). One change: `run_measure(atol=…)`. A circular orbit starts with exact zeros in x₂, p₁ and p_χ, and with the original atol = 10⁻³⁰⁰ the pure relative error control stalls the step (an infinite loop at ε ≤ 3×10⁻⁶). The runs here use atol = 10⁻¹²; controls at atol = 10⁻¹⁰ and rtol = 10⁻¹² are in group A |
| `code/kin/` | item 4: case-L moment network at M1 | `ks_common.py`, `thermo_state.py`, `stage2_moments.py`, `moment_controls.py` copied unmodified from Codex 34f23f8 (its port of our kin_solve, reviewed and reproduced in REVIEW_KINETIC). `bootstrap.py` is replaced so that it reads the review tree in place (the interface SHA-256 sums equal SOURCES.json). `m1_kin.py` patches the M1 inputs at run time. `cubic_m1.py` recomputes γ(ss → XX̄) with the M1 cubic and imports `tracks/kin_rates/code/xsec22.py` read-only. `assemble.py` |

Results: `results/heat_onset.json`, `frag_m1.json`, `circ_transition_{A,B,C}.json`, `cubic_m1.json`, `kin_m1.json`,
`fx_chain_summary.json` (the key numbers of every item).

**Reproduction checks.** The relaxed equilibria reproduce pt_relax and its skeptic: pilot T = 17.5561 with
X̄/net = 0.06914 (φ in equilibrium), and heavy release R = 0 at 13.78 / 0.099 [C]. The M1 T = 0 spinodal is
233.221 m_X³ (F = 1.8566, Ω = 0.67658) and the bulk 449.313 m_X³ (S = 0.15117, Ω = 0.775865), as in m1_vacuum [C].
The network reproduces Codex's gas-only f = 0.0045855 (R = 4.9665) and 0.0047584 (R = 5) to the quoted digits [C].
Our γ(ss → XX̄) equals the kin_rates table to 10⁻⁶ at all 45 temperatures [C].

## 1. Heat budget at e = 10⁻⁶ and the relaxed plasma at onset [C]

**Sources** (kin_source, R in units of the release energy E_X; the X-sector energy per charge before onset is (1 + R) E_X/N_X, i.e. the last column of the table below in m_X at a_on):
- Vectors (R ∝ 1/e): band 4.0538, broad 4.7905. The ratio R(10⁻⁷)/R(10⁻⁶) = 10.006.
- Mediators at e = 10⁻⁶: Y_L 0.80557 plus cascade vectors 0.10716 (∝ 1/e) plus σ 4×10⁻⁸, which is 0.91273. The
  brief's 1.144 is the e = 10⁻⁷ value.
- **Total R(10⁻⁶) = 4.9665 (band) / 5.7033 (broad).**

**Clocks.**
- 1 − 10⁻⁶ of the vector parents has decayed by a = 4.17×10⁹ (e = 10⁻⁶) and by 1.32×10¹¹ (e = 10⁻⁷). Y_L has
  decayed by 1.66×10⁴.
- The relaxed onset is at a_on = 3.81×10¹² (R = 4.97), so every source finishes 915× (vectors) to 2×10⁸× (Y_L)
  before the onset [C].
- The hard daughters (q = 1.19×10⁵ H_i) undergo about 1.5×10³ contact 2→2 collisions on the plasma before onset.
  Only about 11 halvings are needed, so they thermalise [E].
- Number changing then runs through the φ channels (10–10⁴ e-folds, REVIEW_TRANSPORT [Q]).

**Relaxed plasma.** Circular condensate plus Hartree X, X̄, φ gas, at fixed comoving N_X = 2788.18 and energy
(1 + R) E_X (heavy release). The onset is where g²(r_c² + ⟨|δψ|²⟩) + 3λ_φ⟨φ²⟩ = λ_φ v² a²/H_i².

| R | T_conf | X̄/charge before onset (φ eq / φ absent) | charge in condensate | a_on | T_h/m_X at a_on | n_Q / n_cond at a_on (m_X³) | E/Q at a_on (m_X) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 13.78 | 0.099 / 0.136 | 0.969 | 3.74×10¹² | 5.31 | 159.8 / 154.8 | 5.31 |
| 0.913 (mediators, 10⁻⁶) | 16.36 | 0.169 / 0.232 | 0.955 | 3.76×10¹² | 6.28 | 157.5 / 150.5 | 10.1 |
| 1.144 | 16.85 | 0.186 / 0.255 | 0.953 | 3.77×10¹² | 6.46 | 157.0 / 149.6 | 11.3 |
| **4.97 (10⁻⁶ band)** | **21.89** | **0.417 / 0.572** | **0.917** | **3.81×10¹²** | **8.28** | **151.0 / 138.5** | **31.1** |
| 5.70 (10⁻⁶ broad) | 22.54 | 0.457 / 0.627 | 0.912 | 3.82×10¹² | 8.51 | 150.1 / 136.9 | 34.9 |
| 41.7 (10⁻⁷ band) | 35.89 | 1.898 / 2.606 | 0.753 | 4.04×10¹² | 12.8 | 127.1 / 95.7 | 210 |

- The pilot state (n_Q = 8112.7, ρ(1 + R)) gives X̄/charge 0.069 (R = 0), 0.133 (1.14), 0.302 (4.97) and 1.385
  (41.7), with T_h/m_X = 4.75–11.8 at its onset (a = 5.3–5.7×10¹²).
- The relaxed onset of the heavy release is at 0.53–0.57 a_ref, earlier than the eccentric-orbit coordinate
  a_ref = 7.05×10¹². V(0)/n_Q there is 0.63 m_X (M1, 95 m_X⁴), not 3.1–4.0 m_X.
- Just after onset (massive phase, fast s, network initial state), R = 4.97 gives T/m_X = 8.28, X̄ 0.390 and s 0.438
  per charge (P/n_Q = 0.83).

## 2. Onset condensate density and the M1 crossover [D/C/E]

- **[D]** For a pure circular massless condensate, n = √λ_X r³ with g²r² = λ_φv², which gives
  **n_onset = 164.3 m_X³**. In the untilted model this is exactly the homogeneous branch's spinodal: φ = 0 is exact
  for F² ≥ 3, and Ω(n) turns down where φ turns on.
- **[D]** With the M1 tilt, φ_eq has exactly one root for every F and every thermal mass: F² + t_φ − t_lin/S is
  strictly decreasing on (0, 1). So the transition is a **crossover**, with no metastable φ ≈ 0 state. The curvature
  criterion (λ_φ − 6ε) would put it at 76.4 m_X³, but φ is already 0.5–0.8 v at the untilted onset.
- **[C]** Along the relaxed history (R = 4.97; thermal masses from the Hartree state):

  | a/a_on | 0.50 | 0.68 (n_c = 449, bulk) | 0.84 (n_c = 233, spinodal) | 1.0 | 1.5 | 2.0 |
  |---|---:|---:|---:|---:|---:|---:|
  | φ_eq/v | 0.07 | 0.15 | 0.34 | 0.54 (thermal) / 0.67 (T = 0) | 0.88 | 0.95 |
  | T_h/m_X | 16.6 | 12.3 | 9.85 | 8.28 | 5.5 | 4.1 |

- The condensate crosses the M1 bulk density at a/a_on = 0.70 / 0.68 / 0.60 (R = 0 / 4.97 / 41.7) and the spinodal
  at 0.87 / 0.84 / 0.74, with T_h/m_X = 6.1 / 9.9 / 17.3 there.
- Thermal masses frozen at onset move the spinodal to 223 / 206 / 159 m_X³ and the bulk to 465 / 490 / 549 [E].
- **So the relaxed condensate enters the fragmentation regime before the untilted onset, slowly and from the stable
  side: objects form by spinodal decomposition at n ≈ 230 m_X³, T_h ≈ 6–10 m_X (R ≤ 5).**
- **[E] First lumps** (T = 0 dispersion; no thermal damping, noise model or nucleation in the 233–449 band). Near the
  spinodal, Γ_peak = 0.774 δ m_X and k_peak = 0.958 √δ m_X. With H/m_X = 1.3×10⁻¹³, an amplification of 20–60
  e-folds is reached at δ* = 4.5–8×10⁻⁶, where k = 2.0–2.8×10⁻³ m_X. This gives Q ~ n(2π/k)³ = 3×10¹²–7×10¹²
  ((π/k)³: 4–9×10¹¹).
- A sudden quench deep in the band gives 3–5×10⁵ instead.
- Either way the lumps are above Q_bind = 7.13×10⁴. The old symmetric model's fastest lump (2.5×10⁴) was below
  Q_min. The supplied Q = 2×10⁵ is not what the slow passage selects.

## 3. Does the M1 transition preserve the circular condensate? [C]

Copied pt_m1tilt engine, M1 potential, τ-20 charge q = 8117.88 on the circular orbit r₀ = (q/√λ_X)^{1/3} = 63.555,
drift κ = ε_rate √(6λ_X) r₀, same start, window and measures as pt_m1tilt.

| Orbit | ε_rate | n₋/q | N_φ/q |
|---|---|---:|---:|
| circular (J_r = 0) | 10⁻³, 10⁻⁴, 3×10⁻⁵, 10⁻⁵, 3×10⁻⁶, 10⁻⁶ | 2.50–2.84×10⁻⁷ (flat: the measure's floor) | 8.6–8.8×10⁻⁷ |
| circular, ε_t = 0.0018 / 0.0022 | 10⁻³, 10⁻⁵ | 1.9×10⁻⁷ / 4.0×10⁻⁷ | 8.6–8.8×10⁻⁷ |
| δ = r(0)/r₀ − 1 = 10⁻³ (J_r/q = 1.2×10⁻⁶) | 10⁻³–10⁻⁵ | 0.3–2.6×10⁻⁶ | 9×10⁻⁷ |
| δ = 0.01 (J_r/q = 1.22×10⁻⁴) | 10⁻³ … 10⁻⁶ | 0.91–1.07 J_r … 0.44–0.88 J_r | ≤ 3×10⁻⁵ |
| δ = 0.03 (1.09×10⁻³) | 10⁻³, 10⁻⁴ / 10⁻⁵ | 0.94–1.02 J_r / 0.33–0.99 J_r (2 of 4 phases captured) | ≤ 3.9×10⁻⁴ |
| δ = 0.1 (1.19×10⁻²) | 10⁻³ … 10⁻⁶ | 0.35–0.92 J_r | 0.5–3.9×10⁻³ |

- Work identity ≤ 7×10⁻⁸ and charge drift ≤ 7×10⁻⁸ on every run. The rtol 10⁻¹² and atol 10⁻¹⁰ controls change n₋
  by ≤ 2×10⁻⁷.
- **The transition generates no antiparticles and no φ quanta from a circular condensate, at every rate down to
  10⁻⁶ and for ε_t ± 10%.** Resonance capture needs radial action.
- For nearly circular orbits, n₋ + (3/2)N_φ ≲ J_r ≈ 1.2 δ² per charge.
- At thermal-fluctuation-level eccentricity, the zero mode of a fragmentation-scale patch has δ ~ 10⁻³, so J_r/q ~ 10⁻⁶
  [E]. The coherent orbit's post-transition X̄ (1.5 per charge in pt_m1tilt) is absent in the relaxed history. The
  X̄ handed to case L is the thermal 0.4 per charge, not a coherent inventory.

## 4. p for M1 from the relaxed state (case-L moment network) [C]

- M1 inputs: object ω = 0.88172 m_X with E/Q = 0.9374; binding object Q = 71263 with E/Q = 1, ω = 0.94017;
  Y_Q = 5.0794×10⁻⁹; V(0) = 95 m_X⁴; capture ×1.041 (r₉₀²).
- Start: forced chemical equilibrium at the relaxed onset a_on, with ρ = ρ_c(1 + R)(a_ref/a_on)⁴ + V(0).
- **M1 cubic [C].** Our channel budget shows that ss → XX̄ carries more than 99.9% of the pair change at every
  temperature. With V‴(v) = 6(λ_φ − ε)v its rate falls by 9.6% at freeze-out (ratio 0.904–0.93 for T ≤ m_X; 0.9025
  at threshold, where the amplitude is −0.0267 → −0.0253). This is applied as γ_ss × r(T). The s self-reactions
  (≤ 10⁻³ of the pair change) are left unscaled.

| Case (M1, relaxed start unless noted) | R | p = P/n_Q | f = 2.02855 p | 2.13356 p | η ≥ |
|---|---:|---:|---:|---:|---:|
| **gas only (no capture), e = 10⁻⁶ band** | 4.9665 | **9.10×10⁻⁴** | 1.85×10⁻³ | **1.94×10⁻³** | **0.879** |
| gas only, e = 10⁻⁶ broad | 5.7033 | 2.16×10⁻³ | 4.39×10⁻³ | 4.61×10⁻³ | 0.919 |
| gas only, mediators only (10⁻⁶ / 10⁻⁷ value), no vectors | 0.913 / 1.144 | 8.1 / 9.7×10⁻⁶ | | 1.7 / 2.1×10⁻⁵ | 0.8505 |
| gas only, R = 0 | 0 | 5.3×10⁻⁶ | | 1.1×10⁻⁵ | 0.8504 |
| gas only, e = 10⁻⁷ band | 41.71 | 0.723 | 1.47 | 1.54 | fails (by ×154) |
| instant objects Q = 2×10⁵, or binding objects Q = 7.13×10⁴ (R = 1.14, 4.97, 5.70) | | < 10⁻¹⁸ (floor) | | 0 | 0.8503 |
| instant objects, capture × 10⁻³ / 10⁻⁵ / 10⁻⁷ (R = 4.97) | 4.9665 | floor / floor / 5.5×10⁻¹³ | | 0 | 0.8503 |
| instant objects, e = 10⁻⁷ (stop at the kinetic-closure diagnostic) | 41.71 | 0.49 | 0.997 | 1.05 | fails |

Variants at R = 4.97, gas-only p:
- no cubic: 8.61×10⁻⁴ (the cubic raises p by 5.8%);
- start at the M1 spinodal crossing: 8.81×10⁻⁴;
- old Y_Q: 9.28×10⁻⁴;
- old V(0) = 75: 8.82×10⁻⁴;
- rtol 2×10⁻⁹: 9.1023×10⁻⁴ (identical);
- **old placement** (equilibrium start at a_ref with V(0)/n_Q = 4 m_X, M1 inputs): 3.22×10⁻³, which gives
  η ≥ 0.953. The old recipe gives 2.26×10⁻³.

Charge drift ≤ 6×10⁻⁷ and energy-work drift ≤ 1×10⁻⁷ throughout.

**R scan** (gas only, M1, relaxed): p = 7.5×10⁻⁵ (R = 3), 2.7×10⁻⁴ (4), 9.1×10⁻⁴ (4.97), 1.72×10⁻³ (5.5), 2.99×10⁻³
(6), 4.97×10⁻³ (6.5), 7.8×10⁻³ (7), 1.66×10⁻² (8). The test term alone reaches 0.01 at **R = 6.44**, which is
e ≥ 7.4×10⁻⁷ (band) / 8.7×10⁻⁷ (broad) for η = 1.

## 5. The condition on η that remains [C, conditional]

With η the charge fraction in reference objects, the test is 0.0667805(1 − η) + 2.13356 p ≤ 0.01:

- **η ≥ 0.879** at e = 10⁻⁶ (band), if X̄ annihilate only in the gas (no capture by objects; the upper value of p).
- **η ≥ 0.919** with the broad vector source.
- **η ≥ 0.8503** if the objects, which in the relaxed M1 history form at the spinodal (T_h ≈ 6–10 m_X) long before
  X̄ freeze-out (T ≈ 0.1–0.2 m_X), capture X̄. Then p < 10⁻¹⁸, robust to a 10⁷ reduction of the capture kernel, which
  covers the Q^{−1/3} dilution for Q ~ 10¹¹–10¹² lumps.
- e = 10⁻⁷ fails under every rule (p = 0.49–0.72).

**Assumptions.**
- (i) The relaxed state is chemically relaxed through the φ channels before onset [Q REVIEW_TRANSPORT], and the
  injected heat thermalises [E].
- (ii) Moment network: one temperature, fast s, MB rates with μ = 0 Bose factors, vacuum masses up to T/m_X ≈ 8–13
  (one-loop T_c = 18 m_X).
- (iii) Formation is supplied as a rule (gas-only = no capture; instant objects at fixed Q). Neither is a formation
  calculation.
- (iv) The capture kernel is the old object's, scaled geometrically [E].
- (v) The crossover between the bulk crossing and a_on is not evolved kinetically. The two starts that bracket it
  give 8.8–9.1×10⁻⁴.

**What remains is η itself.** The same slow spinodal passage suggests first lumps far larger than the supplied
Q = 2×10⁵ [E]. Objects of different Q have a different E/Q (0.94 → 0.78), and that changes the calibration
(Y_Q ε_Q = D). It is not an η.

## Corrections to the physical picture (this track)

1. **Onset and V(0) per charge.** In the relaxed history, onset is set by the circular condensate plus the thermal
   ⟨r²⟩, at a_on = 0.53–0.57 a_ref and n_Q ≈ 150–160 m_X³. V(0) per charge is therefore 0.6 m_X, not 3.1–4.0. The
   corrected start lowers p at R = 4.97 by 3.5× relative to the a_ref equilibrium start [C].
2. **No sharp transition.** Under M1 the φ transition is a crossover; the instability that matters is the condensate
   spinodal at 233 m_X³, reached before the untilted onset [D/C].
3. **No coherent X̄.** The transition does not create antiparticles from the circular condensate [C]; the X̄ budget at
   onset is thermal (0.42 per charge at R = 4.97).
4. **The M1 cubic is a small upward correction.** It raises p by about 6%, through ss → XX̄ (−9.6%) [C].

## Reproduce

```bash
T=research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_chain/code
python3 -B $T/heat_onset.py                       # items 1-2, 25 s  -> results/heat_onset.json
python3 -B $T/frag_m1.py                          # item 2b, 1 min   -> results/frag_m1.json
cd $T/transition && python3 -B circ_scan.py A && python3 -B circ_scan.py B && python3 -B circ_scan.py C   # 4 + 3 + 9 min (2 workers)
cd $T/kin && python3 -B cubic_m1.py               # 2 min -> results/cubic_m1.json (needed by m1_kin.py)
python3 -B m1_kin.py m1_gas_relaxed_R497 --budget # ~35 s per gas case; object rules: FX_RTOL_OBJ=2e-6 (default), 5-20 s
python3 -B m1_kin.py repro_gas_R497 repro_gas_R5 m1_gas_relaxed_R114 ... m1_gas_relaxed_Rx<R> m1_inst_relaxed_R497_capx<k>
python3 -B assemble.py                            # -> results/kin_m1.json
```
Case outputs go to the scratchpad (`wt2/fx_chain/kin_<case>.json`; caches in `wt2/fx_chain/kwork`); all inputs are read
in place from `grut_formation_review_20261025/tracks` (kin_source, kin_rates, kin_thermo, m1_vacuum, pt_m1tilt MATCH file
via the scratchpad copy of 5f5513b used by pt_m1tilt).

Credits: Q-ball formation by fragmentation of a condensate (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki);
Bose thermodynamics by Bessel series (Kolb–Turner); Vakhitov–Kolokolov stability as used in REVIEW_M1.
