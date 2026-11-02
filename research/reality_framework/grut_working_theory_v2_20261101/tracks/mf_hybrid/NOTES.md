# Track mf_hybrid: a repaired multi-field inflaton sector for M1r, with exported benchmark histories

2 November 2026 · Claude, multi-field inflaton pass (BRIEF_MULTIFIELD.md; GitHub issue #13 context). This is same-session
AI work on one supplied model, not outside replication. Nothing here is evidence for an inflaton, a waterfall, the hidden
sector or any dark-matter branch. Nothing is a validation, an exclusion of the real world, a new theory or a GRUT
prediction. **Every potential, coupling, renormalisation condition, slice convention, decay channel and criterion below
is a model choice; none follows from closure.** Data comparisons (Planck 2018, BICEP/Keck 2021, ACT DR6) are
retrospective.

Labels: [D] derived, [C] computed with this track's code (`code/`, output in `results/`), [E] estimate, [Q] quoted from a
reviewed track, [A] argued, [R] recalled or literature.

Git was read-only and nothing was committed. All runs were in the foreground with one process. No compiled files were
left in the tree (Python run with `-B`). Scratch: `$SP/wt2/mf_hybrid/`. No lattice was run (allowed, but not needed for
the homogeneous stage, see §6).

## 0. Result in brief

1. **The brief's architecture fails as stated: a heavy waterfall cannot hold the energy [D/C].** With the full χ one-loop
   potential along the inflationary valley, the clock's post-waterfall mass obeys m_σ m_χ² = 4√3 X M_P H_i², where
   X = g²σ_c must stay O(1) H_i (§2). For a waterfall heavy enough to reheat perturbatively (m_χ ≳ 3×10³ H_i, λ ≳ 6×10⁻³),
   this forces m_σ = X/√λ ≲ 13 H_i. The lighter field holds V0 on the adiabatic valley. So the clock σ ends up quiet,
   with amplitude ~0.4–1.2 M_P: the single-field obstruction (G1) returns on σ. **G3 is refuted as an architecture; G4 is
   confirmed and sharpened.**
2. **No field holding V0 can reheat perturbatively within a Hubble time unless m ≳ 3.4×10³ H_i [D]**
   (amplitude √6 H M_P/m and y A < m/2). No CW-flat hybrid of this class reaches that.
3. **What works instead [C/E]: a sub-Planckian hybrid with comparable masses (m_χ = 250, m_σ = 600 H_i) that reheats
   non-perturbatively through axial (mass-phase) couplings to singlet neutrinos.** The coupling is
   −(M/2) N̄ e^{iγ5φ/f} N. It has no field-dependent one-loop potential and no kinematic blocking. Its half-period production
   is computed exactly from the fermion mode equations, and agrees with Landau–Zener within ×1.15–2.3 (§4.4). With four
   singlets (f = 1.2–2×10³ H_i, μ5 ≤ 0.5 × 4πf) the radiation fraction is 54% at t = 0.3 H_i⁻¹ and 77% at 0.6. A 13–19% χ
   remnant then decays perturbatively by t ≈ 900 H_i⁻¹, which gives **K = 2.05** (M1r: 1.99). With no
   non-perturbative transfer after t = 0.6 H_i⁻¹, K would be 3.26.
4. **Two exported benchmarks, final for their model choices** (`results/bench_HFa.json`, `results/bench_HFb.json`):

| | HFa | HFb |
|---|---|---|
| clock / waterfall | tilted-hilltop clock + full χ loop; biased real waterfall | same |
| spectator mass | waterfall, m_s = g_s χ, g_s = 8.17×10⁻⁶ | curvature, m_s² = 9H_i² − (3/4)R |
| H_i (slice σ = σ_c); V0^{1/4} | 1.30×10¹¹ GeV; 7.40×10¹⁴ GeV | same |
| n_s; A_s; r; α_s (exact MS) | 0.9649; 2.10×10⁻⁹; 2.75×10⁻⁷; +1.1×10⁻⁴ | 0.9743; 2.10×10⁻⁹; 2.75×10⁻⁷; +1.4×10⁻⁴ |
| pull vs Planck 2018; vs ACT DR6 (retrospective) | 0.0σ; −2.8σ | +2.2σ; 0.0σ |
| N_* (from the computed K) | 52.64 | 52.64 |
| waterfall | fast: ε_H = 1 at t = 0.19 H_i⁻¹ (H = 0.997 H_i), χ = v/2 at 0.19 | same |
| defects | biased: classical/quantum ≥ 33; P(wrong vacuum) ~ e⁻⁵⁵⁰ per volume | e⁻⁶¹⁸ |
| K; T at radiation–matter equality; T_K | 2.048; 2.2×10¹⁴ GeV; 4.35×10¹⁴ GeV | 2.047; same |
| SM (Higgs-coupled) plasma | ρ_SM/ρ ≥ 10⁻³ at τ = 0.64; ≥ 50% at τ = 7.3 [E] | same |
| s_*; isocurvature; energy fraction | 7.3×10¹⁵ GeV (3×10⁻³ M_P); β_iso = 3.8×10⁻³ p²; 1.35×10⁻⁵ | same |

   - p = d ln Ω_DM/d ln s_* comes from mf_ideal/mf_chain.
   - Every field value is ≤ 0.028 M_P.
   - The holder coupling is |g_h| = 2.19×10⁻⁴, so that μ_eff = 10i H_i at the slice.
5. **What this does not settle.**
   - The relic branch: S/Q needs the chain (mf_chain), on these R(t), K and m_s(t).
   - The reheating rates are [E]. They are built on [C] half-period production, without back-reaction, Pauli blocking or
     fragmentation.
   - The axial couplings sit inside, but within ×2 of, their EFT cutoff.
   - The clock and the waterfall carry periodic couplings that a UV completion must explain (§6).

## 1. The estimates G1–G7: verified, refuted or sharpened

| | Claude's estimate | Verdict | Label |
|---|---|---|---|
| G1 | A field with V0 has amplitude √6 H M_P/m (0.82 M_P at m = 3H); fast perturbative decay blocked; the holder must sit on a subdominant field | **stands, generalised**: perturbative Γ ≥ H needs m ≥ (384π)^{1/5}(M_P/H)^{2/5} H = 3.4×10³ H_i at H_i = 1.3×10¹¹ GeV (`structure.json`). The holder is put on the spectator s | D |
| G2 | proposal as written: λ_χ = 4.8×10⁻¹⁵, v = 1.63 M_P, CW 2×10¹³ × λ, y v/(m/2) = 4×10⁷, m² = −4.5H² at χ = 0, walls, holder off | **algebra stands** (v = √24 M_P/m_χ, λ = m_χ⁴/48M_P²H² [D]). The mild-waterfall and δN audit is mf_gemini's | D |
| G3 | a heavy waterfall reheats perturbatively: Γ_χ ≤ λm_χ/32π, m_χ ≳ 3–5×10³ H_i, loop safe | **algebra stands, architecture refuted.** The χ loop on the clock forces m_σ ≲ 13 H_i at such λ. σ, not χ, holds V0, and σ is quiet (§2) | D/C |
| G4 | σ holds ≈ V0 on the valley; σ → −σ makes σ stable; Yukawa on σ capped ~2×10⁻³ (Γ_σ ~ 10⁻⁵ H) | **stands and sharpened.** It holds whenever m_χ ≫ m_σ, which one-loop flatness enforces for heavy χ. For m_σ > m_χ (HF) σ carries 6% of the energy at ε_H = 1 and ≈ 15% after it. The tilt and bias terms break σ → −σ, but the vacuum mixing is negligible (θ ~ 10⁻⁸). The flatness cap y_σ ≲ 10⁻³ (Γ_σ ~ 10⁻⁵ H) is confirmed. **New route:** axial σ → N N, with Γ_pert = 0.29 H_i and a non-perturbative 9.6 H_i at the slice amplitude | C/E |
| G5 | spectator holder; (a) ½g_s²χ²s², (b) −½ξ_sRs² with ξ_s ≈ −3/4; CP sign α + sign(s_*); isocurvature s_* ≳ 2×10¹⁵ GeV × p; energy 1.5(s_*/M_P)² | **both mechanisms work** on the benchmarks. m_s rises from ~0 to 3 H_i by t ≈ 0.2 (t_on ≈ 0.2 H_i⁻¹). β_iso = (H/2πs_*)²p²/A_s = 3.8×10⁻³ p² at s_* = 7.3×10¹⁵ GeV. (b) needs \|ξ_s + 3/4\| ≲ 10⁻³ (a tuning), and works because ε_H ~ 10⁻⁸ keeps R = 12H² | C/D |
| G6 | portals linear in a vacuum fluctuation oscillate in sign; positive χ²\|H\|² or (v² − χ²)\|H\|² keep the Higgs heavy | **sharpened.** On the adiabatic valley, λ_χH(χ² − v²)\|H\|² = −(λ_χH g²/λ)σ²\|H\|². That is positive-definite for λ_χH < 0, but sign-oscillating once χ oscillates about v off the valley (always in HF). The benchmarks set every Higgs portal to zero. The only Higgs tachyonic window is ξ_H R ≥ −0.45 H_i² for ≲ 0.05 H_i⁻¹ in the kinetic-dominated waterfall (exponent ≈ 0.03) | D/C |
| G7 | R_vec ∝ K², K ≥ 1; fast reheating K → 1–1.3 could cut vector heat ×4 | K = 2.05, not 1–1.3. The radiation fraction is 85% by t = 3, where K would read 1.64, but the 13–19% χ remnant restores M1r-like K. The branch is not prejudged | C |

## 2. Structural result: one-loop flatness makes the clock quiet [D + C]

**Relations [D].** Take V0 = λv⁴/4 = 3M_P²H², m_χ² = 2λv², m_σ = gv and σ_c² = λv²/g². Then:
- v = √24 M_P/m_χ;
- σ_c = 2√3 M_P/m_σ;
- X ≡ g²σ_c = m_σ m_χ²/(4√3 M_P). X is half the slope of M_χ²(σ) = g²σ² − λv² at σ_c.

**One-loop flatness [C, `structure.json`].** The χ loop, (M_χ⁴/64π²)(ln M_χ²/μ² − 3/2), enters the clock potential. At
each trial pivot, the tree clock (c1, m2) is retuned so that A_s and n_s are met. N_max is the largest number of e-folds
before σ_c for which this gives a monotonic potential.

| (m_χ, m_σ)/H_i | X | σ_c/M_P | N_max at μ² = (1, 10, 100) × M_χ²(2σ_c) |
|---|---|---|---|
| (250, 600) | 0.29 | 0.0058 | 56, 73, 73 |
| (300, 600) | 0.42 | 0.0058 | 47, 71, 54 |
| (300, 750) | 0.52 | 0.0046 | 38, 62, 35 |
| (1000, 100) | 0.77 | 0.035 | 89, 149, 149 |
| (4000, 3) | 0.37 | 1.15 | 170, 171, 172 |
| (4000, 10) | 1.23 | 0.35 | 170, 327, 414 |
| (4000, 30) | 3.7 | 0.12 | 108, 85, 54 |

- **Reading the table.** The benchmarks need N_* = 52.6. The failure sets in once g²σ_* ≳ 1–2 H_i or X ≳ 0.5–4,
  depending on the renormalisation scale (the tree quartic).
- **Consequence for heavy χ.** At m_χ = 4000 H_i (λ = 0.015, v = 1.2×10⁻³ M_P), a flat clock has m_σ ≲ 10–30 H_i and
  σ_c ≥ 0.1 M_P. χ follows the valley adiabatically, and the energy V0 ends up as σ's oscillation [A on the valley
  argument; the depth released at the onset is ≈ 4(δσ/σ_c)² V0 ≪ V0].
- **What that costs.** The realisable σ channels are then:
  - Yukawa, y ≲ 4×10⁻⁴ at σ ~ M_P: Γ ~ 10⁻⁸ H_i;
  - axial within its EFT: Γ ≲ 10⁻⁶ H_i perturbatively.
  So K ~ 10²–10³, as in TB1/TB2.
- **The "HQ" quiet-clock hybrid.** It does realise "quiet + low scale" (m_σ ≈ 3 H_i, r ≈ 3×10⁻⁷), but not "fast". It is
  not exported. Whether an axial coupling at the EFT edge could empty it (estimated O(1) transfer per half-period, never
  computed) is left open [E].

## 3. The repaired architecture (benchmarks HFa, HFb)

V = (λ/4)(χ² − v²)² + ½g²σ²χ² − cχ + c₁(σ − σ_c) − ½m₂(σ − σ_c)² + V_1-loop + V_s + V_N + const

**The inflaton sector**
- **Clock σ.** A tilted hilltop: the linear term sets the slope, and the negative curvature together with the χ loop gives
  the red tilt.
  - HFa: c₁ = 2.26×10⁴ H_i³, m₂ = 0.0220 H_i².
  - HFb: c₁ = 1.75×10⁴ H_i³, m₂ = 0.0098 H_i².
  - The pivot sits at σ_* = 0.021 (0.019) M_P, with Δσ = 0.015 (0.013) M_P to σ_c.
- **One-loop renormalisation condition.** μ² = 10 M_χ²(2σ_c), equivalently a tree quartic of zero at that scale.
  - At the pivot the loop supplies 58% (HFa) and 74% (HFb) of η, and −0.35 / −0.30 of the slope.
  - So the red tilt is partly radiative, and n_s depends on this renormalisation choice.
- **Waterfall χ (real).**
  - λ = 2.32×10⁻⁷, v = 0.0196 M_P = 4.8×10¹⁶ GeV, g = 1.63×10⁻³, σ_c = 0.0058 M_P, X = 0.289.
  - Vacuum masses m_χ = 250 H_i and m_σ = 600 H_i. Since m_σ > m_χ, χ carries most of the energy.
- **Bias.** c = 10⁵ H_i³, a soft Z2 breaking.
  - During inflation χ_min = c/M_χ², positive everywhere, so the waterfall falls to +v in every patch.
  - Vacuum shift 1.6 H_i; vacuum splitting 2cv ≪ V0; inflationary effect ≤ 10⁻⁶ of V′ at the pivot.

**Spectator and holder**
- **Spectator s.** It is light during inflation. Its holder term, −(g_h e^{iα} s Σ*Y² + h.c.), replaces M1r's
  inflaton holder.
  - Mechanism (a), HFa: ½g_s²χ²s² with g_s = 3H_i/v.
  - Mechanism (b), HFb: ½(m_s0² + ξ_sR)s² with m_s0 = 3H_i and ξ_s = −3/4.
  - s_* = 3×10⁻³ M_P and |g_h| s_* = |10i − μ| = 12.30 H_i.
- **Spectator radiative stability [C/E].**
  - Holder-loop mass² ≈ 2.5×10⁻⁷ H_i².
  - Tree tadpole force from the prepared state 0.6 H_i³, i.e. a drift of 10 H_i over 50 e-folds, against s_* = 5.6×10⁴ H_i.
  - Waterfall-induced inflationary mass g_s²χ_min² ≤ 3×10⁻⁵ H_i².

**Reheating channels (V_N)**
- Four Majorana singlets: −(M_j/2) N̄_j e^{iγ5φ_j/f_j} N_j.
  - σ → N₁: M₁ = 250, f₁ = 1200 H_i.
  - χ → N₂,₃,₄: M = 400, 200, 120 H_i and f = 2000 H_i.
- Each N decays to LH by the type-I seesaw (m_ν = 0.05 eV): M_N = 1.6–5.2×10¹³ GeV, y_ν² = 0.026–0.086, rest width
  0.12–1.37 H_i.
- |m_N| does not depend on φ, so V_N has no field-dependent one-loop term [D]. It does not touch the CMB.
- The two-loop potentials from the lepton-phase mismatch are ~ y_ν⁴M⁴/(16π²)² cos(…) ≈ 400 H_i⁴. Their slope is
  ≤ 10⁻⁵ of V′ and their curvature ≤ 10⁻² of η H² [E].

**Higgs and hidden sector**
- **Higgs.** It has no tree-level portal to σ, χ or s (a model choice), and ξ_H = 0.1 as in M1r.
- **Hidden sector.** M1r at λ_Y = λ_Σ = 0.5. There are no σ- or χ-hidden portals (sequestered by choice). The
  hidden-sector loops do not feed back on σ or χ, because no portal exists [A].

## 4. Benchmark physics, item by item

### 4.1 CMB [C] (`inflation.py`: slow-roll calibration, then exact Mukhanov–Sasaki)

- **Calibration.** The slow-roll fit is retargeted once so that the MS spectrum gives A_s = 2.100×10⁻⁹ and the target n_s
  exactly.
- **N_* bookkeeping.** N_* = 52.64 comes from e-fold bookkeeping with each benchmark's computed K: ρ_r a⁴ = 3K²M_P²H_i²,
  g = 106.75, g_s0 = 3.931, k_* = 0.05 Mpc⁻¹. It converged in one step.
- **Results.** r = 2.75×10⁻⁷ is far below BICEP/Keck's 0.036 (arXiv:2110.00483). α_s = +1.1 / +1.4×10⁻⁴: the χ loop's
  running is tiny.
- **Two n_s targets.** HFa is centred on Planck 2018 (0.9649 ± 0.0042, arXiv:1807.06211) and HFb on ACT DR6 P-ACT-LB
  (0.9743 ± 0.0034, arXiv:2503.14452). The hybrid family can reach either; the tilt is a model choice, not a prediction.
  The pulls are in the table of §0.
- **The other fields during inflation.** χ is heavy during inflation (M_χ ≈ 800 H_i at the pivot), so it has no
  isocurvature or ζ contribution on CMB scales.

### 4.2 Waterfall, small-scale curvature, defects [C] (`waterfall.py`, `diagnostics.py`)

- **The run.** The homogeneous two-field run (with s and linear δχ_k, k ≤ 1 m_χ) starts 6.1 e-folds before σ_c on the
  slow-roll attractor.
- **Speed.** After σ_c, χ reaches v/2 at t = 0.19 H_i⁻¹. ε_H = 1 at t = 0.19 with H = 0.997 H_i: an abrupt end. R swings
  from 12 to −4.5 and then oscillates at ω ≈ 470–490 H_i (≈ 2m_χ).
- **Defects.** The classical (bias-driven) χ exceeds the amplified vacuum spread δχ_rms by ≥ 33 throughout the linear
  stage.
  - A wrong-vacuum patch has ln P ≈ −550 (HFa) / −618 (HFb) per correlation volume, against ln N_vol ≈ 174 for today's
    horizon. **No walls form** [C on a Gaussian tail].
  - Inflating defects away is not needed. Lazarides–Panagiotakopoulos smooth hybrid (hep-ph/9506325) is the alternative
    route [R].
- **Small-scale curvature.**
  - The tachyonic phase lasts ≈ 0.2 H_i⁻¹ (fast, not mild). The time-delay curvature is ζ_wf ≲ (δχ/χ_cl) × H Δt ≈ 6×10⁻³,
    on sub-Hubble scales k ~ m_χ a [E].
  - That gives no PBH. PBH formation needs a mild waterfall lasting several e-folds (Clesse–García-Bellido 1501.07565;
    the waterfall ζ in Lyth 1012.4617) [R].
- **Fragmentation.** By t = 0.6 the linear modes grow ×10³ (k = 0) and ×10² (k = 0.3 m_χ), with δχ_rms/v ≈ 0.005.
  - The field is still coherent when the history switches to fluids. This is because the decays damp the oscillation, so
    χ no longer swings back through 0.
  - Without decay friction, growth reaches ×10⁴–10⁵ by t = 0.6: tachyonic self-resonance (Felder et al. hep-ph/0012142,
    hep-th/0106179 [R]).

### 4.3 Radiative stability along the trajectory [C/D/E]

The full χ loop is inside the CMB calibration (§2, §4.1). The other one-loop terms are §3's list:
- N loops: zero for the axial couplings [D];
- holder, s and tadpole terms on s: negligible [C/E];
- Higgs: σ- and χ-independent (no portal) [D];
- σ loop on χ's quartic: g⁴/16π² = 4.5×10⁻¹⁴ against λ = 2.3×10⁻⁷ [C];
- the χ self-coupling's running: negligible.

The tuning is the clock's tree curvature, which cancels or adds to the loop at the pivot (ratio 1.4–2.9). One more tuning
is the cosmological constant.

### 4.4 Decay channels and the CP-sign bookkeeping [C/E] (`axial.py`, `reheat.py`, `diagnostics.json`)

- **Axial production [C].** Rotating the mass phase gives μ5 = φ̇/2f. Per helicity the Hamiltonian is [[p, M], [M, −p]],
  with p = hk − μ5(t).
  - One half-period of φ = Σ cos mt is integrated exactly (RK4; norm conserved to 10⁻¹⁰).
  - The produced energy agrees with the two-crossing Landau–Zener form n_k = 2P(1 − P), P = exp(−πM²/|μ̇5|), within
    ×1.15–2.3 over A = μ5,max/m = 5–40 and M/m = 0.3–2.
  - At small amplitude it tends to the perturbative pseudoscalar width (M/f)²mβ/16π, within ×2.
  - Cf. Adshead–Sfakianakis 1508.00891 for derivatively coupled fermions after axion inflation [R]. The non-perturbative
    rate per condensate is Γ = (m/π) ε(A)/(S²/2), with ε from the tables.
- **Rates at the slice amplitude [E].**

  | channel | Γ non-perturbative | Γ perturbative | μ5,max/4πf | f / f_edge |
  |---|---|---|---|---|
  | σ → N₁N₁ | 9.6 H_i | 0.29 H_i | 0.49 | 1.4 |
  | χ → N₂N₂ (M = 400) | 7.2 H_i | 0 (2M > m_χ) | 0.42 | 1.5 |
  | χ → N₃N₃ (M = 200) | 3.0 H_i | 0 | 0.42 | 1.5 |
  | χ → N₄N₄ (M = 120) | 1.2 H_i | 5.0×10⁻³ H_i | 0.42 | 1.5 |

  - f_edge = √(mΣ/8π) is where μ5,max reaches 4πf.
  - The species cover successive amplitude ranges. The production of each switches off when μ5,max < πM²/m.
  - Perturbative couplings: y_p = M/f ≤ 0.21.
- **χ → σσ and χ → ss.** χ → σσ is closed (m_σ > m_χ/2). χ → ss has Γ ≈ 10⁻¹³ H_i.
- **The σ tail.** σσ → χχ annihilation (g²) is negligible: Γ/H ~ 10⁻⁷ f_σ. σ's tail is the N₁ perturbative width,
  0.29 H_i.
- **Spectator decay [E].** y_s s N̄_sN_s with y_s = 4.6×10⁻³. This is the largest y_s that keeps s light, |δm_s²| ≤ 0.01 H².
  - The channel unblocks at a ≈ 31, with Γ_s = 1.3×10⁻⁶ H_i.
  - s decays at a ≈ 1.8×10³ with energy fraction 2.4%. That is 1.8% entropy injection into the SM, and **no hidden heat**:
    s → YY is closed (m_s = 3 < 2m_YL) and s → X is loop-only.
  - The release itself sees s only redshifting (mf_ideal's "no decay" option).
- **CP-sign bookkeeping [D/A].** The holder is the only s-odd term.
  - g_s²χ²s², ξ_sRs², m_s0²s² and the s decay coupling are all even under s → −s. The decay coupling is even with
    N_s → iγ5N_s, as inf_reheat's S1.
  - So, by the spurion argument (inf_reheat R1), the CP sign is α + sign(s_*). The mirror history is not the CP image,
    because μ breaks it softly: u = 2μ − 10i [Q].
  - The σ/χ sector violates CP by itself (V_SR odd terms, axial phases). It is sequestered from the hidden sector and does
    not enter the DM sign.
  - The rotated axial couplings give a lepton chemical potential φ̇/2f during reheating. A spontaneous lepton asymmetry is
    a by-product with its own sign; it is not estimated.

### 4.5 Post-inflation history [C on [E] rates]

- **Radiation fraction.** 54% at t = 0.3, 77% at 0.6, 85% at 3. A χ remnant then sits at 13–19% (ε_H ≈ 1.9) and
  decays by t_99 ≈ 890 H_i⁻¹.
- **K.** K = a²H/H_i is 1.51 at t = 0.3, 1.64 at t = 3, and converges to **2.048 (HFa) / 2.047 (HFb)**.
- **Temperatures.** T at radiation–matter equality ≈ 2.2×10¹⁴ GeV. The radiation-extrapolated T at a = 1 is
  T_K = 4.35×10¹⁴ GeV.
- **Equation of state.** In the fluid stage, w = ⅓ f_rad (condensates w = 0). In the homogeneous stage it is exact.
- **Higgs-coupled plasma [E].** The N decays give ρ_SM/ρ ≥ 10⁻³ at τ = 0.64 and ≥ 50% at τ = 7.3. This uses the lab rate
  Γ_N M/E_N, with E_N ≈ ¾μ5,max redshifting.
  - That is earlier than the τ ≈ 8–10 deadline [Q higgs_reheat].
  - N–L–H scatterings, not included, would make it faster.
- **Comparison.** The engine's M1r (3, 0.5) history has K = 1.99 [Q]. The quiet single-field TB benchmarks have
  K ≈ 2000 [Q].

### 4.6 The spectator s [C/E]

- **Turn-on.** s is frozen during inflation: s/s_* = 1.000 to t = 0.15. It turns on over t_on ≈ 0.2 H_i⁻¹:
  - (a) m_s = g_s|χ| oscillates 2.3–3.5 H_i, then settles at 3.00;
  - (b) m_s² = 9 − ¾R, oscillating 2.8–3.3 H_i.
- **Evolution.** s/s_* = 0.54 at t = 0.6 and −0.25 at t = 1.06. It is exact to t = 60, then the WKB envelope (flag
  `s_is_envelope`).
- **Energy fraction.** 1.35×10⁻⁵ after the turn-on, rising ∝ a in radiation to 1.6×10⁻⁴ at t = 100 and 2.4% at its decay.
  There is no early s domination.
- **Isocurvature.** H/2πs_* = 2.83×10⁻⁶, so β_iso = 3.8×10⁻³ p². Planck 2018's uncorrelated CDM-isocurvature bound
  (Table 14 of 1807.06211; the exact entry was not read here, ≈ 10⁻²–4×10⁻² [R]) then needs p ≲ 1.6–3.

### 4.7 Gravitational effects on the hidden sector at the abrupt end [C/E]

- **R(t).** It is exported exactly through the homogeneous stage: 12 → −4.5, then oscillating at ≈ 2m_χ. In the fluid stage
  the export gives the cycle mean ρ_m/M_P² and the oscillation amplitudes 3ρ_σ/M_P² (at 2m_σ = 1200 H_i) and
  3ρ_χ/M_P² (at 2m_χ = 500 H_i).
  - The ξ_X = −5/6 jump of R is part of the release, which mf_chain runs on this R(t).
- **X pairs from the 2m_χ oscillation.** For ξ_X R X², the production is ≈ 10² H_i⁴ per Hubble time, i.e. 5×10⁻⁴ of
  M1r's E_X [E].
- **Minimally coupled hidden fields** (Y, Σ, A): (m/M_P)⁴-suppressed.
- **The Y resonance.** None of the background frequencies (2m_χ = 500, 2m_σ = 1200, m_s = 3 H_i) is near ω_YL = 22.25 or
  44.5 H_i [C]. The R spectrum has 6.6% of its power in 10–60 H_i, which is the 0.4 H_i⁻¹ window trend at 15.7 H_i, not a
  line.

## 5. Exports (`results/bench_HFa.json`, `results/bench_HFb.json`; 0.57 MB each)

- **Slice convention.** a = 1 and H = H_i at σ = σ_c: the end of slow roll, with ε_H ≈ 10⁻⁷ and R = 12H_i². This is the
  static-slice analogue of M1r's end-of-inflation slice. ε_H = 1 follows at t = 0.19 with H = 0.997 H_i, so the ε_H = 1
  slice differs by a = e^{0.19} and H by 0.3%.
- **K.** K = a²H/H_i in radiation, converged.
- **Columns** (3249 rows, t from −6.1 to 3×10⁴ H_i⁻¹):
  - t (H_i⁻¹), N = ln a, a, H/H_i, ε_H, R/H_i²;
  - R_osc_sig and R_osc_chi;
  - σ/M_P and χ/M_P (homogeneous stage);
  - s/s_* (with the `s_is_envelope` flag) and m_s/H_i;
  - the energy fractions f_sig, f_chi, f_rad, f_SM, f_s;
  - K_running.
- **The two stages.** Rows with t ≤ 0.6 are the homogeneous two-field solution with decay friction. Rows after that are
  the fluid continuation.
- **Parameters.** All parameters are in `params` and `extra_params`. The CMB, waterfall, reheating and spectator summaries
  are in the same file.

## 6. Limits

- **Homogeneous stage.** The waterfall stage is homogeneous, with linear δχ_k. There is no lattice, no metric
  perturbations, and no back-reaction of the produced N on the fields beyond the friction term.
- **Pauli blocking [E].** n_k ≤ 1 caps the energy in the μ5 band at ≈ μ5⁴/8π² per helicity, ≈ 30–40% of ρ, until the N
  thermalise. The y_ν scattering rates are estimated at several H_i, which is not computed.
- **Fluid stage.** The rates assume a harmonic condensate of amplitude √(2ρ)/m. Fragmentation into quanta (no zero
  crossings of a coherent mode) would lower the non-perturbative rates. The remnant would then decay perturbatively, so
  K could rise. If every non-perturbative channel stopped at t = 0.6 (σ through N₁ at 0.29 H_i, χ through N₄ at
  5×10⁻³ H_i), K = 3.26 instead of 2.05 [C on the fluid, `code/k_bound.py`].
- **EFT.** μ5,max/4πf = 0.42–0.49. The axial phases wind φ/f ~ 10²–10³ times. A UV completion with radial modes below
  4πf, or a periodic potential, is not provided [A]. It is a real model-building tension for a non-periodic hybrid
  potential.
- **The one-loop treatment.** The χ loop below σ_c is dropped (it vanishes continuously). The renormalisation condition
  (μ², i.e. the tree quartic) is a model choice: at μ² = M_χ²(2σ_c) the e-fold capacity is 56 instead of 73.
- **Not computed:**
  - the release, vectors and S/Q (mf_chain);
  - Higgs survival with ξ_H = 0.1 at H_i = 1.3×10¹¹ (inf_higgs/mf_higgs);
  - the spontaneous lepton asymmetry;
  - the HQ quiet-clock variant beyond the flatness scan;
  - the hidden-sector heat from the N sector, which is zero by sequestering.
- **DM calibration.** H_i = 1.3×10¹¹ GeV is the design value. mf_chain must rescale V0 (λ, v, σ_c move with H_i at fixed
  m/H_i) if the dark matter needs another H_i.

## Code (`code/`) and reproduction (one core, ≈ 6 min)

| file | purpose |
|---|---|
| mfh_common.py | constants, parameter relations, potential with the χ loop (analytic derivatives), N_* bookkeeping |
| inflation.py | slow-roll calibration on the full one-loop potential, exact background + Mukhanov–Sasaki |
| waterfall.py | homogeneous σ, χ, s + linear δχ_k through the waterfall, with axial decay friction and radiation |
| axial.py | exact half-period fermion production for the axial coupling, LZ check → results/axial_check.json |
| reheat.py | ε(A) tables, decay rates, fluid continuation with s → K, T, plasma time |
| bench.py | benchmark driver: calibrate ↔ K iteration, MS retarget, export → results/bench_*.json, axial_tables.json |
| diagnostics.py | defects, ζ_wf, R spectrum, channel table, gravitational production, spectator → results/diagnostics.json |
| structure.py | relations, N_max flatness scan, perturbative bound → results/structure.json |
| k_bound.py | K if the non-perturbative transfer stops at the switch (perturbative widths only) |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_hybrid/code
python3 -B axial.py; python3 -B bench.py HFa HFb; python3 -B diagnostics.py HFa HFb; python3 -B structure.py
```

## Literature

- **Checked this session by web search:**
  - Adshead–Sfakianakis arXiv:1508.00891;
  - Lazarides–Panagiotakopoulos hep-ph/9506325;
  - Clesse–García-Bellido arXiv:1501.07565 (a 2022 review, arXiv:2211.05767, disputes its PBH yield);
  - Felder et al. hep-ph/0012142 and Felder–Kofman–Linde hep-th/0106179;
  - Lyth arXiv:1012.4617;
  - Podolsky–Felder–Kofman–Peloso hep-ph/0507096.
- **Not confirmed here:** hep-ph/0002151 and 1102.5612 (in the brief's list). They are not used.
- **Recalled [R]:**
  - Linde astro-ph/9307002; Copeland et al. astro-ph/9401011;
  - Greene–Kofman hep-ph/9807339; Felder–Kofman–Linde hep-ph/9812289;
  - Lebedev–Westphal arXiv:1210.6987;
  - Planck 2018 arXiv:1807.06211; BICEP/Keck arXiv:2110.00483; ACT DR6 arXiv:2503.14452.
- **Our tracks, read-only:** inf_critic, inf_chain (+ skeptic), inf_reheat (+ skeptic), BRIEF_DECISIVE.
