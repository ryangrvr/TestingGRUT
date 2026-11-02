# Skeptic of track fx_chain: thermal history at e = 10⁻⁶ and the leftover p for GRUT-M1

1 November 2026 · Claude (fx_chain_skeptic), repair pass of BRIEF.md in this directory (issue #13 context). This is
same-session AI work on one supplied model, not outside replication. Nothing here is observational evidence, a
validation, an exclusion or a GRUT prediction. Every repair (relaxed start, formation rule, e = 10⁻⁶, ε = 0.002,
Q = 2×10⁵) is a supplied model choice, not a consequence of closure. Labels: [D] derived, [C] computed, [Q] quoted,
[E] estimate, [A] argued. Git was read-only and nothing was committed. Code is in `skeptic/code/` and results in
`skeptic/results/` (232 kB in total).

## Bottom line

**Every number fx_chain computed reproduces**, most of them with code written fresh for this check: the Hartree
plasma (to 6 digits), the spinodal dispersion (an analytic expansion agrees to 3×10⁻⁴), the circular transition (an
independent integrator), the ss → XX̄ threshold factor (by hand), and the network values p = 9.1023×10⁻⁴, 2.1622×10⁻³,
3.2162×10⁻³ and f = 0.0045855.

The physical picture behind items 4 and 5 does not hold at e = 10⁻⁶, however.

- The network's own trajectories and an independent ideal-Bose criterion agree: at the entropy per charge of the
  relaxed history, s/n_Q = 4.94, **no condensate and no object can exist between T ≈ 1 m_X and T ≈ 0.05–0.1 m_X**.
- The lumps that form at the spinodal (T_h ≈ 10 m_X) evaporate.
- X̄ freeze out (T ≈ 0.1 m_X) in a pure gas, so there is no "capture before freeze-out".
- Bound objects become allowed again only below 0.05–0.1 m_X, from a non-degenerate gas (fugacity 0.087). Making
  them needs a nucleation work of about 1.7×10⁵ T.
- So the η ≥ 0.85–0.92 conditions describe a population this history does not produce. The expected value is
  η ≈ 0, and the frozen test then fails by ×6.9.
- Objects formed early survive the whole cooling only if R_tot ≲ 1.5–2.8, i.e. **e ≳ 2.2–6×10⁻⁶** (band).

## Verdicts on the track's results

| id | fx_chain | Recomputation | Verdict |
|---|---|---|---|
| R1 | R_tot = 4.9665 / 5.7033; sources done 915–2×10⁸× before onset; about 1.5×10³ collisions | Values read back from kin_source. **The cascade vectors scale as e^{−1/2}, not 1/e**: they have a fixed number N_YH, and R ∝ m_A a_decay ∝ e·e^{−3/2}; the interface gives 0.10716 (10⁻⁶) and 0.33888 (10⁻⁷), a ratio of √10. The collision count is conservative: ⟨1/E⟩ of a Bose gas is 0.684/T, not 1/(2.70T), a factor 1.85 | CONFIRMED (scaling statement CORRECTED) |
| R2 | X̄/charge 0.417/0.572, T_h/m 8.28, a_on 3.81×10¹², n_Q 151, V(0)/n_Q 0.63 | New Hartree code (adaptive quad, joint 4-unknown root): T_conf 21.886915, X̄ 0.417243, condensed fraction 0.917373, a_on 0.541233 a_ref, T/m 8.281162, n_Q 151.022, the same at every R. With the tachyonic −κa² kept inside the φ Hartree mass, a_on moves +1.0% (R = 4.97), V(0)/n_Q becomes 0.648 and T/m 8.20; this is negligible | CONFIRMED |
| R3 | n_onset = 164.3 m_X³ = untilted spinodal; M1 crossover | The untilted hom_scan spinodal is 164.317 = √λ_X(λ_φv²/g²)^{3/2}. h(S) = 3C(1−S)²/S + t_lin/S + B(1−S²) is strictly decreasing on (0, 1), so there is one root | CONFIRMED |
| R4 | bulk crossed at 0.68 a_on, spinodal at 0.84 a_on | Crossings reproduced from the profile. Note: the quoted φ/v = 0.15 and 0.34 at the crossings are the T = 0 values; with the history's thermal masses they are 0.14 and 0.28 | CONFIRMED (wording) |
| R5 | Γ = 0.774δ, k = 0.958√δ; Q ~ 4×10¹¹–7×10¹² | Small-k expansion z(c + 4Ω²D) = −ck² − (D + 2aF²)k⁴ with c = 2F²(aD − 2S²): Γ/δ = 0.77462, k/√δ = 0.95791, and c = 0 coincides with min Ω. The hidden X collision rate (about 7×10⁻⁷ m_X at T_h ≈ 10) is below Γ_peak (3.5–6×10⁻⁶), so the collisionless cold-fluid estimate is self-consistent [E]. **But at e = 10⁻⁶ these lumps evaporate by T ≈ 0.8–1.0 m_X (see below), so the estimate does not set the final population** | CONFIRMED as a transient; relevance REFUTED at R = 4.97 |
| R6 | circular: n₋/q < 3×10⁻⁷ "measure floor" | Independent scipy DOP853 integration with an instantaneous free-field n₋: 3.08×10⁻⁷ (10⁻³) and 2.82×10⁻⁷ (10⁻⁴), against fx_chain's 2.84×10⁻⁷; δ = 0.01 gives 1.327×10⁻⁴ (fx_chain 1.307×10⁻⁴). **The "floor" is the start's own radial action.** The tilted minimum χ₀ = 1.08 at M₀ shifts the radial equilibrium by δ_eff = 4.8×10⁻⁴, so J_r/q = 2.8×10⁻⁷, and the floor scales as ε_t⁴ (1.9 / 2.8 / 4.0×10⁻⁷ at 0.0018 / 0.002 / 0.0022). A start at the true equilibrium gives **8.9×10⁻¹⁰**. Not load-bearing for p: the forced-equilibrium start resets X̄ anyway | CONFIRMED (interpretation CORRECTED) |
| R7 | V‴ = 6(λ−ε)v, threshold factor 0.9025, p +5.8% | By hand: M = −0.02 − 0.02 + 0.01333 = −0.02667 (untilted) and −0.02533 (M1), a ratio of 0.95² = 0.9025. no-cubic rerun: 8.6071×10⁻⁴ | CONFIRMED |
| R8 | p = 9.10×10⁻⁴ (band), 2.16×10⁻³ (broad); a_ref 3.22×10⁻³ | Reran the copied engine: 9.1023×10⁻⁴, 2.1622×10⁻³, 3.2162×10⁻³; repro f = 0.0045855. **Mechanism:** p is a function of the entropy per charge, d ln p / d ln(s/n_Q) = 9.3–9.8. The a_on start conserves the pre-onset Hartree entropy to 4×10⁻⁴; isentropic starts at 1.5 a_on and 2 a_on give 9.070 / 9.073×10⁻⁴. The old a_ref recipe carried +13.8% entropy (V(0) dumped as heat at low density), hence ×3.53. A sudden late release at 1.2 / 1.5 a_on gives 1.01 / 1.44×10⁻³, and dropping V(0) gives 7.85×10⁻⁴ | CONFIRMED (placement ≡ entropy) |
| R9 | rule b: p < 10⁻¹⁸, robust to capture ×10⁻⁷ | Floor reproduced. **But the trajectory has objects absent for T ∈ [0.049, 0.968] m_X (Q = 2×10⁵) and [0.092, 0.793] (Q = 3×10⁷).** The X̄ freeze out at 9.7×10⁻⁴ (T ≈ 0.07–0.1) and are captured only after the rule re-creates objects at T < 0.05–0.09. With capture ×0 or ×10⁻¹² the result is 9.71×10⁻⁴, and binding objects give 9.57×10⁻⁴ | Number CONFIRMED; "capture before freeze-out" REFUTED |
| R10 | 2.13356p = 0.01 at R = 6.44 → e ≥ 7.4×10⁻⁷ / 8.7×10⁻⁷ | With R_casc ∝ e^{−1/2}: 7.354×10⁻⁷ / 8.675×10⁻⁷, the same at the quoted precision. η = 1 is not available there, though (objects cannot survive for R ≳ 2–2.8) | CONFIRMED (conditional) |
| R11 | η ≥ 0.879 / 0.919 / 0.8503 at e = 10⁻⁶ | The arithmetic is exact (0.87934, 0.91934, 0.85026). **The premise fails.** At s/n_Q = 4.94 no reference objects can exist for T ∈ [0.05, 1.0] m_X, and re-forming them from the gas needs W ≈ 7×10⁴(E/Q − μ_gas) ≈ 1.7×10⁵ T. The physically expected η at e = 10⁻⁶ is ≈ 0, which gives a test term of 0.0687 (×6.9) [C for absence, E for nucleation] | REFUTED (as the remaining condition at e = 10⁻⁶) |
| R12 | p subdominant; η open; slow passage → Q ≫ 2×10⁵, "a calibration change rather than an η" | (i) The slow-passage lumps evaporate at R = 4.97, so they do not set the final Q. (ii) If Q ~ 10¹² did survive, E/Q = 0.777 changes the test coefficients, not only Y_Q: A = m_X/ε_Q − 1 rises from 0.0668 to 0.287 and B from 2.134 to 2.574. With the recalibrated Y_Q = 6.128×10⁻⁹, p = 8.46×10⁻⁴, so η ≥ 0.973 (gas) / 0.965 (capture) | REFUTED (both parts) |

## The main correction: objects cannot survive the cooling at e = 10⁻⁶

**Criterion [D/C].**
- Objects (or a BEC) hold the gas at μ_X = ω and carry charge but no entropy. Chemical equilibrium then gives
  μ_X̄ = −ω and μ_s = 0.
- The gas must carry all of the conserved entropy (zero portals), so objects can exist at T only if
  s/n_Q < σ(T) ≡ s_gas(T, ω)/q_gas(T, ω).
- Ideal Bose gases with vacuum masses (X, X̄ at m_X; s at √6 m_X), computed with own code in `survival.py`, give:

| holder | ω | σ_min (at T/m_X) | absent at R = 4.97 for T/m_X in | R_crit (survives below) | e_min band / broad |
|---|---:|---:|---|---:|---|
| homogeneous BEC | 1 | 1.317 (T → 0) | (0, 1.33] | 0.09 | none (mediators alone exceed it) |
| binding object Q = 7.13×10⁴ | 0.940 | 2.553 (0.29) | [0.025, 1.11] | 1.52 | 6.0 / 7.1×10⁻⁶ |
| reference Q = 2×10⁵ | 0.882 | 2.928 (0.32) | [0.049, 0.99] | 2.01 | 3.5 / 4.1×10⁻⁶ |
| Q = 3×10⁷ | 0.794 | 3.384 (0.34) | [0.089, 0.83] | 2.64 | 2.3 / 2.7×10⁻⁶ |
| bulk (Q ≳ 10¹⁰, the slow-passage lumps) | 0.776 | 3.472 (0.34) | [0.097, 0.80] | 2.76 | 2.15 / 2.5×10⁻⁶ |

- s/n_Q(R) of the pre-onset plasma: 1.228 (R = 0), 2.058 (0.91), 2.250 (1.14), 2.918 (2.0), 4.936 (4.97),
  5.393 (5.70), 21.8 (41.7).
- R_crit comes from interpolating this curve. e_min uses R_tot(e) = 0.8056 + 4.0538(10⁻⁶/e) + 0.10716(10⁻⁶/e)^{1/2}.

**The network agrees run by run [C].**

| Run (copied engine) | Condensed fraction is zero for T/m_X in | Criterion above |
|---|---|---|
| gas-only, R = 4.97 (ends with 100% free X) | ≤ 1.328 | 1.328 |
| instant Q = 2×10⁵ | [0.049, 0.968] | [0.049, 0.99] |
| instant Q = 3×10⁷ | [0.092, 0.793] | [0.089, 0.83] |
| binding objects | [0.027, 1.09] | [0.025, 1.11] |

- At R = 1.85 the instant objects survive (minimum fraction 0.022 at T = 0.33).
- At R = 2.15 they vanish for T ∈ [0.22, 0.45], which brackets R_crit = 2.01.
- Kinetics does not save them [E]. A Q = 10¹² lump at T ≈ 0.5 m_X evaporates in about 10⁷–10⁸/m_X, compared with
  1/H ~ 10¹⁵/m_X. The gas mean free path exceeds the object spacing, so evaporation is not diffusion-limited.

**Re-formation [E].**
- Below T ≈ 0.05–0.1 m_X the gas is supersaturated with respect to objects. It is non-degenerate, though: with all
  charge in the gas, μ = m_X − T(s/n − 5/2), the fugacity is e^{−2.44} = 0.087 and stays fixed under adiabatic expansion.
- The smallest bound objects (Q ≥ 5.7–7.1×10⁴) cost W ≈ 7×10⁴(E/Q − μ) ≈ 1.7×10⁵ T at every T from 0.09 down to
  0.005 m_X. Thermal nucleation is therefore absent.
- A mean-field gas–liquid instability is also absent: n|λ_eff| ~ 10⁻⁷ ≪ T.
- "Instant objects" (rule b) at T < 0.05 m_X is therefore a strong supplied assumption, and it is the only reason
  rule b gives p → 0.

## Other new issues

1. **p is set by s/n_Q alone, with exponent ≈ 9.7 [C].** Isentropic starts at any position give the same p
   (to ≤ 0.4%), and the a_on start is isentropic to 4×10⁻⁴. The relaxed placement works only by not adding entropy.
   Any entropy production at the crossover or fragmentation of ≥ 1% raises p by ≥ 10%. The η = 1 limit
   (p = 4.69×10⁻³) is reached at +18% entropy per charge.
2. **Gas-only is not the strict upper value of p [C].** Objects that do not capture give 9.57–9.71×10⁻⁴ (+5–7%,
   from the formation heat of the instant rule).
3. **The gas-only rule's homogeneous BEC is inconsistent with the track's own item 2.** It is spinodally unstable
   from 0.84 a_on. It also evaporates at T = 1.33 m_X anyway (R_crit = 0.09), so p = 9.1×10⁻⁴ is in practice the p of a
   pure gas, which is physically what happens at R = 4.97 [C].
4. **The large-Q calibration tightens the test [C].** At E/Q = 0.777: η ≥ 0.965–0.973. It also raises H_i by
   (6.128/5.079)^{2/3} = 1.133 under the M1 coordinated scaling [D, not followed through].
5. **The heavy-release normalization is the pre-M1, pre-holder point.** The holder (+1.8% N_X, if adopted) at fixed E_X
   would lower s/n_Q by about 1.4%, which means p −13% and R_crit up slightly [E, not computed].
6. **The clock convention.** heat_onset uses pt_relax's early clock a = a_m + Cs (C = 1.986). kin_source's late
   background has H ≈ H_i/a² (C → 1), which would double the available conformal time. Conservative [E].

## Corrected physical picture at e = 10⁻⁶ (supplied model)

1. Relaxed pre-onset plasma: circular condensate (92% of charge) plus a Hartree gas, T_h/m_X = 8.3 at a_on, and
   s/n_Q = 4.94 per charge [C]. This is unchanged from fx_chain.
2. The M1 crossover is adiabatic. It is entropy, not timing, that sets p; the network start is isentropic [C].
3. The condensate fragments near the spinodal (T_h ≈ 6–10 m_X) into large lumps [E]. **These lumps then evaporate
   into the hot gas by T ≈ 0.8–1.0 m_X** [C].
4. X̄ freeze out through ss → XX̄ at T ≈ 0.1 m_X in a pure, non-degenerate X gas: p = 9.1×10⁻⁴ (band), 2.2×10⁻³
   (broad) [C].
5. Reference objects are allowed again below T ≈ 0.05–0.1 m_X, but they cannot nucleate thermally
   (W/T ≈ 1.7×10⁵) [E]. Without a non-thermal formation route the dark charge stays a free X gas (η ≈ 0) and the
   frozen test fails by ×6.9 [C/E].
6. **For an object population to exist, R_tot ≲ 1.5–2.8 is needed, i.e. e ≳ 2.2–6×10⁻⁶ (band) /
   2.5–7×10⁻⁶ (broad)**, depending on the object charge [C].
   - This replaces e ≥ 7.4×10⁻⁷, which assumed η = 1.
   - At such e the gas-only p is about 3–6×10⁻⁵ (log interpolation of fx_chain's R = 1.14 and 3 runs) and capture by the surviving objects removes it.
   - η is then set by re-absorption kinetics after the minimum surviving fraction (2% at R = 1.85), which is not
     computed here.

## Reproduce

```bash
S=research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_chain/skeptic
python3 -B $S/code/hartree_indep.py $S/results/hartree_indep.json          # independent Hartree + self-consistent onset (~3 min)
python3 -B $S/code/spinodal_expand.py $S/results/spinodal_expand.json      # analytic spinodal coefficients
python3 -B $S/code/circ_indep.py <out.json> <delta> <eps_rate>             # independent circular transition (20-220 s each)
python3 -B $S/code/survival.py $S/results/hartree_indep.json $S/results/survival.json
python3 -B $S/code/survival_threshold.py $S/results/survival.json $S/results/survival_threshold.json
cd $S/code/kin && python3 -B sk_entropy.py                                  # start entropy vs pre-onset plasma
python3 -B m1_kin.py m1_gas_relaxed_R497 sk_gas_mult1.5_R497_V0x0.367 sk_gas_mult2.0_R497_V0x0.172 sk_gas_mult1.5_R497 \
   sk_gas_mult1.2_R497 sk_gas_mult1.0_R497_V0x0 m1_inst_relaxed_R497_capx0 m1_bind_relaxed_R497_capx0 sk_inst_Q3e7_R497 \
   sk_inst_Q3e7_R497_capx0 sk_inst_Rx1.85 sk_inst_Rx2.15 m1_gas_relaxed_R497_YQv6.128e-09 repro_gas_R497 m1_gas_aref_R497
python3 -B sk_collect.py && python3 -B ../sk_summary.py                     # -> results/sk_kin.json, skeptic_summary.json
```

Engine files `ks_common.py`, `thermo_state.py`, `stage2_moments.py`, `moment_controls.py` are byte-identical to
fx_chain's copies of Codex 34f23f8. `m1_kin.py` and `bootstrap.py` changed only in their paths and gained variant
cases; it reads fx_chain's `heat_onset.json` and `cubic_m1.json` read-only. Run outputs are in
`$SP/wt2/fx_chain_skeptic/`.

Credits: Q-ball thermodynamics and evaporation (Laine–Shaposhnikov 1998; Cohen–Coleman–Georgi–Manohar 1986; Enqvist
et al.), condensate fragmentation (Kusenko–Shaposhnikov; Kasuya–Kawasaki), Bose thermodynamics (Kolb–Turner).
