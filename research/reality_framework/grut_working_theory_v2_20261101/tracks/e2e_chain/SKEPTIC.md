# Skeptic of e2e_chain and e2e_relic (track key e2e_chain_skeptic)

3 November 2026 · Claude, end-to-end pass "M2" (BRIEF_E2E.md; GitHub issue #13 context). This is an adversarial review
of `tracks/e2e_chain` and `tracks/e2e_relic`.

Same-session AI work on one supplied model. **Every history modification, closure, holder amplitude and phase,
interpolation rule, relic rule and criterion used here is a model choice; none follows from closure.** Nothing here
is a validation, an exclusion of the real world, a new theory or a GRUT prediction. Data comparisons are
retrospective.

Labels: [D] derived, [C] computed here, [Q] quoted from a track, [E] estimate, [A] argued, [R] literature (arXiv numbers
verified by web search this session).

**Code and results.** Everything is in `skeptic/code/` and `skeptic/results/`; scratch is
`.../scratchpad/wt2/e2e_chain_skeptic/`.
- e2e_chain's `chain.py` and `e2c_common.py` were imported read-only and unchanged. Through them so were mf_ideal's
  engine and drivers and e2e_relic's `pmap.py` and `relic.py`. Only the scratch path was redirected.
- e2e_chain's `scan.py` was not used, because it writes into e2e_chain/results. `sk_engine.py` is our own driver.
- The final histories were read read-only from e2e_chain's scratch.
- Python ran with `-B`. No binary and no `__pycache__` is in the tree.
- **Regression [C].** Our engine runs reproduce e2e_chain's digits exactly at 11 nodes. For example, hHFa_GM at
  c0 = 1.55 gives S/Q = 4.109785 and e_crit = 1.88×10⁻⁶; hH84_P0 at c0 = 1.45 gives 2.422570.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict in brief

**What stands [C, recomputed]:**
- the K values of the final histories, their convergence and the H continuity at t_sw;
- N_* = 51.55 / 51.76 and r = 2.05 / 1.59×10⁻⁷, from our own entropy bookkeeping;
- the isocurvature floor s_* ≥ 1.13 / 1.36×10¹⁶ GeV, with D = 5.4–6.3 from finite differences;
- M2-X (GM, c0 = 1.55, δ = 0) and M2-Q (P0, c0 = 1.45, δ = 0), which are genuine single-branch roots;
- the relic arithmetic: object mass, radius, flux, impulse, and the free-X cross sections;
- e2e_relic's pair-fraction gap fill. It is reproduced exactly, and changing the relaxed-start interpolation rule moves
  it by ≤ 0.5%.

**What is refuted or must change [C]:**
1. **"M2-Q exists only on the local-blocking closure P0" is refuted under the chain's own relic rule.** On the GM
   closure with the fragment spread σ = 0.3 (K = 12.56), a finer phase scan finds a self-consistent Q-matter pocket:
   - c0 = 1.485–1.52 and δ = +2…+3.5°;
   - H_i = 1.20–1.35×10¹¹ GeV;
   - S/Q(10⁻⁶) = 3.28–3.47 and e_crit ≈ 0.9–1.0×10⁻⁶.

   The coarse scan had missed it. Its quoted "floor 3.49 with the fragment spread" is 3.30. The pocket is bulk-only:
   e_crit is 1.4–1.5×10⁻⁶ for Q = 2×10⁵ objects (§3).
2. **The GM floors are lower than quoted.** On plain GM they are 3.575 at 1.3×10¹¹ (quoted 3.78) and 3.26 at 8.4×10¹⁰
   (quoted 3.49). The N_YL interference minimum is 0.04–0.08 at arg μ_eff ≈ 106–107°, not 2.7–5.6 near 108°. On plain
   GM, M2-Q still has no root (§3).
3. **The fixpoint rule mixes branches.** At P0, c0 = 1.40, the reported root (S/Q 3.458, "M 0.874") is not a root of
   either branch. Done branch by branch, there are **two** self-consistent solutions there:
   - free X at H_i = 1.04×10¹¹ GeV, S/Q 3.56;
   - Q-matter at 1.17×10¹¹ GeV, S/Q 3.41.

   Near the branch boundary the jump in M (from 1.00 to about 0.81, which raises H_req by about 15%) gives two roots
   or none (§2).
4. **The M2-Q phase window is wider than stated.** Negative δ had never been run. At c0 = 1.45, δ = −4°…+4° are all
   self-consistent Q-matter (§4).
5. **The naturalness analysis omits the dominant term.** The held hidden state makes s feel V_min(|u(s)|), the tree-level
   term whose slope at c0 = 1 mf_hybrid quoted as a "tadpole". At M2-X that term has curvature m²_eff = −6.3×10⁻⁴ H²
   (with s_* at the floor). That is 20 to 14,000 times the residual masses e2e_chain used, and of opposite sign.
   - Locally V_min ∝ −c0^6.3, and asymptotically it goes to −7.48|u|⁴.
   - So no Starobinsky–Yokoyama equilibrium can sit in the window, for any quadratic or quartic s potential.
   - Boundedness needs a new supplied quartic, λ_s ≥ 7.48 g_h⁴ ≈ 10⁻¹⁴.
   - The stated P(window) numbers do not apply. The conclusion that c0 is supplied survives, and is stronger (§5).
6. **Smaller items:**
   - The μ_eff components quoted in the NOTES and the digest are wrong (§1).
   - e2e_kinetic's HFa m_s² formula in the kinetic block has the wrong sign for the anharmonic shift; the effect is
     immaterial, ΔS/Q ≤ 7×10⁻⁴ (§1).
   - The interpolation rule carries a ±2% systematic in H_i (§2).
   - The detector-rate geometry factor is 1–2 (§6).
   - The Q = 2×10⁵ seed question applies to every relic with S/Q between 2.93 and 3.47 (§6).

---------------------------------------------------------------------------------------------------------------------

## 1. History splice and conventions (`sk_static.py` A–C; `sk_engine.py hist … anh`)

**Conventions [C].**
- The code's holder is u₀ = μ + c0 e^{iδ}(10i − μ), with μ = 7.159456. It gives:
  - M2-X: μ_eff = −3.938 + 15.500i (|μ| 15.99, 104.25°);
  - M2-Q: μ_eff = −3.222 + 14.500i (14.85, 102.53°).

  `m2_points.json` agrees.
- The NOTES and digest text instead give "−2.38 + 15.81i" and "−1.66 + 14.85i". Those components correspond to
  98.6° and 96.4°, so they are wrong; the moduli and angles quoted beside them are right.
- H_i is the Hubble rate at σ = σ_c, and the reduced M_P = 2.435×10¹⁸ GeV.
- K = a²H/H_i with a = 1 at the slice, which is the same as ρ_r a⁴ = 3K²M_P²H_i².
- N_* recomputed independently, as ln(a_slice H_i/k_*) from entropy conservation with g = 106.75: 51.554 / 51.765,
  against the track's 51.554 / 51.765.

**Final histories [C].**
- K = a²H at the end reproduces every stored value: 15.279, 8.772, 12.561, 15.255, 16.422, 9.207.
- K changes by 3–4×10⁻⁵ over the last decade in t.
- The tables end at t ≈ 1.1–1.7×10⁵, not 10⁸. K has converged by 10⁴, so this is immaterial.

**The m_s² step at t_sw.**
- In e2e_kinetic's export (and e2e_chain's `kin_block`), the HFa spectator mass in the kinetic block is
  m_s² = g_s²(v² + ρ_χ/m_χ²). That treats χ as oscillating harmonically about v.
- In V = λ/4(χ² − v²)² the cubic term shifts the mean: ⟨χ − v⟩ = −(3/2v)⟨δχ²⟩, so ⟨χ²⟩ = v² − 2ρ_χ/m_χ² [D].
- The shipped value is therefore high by 3x = 0.6% at t_sw, where x = ρ_χ/(m_χ²v²) = 2.1×10⁻³. It is part of the +1.1%
  step from the lattice's 8.920 to 9.019. The step is +1.2% at the 8.4×10¹⁰ GM node and +2.1% at the 8.4×10¹⁰ P0 node.
- Re-running with the corrected kinetic-block m_s² (`*_anh`) [C]:
  - M2-X node: S/Q 4.1098 → 4.1091, e_crit unchanged at 1.88×10⁻⁶;
  - 8.4×10¹⁰ node: 3.8584 → 3.8568;
  - c0 = 1.50, δ = +4°: 3.7790 → 3.7784.
- **Not a refutation:** the release is complete by t ≈ 3, so the post-splice m_s² does not matter. The correction
  should still be applied to the export.
- What does matter is the lattice-stage m_s². Its running-mean residual scatter is 0.4–1.6% (std over t = 2–3), and
  the track already showed that a ±3% shift moves M2-X's e_crit by ×5.

---------------------------------------------------------------------------------------------------------------------

## 2. Self-consistency done branch by branch (`sk_static.py` D)

**The flaw in e2e_chain's rule.**
- `fixpoint.py` interpolates f = ln(H_req/H_i) linearly in ln H between the two nodes. Each node's f uses its own
  branch's M: about 1.000 for free X, about 0.80 for Q-matter.
- The root is then classified by the interpolated S/Q.
- When the nodes sit on different branches, f mixes two relic functions, and the "root" belongs to neither.

**Branch-by-branch rule [C on the model choices].**
- For b = free X and b = Q-matter: f_b(node) = ln[H_req,node (M_b/M_node)^(−2/3)/H_node]. M_b is evaluated at the
  node's S/Q with e2e_relic's own functions (p_selfconsistent, or a calibrated relic()).
- A root counts only if the interpolated S/Q there lies on branch b.

| point | e2e_chain root | free-X root | Q-matter root | count |
|---|---|---|---|---|
| GM 1.55, 0° (M2-X) | 1.123e11, S/Q 4.025 | **1.123e11, 4.025 ✓** | 1.27e11, 4.10 ✗ | 1 |
| P0 1.45, 0° (M2-Q) | 9.886e10, 2.433 | 8.67e10, 2.43 ✗ | **9.886e10, 2.433 ✓** | 1 |
| P0 1.40, 0° | 1.126e11, 3.458, "M 0.874" | **1.037e11, 3.556 ✓** | **1.172e11, 3.410 ✓** | **2** |
| GM+σ0.3 1.50, 0° / +4° / 1.48, +3° / 1.55, +3° | not computed | ✗ | ✗ | **0** |

- **M2-X and M2-Q stand.** M2-Q's free-X counter-root would need S/Q > 3.47 at 8.7×10¹⁰, where S/Q is 2.43.
- **P0 c0 = 1.40 is bistable.** It is not the "marginal" single root that was reported.
- **Bands with no self-consistent abundance.** On GM+σ0.3 the abundance cannot be matched at all at several holder
  values next to the pocket. Both effects come from the step in M at S/Q = 3.4716, which the relic rule makes sharp
  (there is no mixed phase). That sharpness is a model artefact [A] that downstream tracks should state.
- **The interpolation rule.** Linear in H instead of ln H moves M2-X to 1.146×10¹¹ (+2.0%) and M2-Q to 1.012×10¹¹
  (+2.3%). That is a ±2% systematic in H_i (±4–5% in r) on top of the "two nodes only" caveat.

---------------------------------------------------------------------------------------------------------------------

## 3. The GM floor, and an M2-Q pocket on GM + fragment spread (`sk_engine.py`; `results/sk_engine.json`)

63 engine runs in all (17 s each), including 7 c0 = 1 bases and the regressions. The new ones step δ by 0.5–1° on the unmodified final histories.

**N_YL near its interference minimum (c0 = 1.50).**

| node | δ = +2° | +2.5° | +3° | +3.5° | +4° |
|---|---|---|---|---|---|
| GM 1.3×10¹¹, N_YL | 2.77 | 1.44 | 1.51 | 2.93 | 5.63 |
| GM 1.3×10¹¹, S/Q(10⁻⁶) | 3.578 | **3.575** | 3.606 | 3.673 | 3.779 |
| GM 8.4×10¹⁰, N_YL | 0.39 | 0.041 | 0.95 | 3.04 | 6.25 |
| GM 8.4×10¹⁰, S/Q(10⁻⁶) | **3.264** | 3.278 | 3.325 | 3.407 | 3.528 |
| GM+σ0.3 1.3×10¹¹, S/Q(10⁻⁶) | 3.303 | **3.300** | 3.327 | 3.389 | 3.486 |

At c0 = 1.515–1.52, δ = +3° the 1.3×10¹¹ node gives N_YL = 0.08.

**Plain GM.**
- The floors are 3.575 (1.3×10¹¹) and 3.264 (8.4×10¹⁰), not 3.78 and 3.49.
- Done branch by branch, every scanned (c0, δ) still has only a free-X root (S/Q ≥ 3.59).
- **"No M2-Q on plain GM" stands**, with a smaller margin.

**GM + fragment spread σ = 0.3 (the track's own GM bracket, K = 12.56).** There is a self-consistent Q-matter pocket
[C at the 1.3×10¹¹ node; E for the root shift]:

| c0, δ | S/Q at node | H_req/H_i at node | root H_i (GeV) | S/Q at root |
|---|---|---|---|---|
| 1.485, +3° | 3.441 | 1.040 | 1.35e11 | 3.468 |
| 1.50, +2° / +2.5° / +3° / +3.5° | 3.303 / 3.300 / 3.327 / 3.389 | 0.963 / 0.974 / 0.985 / 0.996 | 1.255 / 1.268 / 1.281 / 1.295e11 | 3.28 / 3.28 / 3.32 / 3.39 |
| 1.515 / 1.52, +3° | 3.330 / 3.348 | 0.935 / 0.919 | 1.22 / 1.20e11 | 3.30 / 3.31 |

- No 8.4×10¹⁰ GM+σ0.3 history exists. The root shift therefore uses plain GM's slopes of f and S/Q at the same (c0, δ),
  from our matching 8.4×10¹⁰ runs [E]. The roots lie within 1.5–7% of the computed node.
- **Representative point, c0 = 1.50 and δ = +3°** (μ_eff = 15.06 at 106.8°) [C/E]:
  - H_i = 1.28×10¹¹ GeV, K = 12.6;
  - S/Q(10⁻⁶) = 3.32 and e_crit = 0.91×10⁻⁶ (bulk);
  - M = 0.813 m_X and ⟨Q⟩ ≈ 1.5×10²⁰, from the extrapolated 3.2–3.47 band of the Q-matter tables;
  - D ≈ 8.0, so s_* ≥ 1.8×10¹⁶ GeV; r = 2.7×10⁻⁷.
- **Caveats on the pocket:**
  - It is narrow: Δc0 ≈ ±0.02 and Δδ ≈ ±0.75°.
  - It is bulk-only: e_crit,Q=2×10⁵ is 1.4–1.5×10⁻⁶, so Q = 2×10⁵ seeds are not stable at e = 10⁻⁶ (see §6).
  - It sits next to the no-root band of §2.
- **Consequence for e2e_chain's claims.** "The blocking closure decides whether M2-Q exists" is too strong. Three things
  decide it together:
  - the fragment spread at hand-over;
  - phase resolution at the ≤ 1° level;
  - the relic's seed rule.

---------------------------------------------------------------------------------------------------------------------

## 4. Windows (`sk_engine.py`; `sk_static.py` D)

**M2-Q (P0, c0 = 1.45).**
- Previously δ < 0 had never been run on P0. Both nodes are Q-matter for δ = −4°, −2°, 0, +2° and +4°.
- Root S/Q values are 3.30, 2.72, 2.43, 2.59 and 3.26. Root H_i values are 9.0, 9.5, 9.9, 10.3 and 10.8×10¹⁰ GeV.
- e_crit is 0.43–0.83×10⁻⁶.
- **So the phase window is at least −4°…+4°, not 0…+4°.** The c0 window's lower end at 1.40 is bistable (§2).

**M2-X (GM).**
- On the final histories, δ = −2° and −4° at c0 = 1.55 are new runs.
- δ = −2°: root S/Q = 4.38, e_crit = 3.8×10⁻⁶.
- δ = −4°: S/Q = 4.87, and S/Q does not cross 3.47 by e = 10⁻³ at the 1.3×10¹¹ node.
- **The window "−2°…+6°" survives on the final histories.** Across it, e_crit ranges over 1.3–3.8×10⁻⁶ (root values: −2°: 3.8×10⁻⁶; +3°: 1.28×10⁻⁶;
  +4°: 1.32×10⁻⁶; +6°, 1.3×10¹¹ node only: 1.9×10⁻⁶). "e_crit ≈ 2×10⁻⁶" holds only near δ ≈ 0, c0 1.55–1.6.

---------------------------------------------------------------------------------------------------------------------

## 5. Naturalness of c0: the spectator's real inflationary potential (`sk_natural.py`; `results/sk_natural.json`)

**What e2e_chain did.** It put s in Starobinsky–Yokoyama (SY) equilibrium [R, arXiv:astro-ph/9407016, verified:
P ∝ exp(−8π²V/3H⁴)] with V = ½m²s² only, where m² is mf_hybrid's residual masses (4.5×10⁻⁸ to 3×10⁻⁵ H²). Its numbers
are arithmetically right: s_rms = 35.6–919 H, N_relax = 10⁵–6.7×10⁷, ratios 110–2830.

**What it omitted [C].**
- The heavy hidden fields track their held stationary point, so s also feels V_min(|u|), with u = μ + c(10i − μ) and
  c = |g_h| s/|10i − μ|.
- We computed this with the dq_heat skeptic's exact solver at R = 12. V_min depends on |u| only; four phases at
  |u| = 10 agree to all digits.
- mf_hybrid quoted the slope of this term at c = 1 as a "tadpole" of 0.6 H³. Our slope agrees in order (1.2 H³ secant).

| | M2-X (c0 1.55) | M2-Q (c0 1.45) |
|---|---|---|
| V_min, dV/dc, d²V/dc² (H⁴) | −1.88e5, −7.73e5, −2.66e6 | −1.23e5, −5.32e5, −2.17e6 |
| local power law V_min ∝ −c^n | n = 6.36 | n = 6.23 |
| m²_eff = V''/s_u², s_* at the floor (bench s_* 7.3e15) | **−6.3e-4 H²** (−1.5e-3) | **−2.4e-4 H²** (−8.4e-4) |
| ratio to the residual masses 3e-5 / 2.5e-7 / 4.5e-8 | 21 / 2500 / 14000 | 8 / 970 / 5400 |
| classical drift per e-fold; Δc/c over 50 e-folds | 4.0 H; 0.20% (bench: 6.1 H; 0.47%) | 1.9 H; 0.07% |
| asymptote V_min/\|u\|⁴ | → −7.48 (c = 40) | same |

**Reading.**
1. **The residual masses are not the relevant curvature.** The tree-level hidden-sector term is 8–14,000 times larger
   and tachyonic.
2. **There is no equilibrium in the window.** V_min falls faster than c⁴ there (n ≈ 6.3). Any balance against a
   quadratic or quartic s term is therefore a maximum of V_eff, not a minimum: V'' + V'/c and V'' + 3V'/c are both < 0.
3. **The tree potential is unbounded below along s.** V_min → −7.48 (g_h s)⁴ unless a quartic λ_s s⁴ with
   λ_s ≥ 7.48 g_h⁴ is added. That is 9.5×10⁻¹⁵ (M2-X at the floor), 5.5×10⁻¹⁴ (bench s_*) and 2.1×10⁻¹⁵ (M2-Q).
   - This is a new supplied parameter for the card. A one-loop s⁴ term is of order g_h⁴/64π², far too small to supply it.
   - It is harmless for the release, because 12λ_s s² ≪ m_s² = 9 after inflation.
4. **With λ_s above the bound,** the bounded potential's minimum lies outside the window:
   - near c ≈ −0.1…−0.2 for Λ = λ_s/g_h⁴ ≥ 20, where |u| ≈ μ and the holder is nearly off;
   - at large negative c for Λ ≈ 7.5–10.

   The window sits 2–13×10⁶ H⁴ above the minimum, so P_SY(window) ~ exp(−6×10⁷) or smaller.
5. **So:**
   - e2e_chain's "P(window) ~ exp(−5×10³); best case 3–5% with g_h free" does not apply;
   - the "tadpole drift 10 H in 50 e-folds" is the c = 1 bench value; at M2-X it is 200–300 H;
   - c0 is still frozen to ≤ 0.5% through the observable e-folds, and s_* remains a non-equilibrium initial condition
     **[A, strengthened]**.

   The isocurvature tilt from m²_eff is n_iso − 1 ≈ −4×10⁻⁴: negligible.

---------------------------------------------------------------------------------------------------------------------

## 6. e2e_relic (`sk_pmap.py`, `sk_static.py` E)

**Pair-fraction gap fill [C].** The relaxed rows at S/Q = 7, 10 and 14 were re-run through fx_chain's network
(read-only):

| S/Q | e2e_relic's rule (linear in R) | x0 interpolated in ln(1 + R) | x0 + 0.03 |
|---|---|---|---|
| 7 | 0.0220114 | +0.10% | +0.52% |
| 10 | 0.129946 | +0.06% | +0.16% |
| 14 | 0.321496 | +0.02% | +0.07% |

- e2e_relic's rows reproduce to all digits.
- The relaxed onset's x0(R) is concave across the gap (slopes 3.6, 1.5 and 1.0×10⁻³), so linear interpolation in R is
  a choice. It is immaterial (≤ 0.5%).
- The start convention (relaxed vs a_ref) remains the real model choice. **Not refuted.**

**Relic arithmetic at M2-Q [C].** Own recomputation:
- object mass 0.775866 ⟨Q⟩ m_X = 9.61 μg (0.44 M_Pl), R = 0.730 nm, n = 4.91×10⁻¹⁴ m⁻³;
- impulse 3.6×10⁻²⁷ kg m/s;
- M = 0.80188;
- free-X σ/m = 6.00×10⁻⁸ (XX) and 7.07×10⁻⁷ (XX̄) cm²/g, from our own threshold amplitudes, −4λ + 4g²(1/3) and
  −4λ + 4g²(2/3).

All match. The analytic untruncated SHM mean speed is 358 km/s, against the truncated Monte Carlo's 354.

**Corrections and new issues in e2e_relic:**
1. **Rate geometry.** n⟨v⟩A is the rate through a cross-section A: a sphere, or a plane normal to the wind. A flat
   plane in an isotropic flux counts n⟨v⟩A/2, which is 0.27 rather than 0.55 m⁻² yr⁻¹ at M2-Q. Quote the geometry.
2. **Detector context [R, verified at abstract and search-snippet level].** The Windchime white paper
   (arXiv:2203.07242) describes a purely gravitational Planck-mass search as needing about 80 dB of quantum-noise
   reduction with a 100³ array. Carney et al. (arXiv:1903.00492) state "Planck mass or heavier".
   - The impulse scales with the object mass, so 0.44 M_Pl objects are 0.44 times the Planck-mass signal.
   - "At the lower edge of the stated target" should read "below the stated Planck-mass target".
   - No sensitivity is claimed either way.
3. **The Q = 2×10⁵ seed rule (new).**
   - relic() calls Q-matter for every S/Q < 3.4716, on the argument that ripened objects are bulk-like [E].
   - For 2.93 < S/Q < 3.47, Q = 2×10⁵ fragments are themselves unstable. Ripening to 10²⁰ then needs a nucleation or
     seed path that nobody has computed.
   - Affected: M2-X at 3×10⁻⁶ (S/Q 3.05), the GM+σ0.3 pocket of §3, and the c0 = 1.36, 3×10⁻⁶ row (3.28).
   - M2-Q (2.43) is not affected.
4. **The extrapolated band.** S/Q 3.2–3.47 is a linear extrapolation of the Q-matter tables [E, flagged by e2e_relic].
   It is now relevant, because it contains every GM-closure Q-matter root.

**Not checked:** the β exponents, ⟨Q⟩ ∝ Y_Q, f_v(Y_Q), and the P(k) proxy.

---------------------------------------------------------------------------------------------------------------------

## 7. Isocurvature and the spectator (`sk_static.py` C)

- **The floor recomputes exactly.** s_* ≥ D H/(2π√(βA_s/(1 − β))) gives 1.1327×10¹⁶ / 1.3601×10¹⁶ GeV, with
  β_iso < 0.038 for uncorrelated CDM isocurvature with n = 1, Planck 2018 X [R, arXiv:1807.06211; the 0.038 value is
  verified through search snippets].
- **D cross-check [C].** From the c0 grid, with Ω ∝ (H_req)^(−3/2) at a fixed node:
  - GM 1.3×10¹¹: D = 6.29 (c0 1.50–1.55), 5.41 (1.55–1.60), 5.77 (1.534–1.565);
  - P0 1.3×10¹¹: D = 7.38 (1.45–1.50).

  So D at the root is uncertain at about ±10%, and s_* at about ±10%.
- **Spectator spread.** δs/s_* = 6×10⁻⁴ × dS/Q/d ln c0 reproduces exactly. The "s sees the box mean" argument [A] is
  accepted.
- **The s energy fraction at decay,** 0.16% / 0.60% at the floor, reproduces.

---------------------------------------------------------------------------------------------------------------------

## 8. Corrections downstream tracks should apply

1. **μ_eff components.** M2-X is −3.938 + 15.500i, M2-Q is −3.222 + 14.500i.
2. **GM floors.** 3.575 at 1.3×10¹¹ and 3.26 at 8.4×10¹⁰ (plain GM), and 3.30 with the fragment spread. The N_YL minimum
   is 0.04–0.08 at arg μ_eff ≈ 106–107°.
3. **"M2-Q only on P0" → "M2-Q on P0 robustly; on GM+σ0.3 in a narrow bulk-only pocket"** (c0 1.485–1.52, δ +2…+3.5°,
   H_i 1.20–1.35×10¹¹, S/Q 3.28–3.47). It is absent on plain GM.
4. **Use branch-aware self-consistency.** P0 c0 = 1.40 has two roots (1.04×10¹¹ free X; 1.17×10¹¹ Q-matter). Several
   GM+σ0.3 holders have none. Report the M-discontinuity artefact.
5. **M2-Q phase window: at least −4°…+4°** at c0 = 1.45. M2-X window −2…+6° with e_crit 1.3–3.8×10⁻⁶.
6. **Add a ±2% interpolation-rule systematic** to H_i (±4–5% to r).
7. **Replace the SY naturalness numbers** by the hidden-sector-dominated potential (§5).
   - Add λ_s ≥ 7.48 g_h⁴ (≈ 10⁻¹⁴) to the parameter card as supplied.
   - State that the window is never an equilibrium point, and that the classical drift is 200–300 H over 50 e-folds.
8. **Kinetic-block HFa m_s²:** use g_s²(v² − 2ρ_χ/m_χ²). The effect on S/Q is ≤ 7×10⁻⁴.
9. **e2e_relic:** state the rate geometry; correct "lower edge of target" to "below the Planck-mass target"; flag the
   seed rule for 2.93 < S/Q < 3.47.

## 9. New issues

- **N1.** The spectator potential is unbounded below at tree level without λ_s. It is a required supplied coupling.
- **N2.** The window in c0 is a slope of V_eff during inflation, not an equilibrium. The drift points away from it
  (outward without λ_s, inward with it).
- **N3.** The sharp branch rule makes the self-consistency map discontinuous: bistable or empty bands near
  S/Q = 3.47. A mixed-phase relic function is needed to make M continuous.
- **N4.** Q-matter at 2.93 < S/Q < 3.47 needs a seed or nucleation path.
- **N5.** The branch at e = 10⁻⁶ is sensitive to δ at the 1° level, because N_YL varies from 0.04 to 38 over a few
  degrees. Phase scans coarser than 1° can miss Q-matter pockets.

## 10. Limits

- The pocket roots use plain-GM slopes [E]. No GM+σ0.3 history exists at 8.4×10¹⁰.
- No new histories were built, and no lattice or kinetic re-run was done.
- β, Y_Q exponents and P(k) were not re-checked.
- The literature was verified at abstract or search-snippet level:
  - Starobinsky–Yokoyama, astro-ph/9407016;
  - Carney et al., 1903.00492;
  - Windchime, 2203.07242;
  - Planck 2018 X, 1807.06211 (β_iso < 0.038, n_s);
  - ACT DR6, 2503.14452 (0.9743 ± 0.0034);
  - Baxter et al., 2105.00599 (identity only; halo values as recalled).

**Files.**
- Code: `skeptic/code/sk_common.py`, `sk_engine.py`, `sk_static.py`, `sk_natural.py`, `sk_pmap.py`.
- Results: `skeptic/results/sk_engine.json` (every run), `sk_static.json`, `sk_natural.json`, `sk_pmap.json`,
  `sk_summary.json`.
