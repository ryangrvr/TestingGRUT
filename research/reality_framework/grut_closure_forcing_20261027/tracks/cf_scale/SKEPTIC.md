# Skeptic review of track cf_scale (closure in scale)

27 October 2026 · Claude (cf_scale_skeptic), issue #13 parallel pass · same-session AI review, not outside
replication. Git was used read-only and nothing was committed. Code and results are in `skeptic/`.

Labels: [D] derived, [C] computed with own code, [Q] quoted, [A] argued, [R] recalled. Nothing here is a validation, an
exclusion of the real world or a new theory. Every data comparison is **retrospective**, and every postulate is
**adopted by choice**.

## Bottom line

**None of the four verdicts is refuted.** The headline survives: closure in scale gives GRUT no distinctive, testable
content.

Every central number reproduces under an independent method. Four statements need correction:

1. **P4-M forces less than claimed at one end of its own A_λ range.**
   - The claim "g² = λ_φ = 0 forced" (which kills Candidate 1) holds only at fixed points where the nonminimal
     coupling ξ is irrelevant.
   - GRUT's own action contains ξ_X = −5/6 and allows ξ_φ. At the pure-gravity fixed point (A_λ = 2.83, inside the
     track's scan) ξ is relevant, and gravity regenerates the φ–X portal through G²ξ_φξ_X.
   - Where ξ is irrelevant (the SM-content fixed point, A_λ = 0.024), P4-M has a consequence the track missed: it
     forbids GRUT's curvature mass, ξ_X = −5/6.
   - The heavy-window kill (κ = 0, protected by U(1)_X) is robust in both cases.
2. **P4-S fixed-point counts are wrong.** Homotopy continuation finds 8, 12 and 16 real fixed points at A_λ = 1, 2.83
   and 10. The track reported 7, 6 and 4. The conclusion holds and is in fact a theorem for every A_λ ≥ 0.
3. **The P4-SM pull is not a fixed number.**
   - The comparator λ(M_P) = 0 is the large-A_λ (Shaposhnikov–Wetterich) limit.
   - At the track's own SM-content A_λ = 0.024, the critical-surface value depends on how g₂ and g₃ run above M_P.
     The direct-m_t pull then ranges from −3.0σ to −7.8σ.
   - The quoted −2.95σ is the most favourable case, not the prediction.
4. **Three readings of "closure in scale" were not examined:** total asymptotic freedom, classical scale invariance,
   and the multiple-point principle.
   - None gives distinctive, testable GRUT content.
   - The first two are inconsistent with the heavy-window branch as specified.

| ID | Track verdict | Refuted? | Verdict after review | Main correction |
|---|---|---|---|---|
| P4-M | INCONSISTENT WITH GRUT BRANCH | no | INCONSISTENT WITH GRUT BRANCH | Heavy window: robust (κ = 0). Candidate 1: only where ξ is irrelevant. Missed: ξ_X = −5/6 is forbidden there |
| P4-S | INCONSISTENT WITH GRUT BRANCH | no | INCONSISTENT WITH GRUT BRANCH | Analytic for every A_λ ≥ 0 and e* > 0. Real fixed points are 8/12/16, not 7/6/4. Also e(h) = 0.284, outside the window's e range |
| P4-SM | STANDARD | no | STANDARD | Pull ranges from −3.0σ to −7.8σ (direct m_t) or −1.6σ to −5.3σ (cross-section m_t), depending on A_λ and the transplanckian treatment |
| P4b | STANDARD | no | STANDARD | All numbers reproduce exactly. "Window survives" is shown for the mediator-debris gate only |

## 1. Independent recomputations

**β functions [C], `skeptic/code/mv_beta_check.py`.**
- **Method.** This is a different method from the track's.
  - The quartic tensor λ_abcd comes from exact 4-fold finite differences of V₄.
  - The β functions use the Machacek–Vaughn / Luo–Wang–Xiao one-loop tensor formula
    16π²β_abcd = Λ²_abcd − 3Λ^S_abcd + 3A_abcd.
  - The result is projected by least squares onto the 11 operators.
- **Agreement.** It matches all 11 of the track's closed forms to **1.4×10⁻¹⁴** (relative) at 12 random coupling points,
  with projection residual ≤ 3×10⁻¹⁵.
- **Inhomogeneous terms.** s = (6, 486, 7776, 108, 432, 3888)·e⁴ for (λ_X, λ_Y, λ_Σ, l_XY, l_XS, l_YS). They are zero
  for κ and for the φ sector.

**P4-M critical surface, closed form [D, checked C].** For small e, linearising the transplanckian flow with
e ∝ k^(−f_g) gives

λ_i(M_P) = −s_i e⁴ / [16π²(A_λ + 4f_g)],  λ_i(μ < M_P) = λ_i(M_P) − s_i e⁴ ln(M_P/μ)/16π².

This reproduces every small-e row of the track's `as_flow.json` to ≤ 1.8×10⁻⁴. For example, λ_Σ(M_Y) = −1.437×10⁻²⁵ at
A_λ = 0.0237 and e = 10⁻⁷. The track's numbers are correct, and the closed form shows that their size is set by
s_Σ = 6q_Σ⁴ = 7776.

**Monotonicity [D, C].**
- The envelope-theorem argument is correct: at the constrained minimiser, Σq²n·∇V₄ = 4λ_min Σq²n², so
  β_V(n*) ≥ 0 whenever λ_min ≤ 0.
- An independent check with a different minimiser (SLSQP on the sphere) and the MV β functions passed in **60 of 60**
  cases. The minimum derivative was +9.5×10⁻⁷, and the finite-difference derivative agreed.
- This is a one-line consequence of the gauge-only structure (Tr M_A⁴ ≥ 0, no −y⁴ terms). It is the mechanism of
  de Brito et al. (2312.02086), not a new result.

**P4-S: complete solution count [C], `skeptic/code/p4s_homotopy.py`.**
- Total-degree homotopy with 2⁶ = 64 paths and the γ-trick. All 64 paths end at finite, distinct complex solutions at
  every A_λ. The Bézout count is saturated, so the list is complete.

| A_λ | 0.0237 | 0.1 | 0.5 | 0.78 | 0.79 | 1 | 2.83 | 10 |
|---|---|---|---|---|---|---|---|---|
| real fixed points (homotopy) | 0 | 0 | 0 | 0 | 8 | **8** | **12** | **16** |
| real fixed points (track multistart) | 0 | 0 | 0 | 0 | — | 7 | 6 | 4 |
| real fixed points with λ_Σ > 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| largest λ_Σ* | — | — | — | — | −1.71 | −0.816 | −0.200 | −0.052 |

- The multistart missed the fixed points with |λ| ≈ 8–125.
- The largest λ_Σ* values agree with the track (−0.82, −0.20, −0.05).
- **Stronger, analytic statement [D].** The λ_Σ equation is 20λ² + (16π²A − 432e²)λ + C = 0, where
  C = 7776e⁴ + l_XS² + l_YS² (+2l_φΣ² + λ_HΣ²) > 0. For any A_λ ≥ 0 a real root requires
  16π²A − 432e² ≥ √622080·e², and then both roots are negative.
- So for every A_λ ≥ 0 and every e* > 0, no real fixed point has λ_Σ > 0, and V₄ is unbounded below along pure Σ.
  Only fermionic −y⁴ terms in C can evade this.
- The multistart count is therefore irrelevant to the verdict, but the numbers 7/6/4 should not be cited.

**P4b [D, C], `skeptic/code/p4b_and_missed.py`.**
- Closed form for the common λ_Y = λ_Σ: 1/λ₀ = 20 ln(M_P/M_Y)/16π² + 1/λ_crit. This gives:
  - **0.5661** for λ_crit = 4π;
  - **0.3722** for λ_crit = 1;
  - 0.554 under the unitarity criterion 8π/3. The band is criterion-dependent at the ±2% level.
- The pole for λ = 1 is at M_Y·exp[(1 − 1/4π)·16π²/20] = **5.731×10¹⁵ GeV**.
- The κ floor is (κ/[λ_X/(1.5·2^(−1/3))]^(3/4))⁴ = **0.2415**.
- The own integration reproduces Candidate 1's maxima: λ_X 0.1741, λ_φ 0.1931, g 0.5698, m_s/m_X 6.215, and
  e(M_Y) ≤ 0.2511.
- The release table matches `release_check.json` row by row. The corrected passes ×55.7–65.7 are 0.93 times the table
  values.

**Data [Q].**
- pdglive, 2026 edition, fetched 27 October 2026: m_H = 125.13 ± 0.11 GeV (S = 1.5); m_t (direct) = 172.60 ± 0.27 GeV
  (S = 1.5); m_t (pole, from cross sections) = 172.1 ± 0.6 GeV. These agree with the track.
- α_s = 0.1180 ± 0.0009: the pdglive QCD page returned HTTP 500, so this value was not re-verified at the primary
  source.
- Buttazzo et al. eqs. (55)–(61) (1307.3536, fetched): the λ(M_Pl) fit, its slopes and the absolute-stability bound
  M_t < 171.53 ± 0.42 GeV (their eq. 66) are as the track quoted.

## 2. Postulate by postulate

### P4-M (standard asymptotic safety): not refuted. Forcing claim 1 is overstated.

**(a) What is robustly forced [D].**
- **κ = 0.** β(κ) ∝ κ, and U(1)_X rephasing is broken by no other operator of dimension ≤ 4.
  - The scalar basis in ACTION.md (at 04c53f5) includes the Z2_φ-odd terms a_iφ|z_i|² and μ₁φΣ*Y², and the
    curvature terms ξ_i|z_i|²R. All of them are U(1)_X-invariant.
  - The X⁶Σ* torque (α₆ ∝ μκ²) and any independent Σ*X⁶/Λ³ term vanish on a fixed point that respects global
    symmetries.
- **The heavy window cannot exist under P4-M.** This holds independently of everything below.

**(b) Smuggled assumption: ξ.**
- The postulate card lists the supplied relevant couplings as "masses, μ, G, Λ and e". It omits ξ.
- **GRUT's own content.** THEORY.md §2 at 04c53f5 specifies ξ_X = −5/6 as the covariant representative of the
  inflationary −10H²|X|² mass. ACTION.md's F also allows ξ_φφ² and ξ_H.
- **Eichhorn–Pauly 2009.13543 [Q, ar5iv fetched]:**
  - eq. (41): β_{λ_HP} ⊃ g²ξ_dξ_v[160/(1−2λ)³ + 864/(3−4λ)³] + O(ξ²ξ) + O(ξ²ξ²). A portal is regenerated whenever
    both nonminimal couplings are nonzero.
  - ξ is **relevant** at the pure-gravity fixed point ("two relevant directions … scalar mass and non-minimal
    coupling") and **irrelevant** with SM matter.
- **My recomputation [C], from their eq. (32), gravitational part:** θ_ξ = +1.108 at (1.32, 0.13); −0.0185 at
  (7.63, −9.97); −0.0278 at the U(1)_D point.
- **Consequences.**
  - **Where ξ is relevant (pure gravity, A_λ = 2.83).** ξ_X = −5/6 is allowed, and a supplied ξ_φ ≠ 0 is allowed too.
    Then g², λ_φ and λ_X (and the Higgs portals via ξ_H) are not zero on the critical surface. They are functions
    of (ξ_φ, ξ_X, g*, λ*). So "g² = 0, Candidate 1 impossible" is **not established** there [A].
    - I did not evaluate the size: the overall normalisation of eq. (41) as rendered was not independently verified.
    - Since ξ_d²ξ_v²-type terms also enter the self-couplings, the claim "λ_Σ ≤ 0" is likewise not established there.
  - **Where ξ is irrelevant (SM-content or U(1)_D fixed point).** The track's forcing holds. In addition, ξ_X is driven
    to its fixed-point value. That value is ≈ 0 for a hidden scalar with no Yukawas, up to O(q²e²) corrections.
    So **P4-M also forbids GRUT's ξ_X = −5/6**. This element of the forbidden baseline was missed.
  - The weak irrelevance (θ_ξ ≈ −0.02) makes this exact only as a UV-limit statement, the same caveat the track
    attached to A_λ.

**(c) "Candidate 1: no φ binding" uses an incomplete argument [A].**
- ACTION.md's basis allows the Z2_φ-odd trilinear a_Xφ|X|². It is relevant at the SM-content fixed point
  (θ ≈ 1 − 0.024), so P4-M does not forbid it. The track's Z2_φ-even basis omitted it.
- It does not obviously rescue binding, though. With λ_X forced to ≈ 0⁻, integrating out φ gives an attractive X
  quartic −a_X²/(2m_φ²)|X|⁴ with nothing positive to stop collapse, since |X|⁶ is irrelevant.
- So the conclusion "no stable Candidate-1-type object" plausibly survives, but for this reason rather than only
  because g² = λ_φ = 0. The track should state it.

**(d) What actually kills Σ breaking [C].**
- A quartic λ_Σ ≈ −1.4×10⁻²⁵ is not a stability problem. The Fubini bounce exponent 8π²/(3|λ|) ≈ 1.8×10²⁶.
- What kills Σ breaking is that the RG-improved potential has no stationary point below M_P.
  - The bracket 4λ + β_λ = (s e⁴/16π²)[1 − 4/(A + 4f_g) − 4ln(M_P/μ)] is negative at every μ ≤ M_P whenever
    A_λ + 4f_g < 4.
  - At A_λ = 10 the first stationary point is at 0.86 M_P, never near h.
  - So there is no Coleman–Weinberg rescue at h. This confirms the track's conclusion with a sharper reason.

**(e) Relabelling.**
- P4-M's hidden-sector content is the asymptotic-safety "flattening" of the Planck-scale potential (Eichhorn,
  Hamada, Lumma, Yamada 1712.00319 [R]) together with Eichhorn–Pauly 2009.13543 and de Brito et al. 2312.02086. The
  track credits all three.
- The sign of A_λ is robust within the quoted truncation, because 104λ² − 144λ + 51 has discriminant −480 < 0 [D].
- The magnitude for GRUT's matter content (SM + 7 real scalars + 1 vector) was not computed. The track used the SM
  value.

**(f) Minor.**
- ε_BX = 0 "by hidden C" holds only with exact Z2_φ. With μ₁φΣ*Y² present, arg(μ₁/μ₀) is invariant under rephasing
  and breaks hidden C [D].
- Under P4-M, μ₁ is irrelevant and has no source, so it vanishes and the conclusion holds.

### P4-S (e at its interacting fixed point): not refuted, and strengthened

- The conclusion is a theorem for every A_λ ≥ 0 and every f_g > 0 (§1). The track's multistart counts (7/6/4) are
  incomplete; the true counts are 8/12/16.
- **An independent reason to exclude the heavy window [C].**
  - e(M_P) = e* = 0.318 runs down to e(h) = 0.284. With f_g calibrated at the reduced Planck mass instead, the values
    are 0.313 and 0.280.
  - This lies outside the window's own range: e ≳ 3×10⁻⁶ and the Floquet-clean bound e ≤ 0.1 quoted in the track.
  - Even with no quartic, the window fails.
- **Calibration-scale effect.** Using g_Y at 2.435×10¹⁸ GeV (own 2-loop) instead of at 1.22×10¹⁹ GeV changes:
  - f_g by −3.2% (9.52×10⁻³);
  - e* by −1.6% (0.313);
  - A_crit to 0.758.
  The track's note said "a 3% effect on e*"; the effect is about 3% on f_g. It is immaterial.

### P4-SM (λ_H(M_P) ≈ 0): not refuted. data_status is overstated as a single pull.

**Reproduction [C], `skeptic/code/sm_two_loop.py`.**
- Two-loop SM running (with a 3-loop g₃ term) from Buttazzo's NNLO matching at M_t gives λ(M_Pl) = −0.01506 at their
  reference point, against their NNLO −0.0143. The 0.0008 difference is the known size of the 3-loop running.
- Slopes:
  - ∂λ/∂M_t = −0.00649 (theirs −0.0066);
  - ∂λ/∂M_h = 0.0028 (theirs 0.0029);
  - per 0.0007 in α_s: 0.00185 (theirs 0.0018).
- At PDG 2026, the own 2-loop value is λ(M_Pl) = −0.0114 (direct m_t) and −0.0081 (cross-section m_t). The track's
  NNLO-formula values are −0.0105 ± 0.0036 and −0.0072 ± 0.0046.
- Moving to the actual switch scale 2.435×10¹⁸ GeV shifts λ by only −0.00025.
- The thresholds for absolute stability, λ(M_Pl) = 0 and λ(reduced M_P) = 0 are 170.80, 170.83 and 170.80 GeV
  (2-loop). The track's "same near-criticality relation as P_q" is numerically right.
- Conceptually, though, P_q is a one-sided bound (m_t ≤ 170.8). Asymptotic safety with A_λ > 0 gives an equality.

**The comparator was not fixed by the postulate [C].**
- Under the postulate (β_λ += A_λλ above M_P, sharp switch), the prediction is the critical-surface value
  λ_crit(M_P; A_λ). This is 0 only when A_λ dominates the SM one-loop terms.
- With every other SM coupling frozen at its M_P value, λ_crit ≈ 0 for every A_λ (+3×10⁻⁴ at A_λ = 0.024). The reason
  is the near-criticality coincidence: at 2.4×10¹⁸ GeV, −6y_t⁴ (−0.1376) and the gauge term (+0.1368) almost cancel.
- But g₂ and g₃ are asymptotically free with f_g > 0, so above M_P they decrease, as on any Gaussian branch. The
  cancellation then fails and λ_crit becomes positive:

| A_λ | 0.0237 (track's SM-content value) | 0.056 | 0.1 | ≥ 0.3 |
|---|---|---|---|---|
| λ_crit(M_P), frozen | +0.0003 | +0.0001 | +0.0001 | ≈ 0 |
| λ_crit(M_P), g₂ and g₃ asymptotically free, y_t and g_Y at fixed points | **+0.0171** | +0.0061 | +0.0025 | ≤ 4×10⁻⁴ |
| pull, direct m_t | −3.1σ to **−7.8σ** | −3.0σ to −4.7σ | −3.0σ to −3.7σ | −3.0σ to −3.1σ |
| pull, cross-section m_t | −1.7σ to −5.3σ | −1.6σ to −2.9σ | −1.6σ to −2.2σ | −1.6σ to −1.7σ |

- These pulls are retrospective, measured against the formula value shifted to reduced M_P.
- The y_t-on-Gaussian-branch variants are unstable against the UV starting depth (several branches), so they were not
  used.
- **Conclusion.** For A_λ ≳ 0.3 the track's −2.95σ / −1.56σ is right. At the A_λ the track itself assigns to SM matter,
  the comparator is model-dependent and the tension is **larger**, up to −7.8σ. The verdict STANDARD is unchanged,
  since the relation is Shaposhnikov–Wetterich's.
- data_status should give a range, and should state that λ(M_P) = 0 is the large-A_λ limit.
- de Brito et al.'s λ_H* = 5.56×10⁻³ (track: −4.5σ) lies inside this range.

### P4b (perturbative and bounded below to M_P): not refuted

- All numbers reproduce exactly (§1), and STANDARD (EFT hygiene) is correct.
- Two caveats should be visible in the digest:
  1. "The window survives in the band" is demonstrated for the mediator-debris gate only, as the track's own caveat
     says. The Z6-string and Floquet gates depend on λ_Σ (through m_Σ/m_A) and were not rerun at λ = 0.25–0.57.
  2. The band edges depend on the perturbativity convention: 0.566 (4π), 0.554 (8π/3) and 0.372 (|λ| ≤ 1).

## 3. Readings of "closure in scale" the track did not examine [A, D, C]

**TAF (total asymptotic freedom: every coupling → 0 in the UV, no gravity fixed point needed).**
- Credit: Giudice, Isidori, Salvio, Strumia 1412.2769; Holdom, Ren, Zhang 1412.5540 [R].
- An abelian factor has a Landau pole, so the hidden U(1) must be embedded in a non-abelian group, as must U(1)_Y.
- **Minimal SU(2) embedding of charges (1, 3, 6), with q = 2T₃ [D, C]:**
  - X sits in a doublet, Y in j = 3/2 and Σ in j = 3, plus a real adjoint to break SU(2).
  - b₀ = 22/3 − (½ + 5 + 28)/3 − 2/6 = **−4.17**, so the theory is **not** asymptotically free. q = T₃ gives even
    larger representations.
- **Verdict:** INCONSISTENT WITH GRUT BRANCH (minimal embedding; larger groups not explored); the principle itself is
  STANDARD.

**Classical scale invariance ("no exterior yardstick": no dimensionful input).**
- Credit: Coleman–Weinberg 1973; Gildener–Weinberg 1976; Bardeen 1995; Hempfling hep-ph/9604278; Meissner–Nicolai
  hep-th/0612165; Iso–Okada–Orikasa 0902.4050; Hambye–Strumia 1306.2329 [R].
- Σ breaking would have to be radiative.
- **Landau gauge, MS-bar [D, C]:** at the minimum, λ_Σ(μ = M_A) = (4/3)·3(6e)⁴/64π² and
  m_Σ² = 3(6e)⁴w²/8π².
  - At e = 10⁻⁷ and w = h = 5.33×10¹² GeV this gives λ_Σ = 8×10⁻²⁸, m_Σ = **0.37 GeV** and m_A = 3.2×10⁶ GeV.
  - The heavy spot needs λ_Σ ≥ 0.24, which requires e ≈ 0.4. That is outside the window's e range.
  - Generating the EW scale through a portal needs λ_HΣ ~ (m_h/h)² ≈ 5×10⁻²², which is unobservable.
- **Verdict:** INCONSISTENT WITH GRUT BRANCH (heavy spot); the principle is STANDARD; the portal is UNTESTABLE.

**Multiple-point principle: λ(M_P) = β_λ(M_P) = 0.**
- Credit: Froggatt–Nielsen hep-ph/9511371 [R].
- It belongs to the same SM near-criticality family as P4-SM and P_q (thresholds within 0.03 GeV at 2 loops).
- **Verdict:** STANDARD.

## 4. New issues for the track

1. ξ is a supplied relevant coupling at the pure-gravity fixed point, and GRUT already contains ξ_X = −5/6. P4-M's
   forcing of g² = λ_φ = 0 (and of λ_Σ ≤ 0) is therefore fixed-point-dependent. Only κ = 0 is unconditional.
2. Missed P4-M consequence: where ξ is irrelevant, the covariant inflationary mass ξ_X = −5/6 is forbidden.
3. P4-S real fixed-point counts are wrong: 8/12/16, not 7/6/4. The result is analytic, so cite the theorem instead.
4. P4-S has an additional exclusion: e(h) = 0.284, outside the heavy window's e range.
5. P4-SM: λ(M_P) = 0 is only the large-A_λ limit. At A_λ = 0.024 the pull is −3.1σ to −7.8σ (direct m_t), depending on
   the transplanckian gauge flow. The P_q relation is a bound, not an equality.
6. The Z2_φ-odd basis (a_Xφ|X|², μ₁φΣ*Y²) of ACTION.md was not considered. It does not change the verdicts, but the
   "no φ binding" argument should cite the induced unbounded −a²/(2m_φ²)|X|⁴ rather than only g² = 0.
7. Under P4-M the negative quartics are harmless as a vacuum-decay issue (S_E ~ 10²⁶). What kills Σ breaking at h is
   that no stationary point exists below M_P while A_λ + 4f_g < 4.
8. Three closure-in-scale readings were missed: TAF, classical scale invariance and the multiple-point principle.
   None gives distinctive, testable GRUT content.

## Files

Code, all own:
- `skeptic/code/mv_beta_check.py`
- `skeptic/code/p4s_homotopy.py`
- `skeptic/code/sm_two_loop.py`
- `skeptic/code/p4b_and_missed.py`
- `skeptic/code/summary_numbers.py`

Results:
- `skeptic/results/mv_beta_check.json`
- `skeptic/results/p4s_homotopy.json`
- `skeptic/results/sm_two_loop.json`
- `skeptic/results/p4b_and_missed.json`
- `skeptic/results/summary_numbers.json`
- `skeptic/results/skeptic_summary.json`

Reproduce from `tracks/cf_scale/` with `python3 -B skeptic/code/<script>.py skeptic/results/<name>.json`. Runtimes:
mv_beta_check 3 min, homotopy 15 s, sm_two_loop 2 s, p4b_and_missed 3 s.

## Sources

- Eichhorn–Pauly 2009.13543 (ar5iv fetched: θ formula, eq. 11/12 fixed points, eq. 32, eq. 41, ξ relevance,
  shift-symmetry quote).
- Buttazzo et al. 1307.3536 (ar5iv fetched: eqs. 55–61, 64–66).
- PDG 2026 pdglive (m_t, m_H; fetched 27 October 2026).
- de Brito et al. 2312.02086 (as quoted by the track, not re-fetched).
- Recalled [R]:
  - Eichhorn, Hamada, Lumma, Yamada 1712.00319;
  - Shaposhnikov–Wetterich 0912.0208;
  - Giudice, Isidori, Salvio, Strumia 1412.2769;
  - Holdom, Ren, Zhang 1412.5540;
  - Froggatt–Nielsen hep-ph/9511371;
  - Meissner–Nicolai hep-th/0612165;
  - Iso, Okada, Orikasa 0902.4050;
  - Hambye–Strumia 1306.2329;
  - Hempfling hep-ph/9604278.
