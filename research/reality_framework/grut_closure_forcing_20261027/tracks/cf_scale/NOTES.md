# Track cf_scale: closure in scale (asymptotic-safety boundary conditions; perturbativity to M_P)

27 October 2026 · Claude, issue #13, parallel pass "what does closure force that standard physics does not?" ·
same-session AI work, not outside replication. Git was used read-only. Nothing here was committed.

**Scope.** This track is one-loop running of the GRUT hidden sector:
- the real scalar φ;
- the complex fields X, Y, Σ, with hidden charges 1, 3, 6;
- the hidden gauge coupling e;
- gravity, entered as the standard parametrised contributions from the asymptotic-safety (AS) literature, with a
  sharp switch-off at M_P = 2.435×10¹⁸ GeV.

Nothing here is a validation, an exclusion of the real world, or a new theory. The postulates are **adopted by
choice**. None of them follows from the word "closed". Every comparison with data is **retrospective**.

Labels: [D] derived, [C] computed (own code), [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## Verdicts

| ID | Postulate | Verdict | One line |
|---|---|---|---|
| P4-M | Standard AS: every quartic is irrelevant (A_λ > 0); e sits on the Gaussian branch (free, at most e*); masses and the cubic μ are free | **INCONSISTENT WITH GRUT BRANCH** | Forces g² = λ_φ = κ = 0 at all scales and every charged hidden quartic below zero, of size about −q⁴e⁴. Candidate 1 and the heavy window are both unreachable |
| P4-S | Fully predictive AS: e also at its interacting fixed point e* (calibrated on hypercharge) | **INCONSISTENT WITH GRUT BRANCH** | With charges (1, 3, 6), no real hidden-quartic fixed point exists for A_λ < 0.78. For A_λ = 1–10 every real fixed point has V₄ unbounded below |
| P4-SM | The same postulate applied to GRUT's Standard-Model sector: λ_H(M_P) ≈ 0 (Shaposhnikov–Wetterich) | **STANDARD** (retrospectively disfavoured at 1.6–3.0σ) | Predicts m_t = 171.0 GeV, or m_h = 128.8 ± 1.2 GeV. This is the same near-criticality relation as the q-vacuum P_q stability relation |
| P4b | No Landau pole: perturbative and bounded below up to M_P | **STANDARD** | Candidate 1 passes easily. The heavy spot point as stated (λ_Y = λ_Σ = 1) fails at 5.7×10¹⁵ GeV. Inside the band 0.24 ≤ λ_Y = λ_Σ ≤ 0.57 it is repaired, and the window survives there [C] |

**Bottom line.** Closure in scale either forces GRUT's dark sector to vanish (P4-M, P4-S) or forces nothing beyond
standard EFT hygiene (P4b). The weaker branch of P4 (A_λ < 0, quartics relevant) only gives bounds of P4b type. The
one data-facing consequence is the Standard-Model relation λ_H(M_P) ≈ 0. It is Shaposhnikov–Wetterich's, not GRUT's.

None of the four postulates gives GRUT distinctive, testable content. The GRUT-specific computed content is negative:
- the charge assignment (1, 3, 6) and the absence of hidden fermions make the AS-closed hidden sector flat or unstable;
- the AS-closed hidden sector cannot carry the X–φ binding (Candidate 1) or the κ charge transfer (heavy window).

## 1. Conventions and the gravity terms

Tree potential (THEORY.md §2, at 04c53f5), together with every gauge-invariant, Z2_φ-even quartic (11 operators):

V₄ = λ_X|X|⁴ + λ_Y|Y|⁴ + λ_Σ|Σ|⁴ + 2κRe(Y*X³) + l_XY|X|²|Y|² + l_XS|X|²|Σ|² + l_YS|Y|²|Σ|²
     + g²φ²|X|² + l_φY φ²|Y|² + l_φΣ φ²|Σ|² + (λ_φ/4)φ⁴.

The basis is the one in hw_gates gate 4, without the D field. Candidate 1 is g = 0.1 (g² = 0.01), λ_X = 10⁻³,
λ_φ = 0.03, m_s/m_X = √(2λ_φ)/g = √6. The heavy spot point is κ = 3.4587×10⁻³, λ_Y = λ_Σ = 1, e = 10⁻⁷, matched at
M_Y = 4.0×10¹² GeV.

Gravity above M_P (the standard parametrisation, e.g. 1707.01107 and 2312.02086 eq. 31, sharp switch):
- **Quartics:** β_{λ_i} += A_λ λ_i for every quartic monomial.
  - A_λ ≡ −f_λ in Eichhorn's convention β = −f_λ λ + ….
  - A_λ > 0 means the quartic is irrelevant (IR-attractive towards 0): λ_i(M_P) is predicted.
  - Gravity is blind to internal indices, so A_λ is universal [A, Q 2009.13543: "the gravitational contribution to the
    anomalous dimension for the scalar potential is always towards irrelevance"].
- **Gauge:** β_e += −f_g e, universal over gauge groups [Q 2312.02086 eq. 17ff].
  - With f_g > 0, e = 0 is relevant (e is free).
  - The interacting fixed point e* = √(16π² f_g / b_e) is irrelevant (e predicted). It is also the upper bound on
    e(M_P) [Q 1709.07252].

Shaposhnikov–Wetterich [Q 0912.0208] took the positive anomalous dimension A_λ > 0, "favored by explicit computations",
and predicted m_h ≈ 126 GeV. For A_λ < 0 they gave only a range, 126–174 GeV.

## 2. One-loop β functions (sympy, `code/beta_derive.py`) [D]

**Method.** The Landau-gauge RG of the quartic part of V,

16π² dV₄/dlnμ = ½Tr[(V₄'')²] + (3/2)(M_A²)² − 3e² Σ_i q_i² φ_i ∂_iV₄,  with M_A² = 2e² Σ q_z²|z|².

This is the result projected exactly on the 11-operator basis; the residual is identically 0, so the basis is
closed. In the X–φ–e sector:

    16π² β(λ_X)  = 20λ_X² − 12e²λ_X + 6e⁴ + 2g⁴   (+ 18κ² + l_XY² + l_XS² in the full model)
    16π² β(g²)   = g²(8g² + 8λ_X + 6λ_φ − 6e²)    (+ 2 l_XY l_φY + 2 l_XS l_φΣ)
    16π² β(λ_φ)  = 18λ_φ² + 4g⁴                    (+ 4 l_φY² + 4 l_φΣ²)
    16π² β(e)    = (46/3) e³                       (b_e = (1 + 9 + 36)/3)
    16π² β(κ)    = κ(12λ_X + 6 l_XY − 36e²)
    16π² β(λ_Y)  ⊃ 486e⁴,  16π² β(λ_Σ) ⊃ 7776e⁴ (= 6q⁴e⁴),  β(l_XY) ⊃ 108e⁴, β(l_XS) ⊃ 432e⁴, β(l_YS) ⊃ 3888e⁴

The full list is in `results/beta_derive.json`.

**Controls (all pass).**
- Scalar QED: 20λ² − 12e²λ + 6e⁴.
- A lone real φ⁴/4: 18λ².
- The numerical β's at the spot point match hw_gates gate 4 (`gate4_running.json`) to every printed digit: β(λ_X)
  2.7567×10⁻⁶, β(λ_Y) 0.12665, β(κ) 2.6283×10⁻⁷, β(g²) 1.6971×10⁻⁵.

**Structural facts [D, read off the closed forms].**
1. **The φ sector decouples.** {g² = l_φY = l_φΣ = 0} is RG-invariant, and with λ_φ = 0 as well the whole φ sector is
   free.
   - φ is a gauge singlet with Z2_φ and has no Yukawas, so no inhomogeneous term exists.
   - With the Higgs portal λ_Hφ included, the argument is the same [A].
   - Including gravity, the surface is protected by φ's shift symmetry at a fixed point that respects global
     symmetries [Q 2009.13543: λ_HP "is protected by shift symmetry at all scales once set to zero at high scales"].
2. **κ = 0 is RG-invariant.** β(κ) ∝ κ: with κ = 0, X appears only through |X|², so U(1)_X rephasing is an exact global
   symmetry [D].
3. **Monotonicity [D].** In the fermion-free hidden sector, below M_P,
   dλ_min/dlnμ = β_V(n*)/16π² ≥ 0 whenever λ_min ≡ min_{|n|=1} V₄(n) ≤ 0.
   - At the constrained minimiser, Σ q_i² n_i ∂_iV₄ = 4λ_min Σ q_i² n_i², so
     16π² β_V(n*) = ½Tr(V₄'')² + (3/2)M_A⁴ − 12e²λ_min Σ q²n² ≥ 0.
   - It is strict for e ≠ 0 along any charged direction.
   - At λ = 0, β_V = (3/2)e⁴(x² + 9y² + 36s²)²/16π² ≥ 0.
   - Numerical check (`code/theorem_check.py`): 40 random coupling sets with λ_min ≤ 0 all gave
     dλ_min/dt > 0 (min 7.5×10⁻⁶) [C].
   - **Consequence:** if the charged hidden potential is flat or negative at M_P, it stays at or below zero at every
     lower scale.
   - Positive contributions in the IR need fermion loops (−y⁴), and the GRUT architecture has no hidden fermions.

## 3. Gravity inputs (`code/as_flow.py`, part A)

- **Hypercharge at M_P.** g_Y(M_Pl) = 0.6154·√(3/5) = 0.4767 [Q 1307.3536 eq. 61a].
- **Calibrated f_g [D].** If hypercharge sits at its AS fixed point (the predictive reading of 1711.02949 / 1709.07252),
  then f_g = g_Y² b_Y/16π² = **9.83×10⁻³**. Any larger f_g is also allowed, since it only raises the bound.
- **Hidden-coupling fixed point [D].** e*_hidden = g_Y(M_P)·√(b_Y/b_e) = 0.4767 × √(123/276) = **0.318**; e* = 0.450 if
  f_g is doubled. The ratio e*/g_Y = 0.668 is fixed by GRUT's charges alone; it is the only clean "relation" the
  postulate produces.
- **A_λ at leading order** [Q formula of 2009.13543, C evaluation]:
  - 2.83 at the pure-gravity fixed point (g*, λ*) = (1.32, 0.13);
  - 0.0237 with SM matter content (7.63, −9.97);
  - 0.056 at de Brito et al.'s U(1)_D fixed point (2.78, −3.58).
  - The range spans two orders of magnitude. Everything below was scanned over A_λ ∈ [0.024, 10].
- **The nonminimal coupling ξ is truncation-dependent:** relevant at the pure-gravity fixed point, irrelevant at the
  Yukawa or matter fixed point [Q 2009.13543].

## 4. P4-S: e at its interacting fixed point (part B) [C]

With e = e* = 0.318, the U(1)_X-protected κ = 0 and the decoupled φ sector, the fixed points of the six remaining
quartics were searched with 700 multistarts per A_λ.

| A_λ | 0.0237 | 0.1 | 0.5 | 0.78 | 1 | 2.83 | 10 |
|---|---|---|---|---|---|---|---|
| real fixed points | 0 | 0 | 0 | 0 | 7 | 6 | 4 |
| bounded below | — | — | — | — | 0 | 0 | 0 |
| best V₄ min on unit sphere | — | — | — | — | −0.204 | −0.050 | −0.013 |

**Analytic reason [D].** The Σ equation alone is 20λ² + (16π²A − 432e²)λ + 7776e⁴ = 0, with q_Σ = 6. It has real roots
only for A_λ ≥ A_crit = (432 + √622080)e*²/16π² = **0.783**, and those roots are negative.

The fixed point with every direction irrelevant (all six) has λ_Σ* = −0.82 (A_λ = 1), −0.20 (2.83) and −0.05 (10).

**Result.** Under P4-S, GRUT's charge content has no admissible (real, bounded-below) hidden fixed point. The
mechanism is the one de Brito et al. used to rule out U(1)_D vector dark matter [Q 2312.02086: the g_D⁴ term "drives
the quartic coupling towards negative values in the IR"], made stronger here by q_Σ⁴ = 1296.

## 5. P4-M: the standard AS reading (part C) [C]

Method:
1. Transplanckian flow from ln(k/M_P) = 40/A_λ down to M_P, with the exact e(k) on the Gaussian branch.
2. Two different UV initial conditions: all zero, and an arbitrary nonzero set that includes κ = 0.01 and g² = 0.02.
   They agree at M_P to ≤ 2×10⁻¹¹ absolute, which confirms irrelevance.
3. Gravity off below M_P, then running down to M_Y.

| A_λ | e(M_P) | λ_X(M_Y) | λ_Y(M_Y) | λ_Σ(M_Y) | κ, g², λ_φ | V₄ min on unit sphere at M_Y |
|---|---|---|---|---|---|---|
| 0.0237 | 10⁻⁷ (heavy window) | −1.1×10⁻²⁸ | −9.0×10⁻²⁷ | −1.4×10⁻²⁵ | 0 | −3.4×10⁻²⁶ |
| 0.0237 | 0.01 | −1.1×10⁻⁸ | −9.0×10⁻⁷ | −1.4×10⁻⁵ | 0 | −3.4×10⁻⁶ |
| 0.0237 | 0.1 | −1.2×10⁻⁴ | −1.1×10⁻² | −0.44 | 0 | −0.11 |
| 0.0237, 0.1 | 0.286 (0.9e*) | no real critical surface: the quartics run away at ln(k/M_P) = 156 and 122 | | | | |
| 1, 2.83 | 10⁻⁷ | −5×10⁻²⁹ | −4×10⁻²⁷ | −7×10⁻²⁶ | 0 | −1.6×10⁻²⁶ |
| 1, 2.83 | 0.286 | −1.9×10⁻³ | −0.19 | λ_Σ reaches −50 at 1–3×10¹⁶ GeV | 0 | −12.5 |

What P4-M forces [D + C]:
- **φ sector:** g² ≡ λ_φ ≡ l_φY ≡ l_φΣ ≡ 0 at every scale.
- **Charge transfer:** κ ≡ 0.
- **Charged hidden quartics:** below zero at every scale below M_P, with |λ_i| ∝ q_i⁴e⁴, about 10⁻²⁸–10⁻²⁵ at the heavy
  window's e.
- **Free parameters:** only the relevant ones (masses, μ, e ≤ e*, and ξ where it is relevant).

**Benchmarks are unreachable.** The stated GRUT values sit 22–28 orders of magnitude off the critical surface, or have
the wrong sign. If the transplanckian growth λ ∝ (k/M_P)^{A_λ} were followed instead [E, linear], the Candidate 1
couplings run up to M_P (λ_X = 1.07×10⁻³, g² = 0.0108, λ_φ = 0.0356) would leave perturbation theory only at:

| A_λ | 0.0237 | 0.1 | 1 | 2.83 |
|---|---|---|---|---|
| log₁₀(k/M_P) | 107.5 | 25.5 | 2.5 | 0.9 |

So for a weakly irrelevant A_λ the restriction is exact as a statement about the UV limit, but it is not
operationally probeable [A].

**Consequences for the GRUT branches.**
- **Candidate 1:** impossible. m_X = gv needs g ≠ 0; the Friedberg–Lee–Sirlin-type binding needs λ_φ, g² ≠ 0. There
  is no relation λ_X/g² to test (g² = 0 is forced, not related). With λ_X ≈ 0⁻, X is effectively a free massive
  charged scalar: no Q-ball, and no self-interaction relation.
- **Heavy window:** impossible. κ = 0 removes the X⁶Σ* torque (α₆ ∝ μκ²). λ_Σ ≤ 0 removes the
  λ_Σ(|Σ|² − h²/2)² breaking at h = 5.33×10¹² GeV. With quartics this flat, the cubic μΣ*Y² leaves the tree potential
  unbounded along Σ*Y² [A].
- **Dark-matter chain:** the Higgs portals and ε_BX vanish (ε_BX by hidden C [Q REVIEW_REHEATING]; the λ_H-portals as
  in fact 1, up to the ξ_Hξ_X gravitational term of 2009.13543 eq. 41, which is ∝ e⁴ here). The REVIEW_REHEATING
  Higgs-portal obligation becomes moot. Any surviving hidden relic is only gravitationally coupled and gravitationally
  produced, which is standard Planckian dark matter [A].

## 6. P4-SM: the same postulate in GRUT's SM sector (part E) [Q formula, C evaluation; retrospective]

Under P4-M the hidden sector decouples from the Higgs, so GRUT inherits SW exactly. The inputs are 1307.3536 eq. (61e):

λ(M_Pl) = −0.0143 − 0.0066(M_t − 173.34) + 0.0018(α_s − 0.1184)/0.0007 + 0.0029(M_h − 125.15),

together with PDG 2026 (pdglive, fetched 27 October 2026):
- M_h = 125.13 ± 0.11 GeV;
- M_t = 172.60 ± 0.27 GeV (direct; ±0.3 GeV for the Monte Carlo–pole ambiguity added, as in 1307.3536);
- M_t = 172.1 ± 0.6 GeV (pole mass from cross sections);
- α_s = 0.1180 ± 0.0009 (as cited for PDG 2024/2025).

| m_t input | λ(M_P) | pull against 0 | m_h needed for λ(M_P) = 0 | m_t needed for λ(M_P) = 0 |
|---|---|---|---|---|
| direct, 172.60 | −0.0105 ± 0.0036 | **−2.95σ** | 128.8 ± 1.2 GeV | 171.0 GeV |
| cross-section pole, 172.1 | −0.0072 ± 0.0046 | **−1.56σ** | 127.6 ± 1.6 GeV | 171.0 GeV |

This is Shaposhnikov–Wetterich's and Eichhorn–Held's prediction (m_t ≈ 171 GeV [Q 1707.01107]). It is not GRUT's.
A positive fixed-point value such as λ_H* = 5.56×10⁻³ [Q 2312.02086, its FP (iii)] raises the pull to −4.5σ (direct).

It is the same near-criticality relation the unreviewed q-vacuum HANDOFF reached through P_q (absolute stability,
3.1σ). Two different closure strengthenings land on one SM relation. The relation is not GRUT-native, and a dark
portal in a symmetry-broken dark sector can move it [Q 2107.07949].

## 7. P4b: perturbative and bounded below to M_P (`code/p4b_bounds.py`, `code/release_check.py`) [C]

**Candidate 1.** It is matched at m_X = 92.4 MeV; e = 10⁻⁷ is switched on above h. At M_P its couplings are
λ_X = 1.07×10⁻³, g² = 0.0108, λ_φ = 0.0356, with V₄ bounded below. One-at-a-time maxima at m_X:

| Criterion | λ_X | λ_φ | g | m_s/m_X at g = 0.1 |
|---|---|---|---|---|
| no 4π below M_P | 0.174 | 0.193 | 0.570 | ≤ 6.2 |
| \|λ(M_P)\| ≤ 1 | 0.150 | 0.164 | 0.508 | ≤ 5.7 |

The benchmark (√6 = 2.45) is far inside. These are bounds, not relations.

**Heavy spot point.**
- **As stated, it fails.** λ_Y = λ_Σ = 1 at M_Y reaches 4π at **5.73×10¹⁵ GeV**, the same value as hw_gates gate 4.
  No value of the other quartic rescues it while one of them stays at 1.
- **Perturbativity ceiling.** For a common value λ_Y = λ_Σ the maximum is **0.566** (no 4π) or **0.372** (|λ(M_P)| ≤ 1).
- **Boundedness floor.** κ/κ_max = 0.701·λ_Y^{−1/4} < 1 needs λ_Y ≥ **0.241** at M_Y.
  - κ/κ_max is largest at M_Y, because λ_Y grows upward: 0.81 at λ = 0.566, 0.83 at 0.5, 0.95 at 0.3, and 1.05 at 0.2,
    where it fails.
  - V₄ min > 0 all the way to M_P for λ ≥ 0.3.
- **P4b band:** 0.24 ≤ λ_Y(M_Y) = λ_Σ(M_Y) ≤ 0.57, or ≤ 0.37 under the stricter criterion.

**Does the window survive in the band?** The spot release was rerun with the hw_gates C integrator (read-only reuse,
binary built in scratch): Y-relaxed preparation at δ = π/2, window a = 30–60. The control reproduces the spot point
exactly: ε = −0.14568, N_X = −2788.2, 1.000 Y_Q, R = 1.144.

| λ_Y = λ_Σ | ε | yield (Y_Q) | mediator R | R_L (table) | pass |
|---|---|---|---|---|---|
| 1.0 (control) | −0.14568 | 1.000 | 1.144 | 87.9 | ×77 |
| 0.566 | −0.14682 | 1.010 | 1.259 | 88.9 | ×71 |
| 0.45 | −0.14715 | 1.013 | 1.296 | 89.2 | ×69 |
| 0.37 | −0.14738 | 1.015 | 1.331 | 89.4 | ×67 |
| 0.30 | −0.14760 | 1.017 | 1.386 | 89.6 | ×65 |
| 0.25 | −0.14777 | 1.018 | 1.499 | 89.8 | ×60 |

The table R_L is 7–8% high (hw_critic). With it corrected, the passes become ×55–66. The vector debris (R_A = 48.5 at
e = 10⁻⁷) does not depend on λ_Y or λ_Σ. **P4b is satisfiable by a modest coupling change.** It is an obligation, a
correction to the stated spot couplings, and not a closure.

**Hidden gauge coupling.**
- Pure-gauge bounds: e(h) ≤ 0.62 for α(M_P) < 1, and ≤ 0.53 for e(M_P) < 1.
- Including the q_Σ⁴e⁴ drive of λ_Σ (at λ_Y = λ_Σ = 0.3): e(M_Y) ≤ **0.25** (no 4π) or 0.18 (≤ 1).
- These are compatible with the window's e ≳ 3×10⁻⁶ and Floquet-clean e ≤ 0.1.
- The AS bound e(M_P) ≤ e* = 0.318 is similar. Neither is testable at low energy, because the U(1) is Higgsed at
  5×10¹² GeV.

## 8. Postulate cards (Codex foundation gate, items 1–5)

**P4-M (standard AS closure in scale).**
1. **Postulate.** The full GRUT action lies on the UV critical surface of an AS gravity–matter fixed point, with
   - A_λ > 0 for every quartic;
   - f_g > 0;
   - global symmetries (shift symmetry, U(1)_X) respected;
   - the relevant couplings (masses, μ, G, Λ, Gaussian-branch e) supplied.
2. **Compatibility.** It is consistent with the SM sector. It is compatible with "no exterior" only as an added choice.
   In the hidden sector it admits only a flat or negative charged potential (§2 fact 3).
3. **Forbidden baseline [D].** Any matching with g² ≠ 0, λ_φ ≠ 0, κ ≠ 0 or λ_{X,Y,Σ} > 0. Explicitly, Candidate 1
   (g = 0.1, λ_φ = 0.03, λ_X = 10⁻³) and the heavy spot (κ = 3.46×10⁻³, λ_Y = λ_Σ = 1). This is a strict
   S_C ⊊ S_base, from facts 1–3 plus irrelevance.
4. **Discriminator.** In the hidden sector there is none beyond GRUT's own branches. In the SM sector it is λ_H(M_P) = 0,
   i.e. m_t = 171.0 GeV, against PDG 2026 (−2.95σ / −1.56σ, retrospective).
5. **Inputs and failure.**
   - Supplied: the sign of A_λ (truncation-dependent), f_g, the sharp matching, one loop.
   - It fails for GRUT unless hidden fermions with Yukawas are added; that is the escape of Eichhorn–Pauly
     2005.03661 / 2009.13543, where a dark fermion breaks the shift symmetry and gives a nonzero portal.
   - The escape is new field content that GRUT does not have. It would need its own fixed-point analysis.

**P4-S.**
- Items 1–2 are as for P4-M, plus e = e*.
- Item 3: P4-S forbids every GRUT hidden potential, because no admissible fixed point exists (§4).

**P4b.**
1. **Postulate.** Every coupling stays below 4π (or 1), and V₄ stays bounded below, up to M_P.
2. **Compatibility.** Standard.
3. **Forbidden baseline.** The stated heavy spot (pole at 5.7×10¹⁵ GeV).
4. **Discriminator.** None observable.
5. **Inputs and failure.** Fixed by λ_Y = λ_Σ ∈ [0.24, 0.57], and the window survives (§7).

## 9. Escape routes (each is a supplied input, so P4 then forces nothing GRUT-specific) [A]

1. **Hidden fermions with Yukawas to X, φ (or Y, Σ).** The −y⁴ terms make λ positive in the IR, and the broken shift
   symmetry regenerates portals at the fixed point [Q 2009.13543, 2107.07949]. It adds fields and couplings, and the
   portals become fixed-point functions of y, G*, Λ*.
2. **A_λ ≤ 0.** The quartics are relevant (SW's second case). P4 then reduces to bounds of P4b type.
3. **Relevant nonminimal couplings ξ_φ, ξ_X ≠ 0** (pure-gravity-like fixed point). The portal is then regenerated by
   G²ξ_φξ_X [Q 2009.13543 eq. 41]. ξ_X = −5/6 is GRUT's curvature-mass representative; whether ξ is relevant depends on
   the fixed point.
4. **New thresholds between M_P and the hidden scale.**

## 10. Caveats

- One loop, sharp Planck-scale switch, Euclidean FRG-derived signs, no gravitational threshold functions or mass
  dependence in f_λ [Q 2312.02086 eq. 24 has F_λ(m²)].
- Global symmetries are assumed respected by the gravitational fixed point. That is standard in AS truncations, but not
  in swampland-type arguments, which would break U(1)_X and the shift symmetry by Planck-suppressed effects [A].
- f_g calibration uses g_Y at M_Pl = 1.2×10¹⁹ GeV (Buttazzo's), not 2.4×10¹⁸ GeV; that is a 3% effect on e*.
- The q-vacuum HANDOFF is unreviewed. The overlap in §6 is noted, not relied on.
- The release recheck covers the mediator-debris gate only, at fixed δ, H and every other coupling.

## Literature (arXiv)

- Shaposhnikov–Wetterich 0912.0208 (abstract fetched).
- Eichhorn–Pauly 2009.13543 (ar5iv text fetched: irrelevance, shift-symmetry protection, eq. 41, θ formula).
- Eichhorn–Pauly 2005.03661 (Safety in darkness).
- Eichhorn–Pauly–Ray 2107.07949 (abstract).
- Reichert–Smirnov 1911.00012 (abstract via search).
- de Brito, Eichhorn, Frandsen, Rosenlyst, Thing, Vieira 2312.02086 (HTML fetched: f_g, f_λ forms, sharp switch,
  λ_S → negative).
- Eichhorn–Hamada–Lumma–Yamada 1712.00319 (abstract).
- Eichhorn–Versteegen 1709.07252 (abstract).
- Eichhorn–Held–Wetterich 1711.02949 (abstract).
- Eichhorn–Held 1707.01107 (abstract).
- Buttazzo et al. 1307.3536 (PDF text: eqs. 56–61, 64–66).
- Eichhorn review 2606.21522 (abstract only, no numbers used).
- PDG 2026 pdglive (m_t, m_h).

## Files

Code, all own:
- `code/beta_derive.py`: sympy β's and controls.
- `code/rg_common.py`
- `code/as_flow.py`: parts A–E.
- `code/p4b_bounds.py`
- `code/release_check.py`: reuses the hw_gates spot model read-only.
- `code/theorem_check.py`

Results: `results/beta_derive.json`, `as_flow.json`, `p4b_bounds.json`, `release_check.json`, `theorem_check.json`,
`summary.json`.

Reproduce with `python3 -B code/<script>.py results/<name>.json`. Runtimes: as_flow 18 s, p4b 1 min,
release_check 5 min, theorem_check 1 min.
