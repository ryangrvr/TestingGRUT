# Synthesis critic: what does closure force that standard physics does not? (track cf_critic)

27 October 2026. Claude, key `cf_critic`, parallel pass for issue #13 (brief: `../../BRIEF.md`). This consolidates
four gate tracks (`cf_space`, `cf_vacuum`, `cf_scale`, `cf_state`) and their four skeptic reviews. All of it is
same-session AI work, not outside replication. Git was used read-only and nothing was committed.

- Own code: `code/critic_checks.py` (about 80 s; reads sibling results read-only and reuses the gp_stability
  Nielsen–Olesen solver). Results: `results/critic_checks.json`, `results/summary.json`.
- Labels: [D] derived, [C] computed, [Q] quoted (source named), [E] estimate, [A] argued, [R] recalled.
- Every data comparison is **retrospective**. Nothing here validates or excludes a model of the real world, and no new
  theory is claimed. Every postulate below is an added choice, not a consequence of the word "closed".
- Where a skeptic corrected a track and I agree, the corrected statement is used and the track is not re-argued.

## 0. Answer to the owner's question

**(a) Closure as defined (FORMAL_CORE DEF 1.1).** It forces nothing beyond standard unitary quantum theory:
unitarity and its invariants, exact reduced dynamics, entropy relations, the athermality budget, protect-or-form for
records, no-signalling and monogamy (FORMAL_CORE §3). Of 14 "closure forces X" claims, 6 are refuted by legal closed
countermodels, 4 hold by definition and 3 are out of scope (§4 there). Codex's foundation gate (2908e91) states the
same: compactness, finite capacity, an observer dictionary and vacuum equilibration are extra assumptions. This pass
found nothing that changes that.

**(b) Closure plus one extra postulate.** Every strengthening tested either
1. reproduces a published hypothesis with no GRUT-specific number (STANDARD),
2. forbids something no measurement can reach (UNTESTABLE), or
3. removes GRUT's own dark-matter branches (INCONSISTENT WITH GRUT BRANCH).

None is DISTINCTIVE-TESTABLE. What the strengthenings forbid that standard physics allows:

| Strengthening | What standard physics allows and it forbids | Whose result it is |
|---|---|---|
| Compact S³ slices | Ω_K ≥ 0 | closed FLRW |
| Compact T³ slices | Ω_K ≠ 0; net Z6 residue mod 6 (unobservable) | topology searches; Krauss–Wilczek |
| Finite capacity, Banks–Fischler | eternal Λ ≤ 0 | Banks–Fischler; Bousso N-bound |
| Finite capacity, TCC | high-scale inflation (V^{1/4} ≥ 6.9×10⁸ GeV, r ≳ 10⁻³⁰) | Bedroya–Vafa |
| Finite four-volume | eternal expansion (forces Λ_∞ < 0) | Andrei–Ijjas–Steinhardt; Luu–Qiu–Tye |
| Vacuum equilibrium P_q (global minimum) | a metastable electroweak vacuum: M_t ≤ 171.00 GeV + … | absolute-stability bound; Klinkhamer–Volovik reading |
| Asymptotic safety in the SM sector | λ_H(M_P) far from its fixed-point value | Shaposhnikov–Wetterich |
| Purity plus a unitary whole | collapse models (and they are standard QM's prediction too) | standard QM |
| Fundamental clock bound | indefinitely coherent energy superpositions (effect 10²⁵⁺ below reach) | Gambini–Porto–Pullin |

The single **testable** restriction on measured parameters is the electroweak near-criticality line (P_q, and the
equality forms MPP and asymptotic safety). It is published Standard Model physics. Current masses sit on the side it
forbids at 2.5–3.5σ with the direct top mass, 1.6σ with the cross-section pole mass and about 0σ with the
electroweak-subfield α_s [C, §1.1]. Three of the closure-motivated strengthenings (Banks–Fischler, P_q and finite
four-volume) require **mutually incompatible** signs of the final vacuum (§1.2), which is itself evidence that
"closed" does not select them.

## 1. Contradictions resolved

### 1.1 One near-criticality line, three error budgets [C]

cf_vacuum (P3a, SMDR line), cf_scale (P4-SM, Buttazzo et al. 1307.3536 eq. 61e fit) and the q-vacuum HANDOFF quote
3.45σ, −2.95σ and 3.1σ. They are the same line.
- At PDG 2026 inputs the Buttazzo fit gives λ(M_Pl) = −0.01050, zero at M_t = 171.009 GeV. The SMDR line is
  170.997 GeV. They differ by **0.012 GeV**.
- The z values differ only in the error budget (`critic_checks.json` A):

| Budget (direct M_t = 172.60 ± 0.27) | σ_tot (GeV) | z |
|---|---|---|
| no MC-to-pole term, SMDR theory ±0.117 (cf_vacuum) | 0.464 | 3.45 |
| 0.3 GeV MC-to-pole term (cf_scale's choice) | 0.553 | 2.90 (2.96 with the fit's theory error) |
| 0.3 GeV MC-to-pole, Bednyakov's upward theory error 0.36 | 0.649 | 2.47 |
| 0.6 GeV MC-to-pole | 0.759 | 2.11 |
| cross-section pole mass 172.1 ± 0.6 | – | 1.56 |
| electroweak-subfield α_s = 0.1221 ± 0.0027 (PDG 2026 QCD review, via cf_vacuum skeptic) | – | −0.01 |

- The other matchings (cf_vacuum skeptic: B15 170.964, B13 170.902) move z up by 0.1–0.2.
- **Consolidated statement:** the line is M_t,crit ≈ 170.9–171.0 GeV at PDG 2026 inputs; direct-mass tension
  2.5–3.7σ depending on budget and matching; 1.6σ with the cross-section mass; none with the electroweak α_s.
  Retrospective.
- **They are not the same postulate.**
  - P3a is one-sided: M_t ≤ line, and it also needs c₆ ≥ 0 for the |H|⁶/M̄_Pl² operator (cf_vacuum skeptic).
  - P4-SM and MPP are equalities.
  - P4-SM's comparator is not fixed: λ_crit(M_P) = +0.017 at A_λ = 0.0237 with asymptotically free g₂ and g₃. That
    puts the direct-mass pull anywhere from −3σ to −8σ (cf_scale skeptic).

### 1.2 Final-vacuum sign: three strengthenings disagree [D]

| Strengthening | Required final vacuum |
|---|---|
| P2-BF (Banks–Fischler, astro-ph/0307459) | Λ_∞ > 0, eternal de Sitter |
| P_q (q-theory final state, HANDOFF; cf_vacuum P3a) | Λ_∞ = 0 (Minkowski end state) |
| Late-time TCC (1909.11063) | no eternal de Sitter: Λ_∞ ≤ 0 or a rolling field |
| P1-time, finite four-volume (cf_space skeptic) | Λ_∞ < 0, recollapse |

- Banks–Fischler is incompatible with each of the other three.
- P_q is incompatible with finite four-volume.
- P_q with TCC, and TCC with finite four-volume, are compatible.

**Consequence for Candidate 1.** The conditional H₀ exclusion of the 100 GeV, r = 1 preparation needs a total
vacuum ≥ −(2.28 meV)⁴ = −2.70×10⁻⁴⁷ GeV⁴ [C].
- Banks–Fischler supplies it, under the [A] monotone-roll link.
- P_q supplies it only through the calibrated residual, which is itself supplied.
- Finite four-volume reverses it: Luu–Qiu–Tye's (2506.24011) best fit is −6.0×10⁻⁴⁷ GeV⁴ (via cf_space skeptic).

So "Banks–Fischler supplies the vacuum-sign premise" (cf_space) is true but not distinctive. The premise flips
between sibling strengthenings, and present acceleration data already supply it retrospectively.

### 1.3 Verdict labels and track claims corrected by the skeptics (adopted)

| Item | Track said | Final |
|---|---|---|
| P1-S3, P2-BF, P3a-SM, P3d-MPP | MIXED (not one of the brief's five) | STANDARD |
| P1-T3 | UNTESTABLE | STANDARD: Ω_K ≠ 0 and matched circles test it, as they test flat ΛCDM |
| P3a + GRUT reference | INCONSISTENT WITH GRUT BRANCH | STANDARD (conditional). GRUT does not fix M_t, so the conflict is with data, not with the branch. Inherits the SM tension; no internal escape at the reference |
| P2-TCC | heavy window excluded; two-stage rescue dilutes by ≥ 10^54.3 | INCONSISTENT WITH GRUT BRANCH, **conditional** on inflation at H_i being what makes our Hubble volume homogeneous. Under that premise composite TCC forbids every multi-stage history (N_min of the last stage ≥ 19.38 > ln(M_P/H_i) = 16.72 [C reproduced]). The 10^54.3 is not a bound. Credit Tenkanen 1910.00521 |
| P2-BF e-fold cap [D] | forced by finite capacity | needs purity plus field theory valid over the whole patch. Not implied by the covariant bound (Banks–Fischler concede this; Lowe–Marolf hep-th/0402162) |
| P1-S3 data | PR4 removes the DESI deviation | PR4 reduces it: 10³Ω_K = 2.2 ± 1.1, 1.9σ (2511.20757). An exponential-quintessence fit (2504.04226) gives 3σ open, so the S³ failure condition must name the dark-energy model in advance |
| P4-M | g² = λ_φ = 0 and Candidate 1 impossible at every fixed point | only κ = 0 is unconditional (U(1)_X). That kills the heavy window. The Candidate 1 kill holds only where ξ is irrelevant (θ_ξ = +1.108 at the pure-gravity fixed point, recomputed from 2009.13543); there ξ_X = −5/6 is also forbidden |
| P4-S | 7/6/4 real fixed points | 8/12/16 by complete homotopy. Analytic theorem: no real fixed point with λ_Σ > 0 for any A_λ ≥ 0. Also e* gives e(h) = 0.284, above the window's e ≤ 0.1 |
| P5 | [D] energy closure forbids white CSL; X-ray bounds do not apply | withdrawn: out of scope (no conserved total energy with gravity, COND 5.4), and only 1.4×10⁻⁷ eV per nucleon in 13.8 Gyr. Its explicit forbidden baseline (white CSL at the GRW point) is already excluded by data: XENONnT 2506.05507 abstract, "the original values … are excluded experimentally for the first time" [Q, fetched this session]; 9.1σ per the cf_state skeptic's full-text read |
| P6b | Sr clock gap 2.8×10²⁵; Milburn τ < 9.7×10¹⁰ t_P | comparator invalid: the Bothwell observable is common-mode-blind. Still UNTESTABLE; GPP themselves: "too small to be observed in the lab" |
| P6c | forbids the externally switched model | forbids nothing (Θ(H − H_end) is an internal representative). STANDARD |
| P7a | spontaneous-CP window at c/b = 13; "exactly C-invariant" | the window is from the gauge point, not the heavy window (a/b = 2.56 there; it would need a tuned dimension-5 Λ ≈ 2–4×10¹⁵ GeV). C-invariance is a supplied truncation (χΣ*Y² and L_higher add phases). The correction belongs to rp_preparation §4, not REVIEW_REHEATING item 3 |
| P7c | 1.74 X̄ per net charge | definition-dependent (1.80 with an adiabatic split). Strengthens "not forced" |
| cf_vacuum Σ corner | UNTESTABLE, heavy-window results "not established" | UNTESTABLE stands. Sharpened in §4: the visible bath dominates the corner's Σ potential |

The skeptics' extra rows are added to the table: P1-time, total asymptotic freedom, classical scale invariance,
no-boundary, two-time boundary and CPT-symmetric whole.

## 2. Consolidated table (final verdicts)

Credit column: **R** = relabelled published physics, **N** = GRUT-native. **No row is GRUT-native.**

| Postulate | What it forces | Forbidden baseline | Discriminator and status (retrospective) | GRUT consequence | Final verdict | Credit |
|---|---|---|---|---|---|---|
| **DEF 1.1 closure** | standard unitary QM rules only | none observable | none | none | STANDARD | R (FORMAL_CORE §3–4) |
| P1-S3 compact k = +1 | Ω_K < 0; zero total gauge charges; SO(4) charges = 0 | flat/open FLRW | sign of Ω_K. DESI DR2+CMB 0.0023 ± 0.0011 (2.1σ against); PR4 1.9σ; Planck+lens+BAO 0.4σ; exp-quintessence 3σ against | none (homogeneous Gauss already imposed; flat premise of Candidate 1 holds for Ω_K > −0.56) | STANDARD | R closed FLRW, Moncrief |
| P1-T3 (Codex's branch) | Ω_K = 0; zero momentum; 216 Z6 holonomy sectors; flux quantised | curved slices | Ω_K ≠ 0; matched circles L > 6.04–6.14 c/H₀ null (1502.01593, 2211.02603) | no observable depends on L | STANDARD | R topology searches |
| P1-Gauss compact Gauss law | (N_X + 3N_Y) ≡ 0 mod 6 only (Q = 0 already holds on R³ for homogeneous fields and for screened localized objects) | net Z6 residue | none | runs already hold Q = 0 to ≤ 7×10⁻¹⁰ | STANDARD | R Gauss; Krauss–Wilczek |
| P1-QG compact + gravitational path integral | 1-dim Hilbert space without observer | large observer-free Hilbert space | none (e^{−S_obs}) | supports observer-relative dictionary; no number changes | UNTESTABLE | R HUZ 2501.02359; Marolf–Maxfield 2002.08950 |
| P2-BF finite Hilbert space e^{S_Λ} | Λ_∞ > 0 (with no crunch); e-fold cap only with extra premises | eternal Λ ≤ 0; P_q end state | sign of Λ: retrodiction | heavy window compatible (84.5 vs 58.4 e-folds); contradicts P_q | STANDARD | R Banks–Fischler, Bousso |
| P2-TCC | V^{1/4} < 6.9×10⁸ GeV, r < 2×10⁻³¹; acceleration ends within 2.4×10¹² yr | high-scale inflation; eternal dS | any r ≳ 10⁻⁸ falsifies TCC; BK18 r < 0.036 not constraining | heavy window and every inflationary X-charge branch (H = 1.3×10¹¹–4.7×10¹³) excluded, given the homogeneity premise; Candidate 1 compatible but needs a non-inflationary relic | INCONSISTENT WITH GRUT BRANCH (conditional) | R Bedroya–Vafa, BBLV, Tenkanen |
| P1-time finite four-volume | Λ_∞ < 0 | eternal expansion | ρ_DE evolution; fits exist (LQT best fit Ω_Λ = −1.61) | reverses Candidate 1's vacuum premise | STANDARD | R AIS 2201.07704, LQT 2506.24011 |
| P3a P_q global minimum + SM to M_Pl | absolute EW stability: M_t ≤ 171.00 + 0.494ΔM_h + 0.355Δα_s/0.0009 GeV (needs c₆ ≥ 0) | metastable SM (the measured central point) | z = 2.5–3.7 (direct), 1.6 (xsec), ≈0 (EW α_s); decisive: FCC-ee threshold mass (illustrative z ≈ 9–12) | none for the dark branches; no GRUT field escapes at the reference | STANDARD (testable, disfavoured) | R Degrassi, Buttazzo, KV, AHHKM |
| P3a + GRUT Σ corner | λ_Σ ≲ 3×10⁻⁵, λ_ΣY ≳ 0.05, portal in a band | reference λ_Σ = 1 with any portal | sinθ·M_Σ ∈ [30, 210] GeV; unmeasurable for M_Σ ≥ 10⁵ GeV | heavy-window results not established; the visible bath dominates Σ after reheating (§4) | UNTESTABLE | R EEGLS 1203.0237 |
| P3d MPP / P4-SM (equality forms) | M_t = line, or λ(M_P) = λ_crit(A_λ) | every point off the line | two-sided z; P4-SM pull −3σ to −8σ, A_λ-dependent | none | STANDARD | R Froggatt–Nielsen, Shaposhnikov–Wetterich |
| P4-M asymptotic safety, quartics irrelevant | κ = 0 (unconditional); g² = λ_φ = 0 and λ_Σ ≤ 0 where ξ is irrelevant | GRUT hidden potentials | none in the hidden sector | heavy window impossible; Candidate 1 impossible at fixed points where ξ is irrelevant | INCONSISTENT WITH GRUT BRANCH | R Eichhorn–Pauly, de Brito et al. |
| P4-S fully predictive AS | e* = 0.318; no real fixed point with λ_Σ > 0 | every stable Σ vev | none | no dark branch survives | INCONSISTENT WITH GRUT BRANCH | R Eichhorn et al. |
| TAF / classical scale invariance | not asymptotically free (b₀ = −4.17) / radiative Σ breaking (m_Σ = 0.37 GeV at e = 10⁻⁷) | heavy spot | none | heavy spot impossible | INCONSISTENT WITH GRUT BRANCH | R |
| P4b perturbative + bounded below to M_P | λ_Y = λ_Σ ∈ [0.24, 0.57]; e(M_Y) ≤ 0.25 | spot as stated (Landau pole about 10¹⁶ GeV) | none (consistency) | window survives the mediator gate (1.01–1.02 Y_Q); string and Floquet gates not rerun | STANDARD | R EFT practice |
| P5 purity + unitary whole | no fundamental collapse | collapse and classical-quantum gravity models (not standard QM) | XENONnT already excludes white CSL at GRW (data, not closure) | none | STANDARD | R |
| P6a / P6b / P6c relational time | nothing / GPP decoherence / nothing | none / long-coherent energy superpositions / none | GPP effect 10²⁵–10⁴⁶ below reach | none | STANDARD / UNTESTABLE / STANDARD | R Page–Wootters, GPP |
| P7a–d DM history | Sakharov conditions in finite form; Q = 0; 1.74 not forced | none beyond Sakharov | none | coherent asymmetry needs a supplied C/CP-odd input | STANDARD | R Sakharov |
| No-boundary / two-time / CPT-symmetric whole (missed) | measure-dependent scale preference / not gated / m_ν1 = 0, Majorana ν, no long-wavelength GW, DM = 4.8×10⁸ GeV ν_R | – | CPT-symmetric is testable | CPT-symmetric DM ≠ X; bare no-boundary suppresses the heavy-window scale | UNTESTABLE / not scored / INCONSISTENT WITH GRUT BRANCH | R Hartle–Hawking; ABL; Boyle–Finn–Turok 1803.08928 |

## 3. Ranking by promise for distinctive, testable content (strict: forced vs chosen)

Nothing in this list is forced by closure; each is chosen. The ranking is by what the choice would buy GRUT.

1. **P_q / near-criticality family (P3a, MPP, P4-SM).**
   - The only restriction on parameters people will measure, with a decisive future test (FCC-ee threshold scan with
     α_s; HL-LHC with lattice α_s).
   - Not distinctive: it is the published stability line, and GRUT's fields add nothing measurable to it.
   - Retrospectively on the forbidden side.
   - Its strictness rests on two supplied clauses: the residual offset δμ may not absorb ΔV, and c₆ ≥ 0.
   - The only GRUT-specific handle is whether GRUT's heavy window survives in the Σ corner (§5).
2. **P2-TCC.** Sharp and falsifiable (any r ≳ 10⁻⁸), but a liability: it removes every inflationary GRUT
   dark-matter branch, conditional on the homogeneity premise. It gives GRUT content only by elimination.
3. **P1-S3.** Testable through the sign of Ω_K and already mildly disfavoured in some model classes. It changes no
   GRUT number.
4. **P5 and the CPT-symmetric whole.** Testable only against other modified theories (collapse) or with a dark-matter
   identity that is not GRUT's.
5. **Banks–Fischler, finite four-volume, P1-T3/Gauss/QG, P4b, P6, P7.** These are retrodictions, consistency
   conditions or unobservable statements.
6. **P4-M, P4-S, TAF, classical scale invariance.** They remove GRUT's dark sector outright.

## 4. Quick check: the P3a Σ corner after reheating [C, E] (new)

The corner is the only GRUT-internal way to satisfy P_q if the SM is metastable. It needs λ_Σ ≲ 3×10⁻⁵ and a Higgs
portal λ_HΣ = ±2√(δλ_Σ), with δ inside the window in `cf_vacuum/results/threshold.json`, which exists for both signs.
Hence **λ_HΣ/λ_Σ = 2√(δ/λ_Σ) ≈ 30 to 4×10⁹ ≫ 1**. hw_gates gate 3's "no strings" argument rested on three
premises:
- Σ is frozen at its vev (m_Σ = 56.6 H);
- the hidden temperature is ≤ 0.25 h;
- **Σ has no visible coupling.**

The corner removes the first and the third. With the visible reheating temperature T_RH = 3.08×10¹⁴ GeV = 58 h and a
thermal Higgs (⟨|H|²⟩ = T²/6):

- **Positive portal.**
  - The effective potential restores the hidden U(1) above T_c = h(9λ_Σ/δ)^{1/4}. Across the window
    **T_c = 1×10⁸ – 1.1×10¹² GeV**, which is 3×10⁻⁷ to 3.5×10⁻³ of T_RH.
  - The bath reaches T_c when radiation is only 10⁻²⁶–10⁻¹⁰ of the reheating density. That is essentially at once,
    during the epoch in which the charge transfer is fixed (90% by τ ≲ 3; REVIEW_REHEATING).
  - The condensate tracks: m_T/H = 4×10³ – 3×10⁴ at T_c.
  - Whether Z6 strings then form by the Kibble mechanism is **undecided**, and the lean is no:
    - Σ quanta thermalise only below T_c (Γ/H = 1 at ≤ 0.16 T_c [E]);
    - the condensate is coherent and carries hidden charge;
    - vacuum fluctuations at T_c are ≤ 0.42 of the redshifted coherent amplitude [E].
  - If a network formed, Gμ = 7.6×10⁻¹⁴ – 5.7×10⁻¹² [C]. ε(β) comes from gate 3's large-β law, and from the
    gp_stability solver for type I (ε(1) = 1.000 BPS check, ε(10⁻⁸) = 0.128). That is 19–1400× above the LVK O3
    model-B edge (2101.12248) and 0.001–0.07 of the NANOGrav 15-yr limits (2306.16219), and LISA-visible
    [retrospective/illustrative].
- **Negative portal.**
  - There is no restoration. The thermal minimum moves to |Σ|_T = (δ/λ_Σ)^{1/4}T/√6, above the vev for
    T > T_* ≈ 1×10⁸ – 1.1×10¹² GeV.
  - At T_RH that is 3×10² – 3×10⁶ vev. The condensate follows only for M_Σ ≳ 10⁷ GeV (m_T/H at T_RH = 4×10⁻³ … 75).
  - The Higgs thermal mass² falls to (c_H − δ/3)T², which is about 0 near δ_max for M_Σ ≥ 2×10¹⁰ GeV. That is a
    high-T Higgs-stability question [E], with c_H ≈ 0.11 estimated.

**Reading [A].** In either sign the corner's Σ is controlled by the visible bath through the release, not frozen at
its vev. cf_vacuum's "heavy-window results not established in the corner" is therefore an expectation of change, not
only an absence of proof. The corner is internally consistent as a vacuum statement. As a dark-matter branch it is
unbuilt.

## 5. The single most valuable next calculation

**Re-solve the heavy-window release in the P_q Σ corner with the visible bath switched on.**

**Inputs:**
- λ_Σ ∈ {1.8×10⁻¹⁰, 1.8×10⁻⁶, 2.8×10⁻⁵};
- λ_ΣY = 0.06;
- λ_HΣ = ±2√(δλ_Σ) at δ_min and δ_max (threshold.json);
- κ, μ, M_Y, h and e unchanged;
- the REVIEW_REHEATING quiet inflaton corner (m_I = 3–5 H_i, Γ_I = 0.2–0.3 m_I);
- the Higgs-bath thermal mass λ_HΣ(T²/6)|Σ|² and a dynamical Σ radial mode.

**Outputs:**
- Y/Y_Q and the mediator debris R;
- the |Σ| and phase history;
- for the positive portal, whether the phase is randomised at T_c (strings);
- for the negative portal, how far |Σ| is pushed and whether the thermal Higgs mass stays positive.

**Why this one.** It is the only calculation found that can attach a GRUT-specific branch to a future measurement.

**Failure condition.** The corner is closed if, for both portal signs:
- no corner point reaches Y_Q for any Γ_I/m_I in the perturbative range with m_I ≤ 5 H_i; or
- every point that does reach it fails the debris gate, or (positive portal) forms strings with Gμ above the LVK
  model-B/C edge.

Then P_q + GRUT heavy window ⇒ absolute SM stability. That combination would be falsified by a ≥ 5σ short-distance
top-mass determination above the frozen line, with α_s and c₆ ≥ 0 stated in advance. If a point survives, the
combination has no measurable consequence beyond the SM line, and P_q gives GRUT nothing distinctive.

**Runner-up (cheap bookkeeping, to do before any FCC-ee data).** Refreeze the line in m_t(m_t) or y_t(M_t). Use the
asymmetric theory error (+0.17/−0.36) and state c₆ ≥ 0 and the cutoff convention.

## 6. Tensions between closure postulates and the GRUT dark-matter branch

1. **Capacity versus the heavy-window scale.**
   - TCC allows ln(M_P/H_i) = 16.72 e-folds at H_i = 1.33×10¹¹ GeV against 58.4 needed [C].
   - Composite TCC forbids every multi-stage history descending from one H_i patch.
   - It excludes all inflationary X-charge branches, conditional on the homogeneity premise.
   - Banks–Fischler is compatible (84.5 e-folds), but it contradicts the P_q end state used by the vacuum branch.
2. **Vacuum stability versus portal obligations.**
   - If metastability holds, P_q needs Higgs-coupled physics below Λ_λ ≈ 5×10¹⁰ GeV.
   - Through φ it would need b ≥ 0.117. That is 23× the invisible-width bound and **≈10⁹× the detector-level portal
     b ≈ 1–1.5×10⁻¹⁰**, at which the dark-matter chain already needs ≥ 97–99% erasure of pairs injected at onset
     (REVIEW_REHEATING) [C ratio].
   - The only internal route is the Σ corner. It needs a Σ–Higgs portal that REVIEW_REHEATING's ledger never
     considered, re-introduces a Higgs-mass tuning of 10⁻¹¹–10⁻¹⁹, and puts Σ under visible-bath control during the
     release (§4).
3. **Vacuum stability versus the inflation scale (standard physics, independent of P_q).**
   - H_i ≈ 0.5 h_max at central M_t. The heavy window needs 0.05–0.07 ≲ ξ_H ≲ O(1).
   - Under P4-M at the SM-content fixed point ξ is irrelevant (driven to 0). That removes this remedy and also forbids
     ξ_X = −5/6.
   - Under P3a with a stable SM the obligation disappears.
4. **Scale closure versus the hidden sector.** κ = 0 is protected by U(1)_X under P4-M, which kills the X⁶Σ* charge
   transfer. P4-S, TAF and classical scale invariance kill the heavy spot.
5. **P4b versus the reference couplings.**
   - λ_Y = λ_Σ = 1 hits a Landau pole at about 1.1–1.4×10¹⁶ GeV, so "GRUT valid to M_Pl" (a P_q premise) fails at the
     reference.
   - In the P4b band the hidden P_q bound passes with margin 64 (λ_Σλ_Y ≥ 0.058 against 9.05×10⁻⁴) [C].
   - There M_Σ = 3.7–5.7×10¹² GeV = 75–115 Λ_λ, so there is no Higgs escape [C].
6. **Combined strengthenings.**
   - P_q ∧ P4-M: portals → 0 and λ_Σ ≤ 0, so no escape. The SM must be stable, and λ(M_P) = λ_crit.
   - Banks–Fischler ∧ P_q: contradictory.
   - P5 ∧ P6b: reversibility needs the external time that P6 removes (cf_state skeptic).
7. **Preparation.** No closure-type condition selects the sign of the X asymmetry. It needs a supplied C/CP-odd input,
   or a tuned dimension-5 coefficient (Λ ≈ 2–4×10¹⁵ GeV) for a spontaneous canted minimum at the heavy point,
   untested dynamically.

## 7. Limits

- My checks reuse sibling numbers (threshold windows, T_RH, ε(β) law, string limits) as [Q] from the files named.
- The thermal analysis in §4 uses:
  - leading-order thermal masses and a high-scale c_H estimate;
  - a 2→2 rate estimate for Σ thermalisation;
  - an adiabatic-redshift estimate of the coherent amplitude.

  It is not a release simulation.
- String limits are loop-model dependent and retrospective.
- The XENONnT statement is quoted from the abstract fetched this session. The 9.1σ figure is the cf_state skeptic's
  full-text read.
- Skeptic corrections were adopted where I found them sound. I did not independently recompute the skeptics' 2-loop
  SM runs, homotopy counts or CSL rates.
