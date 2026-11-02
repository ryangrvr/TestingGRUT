# GRUT working theory, version 2: the corrected physical picture

2 November 2026; **updated 3 November 2026 with the decisive pass** (`tracks/dq_heat`, `dq_entropy`, `qm_micro`,
`qm_transfer`, `higgs_reheat`, each with a skeptic, and the combiner `tracks/dq_critic`; brief BRIEF_DECISIVE.md).
Where the update changes a statement, the section says so. Claude, for the owner (issue #13). This replaces the physical picture in
`research/reality_framework/grut_closed_whole_theory_20261019/WORKING_THEORY.md` and the dark-matter picture assembled
in Codex's M1 packets (d4747cb → f8186fa). It rests on:
- the reviews in `research/reality_framework/grut_formation_review_20261025/`: REVIEW_M1.md, REVIEW_TRANSPORT.md,
  REVIEW_EARLY.md, REVIEW_TIMING.md, REVIEW_KINETIC.md and earlier;
- the closure pass: `research/reality_framework/grut_closure_forcing_20261027/SYNTHESIS.md`;
- the repair pass in this directory: `tracks/fx_holder`, `fx_chain`, `fx_formation`, `fx_observables`, each with an
  adversarial SKEPTIC.md, and the consistency critic `tracks/wt_critic`.

Everything here is same-session AI work (Claude and Codex) on one supplied model. It is not outside replication.
Nothing here is observational evidence, a validation, an exclusion or a confirmed GRUT prediction. Labels: [D]
derived, [C] computed, [E] estimate, [A] argued, [R] recalled. All data comparisons are retrospective.

## 0. In one paragraph

GRUT currently has two layers.

**The foundation, a closed-whole framework,** is a consistent way of organising standard quantum physics. It has no
law of its own. Closure as defined forces nothing beyond unitary quantum theory, and every strengthening we tested is
either published physics, untestable, or fatal to GRUT's own branches.

**The matter sector, GRUT-M1 repaired ("M1r"),** is an explicit hidden dark sector added by choice. It contains:
- a complex scalar X (m_X ≈ 92 MeV) whose mass comes from a binding scalar φ;
- heavy mediators Y and Σ;
- a hidden U(1);
- a renormalizable coupling to the inflaton that sets the matter–antimatter sign during inflation.

Inflation prepares a charged X condensate. It fragments, relaxes chemically, and is heated by mediator and vector
decays. A single number then decides what the dark matter is: the hidden **entropy per charge, S/Q**, at the time the
hidden sector cools.
- **S/Q > 3.47.** Charge-stabilised objects (Q-balls) cannot survive. The dark matter is a cold, asymmetric gas of
  free X particles.
- **S/Q < 3.47.** The dark matter is macroscopic Q-matter (Q ~ 10¹²–10²⁰) with a few-percent warm vapour.

The decisive pass computed S/Q for the actual repaired release, with nothing borrowed. At the declared hidden gauge
coupling e = 10⁻⁶ it is 4.955 on the static end-of-inflation slice M1r specifies, so the dark matter is free X. The
critical coupling is e_crit = 2.06×10⁻⁶.

But the branch at e = 10⁻⁶ is **not robust**. It depends on the inflaton's history before the slice, which M1r does
not supply: every rolling-inflaton embedding tried gives S/Q = 3.1–3.4, which is Q-matter. The Q-matter branch, where
it holds, is distinctive in kind: the objects ripen to Q ~ 10²⁰. That is about 0.2–7 Planck masses, roughly 1 nm
across.

Either way the abundance is fitted (H_i is calibrated to it). Small-scale clustering is CDM-like in both branches, to
4–6 decades below current bounds. The conservative model's only distinctive observable is the Planck-mass objects of
the Q-matter branch, visible only to proposed gravitational direct detection.

The open inputs now sit mainly in the inflaton sector: its history decides the branch, and its decay channel to the
Standard Model is not realisable as written (section 7).

## 1. The foundation: what closure does and does not do

| Statement | Status |
|---|---|
| Closure (a closed whole, no exterior, every bath is the rest) forces unitarity, the reduced dynamics, the entropy inequalities, no-signalling, quasi-recurrence | standard quantum theory [D]; FORMAL_CORE.md |
| Read as purification, closure singles out quantum theory among probabilistic theories | published (Chiribella–D'Ariano–Perinotti) [R] |
| Closure fixes the state, the couplings, the vacuum density, the arrow's direction or a clock | **no**: 6 of 14 claims refuted by countermodels, 4 true by definition, 3 out of scope |
| Stronger closure postulates (compact space, finite capacity, TCC, q-theory equilibrium, asymptotic safety, purity, internal time) | each is published physics, untestable, or inconsistent with GRUT's branches (SYNTHESIS.md). TCC would exclude the inflationary dark-matter branch. The postulates disagree on the sign of the final vacuum energy |
| Poincaré recurrence plus typical observers (Boltzmann-brain argument) | prefers a de Sitter phase that ends, i.e. a metastable vacuum, the side the measured Higgs and top masses favour; it opposes the q-theory "absolute stability" reading. Measure dependent, weak, published (Dyson–Kleban–Susskind) [A] |
| The matter–antimatter sign | **cannot** come from closure: every closure condition is conjugation-even (preparation theorem, REVIEW_REHEATING.md). It must be a supplied CP-odd input |

**Consequence.** GRUT's foundation is a framework, not yet a theory with independent content. A genuinely new,
explicitly stated dynamical law would be needed for that, and none is on the table. Everything below is a model choice
made inside the framework.

## 2. The matter-sector model M1r

### Fields and potential

- **Light hidden sector.** A complex scalar X with hidden charge 1 and a real binding scalar φ, with:
  - U = λ_X|X|⁴ (λ_X = 10⁻³);
  - g²φ²|X|² (g = 0.1);
  - V_φ = λ_φ/4 (φ² − v²)² − εv(φ − v)³, with λ_φ = 0.03, v = 10 m_X and ε = 0.002.
- **Masses.** m_X = g v = 0.0924 GeV and m_φ = √6 m_X. The ε tilt leaves exactly one vacuum, φ = +v [D], so there
  are no domain walls.
- **Heavy hidden sector.** Complex Y (charge 3) and Σ (charge 6), with:
  - M_Y = 30 H_i, ⟨Σ⟩ = h/√2 with h = 40 H_i, μΣ*Y² (μ = 7.16 H_i) and κY*X³ (κ = 3.46×10⁻³);
  - **λ_Y = λ_Σ = 0.5**, inside the band that is perturbative to M_P and keeps the window. This is changed from
    M1's 1, which had a Landau pole at 5.7×10¹⁵ GeV.
- **Hidden U(1).** Gauge coupling e, Higgsed by Σ (m_A = 6eh). Declared value: e = 10⁻⁶.
- **Phase holder (repaired).** −(g_I e^{iα_I} φ_I Σ*Y² + h.c.), with φ_I the inflaton:
  - dimension 4; g_I ≈ 8×10⁻⁷; α_I = 2.192;
  - it sets μ_eff = 10i H_i at the end of inflation;
  - it switches off as the inflaton decays [C, fx_holder].
  - It replaces Codex's dimension-5 curvature operator, whose own scale was ≈ H_i, below every field amplitude.
- **Gravity and the inflaton.** Nonminimal ξ_X = −5/6; Higgs ξ_H = 0.1 (needed for Higgs survival during
  inflation); an oscillating inflaton with (m_I, Γ_I) = (3, 0.5) H_i; zero Standard Model portals.

### Inputs, with their status

| Input | Status |
|---|---|
| g, λ_X, λ_φ, v, ε | supplied. ε is sensitive: the one-vacuum window is ε ∈ (0.0017, 0.058), and the coherent-orbit transition depends on ε within ±5% |
| M_Y/H_i, h/H_i, μ/H_i, κ, λ_Y = λ_Σ = 0.5 | supplied (perturbative band) |
| e | supplied. **It selects the dark-matter branch** (section 4) |
| α_I and the sign of φ_I during inflation | supplied: this is the CP sign. α_I ≡ α_I + π only if the inflaton sector is reflection-even (a fermionic decay channel is needed; φ_I\|H\|² would break it) |
| Inflaton history (m_I, Γ_I), sequestering of all other inflaton–hidden portals (≲ 5×10⁻¹⁵), inflaton quartic ≥ 1.4×10⁻²³ | supplied (an approximate shift symmetry) |
| Renormalisation condition of the induced portals (λ_IΣ, A_IΣ, λ_IY, A_IY) | supplied. Vanishing at M_Y versus at M_P moves H_i by 17% |
| H_i | **fitted to the dark-matter abundance** (section 5) |
| ξ_X, ξ_H | supplied |
| Vacuum density ρ_vac = 2.52×10⁻⁴⁷ GeV⁴ | measured input, underived |
| Target D = 4.4×10⁻¹⁰ GeV per unit entropy | supplied (rounded; Planck 2018 gives 4.37×10⁻¹⁰ [R]) |

**Retired from M1:**
- Q = 2×10⁵ as "the" object. In M1r no object of that size survives in either branch.
- The supplied 10% excitation pulse and the 20 MeV clock behind the gas forecast.
- η = 1 and p = 0 as inputs.
- The dimension-5 holder.

## 3. The physical history, stage by stage (M1r at e = 10⁻⁶)

| # | Stage | Physics | Numbers | Source |
|---|---|---|---|---|
| 1 | Inflation | The inflaton coupling holds μ_eff = 10i H_i, so the X–Σ relative phase sits at −π/2 and the CP torque is maximal. The hidden fields are heavy spectators (m² ≥ 10 H_i² > 9H_i²/4), so they track the local clock and their isocurvature is suppressed as e^{−3N} | prepared state (90.80, 6.89, 28.58) H_i; the charge isocurvature from the five hidden modes is ~10⁻⁷⁰ after 50 e-folds | [C] fx_holder; m1_holder |
| 2 | Release and reheating | As the inflaton oscillates and decays, the phase is released and a net X charge is kicked; reheating completes within about a Hubble time | Y/Y_Q = 0.9907 at fixed H_i; ε (charge/energy) = 0.127; T_RH ≈ 2–5×10¹⁴ GeV | [C] fx_holder, REVIEW_REHEATING |
| 3 | Early breakup | The eccentric charged orbit is linearly unstable (quartic self-resonance; vacuum seed 8–25× the 0.03 rad control) and fragments | strong charge contrast at τ ≈ 60–70; φ pumped but held to ~0.3% by depletion and detuning | [C] REVIEW_TIMING, REVIEW_EARLY |
| 4 | Prethermal plateau | Number-conserving 2↔2 scattering moves X̄ from the coherent orbit into waves | about 1.75 X̄ per net charge, constant over the pilot's span | [C] REVIEW_EARLY, 5f5513b |
| 5 | Chemical relaxation | The hidden sector is scale-free before φ condenses, so rates are constant in conformal time, with about 10¹² conformal units available. The φ channels (mixed g²λ_φ 2↔4, φ self-interactions, XX̄ ↔ φφ) relax the plateau | 10–10⁴ e-folds available; the plateau is not protected | [C/E] pt_relax and its skeptic |
| 6 | Heat injection | Mediator (Y_L) and longitudinal-vector decays heat the hidden plasma, all well before φ onset | total R = 4.97 (band) to 5.70 (broad) in units of the X energy; all sources finish 10³–10⁸ times before onset | [C] fx_chain, kin_source |
| 7 | Relaxed plasma before onset | Circular condensate holding 92% of the charge, plus a Hartree X, X̄, φ gas | heat R_total = 4.006 (vectors 3.608, Y_L 0.379, finite-k Y ≤ 4×10⁻⁴, cascade 0.019); **S/Q = 4.955** on the static slice (band 4.92–4.99); frozen yield Y/Y_Q = 0.99102 after the X ↔ Σ exchange freezes (decays as a⁻²) | [C] dq_heat, dq_entropy, dq_critic |
| 8 | φ crossover | The tilt makes φ condensation a crossover: no metastable φ ≈ 0 state and no walls. A circular condensate makes no X̄ and no φ quanta in the transition (resonance capture needs radial action) | onset a ≈ 3.8×10¹² (T_vis ≈ 120 → 80 GeV) | [C] fx_chain; pt_m1tilt |
| 9 | Spinodal decomposition | The condensate density passes the spinodal (233 m_X³) at T_h ≈ 10 m_X and decomposes quasi-statically into lumps | first lumps Q ≈ 4×10¹¹–10¹³, Γ/H ≈ 10⁷ | [C/E] fx_formation, fx_chain |
| 10 | **Evaporation** | Objects carry charge but no entropy, so they can coexist with the cooling vapour only if S/Q < min_T(s_gas/q_gas) = **3.47** (bulk; 2.93 for Q = 2×10⁵). At S/Q ≈ 4.9 every lump evaporates | evaporated by T_h ≈ 0.8 m_X (T_vis ≈ 7.5 GeV); re-formation impossible (fugacity 0.09–0.11 fixed by adiabatic expansion; nucleation barrier ≥ 10⁵ T; no small bound states below Q_min = 5.7×10⁴) | [C] fx_formation and both skeptics |
| 11 | X̄ freeze-out | Annihilation in a pure gas | p ≈ 6×10⁻⁴–1.3×10⁻³ pairs per charge (d ln p/d ln(S/Q) ≈ 9.5, so any entropy production raises it) | [C/E] fx_chain, wt_critic |
| 12 | Today | A cold, asymmetric gas of free X (92 MeV complex scalar) with a 0.1% X̄ admixture. Both are stable, and every annihilation channel at rest is closed: no lighter hidden states, and m_A = 6eh ≈ 3×10⁷ GeV | free-streaming 0.3–0.5 pc (cold); σ/m ≈ 6×10⁻⁸ cm²/g; ΔN_eff ≈ 10⁻⁸ | [C] wt_critic, fx_observables |

## 4. The branch point: what decides what the dark matter is

**Decisive-pass result (supersedes the transferred estimates below).** S/Q(e) for the actual M1r release (portals
vanish at M_Y):

| e | 10⁻⁷ | **10⁻⁶** | 1.5×10⁻⁶ | 2×10⁻⁶ | 3×10⁻⁶ | 5×10⁻⁶ | 10⁻⁵ | → ∞ |
|---|---|---|---|---|---|---|---|---|
| S/Q | 22.7 | **4.955** | 4.015 | 3.516 | 2.993 | 2.552 | 2.204 | 1.84 |
| Branch | free X | **free X (×1.43)** | free X | free X (0.4–2.2% margin) | Q-matter | Q-matter | Q-matter | Q-matter |

- **Critical couplings:** e_crit = 2.06×10⁻⁶ (bulk), 3.19×10⁻⁶ (Q = 2×10⁵), 4.99×10⁻⁶ (binding objects). With the
  induced portals vanishing at M_P instead: S/Q(10⁻⁶) = 4.09 and e_crit = 1.46×10⁻⁶.
- **Vector heat scaling, settled:** R_vec ∝ 1/e (computed exponent −1.007). N_A ∝ e^{−1/2} and the deposit per
  vector ∝ e^{−1/2}; an old code comment claiming otherwise was wrong.
- **The branch at e = 10⁻⁶ is not robust:**
  - **Inflaton history before the slice.** The vector number is set when the relevant modes left the horizon, 3–4
    e-folds before the slice, through the Higgs norm F there. Free X at 10⁻⁶ needs a horizon-exit factor
    f_eff ≥ 0.47. Every rolling-inflaton embedding tried gives 0.35–0.44, hence S/Q(10⁻⁶) = 3.1–3.4, which is
    Q-matter, with e_crit = 0.75–0.93×10⁻⁶. Any history obeys S/Q(10⁻⁶) ≤ 5.14.
  - **Inflaton decay history.** The vector heat scales as K² (K = a²H in radiation). Over (m_I, Γ_I) = (2–4, 0.25–1)
    H_i, S/Q(10⁻⁶) = 4.3–8.0 [E].
  - **Consequence.** The statement below that "everything uncomputed can only raise S/Q" is **false**: the inflaton
    history can lower it.
- **Post-relaxation entropy production** is +0.1–0.2% of S/Q (crossover adiabatic; spinodal latent heat
  ΔS/Q ≈ +0.005; evaporation reversible; X̄ burn-out zero).

*Original section (transferred estimates), kept for the record:*

The decisive number is the hidden entropy per charge after relaxation.

- **Free-X branch, S/Q > 3.47.** This is where M1r sits at e = 10⁻⁶: S/Q = 4.68–5.13, a margin of ×1.35–1.48. It is
  also where Codex's original M1 sits (S/Q ≈ 5.9).
  - The dark matter is a cold asymmetric scalar gas.
  - The mass per charge is 1.001–1.003 m_X.
  - This is robust in direction: everything not yet computed (crossover entropy, fragmentation, finite-k heating) can
    only raise S/Q.
- **Q-matter branch, S/Q < 3.47.** This needs less vector heat: e ≳ 2×10⁻⁶ for bulk objects (≈ 2.9×10⁻⁶ for
  Q = 2×10⁵), since the vector heat falls with e. Alternatively it needs a light hidden neutral state to carry the
  latent heat, which lies outside M1.
  - The dark matter is macroscopic Q-matter. ⟨Q⟩ = 10¹²–10²⁰ is set by mergers and Ostwald ripening, not by
    equilibrium.
  - E/Q ≈ 0.776 m_X; mass per charge 0.78–0.80 m_X.
  - A latent-heat vapour of free X holds 3–11% of the charge. It is warm (v ≈ 0.1c at capture freeze-out).
- **Heat threshold.** Objects of any kind survive only for total heat R < 1.5 (binding objects), 2.0 (Q = 2×10⁵) or
  2.8 (bulk).

**What is not yet computed:** S/Q of the actual M1r release end to end. The vector abundance on the repaired
background, finite-k Y production and the X ↔ Σ charge exchange to freeze-out are currently transferred from the
spot point. This is the decisive calculation (section 8).

## 5. The abundance and the mass test

- **The abundance is fitted, not predicted.** At a fixed dimensionless history Y_Q ∝ H_i^{3/2} (exponent measured
  1.5000000002), so H_i absorbs the mass per charge. This was already true in M1.
- **The frozen test** 0.0667805(1−η) + 2.13356p ≤ 0.01 assumed Q = 2×10⁵ objects (η = 1).
  - **Free-X branch:** η ≈ 0, so it reads 0.068–0.070 and fails by ×7.
  - **Q-matter branch:** the mass per charge is 15–17% below the reference.
  - **Diagnosis:** both are calibration failures, not model failures. The reference object was the wrong relic.
- **Recalibrated on the computed relic:**
  - **Free-X branch:** H_i = 1.29×10¹¹ GeV, or 1.08×10¹¹ if the induced portals vanish at M_P. The mass side closes
    to 0.15%; the yield side keeps a ±0.5% window spread at fixed inputs.
  - **Q-matter branch:** H_i = 1.50–1.52×10¹¹ GeV, with a ±1.2% spread from the vapour fraction until sticking and
    ⟨Q⟩ are computed.
- **Supplied-input sensitivity:** H_i ∝ φ_I(slice)^{−2.23}; the inflaton history moves it by −16% or +21%; the portal
  renormalisation condition by −17%. The 1% program tolerance cannot be met until those inputs are fixed.

## 6. What an observer could see

**Decisive-pass corrections (supersede the table below where they conflict).**
- **The Q-matter relic, benchmark e = 3×10⁻⁶, S/Q = 2.99:**
  - X sticking on Q-matter is high: s̄ = 0.86–0.97 at T_h = 0.1–0.5 m_X, 0.47 at freeze-out.
  - Ripening sets ⟨Q⟩ ≈ 1.2×10²⁰ (0.4–1.5×10²⁰ across the whole branch).
  - The vapour holds 12–17.5% of the charge.
  - The mass per charge is M = 0.80–0.82 m_X. An energy-conservation correction rules out net capture after
    freeze-out, because the objects have no heat capacity and no light neutral exists.
- **Small-scale power is CDM-like in every e range.** A new kinetic Boltzmann code (validated against CLASS fits)
  gives 1 − P/P_ΛCDM ≈ 10⁻⁸–10⁻⁵ at k = 10–100/Mpc. That is 4–6 decades below current Lyman-α, satellite and lensing
  bounds, and 2 decades below a 20 keV-class future probe.
  - The "O(10%) deficit near 10–50/Mpc" scenario does not occur.
  - A testable warm-vapour signal would need ⟨Q⟩ ≈ 4×10⁷–3×10¹² at s ≈ 1, which ripening never produces. In practice
    the effective thermalising capture is ζ_eff ≈ 4–6×10⁻⁹, limited by phonon trapping.
- **What is distinctive instead: the objects themselves.** At Q ~ 10²⁰ each object weighs 0.2–7 Planck masses (≈ 0.7
  at the benchmark, ≈ 1.5×10⁻⁸ kg) and is about 1 nm across. They arrive at roughly 0.3 per m² per year at the local
  dark-matter density.
  - Only proposed gravitational direct detection reaches such objects (arXiv:1903.00492; Windchime, 2203.07242), and
    no scheduled experiment does.
  - This is the conservative model's one distinctive prediction, and it exists only in the Q-matter branch.
- **Free-X branch:** unchanged; CDM-like, with no non-gravitational signature.
- **Tensor ratio:**
  - free X, H_i = 1.290–1.292×10¹¹ GeV: r = 2.7×10⁻⁷;
  - Q-matter, H_i = 1.483–1.499×10¹¹: r = 3.6×10⁻⁷;
  - with portals vanishing at M_P: 1.9×10⁻⁷ / 2.5×10⁻⁷.
  - These are floors that assume H at CMB exit equals H_i.

*Original section, kept for the record:*

| Channel | Free-X branch (e = 10⁻⁶) | Q-matter branch (e ≳ 2×10⁻⁶) |
|---|---|---|
| Large-scale structure | CDM on every accessible scale (free streaming 0.3–0.5 pc) | CDM-like only if ζ = s·⟨Q⟩^{−1/3} ≲ 10⁻⁵ (sticking s). For ⟨Q⟩ ~ 10¹² and s ~ 1 the warm vapour gives an O(10%) deficit near k ≈ 10–50/Mpc, near current Lyman-α and satellite reach. The vapour is not tightly coupled (its diffusion length exceeds the acoustic scale), so a kinetic or two-fluid transfer calculation is needed [C/E] |
| Self-interaction | σ/m ≈ 6×10⁻⁸ cm²/g | 10⁻³–10⁻⁹ cm²/g; 2–6 decades below the SIDM range |
| Direct detection, conservative (zero portals) | none | none |
| Direct detection, with a φ²\|H\|² portal | a single 92 MeV X scatters with σ ~ 10⁻⁵⁸ cm² at the costed coupling: invisible | composite form-factor zero (Xe 27.6 keV at Q = 2×10⁵, ∝ Q^{−2/3}) needs ⟨Q⟩ ≲ 10⁶; a detectable b ≈ 2×10⁻¹⁰ adds freeze-in entropy ΔS/Q ~ 10³–10⁴ and destroys the objects. The trilinear Aφ\|H\|² (allowed by the tilt) is marginal |
| Indirect, ΔN_eff | none (all hidden states massive; ΔN_eff ≈ 10⁻⁸) | none |
| Inflation | r = 2H_i²/(π²A_s M_P²) ≈ 2.5×10⁻⁷ at H_i = 1.29×10¹¹ GeV | r ≈ 3.5×10⁻⁷ |

**The one link to a measurable quantity is the inflation scale.** Within the benchmark's fixed dimensionless ratios,
the dark-matter abundance fixes H_i ≈ 1.1–1.5×10¹¹ GeV, so r ≈ 2–4×10⁻⁷. That is far below planned CMB reach
(r ~ 10⁻³). A detection of primordial tensor modes at r ≳ 10⁻³ would need H_i ~ 10¹³ GeV and would overproduce the
dark matter about 10³-fold unless the supplied ratios change [D/E]. It is a consistency link, not a distinctive
signature.

**Honest summary.** The conservative model is secluded dark matter that looks like CDM. Its only potentially
distinctive observable is in the Q-matter branch: the warm-vapour suppression of small-scale power. Its size depends
on two uncomputed quantities, the sticking probability and ⟨Q⟩, and its shape is generic.

## 7. What the repairs fixed and what they did not

| Problem (where found) | Status in M1r |
|---|---|
| φ domain walls (wall audits) | fixed by the ε tilt: one vacuum, no static wall [D] |
| Holder EFT control: dimension-5 operator with scale ≈ H_i (REVIEW_M1) | **fixed**: renormalizable inflaton trilinear, technically natural, loop corrections ≤ 10⁻⁹ H_i², built-in switch-off [C]. Costs: inflaton sequestering, the CP sign now in α_I plus sign(φ_I), strong history and clock sensitivity |
| Landau pole at 5.7×10¹⁵ GeV (closure pass) | **fixed**: λ_Y = λ_Σ = 0.5; κ and tree bounds hold |
| "Slow matching returns 1.76 X̄" (5f5513b) | resolved: the plateau relaxes chemically before onset, and a relaxed circular condensate passes the crossover without making X̄ [C] |
| Coherent-orbit resonance capture under the tilt (pt_m1tilt) | made irrelevant: relaxation removes the radial action before onset [C] |
| The forced-equilibrium start of the kinetic network (REVIEW_KINETIC) | justified by relaxation; p now depends only on S/Q, and placement no longer matters [C] |
| Supplied object size Q = 2×10⁵ and η = 1 | replaced by a computed relic (free X, or Q-matter with ⟨Q⟩ from ripening) |
| The 1% mass test | recalibration closes the mass side; the yield side depends on supplied inflaton and portal inputs |
| CP sign | **not fixable by closure**; supplied (α_I and the sign of φ_I) |
| Vacuum density | **not fixed**; measured input |
| Distinctive prediction | **none yet** |
| Higgs survival during inflation | **computed** (higgs_reheat + skeptic). ξ_H = 0.1 is safe during inflation: ξ_min = 0.044–0.055, a margin of about ×2. With gravitational-only reheating it is safe only below a critical top mass of 172.79–172.89 GeV (PDG 172.57 ± 0.29, so +0.8 to +1.1σ). At M_t + 1σ it fails at both H_i once the field normalisation is fixed (D = +13, +27). Above the critical mass, survival needs a Higgs-coupled plasma within τ ≈ 8–10 H_i⁻¹ |
| **The inflaton's decay channel to the Standard Model** (new) | **not realisable as written.** The reflection-even inflaton needed for the CP-sign bookkeeping forbids φ_I\|H\|²; M1r has zero SM portals; and a gauge-singlet inflaton has no renormalizable Yukawa to chiral SM fermions, so the decay must go through singlets (e.g. right-handed neutrinos). The supplied fast decay Γ_I/m_I = 1/6 needs N_f y² ≈ 4 and is kinematically blocked at φ_I ≈ 0.8 M_P (y φ_I ≫ m_I) except near zero crossings, i.e. it is non-perturbative. The reheating history, and with it K, the vector heat and the Higgs plasma timing, must be constructed |
| **The inflaton's history before the slice** (new) | not supplied. It decides the branch at e = 10⁻⁶ and sets H_i and r. The engine's post-inflation oscillator (m_I = 3 H_i, amplitude 0.82 M_P) is phenomenological; a potential that inflates at H ≈ 1.3–1.5×10¹¹ GeV with r ~ 10⁻⁷ must be a plateau or hilltop shape, not this quadratic |
| TCC tension (closure pass) | unresolved by choice: TCC is not adopted |

## 8. The decisive calculations

**Status after the decisive pass.** Items 1–3 below have been computed:
- **S/Q end to end:** free X at e = 10⁻⁶ on the static slice; Q-matter above 2.06×10⁻⁶.
- **The Q-matter microphysics and observable:** CDM-like; Planck-mass objects.
- **X ↔ Σ freeze-out:** Y/Y_Q = 0.99102.
- **Higgs survival:** above.

The free-X failure condition ("vector heat below about 0.5 of the transfer") is not met on the static slice, but it
**is met** on every rolling-inflaton path tried.

**The open inputs that now decide everything, ranked:**
1. **The inflaton's potential and history before the slice.** It decides the branch at e = 10⁻⁶ and sets H_i and r.
2. **The inflaton's decay channel and history.** It sets K and hence the vector heat, the Higgs plasma timing, and
   whether reheating is perturbative at all.
3. **The induced-portal renormalisation condition** (M_Y vs M_P): 17% in H_i, and e_crit 2.06 vs 1.46×10⁻⁶.
4. **The top mass, together with the Higgs decay criterion:** survival without a Higgs-coupled plasma fails at
   M_t + 1σ.
5. **Dephasing of the X ↔ Σ exchange at fragmentation:** a lattice with Y and Σ; H_i +1.5%/−0.4%.
6. **Completeness of chemical relaxation before onset:** the one omitted process that could lower S/Q.
7. **Interior trapping in Q-matter** (perfect vs Beliaev): f_v 0.12 vs 0.175, M up to 0.818.
8. **The pair-fraction map p(S/Q).**

The next construction is therefore an **explicit inflaton sector**, consistent with H ≈ 1.3–1.5×10¹¹ GeV,
n_s ≈ 0.965 and A_s. It needs a realisable reheating channel to the Standard Model (for example through singlet
fermions that decay to SM states) that preserves the CP-sign bookkeeping and gives the Higgs a plasma mass in time.
That construction decides the branch, H_i and r together.

*Original section, kept for the record:*

1. **S/Q of the actual M1r release, end to end.** Compute the inflationary longitudinal-vector abundance on the
   λ = 0.5, trilinear-held background through the (3, 0.5) era, the finite-k Y_L production, and the relaxed Hartree
   S/Q under the adopted portal renormalisation condition.
   - **Failure condition for the free-X branch:** S/Q ≤ 3.47, i.e. vector heat below about 0.5 of the transferred
     estimate. The relic would then be Q-matter and H_i would rise about 17%.
2. **If the Q-matter branch is wanted:** compute the sticking probability and the saturated vapour density next to a
   finite object (a grand-canonical, thermostatted measurement or a quantum-kinetic vapour coupled to the classical
   object; the 48³ boxes cannot hold S/Q ≥ 3). Then compute ⟨Q⟩ from coarsening, and a kinetic or two-fluid transfer
   function for the vapour against Lyman-α and satellite bounds.
3. **Two smaller items:** the X ↔ Σ charge exchange to freeze-out (it sets N_X/s to better than 0.2%), and the Higgs
   field's evolution through reheating with ξ_H = 0.1.

## 9. Bottom line for the owner

**After the decisive pass (3 November).** The dark-matter sector is now worked out end to end at the level of
homogeneous, linear and kinetic calculations, every one adversarially checked. Its outcome is a clean either/or set
by one number, the entropy per charge:
- **A cold asymmetric gas of 92 MeV scalars.** Ordinary-looking dark matter.
- **Planck-mass, nanometre Q-matter nuggets.** Distinctive in kind, detectable only by proposed gravitational
  sensors.

Which one nature would pick in this model depends on the inflaton, which the model does not yet describe properly.
In particular, it has no realisable way to reheat the Standard Model, and that same gap makes the Higgs's survival
marginal. Small-scale clustering gives no signal in either branch.

So the honest state of the theory:
- the foundation adds no law;
- the matter sector is a consistent hidden-sector model whose remaining freedom has been pushed into one place, the
  inflaton sector;
- its one distinctive possible signature (Planck-mass objects) is real but far from testable.

*Original section, kept for the record:*

The foundation organises physics but does not yet add any. The matter sector is now a coherent, internally
consistent, renormalizable model with a definite history, at the level of the calculations done here. Its dark matter
at the declared coupling is a cold asymmetric scalar gas: an ordinary-looking dark-matter candidate, not the
charge-stabilised objects GRUT set out to describe. Those survive only if the hidden sector stays colder (e ≳ 2×10⁻⁶),
and then they are huge (Q ~ 10¹²–10²⁰), with a possible small-scale clustering signature that is not yet computed.

The model fits the dark-matter abundance; it does not predict it. It has no distinctive signal in its conservative
form. Getting one requires choosing the Q-matter branch and computing its vapour, or adding a portal at a cost the
calculations quantify.

Candidate 1 is unchanged. Closure supplies no new law. The vacuum density is underived.
