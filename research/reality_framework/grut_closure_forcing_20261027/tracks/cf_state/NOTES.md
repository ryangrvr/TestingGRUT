# Track cf_state: closure of the whole state and of time, and what closure forces in the dark-matter history

27 October 2026. Claude, closure-forcing parallel pass for issue #13 (brief: `../../BRIEF.md`). This is same-session AI
work, not an outside replication. Git was used read-only and nothing was committed. Codex's packets were read with
`git show` only. Downloaded papers are in the scratchpad (`.../scratchpad/closure/cf_state/dl/`), not in the repository.

**Labels.**
- [D] derived here (symbolic checks in `code/derivations.py`)
- [C] computed by `code/`, with outputs in `results/`
- [Q] quoted from a source read in this pass (read level in §7)
- [E] estimate
- [A] argued
- [R] recalled, not re-read in this pass

**Scope.** Nothing here is a validation, an exclusion of the real world or a new theory. Every comparison with an
existing measurement is **retrospective**. Every postulate is an **added choice**, not a consequence of the word
"closed": FORMAL_CORE §4 lists "a pure global state" as BY DEFINITION, and DEF 1.6 lists the clock as SUPPLIED.

**Files.**
- `code/common.py`: constants, and every quoted number with its source.
- `code/p5_collapse.py`: computes `results/p5_collapse.json`.
- `code/p6_clocks.py`: computes `results/p6_clocks.json`.
- `code/p7_dm.py`: computes `results/p7_dm.json`.
- `code/derivations.py`: computes `results/derivations.json`.
- `results/summary.json`: verdicts and key numbers.

---------------------------------------------------------------------------------------------------------------------

## 0. Answer

**None of the three families gives GRUT content that is both distinctive and testable.**

**Closure of the state (P5).** P5 forbids objective collapse and fundamental information loss, but only as against
those modified theories. Standard unitary quantum mechanics forbids exactly the same things.
- Purity by itself forbids nothing observable. A collapse trajectory also stays pure.
- The content that forbids collapse is narrower: a linear unitary law for the whole, with no classical noise field
  outside it.
- At the level of one-time statistics, even that content is not a strict restriction. Every collapse master equation
  equals the reduced dynamics of a closed unitary whole with a hidden bath. We checked this to machine precision
  [C].
- One genuinely closure-flavoured consequence was found: a finite closed whole with conserved energy cannot sustain
  white-noise CSL heating forever. Any localising noise must therefore be a finite-capacity internal bath, which
  places it in the dissipative or coloured class [D].
  - This consequence is a relabelling of the dissipative CSL model of Smirne–Bassi (arXiv:1408.6446).
  - It removes the X-ray bounds, since those apply to white noise.

**Closure of time (P6).**
- With ideal internal clocks, P6 forces nothing. This is Codex's RECORDS.md result, reproduced here by a different
  route [C].
- With a fundamental clock bound (the Gambini–Porto–Pullin and Ng–van Dam reading), P6 forbids indefinitely coherent
  energy superpositions. The predicted decoherence is far below anything measurable:
  - a decoherence exponent of 2.5×10⁻²⁶ for a Sr optical-clock superposition over 37 s, where the measured
    sensitivity is about 0.7, a gap of about 10²⁵;
  - about 10⁻⁴⁷ for the 25 kDa molecule;
  - about 10⁻²³ for a proposed 10⁻¹⁴ kg superposition.
- In the dark-matter history, P6 forces only that the Hubble-induced X mass be a function of internal fields. The
  internal representatives disagree after inflation, so P6 does not fix the release history [C].

**The dark-matter history (P7).**
- (a) The forced negative is true for expectation values, but it rests on hidden charge conjugation C, not on the
  toy's antiunitary K. In the toy, a K-invariant state transfers light charge at t ≠ 0 [C].
  - With C exact, as in the stated Lagrangian, the sign of the X asymmetry is a convention.
  - Branch-wise or domain-wise asymmetry is not excluded, so "closure cannot generate the asymmetry" is too strong as
    worded.
  - What is forced is narrower: either a C/CP-odd input, or a spontaneous breaking that is coherent across the
    observable universe and survives the release.
- (b) Neutrality is STANDARD.
- (c) The 1.74 X̄ per net charge is not forced. It is a function of the K-even datum E/E_min(Q) = 8.600 of the matched
  orbit. Closure-type selection rules give 0, about 0.1 or any value [C].
- (d) Nothing else is forced beyond Sakharov's conditions in finite form.

| Postulate | Verdict |
|---|---|
| P5 global purity + unitary whole | **STANDARD** (distinctive only against collapse and classical-gravity models) |
| P6a relational time, ideal clocks | **STANDARD** (no strict restriction) |
| P6b relational time + fundamental clock bound | **UNTESTABLE** (forbidden effect 10²³–10⁴⁷ below sensitivity) |
| P6c no externally timed switches in the dark-matter history | not forced beyond the existence of an internal representative (**STANDARD**) |
| P7a C/CP-odd input forced | **STANDARD** (Sakharov), with the scope corrections of §3.1 |
| P7b gauge neutrality | **STANDARD** |
| P7c antiparticle budget | not forced (**STANDARD** dynamics, supplied datum) |
| P7d other dark-matter observables | none found |
| Family verdict | **MIXED** (STANDARD + UNTESTABLE) |

---------------------------------------------------------------------------------------------------------------------

## 1. P5: global purity and a unitary whole (the purification reading)

### 1.1 The exact postulate

There are two separable parts.

- **P5-state.** The whole Ω is in a pure state.
- **P5-law.** The whole evolves by a linear unitary group U(t) = e^{−iHt}, where H is time independent and is defined
  on the whole. The whole includes gravity as a quantum subsystem. There is no classical stochastic field, classical
  sector or measurement device outside Ω.

The operational reading is the purification postulate of Chiribella–D'Ariano–Perinotti (Phys. Rev. A 84, 012311 (2011))
[Q via FORMAL_CORE §0].

### 1.2 What it forces (derivations)

1. **P5-state alone forbids nothing observable [D].**
   - Born probabilities are linear in ρ_Ω, so a mixed ρ_Ω reads as an ensemble of pure wholes (FORMAL_CORE §4, row 1).
   - Individual trajectories of a collapse model also stay pure. In our run, the nonlinear CSL/QSD trajectories are
     normalised pure states at every step [C].
   - Therefore "the whole is pure" does not distinguish collapse from no collapse. The forbidding content is P5-law.
2. **One-time statistics are not strictly restricted [C].** We took a two-site qubit with tunnelling Δ = 1 and
   localisation L = σ_z at γ = 0.5, run to t = 3. Four dynamics give the same reduced state:

   | Dynamics | Max \|Δρ\| from the Lindblad solution | Reference |
   |---|---:|---|
   | (2) nonlinear QSD collapse average, 20,000 trajectories | 1.8×10⁻³ | 1/√N = 7×10⁻³ |
   | (3) linear random-unitary average (classical noise, no collapse) | 3.3×10⁻³ | same |
   | (4) closed unitary whole: system + 12 ancilla qubits | 8.9×10⁻¹⁶ | against the same Trotterised Lindblad map; 6×10⁻³ against the exact solution is step error |

   The general statement is Gisin's: no-signalling forces the averaged dynamics to be linear (Helv. Phys. Acta 62, 363
   (1989)) [R]. Two consequences follow.
   - What collapse experiments bound (heating, X-ray emission, loss of contrast) is a Lindblad coefficient.
   - Any such term is reproduced inside S_base = "QM with a supplied hidden bath".
   - P5 is therefore a strict restriction only once the field content Θ is fixed.
3. **What does differ [C]: individual-outcome ontology and reversibility.**
   - At H = 0, the QSD trajectories end with |⟨L⟩| > 0.99 in 97.4% of runs. They collapse "up" with frequency
     0.306 ± 0.007, against the Born weight 0.3.
   - Random-unitary trajectories do not localise.
   - Only the closed whole (4) can be undone: applying the global inverse returns the system to its pure initial
     state with fidelity 1 − 1×10⁻¹⁵.
   - P5-law makes every decoherence reversible in principle. Collapse makes it fundamental. No experiment can access
     a hidden bath to test the difference.
4. **Energy closure restricts the collapse-noise class [D].**
   - THM 2.1 conserves ⟨H⟩. If H is bounded below and the total energy is finite, no subsystem's energy can grow
     without bound.
   - White-noise mass-proportional CSL heats every nucleon at dE/dt = 3λħ²/(4 m₀ r_C²). At GRW parameters that is
     5.0×10⁻⁴⁴ W per nucleon, or 1.8×10⁸ W for the Earth [C]. It continues forever.
   - Under P5, any localising noise must be an internal bath that pays for the heating, with bounded capacity. It
     therefore saturates, as in the dissipative or coloured CSL class.
     - Credit: Smirne–Vacchini, dissipative GRW, PRA 90, 062135 (2014), arXiv:1408.6115 [Q abstract].
     - Credit: Smirne–Bassi, dissipative CSL, Sci. Rep. 5, 12518 (2015), arXiv:1408.6446 [Q title].
     - Coloured/dissipative bounds: arXiv:1601.03672 [Q title].
   - The keV X-ray bounds apply to white noise only, which arXiv:2202.01343 itself states [Q]. Under P5, they
     therefore constrain a model class that P5 already forbids.
   - This is the one closure-flavoured consequence found in this family. It is a relabelling of dissipative CSL.
5. **Gravity [A, R].**
   - With gravity a quantum subsystem of Ω, gravitationally induced entanglement (GIE) is allowed:
     - Bose et al., PRL 119, 240401 (2017), arXiv:1707.06050 [Q];
     - Marletto–Vedral, arXiv:1707.06036 [R].
   - Diósi–Penrose (DP) as fundamental collapse is forbidden.
   - Classical–quantum (CQ) gravity is forbidden: its classical metric is an exterior. CQ gravity has a
     decoherence–diffusion trade-off (Oppenheim et al., Nat. Commun. 14, 7910 (2023), arXiv:2203.01982; Layton et
     al., JHEP 2023, 163, arXiv:2307.02557 [Q abstracts]).
   - Fundamental black-hole information loss is forbidden (Hawking 1976 [R]). It is untestable.
   - Gambini–Porto–Pullin's relational information loss (§2) is compatible with P5. Their PRL says explicitly that the
     argument "is based on ordinary (unitary) quantum mechanics" [Q, hep-th/0406260]. The whole stays unitary; the
     clock-conditioned description does not.

### 1.3 Baseline forbidden and strictness

- **Against standard QM with supplied Θ:** S_C = S_base. P5 is DEF 1.1 plus purity, which is BY DEFINITION
  (FORMAL_CORE §4). Gate item 3 fails.
- **Against S_base' = QM ∪ the collapse family:** P5 forbids:
  - GRW/CSL with λ > 0 as fundamental dynamics;
  - DP as fundamental;
  - CQ gravity;
  - the Hawking $-matrix.

  The explicit excluded baseline is white CSL at (λ, r_C) = (10⁻¹⁶ s⁻¹, 10⁻⁷ m). Its averaged dynamics has unbounded
  energy (item 4).
- **Within GRUT's fixed field content** (SM + metric + φ, X, Y, Σ + hidden U(1)), no field produces a CSL-like rate
  for laboratory matter.
  - With ε_BX = 0, the hidden photon does not couple.
  - The φ mass is about √(2λ_φ) v ≈ 0.23 GeV for Candidate 1, so φ is Boltzmann-absent at laboratory temperatures [E].
  - So P5 + Θ_GRUT predicts λ_eff = 0. That is the same prediction as SM + GR.

### 1.4 Discriminator and comparator (retrospective)

Collapse-model numbers at r_C = 10⁻⁷ m; current best bounds where quoted:

| Quantity | GRW / theory value | Current bound | Source |
|---|---|---|---|
| CSL λ, white noise, mass-proportional, coherent nuclei | GRW 10⁻¹⁶ s⁻¹; Adler 4×10⁻⁸±² | λ/r_C² < 0.494 s⁻¹m⁻² ⇒ λ < 4.9×10⁻¹⁵ | Majorana, PRL 129, 080401 (2022), arXiv:2202.01343 [Q] |
| CSL λ, quasi-free electrons only | – | λ/r_C² < 17.4 ⇒ 1.7×10⁻¹³; first full exclusion of Bassi et al.'s suggested m-p value under that assumption | same [Q] |
| CSL λ, Gran Sasso X-ray | – | < 5.2×10⁻¹³ | Donadi et al., EPJC 81, 773 (2021), via arXiv:2203.04231 [Q] |
| CSL λ, non-X-ray (bulk heating / cantilever / LISA Pathfinder) | – | 3.3×10⁻¹¹ / 2.0×10⁻¹⁰ / 3.8×10⁻⁹ | Carlesso et al., Nat. Phys. 18, 243 (2022), arXiv:2203.04231 [Q] |
| CSL, atom interferometry (2026) | – | λ < 1.27×10⁻⁵ at r_C = 10⁻⁵ m | Zhu et al., arXiv:2605.26860 [Q abstract] |
| DP R₀ | Penrose ~5×10⁻¹² m for Ge | R₀ > 2.54×10⁻¹⁰ m (95% CL) | arXiv:2202.01343 and erratum PRL 130, 239902 (2023) [Q]; earlier 0.54×10⁻¹⁰ m (Donadi et al., Nat. Phys. 17, 74 (2021), arXiv:2111.13490 [Q]) |
| DP "effective collapse" | – | R₀ ≲ 10⁻⁴ m (window ≈ 4×10⁻¹⁰–10⁻⁴ m) | Figurato et al., NJP 26, 113004 (2024), arXiv:2406.18494 [Q via search summary] |

Our predictions [C], `results/p5_collapse.json`:

- **Fein et al. molecule.** 25 kDa, v = 301 m/s from λ_dB = 53 fm, flight 6.6 ms, separation 266 nm [R].
  - CSL contrast exponent: 3.4×10⁻¹⁰ (GRW); 0.138 (Adler central value); 1.7×10⁻⁸ at the Majorana bound.
  - Sensitivity: −ln 0.9 = 0.105 for visibility > 90% of expected. That corresponds to λ ≈ 3×10⁻⁸, comparable to
    Adler's central value [E: not the authors' analysis].
  - DP time: 10¹² s.
- **Proposed 10⁻¹⁴ kg sphere** (R = 1.03 μm, Δx = 250 μm, 2.5 s; Bose et al. parameters [Q]). It has not been built.
  - CSL decoherence times: 34 ns at GRW; 0.7 ns at the Majorana bound. The GRW point would be decided immediately.
  - DP: E_G = (6/5 R − 1/d) G m² = 7.8×10⁻³³ J, so **τ_DP = 13.6 ms**, against a 2.5 s protocol. This is independent of
    R₀ at the bound; the nuclear-granular part is 1×10⁻⁴⁰ J.
  - That experiment would discriminate P5 (and standard QM) against DP and GRW-CSL.
  - Status: no GIE experiment has been performed as of October 2026; the proposals remain "beyond current
    technological capability" [Q, seminar abstract via search, low read level].
  - A 2025–26 membrane–qubit experiment states that the predicted DP collapse time is comparable to its measured
    mechanical decoherence time (arXiv:2505.21481 [Q abstract]).

### 1.5 Credit; GRUT-native?

- Credit: standard unitary QM; Everett [R]; Page–Wootters [Q]; the collapse-model literature (Bassi et al., RMP 85, 471
  (2013), arXiv:1204.4325 [R]); Gisin [R]; dissipative CSL [Q].
- **Nothing is GRUT-native.** GRUT's FORMAL_CORE already adopts DEF 1.1 (a unitary closed whole).

### 1.6 Consequences for the GRUT branches

- **Candidate 1:** none.
- **Heavy window:** none.
- **Dark-matter chain:** none.
  - Classical-statistical evolution and decoherence are already unitary-QM tools.
  - P5 forbids a collapse-driven "formation".
  - P5 makes CSL-type heating constraints on early-universe or IGM heating irrelevant (they are bounds on λ; under
    P5, λ = 0).

### 1.7 Verdict and failure condition

**STANDARD.** P5 is distinctive against collapse and CQ-gravity models, but not against standard physics. Its single
closure-flavoured sharpening (no unbounded white noise) is a relabelling of dissipative CSL.

P5 fails if any of the following is confirmed:
- a reproducible anomalous Lindblad term with the CSL fingerprint (mass-proportional, r_C-dependent heating or
  spontaneous emission) that no internal bath in Θ accounts for, including the GRUT hidden sector;
- a white-noise spectrum persisting without saturation;
- DP-type collapse at τ ≈ ħ/E_G in a GIE-class experiment;
- absence of GIE where linear quantum gravity predicts it.

---------------------------------------------------------------------------------------------------------------------

## 2. P6: no external time (observables relative to internal clocks)

### 2.1 Exact postulates

- **P6a.** Every prediction is a conditional probability p(a | T = t) relative to a clock subsystem inside Ω (Page &
  Wootters, PRD 27, 2885 (1983) [Q]). The clocks are ideal.
- **P6b.** P6a, plus the bound that no physical clock does better than the Salecker–Wigner limit optimised against
  black-hole formation (Ng–van Dam, MPLA 9, 335 (1994), hep-th/9406110; Amelino-Camelia, gr-qc/9603014 [Q via GPP
  refs]):

  δT² ≥ (3/2) t_P^{4/3} T^{2/3}.

  This is the reading of Gambini–Porto–Pullin, PRL 93, 240401 (2004). **The arXiv number is hep-th/0406260, not the
  brief's gr-qc/0406260 [Q].**
- **P6c.** In the GRUT history, no control or switch is timed externally. Every time dependence is a function of
  internal fields.

### 2.2 What P6a forces: nothing beyond standard QM [C, D]

Page–Wootters history state with a non-interacting clock (ħ = 1). The conditional coherence is

    |ρ₀₁(T)|/|ρ₀₁(0)| = |∫dt w(t) e^{−iωt}| / ∫dt w(t),   w(t) = |ψ_C(x_T, t)|².

- **Ideal clock** (linear dispersion, no spreading).
  - The coherence is independent of T.
  - It tends to 1 as the packet narrows: 1 − 1.25×10⁻⁷ at width 10⁻³ for ω = 0.5, and 1 − 2×10⁻⁶ for ω = 2, at every
    T.
  - The conditional probabilities are those of standard QM. This is Codex's RECORDS.md result ("the same constraint
    permits every branch probability"), reached here by a different route.
- **Free-particle clock of mass M.**
  - The spreading makes coherence fall with T.
  - At the optimal width s₀² = T/2M [D], the timing variance is δT² = T/(Mv²): the Salecker–Wigner form [D].
  - The computed coherence matches exp(−ω²δT²/2) to 2×10⁻⁴–8×10⁻³ relative (the residual is the non-Gaussian
    t-dependence of the width).
  - An FFT-grid evolution of the packet agrees with the analytic integral to ≤ 8×10⁻⁷.
- **Gate result.** The decoherence of P6 is entirely conditional on clock imperfection. Ideal internal clocks forbid
  nothing (S_C = S_base), so P6a fails gate item 3.

### 2.3 What P6b forces, and how far it is from data [C]

P6b forces a universal Lindblad term in the energy basis:

    ∂ρ/∂T = −i[H,ρ] − σ(T)[H,[H,ρ]],   ln|ρ₁₂(T)/ρ₁₂(0)| = −(3/2) t_P^{4/3} T^{2/3} ω₁₂².

This is GPP eq. (6) [Q]; the integral is checked [D].

Variants:
- GPP's earlier sub-optimal clock: σ ≈ t_P, D = ω² t_P T (NJP 6, 45 (2004), gr-qc/0402118, eq. 55 [Q]).
- Milburn's intrinsic decoherence (PRA 44, 5401 (1991)) [R]: D = ω²τT/2, here at τ = t_P.

The standard GR + QM effects are included for comparison. They are not closure postulates.
- Pikovski et al., Nat. Phys. 11, 668 (2015), arXiv:1311.1095: τ_dec = √(2/N) ħc²/(k_B T g Δx) [Q].
- Castro-Ruiz–Giacomini–Brukner, PNAS 114, E2303 (2017), arXiv:1507.01955: t_d = 2ħc⁴x/(√N G ΔE²) [Q].

| System (T) | ω or ΔE | D, GPP optimal | D, σ = t_P | D, Milburn τ = t_P | Measured sensitivity (retrospective) | Gap (optimal clock) |
|---|---|---:|---:|---:|---|---:|
| ⁸⁷Sr clock superposition, 37 s (Bothwell et al., Nature 602, 420 (2022), arXiv:2109.12238) | 2.70×10¹⁵ s⁻¹ | 2.5×10⁻²⁶ | 1.5×10⁻¹¹ | 7.3×10⁻¹² | all extra decoherence: 1/52.8 s ⇒ D ≤ 0.70 [Q] | 2.8×10²⁵ |
| 25 kDa molecule, 6.6 ms; ΔE = mgd, separation taken vertical (worst case) [E] | 1.03×10⁶ s⁻¹ | 1.1×10⁻⁴⁷ | 3.8×10⁻³⁴ | 1.9×10⁻³⁴ | D ≤ 0.105 (visibility > 90%) [E] | 9×10⁴⁵ |
| Proposed 10⁻¹⁴ kg, 250 μm vertical, 2.5 s | 2.3×10¹⁷ s⁻¹ | 3.0×10⁻²³ | 7.3×10⁻⁹ | 3.6×10⁻⁹ | not built | – |
| ⁸⁷Rb, 0.54 m, ~2 s (Kovachy et al. 2015) [R] | 7.2×10⁹ s⁻¹ | 2.5×10⁻³⁸ | 5.7×10⁻²⁴ | 2.8×10⁻²⁴ | – | – |
| GRUT heavy-point release (ω ~ H_i, T ~ 1/H_i) | ω t_P = 1.1×10⁻⁸ | 3.6×10⁻¹¹ | 1.1×10⁻⁸ | – | – | – |

- **Clock comparison.** The Ng–van Dam fractional timing floor at T = 10⁴ s is 3.8×10⁻³²; the best measured fractional
  frequency uncertainty is 7.6×10⁻²¹ [Q, arXiv:2109.12238]. Clocks would need to improve about 2×10¹¹-fold.
- **Retrospective bound on a Milburn-type τ.** The Sr data give τ < 5.2×10⁻³³ s = 9.7×10¹⁰ t_P. That does not reach
  τ = t_P.
- **What observability would take** (D = 10⁻³ in 1 s):
  - GPP optimal clock: a coherent superposition split by 1.2 TeV. Equivalently, a 1.9 g mass held in a 10 μm vertical
    superposition for 1 s.
  - The σ = t_P variants: about 90–130 keV splittings for about 1 s.
- **Standard time dilation** (not a closure effect): Pikovski gives D = 4.8×10⁻¹⁹ for the molecule at T_int = 500 K
  [E]. For the proposed sphere with a vertical split, D = 3.6×10⁻⁵ at 1 K and 3.2 at 300 K. Castro-Ruiz t_d for two Sr
  atoms 1 mm apart: 2×10⁴⁴ s.
- **GRUT.** The realistic-clock decoherence during the heavy-window release is ≤ 10⁻⁸. It has no effect on the
  dark-matter chain.

### 2.4 P6c in the dark-matter history: internal clocks do not select the release [C]

Codex replaced the externally switched −10H²|X|² by ξ_X R|X|² with ξ_X = −5/6 (04c53f5, ACTION.md). P6c requires some
representative of this kind. It does not say which.

We solved an oscillating m²χ²/2 inflaton era, with R = 4V − 2K [D] in reduced Planck units. Each representative below
is normalised to −10H² in de Sitter:

| Representative | Matter era mean (range), in H² | Radiation era |
|---|---|---|
| ξR, ξ = −5/6 (Codex) | −2.54 (−10 to +5) | 0 |
| c ρ_tot/M², c = −10/3 | −10 (constant) | −10 |
| c V(χ)/M², c = −10/3 | −5.0 (−10 to 0) | 0 |
| c (∂χ)²/2M², c = −10/ε (ε = 0.01) | −1490 (−3000 to 0) | 0 |

- Codex's −2.5 and REVIEW_REHEATING's "between −10 and +5 H²" are reproduced.
- **No external time** removes the externally timed switch. It leaves a one-function family of release histories. The
  representative is a supplied coupling, as Codex says.

### 2.5 Credit, GRUT consequences, verdicts, failure conditions

- **Credit.** Page–Wootters; Gambini–Porto–Pullin; Salecker–Wigner; Ng–van Dam; Amelino-Camelia; Milburn; Pikovski et
  al.; Castro-Ruiz et al.; Codex RECORDS.md (04c53f5) for the finite exact model. Nothing is GRUT-native.
- **P6a: STANDARD.**
- **P6b: UNTESTABLE.**
  - It is a strict restriction: it forbids energy superpositions coherent beyond T_D = ħ³/(ΔE³ t_P²) × O(1).
  - It sits 25 to 46 orders below sensitivity in the systems asked about.
  - It would fail if a coherent energy superposition with (3/2) t_P^{4/3} T^{2/3} ω² ≳ 1 were observed without the
    predicted loss. That needs about TeV splittings for about a second.
- **P6c: not a selection.**

---------------------------------------------------------------------------------------------------------------------

## 3. P7: what closure forces in the GRUT dark-matter history

### 3.1 (a) The preparation theorem: assumptions checked

**Statement (REVIEW_REHEATING item 3; tracks/rp_preparation §4).** Suppose:
- S is a symmetry with S H S⁻¹ = H;
- every selection condition is S-even;
- the asymmetry observable A is S-odd.

Then the admitted set is S-closed, its A-values come in ± pairs, and a unique selection is S-invariant with ⟨A⟩ = 0.

We reproduced the toy with our own code. We worked in the neutral basis (dimension 58 at the (24,12,9) cutoff),
whereas rp_preparation used a Kronecker-product construction. The torque is ∓1.3576447. The exact formula
12·0.04·2√2 = 1.3576450 is checked with sympy [C, D].

Assumption audit:

1. **Which S carries the dark-matter conclusion [C, D].**
   - In the toy, S = K is antiunitary, of time-reversal type, and L = n_a + 3n_b is **K-even**. Only the instantaneous
     torque is K-odd.
   - For a K-invariant state at t = 0, the torque is odd in t and ⟨L⟩ is even in t. Neither vanishes at t ≠ 0.
   - Computed for the real state δ = 0: ⟨L⟩ goes from 11.189 to 11.646 at t = ±0.5, and the torque is +0.56 at t = 5.
   - For the energy-minimising δ = π state ("the torque-free phase"): ⟨L⟩ goes from 11.189 to 10.800, and the torque
     is −0.50 at t = 5.
   - So the toy theorem constrains only the selection instant.
   - The dark-matter asymmetry (net X number) is C-odd and T-even. The statement that carries it is the
     unitary-C version: if [C, U(t)] = 0 and Cψ₀ = ψ₀, then ⟨Q_X⟩(t) = 0 at every t [D]. This is gp_phase §2's
     real-data argument in its quantum form.
   - Conclusion: the forced negative is valid, but through C (or CP), which the three-mode toy does not contain.
2. **Is H S-invariant in GRUT?** We counted rephasing invariants (rank of the charge matrix) [C]:
   - 0 physical phases at dimension ≤ 4 with exact Z₂_φ, for {κ, μ₀};
   - 1 phase if Z₂_φ is broken, because μ₁φΣ*Y² makes arg(μ₁μ₀*) physical;
   - 8 or more at dimension 5 (gp_phase §3 [Q]).

   With exact Z₂_φ, ε_BX = 0 (REVIEW_REHEATING: forbidden to all orders by hidden C) and renormalisable couplings, the
   stated GRUT Lagrangian is **exactly C-invariant** [D].
3. **Are the conditions S-even?**
   - Under K: energy, neutrality, entropies, mutual information, purity and stationarity are even.
   - Under C: the same list is even, provided occupation conditions are particle–antiparticle symmetric. A condition
     such as "n_X = 6" is not C-even.
   - The past hypothesis, as a one-time entropy condition, is even.
4. **Unique selection versus branches [C].**
   - The record-correlated whole (ψ₊|0⟩ + ψ₋|1⟩)/√2 is invariant under K ⊗ X_R.
   - Its global hidden torque is exactly 0. The branch-conditional torques are ∓1.3576. The hidden entropy is 0.689
     nats.
   - The branches have identical K-even data: E = 14.0849, the same occupations, and the same mode entropies.
   - The naive coherent cat (ψ₊ + ψ₋) is K-invariant but is not a clean example. The overlap ⟨ψ₊|ψ₋⟩ = 0.087 changes
     its K-even data (E = 15.40).
   - So the theorem constrains the expectation in an S-symmetric whole, not the asymmetry in a decohered branch or a
     causal domain.
5. **With hidden C exact, the sign is a convention [A].** "X dark matter" and "X̄ dark matter" are the same world
   relabelled. The physical question is whether |η_X| ≠ 0 coherently over the observable universe.

   **Spontaneous route [D, C].** The gauge-invariant angles A = 3θ_X − θ_Y and B = 2θ_Y − θ_Σ, plus a real Σ*YX³
   term of frustrated sign, give the phase potential

       V = −a cos A − b cos B + c cos(A+B).

   - Its global minimum violates CP spontaneously, as a degenerate ± pair, iff 1/a, 1/b and 1/c satisfy all three
     triangle inequalities [D: stationarity gives a sin A = b sin B = c sin(A+B), the law of sines].
   - Validated on 597/597 random triples [C].
   - Every minimum has zero gradient, so zero torque.
   - At gp_phase's natural ratio c/b = 13 (inflationary point, Λ = M_red [Q]), the window is
     a/b ∈ (13/14, 13/12) = (0.929, 1.083) [C]. That is a ±8% coincidence of the κ and μ terms.
6. **Coherence over the observable universe [E].** If the sign is set at release, at about one Hubble length of
   H_i = 1.33×10¹¹ GeV with T_RH = 2–5×10¹⁴ GeV:
   - today's Hubble volume holds 3×10⁷⁵–5×10⁷⁶ independent patches;
   - the net asymmetry is about N^{−1/2} ≈ 10⁻³⁸ of the local value;
   - this agrees with gp_phase's random-sign release (mean ε_X = −0.002 ± 0.011 [Q]).

   A coherent sign must be fixed **before the observable universe left the horizon**, during inflation, and held
   through the release. That is the hw_gates gate-2 holder requirement, which a spontaneous ±δ₀ could meet in place of
   a supplied CP-odd parameter. Its survival through the chaotic release is not computed here.

**Corrected statement of the forced negative.** Under exact hidden C, no C-even closure condition selects a state with
⟨Q_X⟩ ≠ 0. This holds at all times, by C rather than the toy K.

A coherent asymmetry over the observable universe needs one of the following:
- an explicit C/CP-odd input:
  - arg(μ₁μ₀*), which requires broken Z₂_φ;
  - a phase at dimension 5 or higher;
  - a supplied initial phase;
- a spontaneous breaking that is coherent across the present horizon and survives the release.

Closure supplies none of these and selects none of them.

- **Verdict: STANDARD.** This is Sakharov's C/CP condition (JETP Lett. 5, 24 (1967) [R]) in finite form.
- REVIEW_REHEATING's "a holder with a supplied CP phase remains the only route" should read "...a supplied CP phase or
  a spontaneously chosen, horizon-coherent phase...".
- **Failure condition.** A GRUT preparation condition proved C-odd, or a C-violating term found in the stated
  Lagrangian; either would void the hypotheses.

### 3.2 (b) Gauge neutrality [D]

- Gauss's law ∇·E = e J⁰ with homogeneous E gives J⁰ = n_X + 3n_Y + 6n_Σ = 0 pointwise, on R³ or T³. The homogeneous
  GRUT runs therefore already impose it (cf_space NOTES, P1 Gauss items; gp_phase "Gauss residual < 3.3×10⁻¹¹" [Q]).
- The X asymmetry is balanced by opposite Y/Σ charge. In the Higgs phase, the Σ condensate absorbs charge in units of
  6, so X number is an accidental global charge.
- On a compact slice, the only residue is the unobservable mod-6 condition (cf_space).
- **Verdict: STANDARD**, adding nothing to the dark-matter observables.

### 3.3 (c) The conserved antiparticle budget (≈ 1.74 X̄ per net charge) [C, D]

We wrote our own exact line decomposition of the homogeneous quartic orbit H = |π|²/2 + (λ/4)|χ|⁴, with λ = 10⁻³, from
en_particles' τ = 20 mean state:
- X̄ per Q₀ = **1.73778835311**, against en_particles' 1.73778835304: 4×10⁻¹¹ relative;
- line charge and direct charge agree to 2×10⁻¹⁴;
- r_max/r_min = 8.0625;
- **E/E_min(Q) = 8.6001**, where E_min(Q) = (3/4)λ^{1/3}Q^{4/3} [D]. The homogeneous mode alone has E_min/E = 11.63%,
  close to REVIEW_EARLY's lattice-wide bound of 11.594% of the energy.

The quartic orbit is scale invariant, so X̄/Q depends only on E/E_min(Q):

| E/E_min | 1 | 1.24 | 1.98 | 3.34 | 5.50 | 8.67 | 13.1 | 27.1 |
|---|---|---|---|---|---|---|---|---|
| X̄/Q | 0 | 0.071 | 0.277 | 0.621 | 1.109 | 1.752 | 2.565 | 4.764 |

The harmonic check (e−1)²/4e = 1.546 at e = 8.06 [D] reproduces REVIEW_EARLY's 1.55.

The budget is a function of one K- and C-even datum, the energy per charge^{4/3} of the matched mode. What closure-type
rules select:
- **energy minimum at fixed Q:** the circular orbit, so 0;
- **stationary (time-independent) density:** circular again, so 0;
- **maximum entropy at fixed (Q, E):** the equilibrated value of about 0.10–0.12 (REVIEW_TIMING/REVIEW_EARLY [Q]);
- **purity:** any value, because |α⟩_X|β⟩_X̄ is pure for every α, β.

The value 1.74 matches none of these. It is set by the release history: the torque-to-radial ratio when the field
leaves the Hubble-mass phase.

Its conservation through τ 1460 is number-conserving 2↔2 dynamics (Micha–Tkachev [R]), not closure. THM 2.1 conserves
the global spectrum, not particle number, and a pure closed whole can equilibrate (COND 2.10).

**Verdict: not forced**, a supplied datum. The failure condition does not apply.

### 3.4 (d) Other closure-forced observables of the dark-matter objects

- **Sakharov in finite form [D, Q].**
  1. X number must fail to commute with H while the asymmetry is written, and commute afterwards (protect-or-form, THM
     2.15).
  2. A C/CP-odd input is needed (§3.1).
  3. Departure from stationarity is needed: Tr(ρ[H, L]) = 0 for [ρ, H] = 0 (rp_preparation).

  All three are standard and select no rate or number.
- **Athermality budget, THM 2.14 [E].**
  - An asymmetry η_X = n/s = 4.8×10⁻⁹ (all dark matter as m_X = 92 MeV quanta) costs a relative entropy of about
    8×10⁻¹⁶ per entropy unit.
  - The initial condensate's athermality is O(1) per entropy unit.
  - The budget is about 15 orders from binding.
- **P6c** (internal switch representatives, §2.4): no selection.
- **P5:** nothing for the objects.
- We found no closure-forced observable of the dark-matter objects.

---------------------------------------------------------------------------------------------------------------------

## 4. Foundation-gate table (Codex 2908e91, five items)

| Item | P5 | P6a | P6b | P7a | P7b | P7c |
|---|---|---|---|---|---|---|
| 1. Postulate and what it constrains | §1.1: unitary linear law of the whole, no exterior noise | §2.1: clock-conditional probabilities | §2.1: plus the Ng–van Dam clock bound | exact hidden C, C-even selection | Gauss law | none proposed; tested "purity of preparation" |
| 2. Compatibility | consistent; with gravity quantum, GIE allowed | consistent (Codex 04c53f5 model) | consistent; GPP argue unitarity of the whole is kept | consistent | consistent | – |
| 3. Forbidden baseline (strict?) | collapse/CQ/$-matrix; not strict against QM + hidden bath at one-time level; white-CSL heating strictly forbidden by energy closure | none (S_C = S_base) | coherent energy superpositions beyond T_D: strict | ⟨Q_X⟩ ≠ 0 for a unique C-even selection: strict but standard | net hidden charge: strict, already imposed | none |
| 4. Discriminator and comparator | λ, R₀, GIE / τ_DP = 13.6 ms (table §1.4) | – | Sr clock D = 2.5×10⁻²⁶ vs 0.70 | none observable (the sign is a convention) | none | 1.74 vs closure rules 0 / 0.1 / any |
| 5. Supplied inputs and failure | §1.7 | – | §2.5 | §3.1 | – | – |
| Verdict | STANDARD | STANDARD | UNTESTABLE | STANDARD | STANDARD | not forced |

---------------------------------------------------------------------------------------------------------------------

## 5. Corrections to records and to the brief

1. The brief's gr-qc/0406260 is **hep-th/0406260** (GPP, PRL 93, 240401 (2004)). gr-qc/0402118 is GPP, NJP 6, 45
   (2004).
2. REVIEW_REHEATING item 3 and rp_preparation §4 say the theorem "holds for either [C or T]". The finite K-theorem
   holds only at the selection instant, and the toy transfers light charge from K-invariant states at t ≠ 0 (§3.1.1).
   The dark-matter conclusion should cite the C version (gp_phase §2). This is a scope correction; no number changes.
3. "A holder with a supplied CP phase remains the only route" (REVIEW_REHEATING table) omits spontaneous, horizon-coherent
   CP breaking from real couplings. That route exists in a ±8% window at c/b = 13 and is otherwise subject to the same
   holder and release obligations (§3.1.5–6).
4. The brief's "Donadi et al. Nature Phys 2020" is Nature Phys. 17, 74 (2021), online 2020. Its R₀ bound (0.54×10⁻¹⁰ m)
   has been superseded by R₀ > 2.54×10⁻¹⁰ m (Majorana 2022 with the 2023 erratum).

---------------------------------------------------------------------------------------------------------------------

## 6. Limits

- **P5 Part A** is a two-level illustration of a general theorem (Gisin; Stinespring dilation). Multi-time statistics
  are not compared beyond the fresh-ancilla collision model.
- **Collapse-model numbers** use point/sphere CSL formulas at r_C = 10⁻⁷ m with white noise. Bounds are quoted, not
  reanalysed. "Fein visibility > 90% of expected" is used as a crude sensitivity [E].
- **P6 Part B**:
  - The gravitational ΔE assumes vertical separation, a worst case. Bose et al. and Fein et al. separate horizontally,
    in which case ΔE_grav ≈ 0 and the GPP effect is smaller still.
  - The Milburn τ = t_P reading is a choice.
- **P6c representatives** are a minimal set, not exhaustive.
- **P7a**:
  - The C statement is derived, not simulated with antiparticle modes.
  - The spontaneous-CP criterion is for the three-angle tree potential plus one real dimension-5 term.
  - Survival of a spontaneous sign through the chaotic release (gp_phase) is not computed.
- **P7c** uses the homogeneous orbit only. The fluctuations' share, which is conserved in total, is taken from
  REVIEW_EARLY.

## 7. Sources and read levels

- **Read in full text this pass (PDF via arXiv):**
  - GPP, hep-th/0406260 and gr-qc/0402118 (§§IV–VI);
  - Majorana, arXiv:2202.01343 (results section);
  - Donadi et al., arXiv:2111.13490 (bound section);
  - Carlesso et al., arXiv:2203.04231 (bound sections, Table 1);
  - Pikovski et al., arXiv:1311.1095 (eq. 3);
  - Castro-Ruiz et al., arXiv:1507.01955 (eq. 8);
  - Bothwell et al., arXiv:2109.12238 (coherence section);
  - Young et al., arXiv:2004.06095 (abstract and coherence);
  - Bose et al., arXiv:1707.06050 (parameters).
- **Abstracts only:** arXiv:2605.26860, 2505.21481, 2402.13057, 1907.06896, 2203.01982, 2307.02557, 1408.6115,
  1408.6446 (title/journal), 1601.03672 (title).
- **Abstracts or repository pages:** Fein et al. 2019 (edoc.unibas.ch/75500; not on arXiv; the grating period 266 nm is
  [R]).
- **Search summaries (lowest read level):** the erratum PRL 130, 239902; Figurato et al. arXiv:2406.18494; the
  GIE-status seminar abstract.
- **Recalled [R]:**
  - Milburn PRA 44, 5401 (1991);
  - GRW PRD 34, 470 (1986);
  - Ghirardi–Pearle–Rimini PRA 42, 78 (1990);
  - Diósi PRA 40, 1165 (1989);
  - Penrose GRG 28, 581 (1996);
  - Gisin (1989);
  - Hawking (1976);
  - Sakharov (1967);
  - Kovachy et al. Nature 528, 530 (2015);
  - Micha–Tkachev;
  - Marletto–Vedral.
- **Program records:** FORMAL_CORE.md, the vacuum SYNTHESIS.md, the q-vacuum HANDOFF.md, FOUNDATION_GATE.md
  (2908e91), THEORY.md, RECORDS.md, REPORT.md and ACTION.md (04c53f5), REVIEW_REHEATING.md, REVIEW_EARLY.md,
  tracks/rp_preparation, tracks/en_particles, tracks/gp_phase, and ../cf_space/NOTES.md.

## 8. Reproduction (about 2 minutes in total)

```
cd research/reality_framework/grut_closure_forcing_20261027/tracks/cf_state/code
python3 -B p5_collapse.py ../results/p5_collapse.json   # ~60 s (20,000 trajectories)
python3 -B p6_clocks.py   ../results/p6_clocks.json     # ~20 s (scipy quad IntegrationWarnings at tiny coherences are harmless)
python3 -B p7_dm.py       ../results/p7_dm.json         # ~30 s
python3 -B derivations.py ../results/derivations.json   # ~3 s
```

`results/summary.json` was assembled from the four outputs; the assembling snippet is not saved as a script.
