# Track `ledger`: what GRUT closure supplies to q-theory, and the exposed-data ledger

21 October 2026. This is one track of the construction pass in `../../BRIEF.md`. Part A maps GRUT's closure statements
onto Klinkhamer-Volovik (KV) q-theory. Part B lists the exposed data that any q-theory construction can be compared
with.

**Labels.**
- [D] derived here.
- [C] computed by `code/ledger_build.py`. The output is in `results/ledger_numbers.json`, with 8/8 checks passing,
  and in `results/ledger.json`, which has 15 entries. The script is deterministic and runs in under 1 s.
- [Q] quoted from a source at the read level given in §6.
- [E] estimate; [A] argued, not a theorem; [R] recalled, not checked against a source in this track.

**Scope of the claims.**
- No new theory, empirical validation or confirmatory test is claimed.
- Every observational number is already exposed, so every comparison is **retrospective**.
- Candidate 1 stays paused, with Codex's status unchanged: "Candidate 1 is not independently test-ready; universal
  equivalence to cold dark matter remains unproved." §3 is conditional throughout.

**Sibling records used.**
- `../equations/NOTES.md`, the blind derivation from the action. Cited as "equations §n".
- `../literature/NOTES.md`, the KV source extraction. Cited as "literature §n".
- I spot-checked their KV statements that are load-bearing here against the downloaded TeX (§6). Where both tracks
  reach the same statement, I cite them and do not re-derive it.

**Notation collision (read before §3).** Candidate 1's "ρ_Λ(μ₀) = 0" uses μ₀ = g_X v = 0.1 v, which is a
renormalization scale. KV's μ₀ is the equilibrium value of the 3-form integration constant, a chemical potential.
The two are unrelated. Below, μ means KV's μ, and the RG scale is written μ_RG.

---------------------------------------------------------------------------------------------------------------------

## 0. Result

**A1. "No exterior ⇒ the vacuum of the whole does not gravitate" cannot be derived from DEF 1.1–1.8 [D/A].**
- DEF 1.1 supplies exactly one of the eight conditions that Volovik lists for a zero vacuum energy: (vi), "does not
  interact with environment" (gr-qc/0604062 §8, opening paragraph) [Q].
- In the finite class, closure forces none of the other seven: equilibrium, stationarity, homogeneity, zero curvature,
  no excitations, T = 0, no boundaries (FORMAL_CORE §4). It also supplies no vacuum variable, no volume or pressure
  variable, and no gravity.
- "Zero external pressure" is true by DEF 1.1, but only vacuously: nothing outside exerts any. The conclusion that
  matters is zero **vacuum** pressure, and that needs the added structure.
- Countermodel: q-theory solutions with μ ≠ μ₀ (equations §5.2, §5.4 case C) are closed wholes, with no reservoir and
  no drive, whose gravitating vacuum energy is nonzero.
- §1.

**A2. What has to be added.**
- a 3-form vacuum variable with ε(q) and K(q) (algebra and dynamics slots);
- a selection of its conserved Gauss-law constant μ (preparation slot);
- the reading of "equilibrium of the whole" as the **final** (asymptotic, T → 0, post-transition) vacuum. This is a
  final boundary condition, a kind that FORMAL_CORE §7 lists as outside its table.

The other two readings of the GRUT wording fail:
- **Initial equilibrium** leaves −Σ Δ_i. For Candidate 1 plus the SM this is −2.62×10⁹ GeV⁴ = −1.0×10⁵⁶ ρ_Λ, an AdS
  crunch [C].
- **Zero pressure re-imposed after every transition** needs μ to change, which the 4-form dynamics forbids. That is an
  inconsistency, not an option (literature §3.3).

The candidate postulate (§1.4), with ε_∞ ≡ ε(q) plus every q-independent vacuum energy of the final matter vacuum, to
all orders:

    (P_q)   μ_Ω = μ₀[ε_∞],     with  ε_∞'(q₀) = μ₀,   ε_∞(q₀) − μ₀ q₀ = 0,   ε_∞''(q₀) > 0

**A3. Closure contributes one structural statement [D/A].** THM 2.15 (protect or form), applied to the conserved μ̂,
says:
- if μ is conserved, the dynamics can neither write nor erase it, so its value is preparation data;
- if dissipation lets μ relax, it is no longer protected.

This is KV's own concession in GRUT's words ("the fine-tuning moves to μ": 0907.4887 §3, 1601.00601 §1). **Nothing
is GRUT-native.** Every physical element is in KV or Volovik 2006, and GRUT supplies wording and slot bookkeeping
(§2).

**A4. Candidate 1, if q-theory were adopted (conditional; Candidate 1 is not unpaused).**
- The uncancelled −1.2265×10⁶ GeV⁴ (−4.9×10⁵² ρ_Λ) is absorbed by re-selecting μ: δμ₀ = c/q₀, and
  δq₀/q₀ ≈ δK/K ≈ −3.5×10⁻⁶⁸ at natural scales [C].
- P_q is the "stated physical total-vacuum matching rule" that the quantum report asked for. It matches the
  asymptotic total to 0, not to ρ_Λ.
- The H₀-exclusion premise ρ_vac,tot ≥ −(2.28 meV)⁴ then follows from P_q. It no longer has to be assumed.
- The hidden Z2 transition releases ΔV = λ_φv⁴/4 = 2.5×10⁹ GeV⁴, which is 21 times the SM electroweak drop. It
  gravitates as in GR until it is released, and at T_c it is only 1.2×10⁻⁵ of the radiation density.
- That transition leaves no constant remnant. Its transient second-order term is 1.4×10⁻⁸ ρ_Λ [C].
- §3.

**B. Ledger (§4, `results/ledger.json`).** The 15 entries are all exposed. Read levels are abstract or section, plus
one recalled number.
- Without screening, the q-mode's Yukawa strength 1/3 requires m_q > 2.9×10⁻¹⁵ eV (Cassini) and m_q ≳ 4–5 meV
  (Eöt-Wash) [C/E].
- In that allowed regime the construction is GR plus a constant ρ_V on cosmological scales:
  - w = −1 to better than 10⁻³⁸;
  - Ġ/G ≈ 0;
  - G_cos = G_lab;
  - no early dark energy;
  - no late vacuum decay.
- So the DESI DR2 w₀w_a preference (2.8–4.2σ) is a retrospective tension for the construction exactly as it is for
  ΛCDM.
- The only distinctive output, the q-oscillation behaving as dark matter, has no computed amount. Transitions deposit
  none.

---------------------------------------------------------------------------------------------------------------------

## 1. What GRUT closure can and cannot supply to q-theory

### 1.1 The closure statements and what each carries [Q from FORMAL_CORE; A]

| Statement | Content | What it carries into a vacuum-energy argument |
|---|---|---|
| DEF 1.1 | finite ℋ_Ω, time-independent H, ρ_Ω(0) supplied, unitary law, "no reservoir, drive or measurement device outside Ω", t supplied | **no exterior.** No geometry, volume, pressure or gravity. Every ρ_Ω(0) is admitted (§4, "a preferred initial state: BY DEFINITION") |
| DEF 1.2 | access structure supplied | nothing about a global vacuum variable. A 3-form Gauss-law constant is a global (central) label, not a subsystem |
| DEF 1.3 | bath = commutant | for the whole, the commutant of B(ℋ_Ω) is ℂ·1: **the whole has no bath.** This is a one-line restatement of "no environment" |
| DEF 1.6 | t external; "in GR no preferred global t and, in general, no conserved total energy … without a stated asymptotic or gauge structure" | P_q supplies such an asymptotic structure (Minkowski in the matter-free limit). It is added, not derived |
| DEF 1.7 | Θ = (algebra, dynamics, preparation, access, instruments, clock), all supplied | q-theory fills the algebra and dynamics slots, and the μ-selection fills the preparation slot |
| THM 2.1 / 2.15 | operators commuting with H are conserved. Protect or form: a conserved O can be neither written nor erased by the closed dynamics | applied to μ̂: a conserved μ is preparation data (§1.5) |
| COND 2.10 / 2.17 | equilibration on average (non-degenerate gaps); relaxation only for an **infinite complement** | the whole has no complement. Equilibration on average is relative to the supplied preparation (d_eff) and says nothing about a pressure |
| COND 2.18 | a distinct prediction differs by response, state or higher cumulants | q-theory + P_q differs from GR + unimodular/HT only by (b), the state (μ), and by (a), the response at ω ~ m_q (Planckian at natural scales) [A] |
| COND 2.20 | apparatus-inside limits | no bearing on the vacuum value. An internal observer can infer μ from the local Λ_eff, but cannot verify a final-state condition before it holds [A] |
| §4 | not forced: thermal complement, relaxation, preferred state. SCOPE: compactness, "closed system ≠ closed universe", global energy conservation in cosmology | equilibrium of the whole and the asymptotic state are **not forced** |
| COND 5.4 | with gravity: no dressing-free subsystems, no preferred t, no conserved total energy | the q-theory argument lies outside the finite class. Using it extends Θ |
| §7 | slot table. "Known kinds outside it: two-time or final boundary conditions (ABL 1964; Gell-Mann–Hartle 1994; Horowitz–Maldacena 2004)" | P_q in its consistent reading is of this final-boundary kind (§1.3 d) |

### 1.2 Volovik's conditions for a zero vacuum energy against GRUT closure [Q; A]

Volovik (gr-qc/0604062 §8, first paragraph; read in the TeX) lists the assumptions behind his zero result: (i)
complete thermodynamic equilibrium; (ii) static or stationary; (iii) homogeneous; (iv) zero curvature; (v) no
excitations (matter); (vi) no interaction with an environment; (vii) T = 0; (viii) no boundaries. In his §2.2
(eqs (3)-(4)), cooling and the absence of external pressure give T = 0 and p = 0, and so Λ_natural = ε_vac = −p_vac
= 0, "if we believe that our Universe does not interact with the environment".

| Volovik condition | Supplied by GRUT closure? | Where it would come from |
|---|---|---|
| (vi) no environment | **yes**, by DEF 1.1 (and DEF 1.3 for the whole) | GRUT wording of the same premise |
| (viii) no boundaries | partly: no exterior, but FORMAL_CORE §4 says compactness is not implied, and a "closed system ≠ closed universe" | supplied geometry |
| (i) equilibrium of the whole | **no**: §4 (relaxation, thermality and a preferred state are not forced); COND 2.17 needs an infinite complement | a state postulate |
| (ii) stationary, (iii) homogeneous, (iv) flat | **no**: no geometry in the finite class (COND 5.4) | the asymptotic state of P_q |
| (v) no excitations, (vii) T = 0 | **no**: the preparation is free | the asymptotic, diluted state |
| a vacuum variable q with Gibbs-Duhem P = −ε + μq | **no**: field content is an entry of Θ (§4, "BY DEFINITION") | KV's 4-form (algebra and dynamics slots) |
| identification of the gravitating term as ε − μq | **no**: no gravity in the finite class | KV's field equations (equations §1.3; literature §2.3) |

### 1.3 Is "no exterior ⇒ zero external pressure on the vacuum of the whole" derivable? [D/A]

No. There are five separate gaps.

**(a) The vacuous part and the part that does the work.**
- "P_ext = 0" holds trivially under DEF 1.1: there is no reservoir to exert a pressure.
- The KV/Volovik step that does the work is mechanical equilibrium, P_vac = P_ext. It needs an adjustable extensive
  variable (a volume V at fixed charge Q = qV in 0711.3170 §II.A) and relaxation to equilibrium.
- DEF 1.1 has neither. The finite model has no V. Relaxation is not forced (§4), and COND 2.10 gives only equilibration
  on average relative to the supplied preparation.

**(b) The dynamical theory has no adjustable volume.**
- In FRW the 4-form conserves μ, not Q = q a³ (literature §3.3). Nothing in the dynamics pushes the vacuum to P = 0.
- P = 0 holds asymptotically iff the conserved μ already equals μ₀ (equations §5.3; literature §3.1).
- The Gibbs-Duhem condition is therefore a condition on a constant of the motion, i.e. on the state.

**(c) Countermodels.**
- Every μ gives a legal closed whole with no exterior:
  - μ = μ₀[ε_∞]: a Minkowski end state;
  - μ above it: de Sitter (the toy run in equations §5.4 case C ends at ρ_V = 0.020032 = Δ + (4/3)Δ²/(K m_q²));
  - μ below it: AdS and a crunch (case B).
- In GR alone, de Sitter with any Λ > 0 has no exterior either. FORMAL_CORE §4 notes that spatially open FLRW models
  have no exterior.
- So "no exterior" is compatible with every value of the gravitating vacuum energy: **NOT FORCED.**

**(d) Which equilibrium.** The GRUT phrase "the equilibrium vacuum of the whole" has three readings.

| Reading | Formal content | Result |
|---|---|---|
| initial: the whole is prepared in equilibrium, and μ "remembers" it (KV 0806.2805 §V.A, read in the TeX) | μ = μ₀[ε_initial] | ρ_V,∞ = −Σ Δ_i (equations §5.3). For SM + Candidate 1: −2.62×10⁹ GeV⁴ = −1.04×10⁵⁶ ρ_Λ, AdS crunch [C]. **Excluded** (retrospective) |
| every equilibrium vacuum, re-imposed after each transition (Volovik 2006 §7; 0711.3170 §II.C) | μ follows the current vacuum | requires Δμ = Δ/q₀ at each transition, which the conserved-μ dynamics forbids (literature §3.1 item 4, §3.3). **Inconsistent** with the action, not a refinement |
| final: the asymptotic T → 0 vacuum after all transitions (KV 1102.3152 §6, via literature §3.2) | μ = μ₀[ε_∞] | ρ_V,∞ = 0 to all orders (equations §5.3). Intermediate plateaux gravitate as in GR. **The only consistent, non-excluded reading** |

The consistent reading constrains the state through its asymptotic future. It is a final boundary condition, which
FORMAL_CORE §7 lists among the "known kinds outside" its table (credited there to ABL, Gell-Mann–Hartle and
Horowitz–Maldacena).

**(e) GRUT's premise favours the reading the dynamics rejects.** "No reservoir" is the thermodynamic setting of an
isolated droplet: fixed Q, free V, and P adjusting to zero. That is the static picture of 0711.3170 §II and Volovik
§7. The 4-form dynamics instead conserves μ, which in thermodynamic language is the grand-canonical variable. So GRUT
closure, read literally, supports the static picture, and the static picture is the one that is inconsistent with the
equations [A].

### 1.4 What exactly is added, and the candidate postulate [D/A]

Added, with its slot (FORMAL_CORE DEF 1.7, §7):
1. **Algebra and dynamics.** A 3-form A with q = F~/sqrt(−g), and supplied functions ε(q) and K(q). The vacuum
   variable is not supplied by closure. Credit: KV 0711.3170 §V.C, 0806.2805 eqs (2.1)-(2.8).
2. **Gravity.** The action of the BRIEF. This is outside the finite class (COND 5.4).
3. **Preparation.** The value of the conserved Gauss-law constant μ. This is a superselection label. The q itself is
   not conserved in FRW; μ is (literature §3.3; equations §1.2). The BRIEF's phrase "the conserved q" should read
   "the conserved μ".
4. **Equilibrium of the whole, as a final condition.** The asymptotic, matter-free, T → 0 state has zero vacuum pressure
   for the full ε_∞.

**Candidate GRUT postulate (one equation; credited to KV 0905.1919 §II eqs (2.3a,b) for the form, Volovik 2006 §2.2 eq
(4) for the "no environment" justification, and KV 1102.3152 §6 for the final-state reading):**

    (P_q)   μ_Ω = μ₀[ε_∞]:   ε_∞'(q₀) = μ₀,   ε_∞(q₀) − μ₀ q₀ = 0,   ε_∞''(q₀) > 0,
            ε_∞(q) ≡ ε(q) + (every q-independent vacuum energy of the final, T → 0, post-transition matter vacuum,
                             renormalized to all orders)

Equivalently, P_vac,∞ ≡ μ_Ω q₀ − ε_∞(q₀) = 0 with ε_∞'(q₀) = μ_Ω: the Gibbs-Duhem pressure of the final vacuum
vanishes.

**Scope.**
- (i) The classical action of the BRIEF on one connected spacetime. μ is the 3-form Gauss-law constant. There is no
  membrane nucleation (Brown-Teitelboim would change μ; equations §1.5). K(q₀) > 0 is required.
- (ii) A history whose matter dilutes, so that a matter-free asymptotic state exists.
- (iii) Existence: c_tot must lie in the range of P(q) = q ε' − ε on the stable branch. For a convex branch through
  q = 0 this needs ε(0) + c_tot > 0 (equations §5.1).
- (iv) Slot: preparation, imposed as a final condition on a superselection label.
- (v) It fixes ρ_V,∞ = 0 exactly, and covariantly under ε → ε + c: μ₀ moves by c/q₀ + O(c²) (equations §5.1).
- (vi) It does **not** fix today's ρ_V to ρ_Λ,obs. That needs the calibrated offset μ − μ₀ = −ρ_Λ/q₀ + O(ρ_Λ²), of
  relative size 7.2×10⁻¹²¹ at natural scales (equations §6.4).
- (vii) Gate Q3: no Lagrangian parameter is retuned, but the selected label moves with every threshold. The tuning is
  carried by the state (equations §6.3), as in unimodular/HT gravity, and is not removed.

### 1.5 What closure does supply [D/A]

1. **Protect or form, applied to μ̂.** In minisuperspace the momentum conjugate to the 3-form zero mode is
   K'R/2 − ε' = −μ, and it is conserved (equations §1.2, checked in `frw_minisuperspace.py`). Quantum mechanically,
   [H, μ̂] = 0. By THM 2.15 the closed dynamics can neither form nor erase the value of μ. So:
   - **Protected (KV's conserved-μ theory).** The zero is stable against thresholds because μ cannot move, and for the
     same reason μ₀ must be put in as preparation (P_q).
   - **Formed (KV's dissipative variants: 1601.00601, 1601.04676, 2207.03453).** μ_eff relaxes, so it is no longer
     protected. Minkowski is then reached only from a finite domain of initial data (1601.00601), from no attractor
     with the dissipation of 1601.04676, or with a brane realization in which μ is a parameter of the action
     (2207.03453) (literature §3.2).

   This dichotomy is GRUT's sharpest contribution here. Its content is KV's own concession (0907.4887 §3; 1601.00601
   §1; 2207.02826 §I). The theorem is standard, and its application to μ is the bookkeeping [A].
2. **Slot classification.** P_q is a preparation postulate of final-boundary type on a superselection label.
3. **It avoids one GRUT cost.** In the finite class, a state rule F(H) that is equivariant under unitary conjugation is
   stationary (calc 1 S9). P_q fixes only the μ sector, not the state within it, so the evolution of matter inside the
   sector is unaffected [A]. COND 5.3 (Fewster-Verch) concerns natural *states* and does not address a superselection
   label [A].
4. **COND 2.18 reading.** q-theory + P_q differs from GR + an integration-constant Λ (HT) only by:
   - (b) the state, i.e. the selected μ;
   - (a) a response at ω ≈ m_q, plus O(H²/m_q²) corrections (equations §3.4).

   At natural scales m_q ≈ 2×10¹⁸ GeV, so the response channel is unobservable.

---------------------------------------------------------------------------------------------------------------------

## 2. Credit: what is KV/Volovik, what is GRUT wording, what is GRUT-native [Q via literature §§2-3 and own spot-checks; A]

| Element | Already in | GRUT wording | GRUT-native? |
|---|---|---|---|
| Isolated vacuum, no environment ⇒ T = 0, p = 0 ⇒ Λ = 0 | Volovik 2006 §2.2 eqs (3)-(4) and §8 conditions (i)-(viii); KV 0711.3170 §II.A, §III.B eq (31); 0905.1919 §II ("self-sustained … without external pressure, at P = 0") | DEF 1.1 "no exterior"; "the whole has no bath" (DEF 1.3, trivial commutant) | no |
| Gravitating vacuum energy is ε − μq (= −P), not ε | 0711.3170 eqs (11), (28)-(30); 0806.2805 eq (4.4) | none | no |
| Conserved vacuum variable (4-form) with K(q) | 0711.3170 §V.C; 0806.2805 eqs (2.1)-(2.8); 0905.1919 eqs (2.1)-(2.2) | "a preferred vacuum variable" (algebra slot) | no |
| μ a global integration constant; choice μ = μ₀ | 0806.2805 §III and §V.A; 0905.1919 §II ("we choose a value μ₀ … the CC problem would be completely solved if we could explain … μ₀") | preparation slot, superselection label | no |
| Final-state μ₀ (asymptotic T = 0 vacuum); plateaux gravitate before release | 1102.3152 §6 (literature §3.2) | "equilibrium of the whole" read as a final boundary condition (FORMAL_CORE §7 "known kinds outside") | no |
| Shift invariance through compensation in μ | Volovik 2006 §4 (μ → μ + U₀ leaves ε_vac invariant); 0711.3170 eq (14) | the obstruction pass's H → H + cI framing | no |
| "The fine-tuning moves to μ" | 0907.4887 §3; 1601.00601 §1; 1601.04676 §II.B; 2207.02826 §I | THM 2.15 protect-or-form applied to μ̂ (§1.5) | no: the theorem is standard, the application is bookkeeping |
| Lorentz invariance ⇒ vacuum with zero 3-momentum in every frame has zero total ε~ | 0711.3170 §III.B eq (32) | none | no |
| A single selection that covers every heavy threshold, made at the final state | 1102.3152 §6; equations §6 | none beyond the slot label | no |

**Verdict.**
- Physics: nothing GRUT-native. Every assumption that P_q uses is in KV or Volovik 2006.
- Wording: "no exterior", "bath = rest", "equilibrium of the whole".
- Organization: the slot classification, the final-boundary reading under FORMAL_CORE §7, and the protect-or-form
  restatement.
- GRUT's contribution is therefore small and organizational. The BRIEF's "almost verbatim" is accurate.

---------------------------------------------------------------------------------------------------------------------

## 3. What q-theory would change for Candidate 1, if adopted (conditional; Candidate 1 stays paused)

### 3.1 Inputs, all from the repository [Q]

- **Benchmark.** v = 1 TeV, the scale at which the ledger quotes −1.2265×10⁶ GeV⁴ (ledger_s08 §13c).
  λ_φ = 0.01 at μ_RG = 0.1 v, and the polynomial constant C = 0.0025 v⁴ = λ_φv⁴/4. With that C, the tree-level broken
  minimum sits at V = 0 (`grut_issue13_wall_history_20261020/modules/thermal_potential/NOTES.md` §1; check passed).
- **One loop.** The new-scalar one-loop vacuum value is ρ_vac,new = −1.22652632365×10⁻⁶ v⁴, with "No adjustment has
  been applied here", and the quantum asymptotic flatness "requires a stated physical total-vacuum matching rule"
  (`grut_candidate_1_quantum_20261018/REPORT.md` §6).
- **H₀ exclusion.** H₀ ≥ 116.66 km/s/Mpc rejects the 100 GeV, r = 1 preparation. Its premises are flat FRW, standard
  BBN, preserved charge, and ρ_vac,tot ≥ −2.72×10⁻⁴⁷ GeV⁴ = −(2.28 meV)⁴ (`ledger/ledger_s08_candidate1.md` §13b).
- **Hidden transition.** The hidden Z2_φ transition has T_c/v ∈ [1.557140, 1.564147]. It is perturbatively second
  order in every resummed scheme, with φ_c = 0 and L = 0, and is expected to be 3D-Ising continuous (thermal_potential
  §0).

### 3.2 The uncancelled one-loop term [C, D]

- **Size.** −1.2265×10⁶ GeV⁴, which is −4.87×10⁵² ρ_Λ and −3.49×10⁻⁶⁸ M⁴.
- **Under P_q.** The term is q-independent, so it is part of ε_∞. The asymptotic ρ_V is then 0 whatever its value or
  sign. The SM loops, of opposite sign and |ρ| ~ 10⁵²–10⁵³ ρ_Λ (ledger_s08 §13c), are absorbed the same way.
- **The cost is a re-selected label.**
  - δμ₀ = c/q₀ and δq₀ = c/(q₀ ε'') (equations §5.1).
  - At natural scales (q₀²ε'' = M⁴, q₀K' ~ K), δq₀/q₀ ≈ δK/K ≈ c/M⁴ = −3.5×10⁻⁶⁸.
  - This is the same kind of number as Volovik 2006 §7.3 eq (20), δG/G ~ 10⁻⁶⁵ for the electroweak case [Q via
    literature §6]. It is unobservable.
- **Existence.** The term needs an equilibrium for the shifted function. On a convex branch through q = 0 this needs
  ε(0) + c_tot > 0. That holds trivially if ε(0) ~ +M⁴, and it is not a constraint on Candidate 1 unless ε is
  specified [E].

### 3.3 The vacuum-matching obligation [D/A]

The quantum report asked for "a stated physical total-vacuum matching rule". P_q is such a rule: stated, credited,
and not derived.
- It replaces the benchmark value statement ρ_Λ(μ_RG) = 0, an R1-type value statement in the vacuum pass (SYNTHESIS
  §4).
- The value statement constrains a coupling, which has to be readjusted at every order. P_q constrains a state label
  instead, and that label is re-selected at every order.
- It matches the asymptotic total to **0**. Matching ρ_Λ,obs needs the calibrated offset of §1.4 (vi). That is the
  cosmological-constant problem moved into the state, not solved.

### 3.4 The H₀-exclusion premise [C]

Under P_q (final reading) with m_q ≫ H₀, today's vacuum energy is:

| Case | ρ_V,today |
|---|---|
| without the calibrated offset | 0 + ρ_m²/(12 K m_q²) = 4.8×10⁻¹⁶⁹ GeV⁴ |
| with the offset | +ρ_Λ |

Either way ρ_vac,tot ≥ −(2.28 meV)⁴ holds.
- The conditional exclusion of the 100 GeV, r = 1 preparation would then rest on a stated postulate instead of an
  unstated matching.
- It stays conditional on flat FRW, standard BBN and preserved charge.
- It does not identify, validate or rule out any field.
- Under the initial reading the premise fails by a factor 9.6×10⁵⁵, so that reading excludes itself (§1.3 d).

### 3.5 The hidden Z2 transition as a vacuum-energy change [C, D, E]

**Size.**
- ΔV = V₀(0) − V₀(v) = C = λ_φv⁴/4 = 2.5×10⁻³ v⁴ = 2.5×10⁹ GeV⁴ at tree level. The one-loop shift of the broken
  minimum is 4.9×10⁻⁴ of this.
- That is 21.0 times the SM electroweak drop Δ_EW = m_H²v_EW²/8 = 1.19×10⁸ GeV⁴, and 9.9×10⁵⁵ ρ_Λ.
- The transition is continuous, so the vacuum energy is released gradually below T_c as the condensate grows
  (φ_b/T_c ≈ 1.1 √|ε| near T_c, thermal_potential §0). It is not released as a latent-heat step.

**Under P_q (final reading).**
1. **Before and during release, the plateau gravitates as in GR**, with ρ_V,pre = ΔV + (4/3)ΔV²/(K m_q²) (equations
   §5.3).
   - At T_c, ΔV/ρ_r = 1.15–1.17×10⁻⁵ for ξ = T_h/T_SM = 1, and 1.2×10⁻⁹ for ξ = 0.1.
   - So there is no vacuum-dominated episode, which is the same as in GR.
   - H(T_c) = 3.47–3.50×10⁻¹² GeV for ξ = 1.
2. **No constant remnant.**
   - m_q/H(T_c) ≈ 5.7×10²⁹ at natural m_q, so the q-mode follows adiabatically.
   - The oscillation energy deposited carries a factor (x/sinh x)² with x = π m_q τ/2 ~ 10³⁰, which is zero for all
     purposes.
   - The second-order term (4/3)ΔV²/(K m_q²) = 3.6×10⁻⁵⁵ GeV⁴ = 1.4×10⁻⁸ ρ_Λ is a transient, not a constant
     (equations §5.3).
   - Check [D]: KV's kick amplitude (9/8)κ_M² H⁴ K'²/ε'' (literature §5.2) gives the same order when κ_M is taken from
     the plateau's trace, κ_M ≈ 4ΔV/ρ ≈ 4.7×10⁻⁵, with H(T_c)⁴ ≈ 5.8–6.0 ρ_Λ.
3. **A scale coincidence, not a prediction [E/A].**
   - H(T_c)⁴ ≈ 6 ρ_Λ at v = 1 TeV is the AHHKM scaling, ρ ~ (T²/M)⁴, credited to astro-ph/0005111 (literature §5.5).
   - It would become a frozen remnant only with KV's phenomenological relaxation law (0905.1919 eq (4.1)), which is
     not derived from the action, plus κ_M ~ 1. This transition supplies κ_M ≈ 5×10⁻⁵.
   - v is a free input of Candidate 1. Choosing it to make H(T_c)⁴ = ρ_Λ would calibrate v, not predict ρ_Λ.
4. **Under the initial reading**, the drop adds −2.5×10⁹ GeV⁴ to the final residual (§1.3 d). Under the static
   reading it would require Δμ = ΔV/q₀, which the dynamics forbids.

A related exposed-data channel, recorded only:
- Chung-Long (1108.5193; literature §7 item 5) note that self-tuning which removes pre-transition vacuum energy would
  shift relic abundances, but only for first-order transitions.
- The final reading does not remove the plateau. The SM transitions and this Z2 transition are crossover or second
  order.
- So no effect is expected [A].

### 3.6 What q-theory would not change [A]

- The particle physics of Candidate 1 is unchanged: S_m does not depend on q.
- The only new field is the q-mode: a scalar with a gravitational-strength Yukawa (α = 1/3) and m_q ~ 2×10¹⁸ GeV at
  natural scales. It has no effect on solitons, walls or relic yields at any computed order.
- The q-mode does not bear on the CDM-equivalence question. That question stays unproved.
- Wall and network energy is matter (T^m), not ε_∞.
- Candidate 1 stays paused, and nothing here is a test of it.

---------------------------------------------------------------------------------------------------------------------

## 4. Exposed-data ledger

### 4.1 Table

All rows are already exposed, so every comparison is retrospective. "Read" is the read level this track achieved. The
machine-readable version, with every number and secondary source, is `results/ledger.json`.

| # | Observable | What the construction produces | Current constraint | Source | Read |
|---|---|---|---|---|---|
| 1 | Ġ/G today, LLR | Ġ/G = −K̇/K; heavy branch 3Ω_m H₀³/m_q² | (−5.0 ± 9.6)×10⁻¹⁵ /yr | Biskupek, Müller, Torre 2021, arXiv:2012.12032 (Universe 7, 34 [R]) | abstract |
| 1b | same, older LLR | | (4 ± 9)×10⁻¹³ /yr; (7.1 ± 7.6)×10⁻¹⁴ /yr | Williams et al., PRL 93, 261101 (2004), gr-qc/0411113; Hofmann-Müller, CQG 35, 035015 (2018) | abstract; **[R]** (publisher page blocked) |
| 2 | Ġ/G, binary pulsar | same, for strongly self-gravitating bodies | (−0.1 ± 0.9)×10⁻¹² /yr | Zhu et al., MNRAS 482, 3249 (2019), arXiv:1802.09206 | abstract |
| 3 | Ġ/G, helioseismology | same, over the solar age | (1.25 ± 0.30)×10⁻¹³ /yr, claimed nonzero at 4σ | Bonanno, Fröhlich, arXiv:1707.01866 (no journal reference in the arXiv metadata) | abstract |
| 4 | G at BBN | K(q_BBN)/K(q₀) | G_BBN/G₀ = 0.99 (+0.06 −0.05) (2σ) | Alvey et al., EPJC 80, 148 (2020), arXiv:1910.10730 | abstract |
| 5 | G at recombination | cosmological G vs G_lab | λ_G = (G/G_N)^{1/2} = 1.0056 ± 0.0091, i.e. G/G_N = 1.011 ± 0.018 (68%) | Lamine et al., A&A 697, A109 (2025), arXiv:2407.15553, Table 1 (Planck PR4 + DESI DR1 + BBN) | section |
| 5b | G, radiation era to now | | variation < 3% (95%); Δ = −0.022 ± 0.023 | Ballardini, Finelli, Sapone, JCAP 06 (2022) 004, arXiv:2111.09168 | abstract |
| 6 | Yukawa, torsion balance | Yukawa with α = 1/3, λ = 1/m_q | \|α\| = 1 excluded for λ ≥ 38.6 μm (95%) | Lee et al., PRL 124, 101101 (2020), arXiv:2002.11761 | abstract |
| 7 | PPN γ, Cassini | γ = (1 − e^{−mr}/3)/(1 + e^{−mr}/3) | γ − 1 = (2.1 ± 2.3)×10⁻⁵ | Bertotti, Iess, Tortora, Nature 425, 374 (2003), DOI 10.1038/nature01997 | abstract |
| 7b | PPN γ, LLR | | (1.7 ± 1.6)×10⁻⁴ | arXiv:2012.12032 | abstract |
| 8 | w₀, w_a | w of ρ_V(t); heavy branch w = −1 | DESI+CMB: w₀ = −0.42 ± 0.21, w_a = −1.75 ± 0.58 (3.1σ). +Pantheon+: −0.838 ± 0.055, −0.62 (+0.22 −0.19) (2.8σ). +Union3: −0.667 ± 0.088, −1.09 (+0.31 −0.27) (3.8σ). +DESY5: −0.752 ± 0.057, −0.86 (+0.23 −0.20) (4.2σ) | DESI DR2 II, PRD 112, 083515 (2025), arXiv:2503.14738v3, Tables 5-6 | section |
| 9 | Ω_m, flat ΛCDM | background matter with constant ρ_V | DESI 0.2975 ± 0.0086; DESI+CMB 0.3027 ± 0.0036 (H₀ = 68.17 ± 0.28); Planck 2018 0.315 ± 0.007 | arXiv:2503.14738 Table 5; Planck 2018 VI, arXiv:1807.06209 | section; abstract |
| 10 | Ω_K | none (the construction assumes k = 0) | DESI+CMB 0.0023 ± 0.0011; w₀w_a+Ω_K with SNe 0.0011 / 0.0006 / 0.0008 ± 0.0013 | arXiv:2503.14738 Table 5 | section |
| 11 | early DE near recombination | pre-transition plateaux | f_EDE < 0.088 (ACT), < 0.12 (P-ACT, P-ACT-LB), < 0.087 (Planck alone, Hill et al.) (95%), for log₁₀ z_c ∈ [3, 4.3] | ACT DR6, arXiv:2503.14454 eqs (14), (16); Hill et al., PRD 102, 043507 (2020), arXiv:2003.07355 | section; number via ACT text |
| 11b | constant early DE fraction | | Ω_e < 0.0071 (TT+lensing+BSH), < 0.0036 (TT,TE,EE+lowP+BSH) (95%) | Planck 2015 XIV, A&A 594, A14, arXiv:1502.01590, conclusions | section (TeX) |
| 12 | extra energy at BBN | vacuum or q-oscillation energy at 1 MeV | N_ν = 2.898 ± 0.141, N_ν < 3.180 (2σ) (BBN+CMB); N_eff = 2.86 ± 0.13 (ACT), 2.99 ± 0.17 (Planck+BAO) | Yeh et al., JCAP 10 (2022) 046, arXiv:2207.13133; arXiv:2503.14454; arXiv:1807.06209 | abstract |
| 13 | Ω_c h² | dust-like q-oscillation | 0.120 ± 0.001 | Planck 2018 VI, arXiv:1807.06209 | abstract |
| 14 | oscillating ultralight scalar DM | q-oscillation for ultralight m_q | m > 2×10⁻²⁰ eV for all of the DM (95%); percent-level fractions for 10⁻³³–10⁻²⁴ eV | Rogers, Peiris, PRL 126, 071302 (2021), arXiv:2007.12705; Hlozek, Marsh, Grin, arXiv:1708.05681 | abstract |
| 15 | vacuum decay into photons | only in dissipative variants | \|y\| < 1.5×10⁻⁵, \|μ\| < 9×10⁻⁵ (95%); T(z) = T₀(1+z)^{1−β} with β_max ≈ 3×10⁻³ | Fixsen et al., ApJ 473, 576 (1996), astro-ph/9605054; Opher, Pelinson, MNRAS 362, 167 (2005), astro-ph/0409451 | abstract |

### 4.2 What the ledger implies for a q-theory construction in this action [C, D, E]

**Bounds on the q-mode (unscreened Yukawa with α = 1/3; equations §4.3).**

| Bound | Value | Label |
|---|---|---|
| Cassini | m_q > 2.9×10⁻¹⁵ eV for impact parameter R_sun; 1.8×10⁻¹⁵ eV for 1.6 R_sun (needs m b > 10.27) | [C] |
| Eöt-Wash | m_q > 5.1 meV at the α = 1 point; for α = 1/3, about 4 meV | [C]; [E], since the α = 1/3 crossing was not read from the exclusion curve |
| LLR | m_q > 58 H₀ = 8×10⁻³² eV in the quasi-static branch | [C] |
| BBN | m_q ≳ 6×10⁻¹⁶ eV | [E] |

- Natural m_q = 2.0×10¹⁸ GeV passes all of them by about 30 orders of magnitude.
- **Every ultralight-q regime is excluded unless the q-mode is screened.** That covers m_q ~ H₀ (needed for a
  matter-induced ρ_V ~ ρ_Λ, equations §6.4) and the fuzzy-dark-matter window (row 14).
- Screening would need ε(q) and K(q) designed as in chameleon f(R) dark-energy models (Hu-Sawicki, arXiv:0705.1158
  [R]). That is a design KV do not make [A].

**Cosmology in the allowed regime.**
- **Rows 1-5.** Ġ/G ≈ 3.4×10⁻¹³¹ /yr at natural m_q [C], and G_cos = G_lab. The data are met trivially, so they do not
  test the construction.
- **Row 8.**
  - Deviation from w = −1: the heavy-branch ρ_V moves by ρ_m²/(12 K m_q²). That is 1.9×10⁻¹²² ρ_Λ at natural m_q, and
    8.8×10⁻³⁹ ρ_Λ even at the Cassini bound [C], so w = −1.
  - Comparison: the construction sits at the ΛCDM point. DESI DR2 combinations disfavour that point at 2.8–4.2σ (3.1σ
    without SNe).
  - Consequence: this is a **retrospective tension shared with ΛCDM**. The construction cannot produce w₀ > −1,
    w_a < 0 without screening or extra dynamics.
- **Rows 11-12.**
  - Every SM and Candidate 1 plateau is gone by T ≈ 150 MeV. At 1 MeV, ρ_Λ/ρ_r = 7×10⁻³⁶.
  - The BBN bound corresponds to an extra fraction of at most 0.028 at 1 MeV [E: treated as extra radiation].
  - The QCD bag plateau, B/ρ_r = 0.04–0.26 at T = 155 MeV, is present in GR and in q-theory with final μ₀ alike [C/R].
  - Early-DE limits are therefore met trivially. They are not a test.
- **Rows 13-14.**
  - The q-oscillation is the only distinctive output. Its amount is not computed in KV (1612.02326; literature §4.6).
  - It is zero from transitions, because of (x/sinh x)².
  - A primordial amplitude would be an input. Matching it to Ω_c h² = 0.120 would be a **calibration**.
  - A Planck-mass oscillation decays within roughly 3×10⁻⁴⁰ s if any non-conformal matter couples to it. The estimate
    uses the scalaron rate Γ ≈ m³/(192π M²) per minimally coupled scalar [E/R], is not computed for this action, and
    is load-bearing for any "q-dark-matter" claim.
- **Row 15.** The non-dissipative construction has no late vacuum decay. KV's dissipative variants (0905.1919 eq
  (4.1); 1601.00601) produce gravitons and neutrinos, which row 12 bounds as dark radiation.

**Overall [A].**
- In the regime that local gravity allows, the construction's cosmological output is GR with a constant ρ_V.
- P_q sets that ρ_V to 0, or to ρ_Λ with a calibrated offset.
- No ledger row offers a held-out test. The one quantity not already constrained to its ΛCDM value is the q-dark-matter
  amount, and it is uncomputed.

---------------------------------------------------------------------------------------------------------------------

## 5. Limits

- §1 is argued [A] on FORMAL_CORE's finite-class statements. Using q-theory at all extends Θ beyond the finite class
  (COND 5.4), so "closure supplies" means "the closure rules that survive that extension". Quantum 4-form issues are
  not analysed: the operator ordering of μ̂, membrane nucleation, and the self-adjointness questions raised in the
  vacuum pass.
- The protect-or-form application (§1.5) uses the minisuperspace conserved momentum. A field-theoretic statement would
  need the Gauss-law superselection structure of the 3-form, which I recall from Brown-Teitelboim [R] and did not
  derive.
- The KV source statements are the `literature` track's, at its read levels. I re-read in the TeX only these:
  Volovik 2006 §2.2, §4 and §8 (opening); 0905.1919 §II (equilibrium conditions and the self-sustained assumption);
  0806.2805 §V.A ("remember the original perfect equilibrium"); 0711.3170 (Lorentz-invariance passage, by search).
- Candidate 1 numbers are quoted from the repository and were not recomputed, except C = λ_φv⁴/4, the ratios and the
  plateau fractions. The plateau fraction uses g* = 106.75 + 4 hidden bosons and ignores hidden thermal masses (≤ 3%
  of T_c).
- Several bounds in §4.2 are estimates: the α = 1/3 Eöt-Wash crossing, the BBN trace fraction κ ≈ 0.1, the
  extra-radiation equivalence at BBN, and the Planck-mass decay rate. The Cassini bound uses e^{−m b} at the impact
  parameter and not the full line-of-sight integral. Its error is O(1) in the logarithm, i.e. a factor < 2 in m.
- Read levels:
  - Abstract only: most rows. The Hofmann-Müller 2018 number is **recalled** (IOP page blocked).
  - Not read: the full confidence-level conventions of the LLR analyses. LLR groups often quote inflated formal
    errors [R].
  - Section level: DESI Tables 5-6, ACT DR6 EDE section, Lamine et al. Table 1, Planck 2015 XIV conclusions.
- The ledger is not exhaustive. It does not include:
  - MICROSCOPE/Eöt-Wash equivalence-principle tests. They are not needed, because the q-coupling is universal through
    the metric.
  - GW170817 c_T, which is unaffected for ω_BD = 0 scalar-tensor theory [R].
  - Atomic-clock tests of constants.
  - Binary-pulsar dipole radiation for a light scalar. A light scalar is already excluded by Cassini.

## 6. Sources and read levels (this track)

| Source | Identifier | Read level | Used for |
|---|---|---|---|
| FORMAL_CORE, WORKING_THEORY | `research/reality_framework/grut_closed_whole_theory_20261019/` | full (FORMAL_CORE §§0-7); WORKING_THEORY §§0-2 and s08 rows | §1 |
| CAMPAIGN_LEDGER, ledger_s08 | same packet, lines 28, 1092, 1181; `ledger/ledger_s08_candidate1.md` §13 | sections | §3 |
| Vacuum pass | `grut_vacuum_principle_20261021/` SYNTHESIS; tracks grut_native §§0, 5-7, 9; review_scope; obstruction §§0, 1.5, 6-7 | full / sections | §§1-2 |
| Candidate 1 quantum report | `grut_candidate_1_quantum_20261018/REPORT.md` §§1, 6 | section | §3 |
| Wall-history thermal module | `grut_issue13_wall_history_20261020/modules/thermal_potential/NOTES.md` §§0-1 | section | §3.5 |
| Sibling tracks | `../equations/NOTES.md`, `../literature/NOTES.md` | full | throughout |
| Volovik 2006 | IJMPD 15, 1987, gr-qc/0604062 | sections §2.2, §4, §8 opening (TeX, scratchpad) | §1.2, §2 |
| Klinkhamer, Volovik | 0905.1919 §II; 0806.2805 §V.A; 0711.3170 (search); others via literature track | sections (TeX) | §§1-2 |
| DESI DR2 II | PRD 112, 083515 (2025), arXiv:2503.14738v3 | section (abstract; Tables 5-6; HTML) | rows 8-10 |
| Planck 2018 VI | A&A 641, A6 (2020), arXiv:1807.06209 | abstract | rows 9, 12, 13 |
| Planck 2015 XIV | A&A 594, A14 (2016), arXiv:1502.01590 | section (conclusions, TeX) | row 11b |
| ACT DR6 extended models | Calabrese et al., arXiv:2503.14454v2 | abstract; EDE section (HTML) | rows 11, 12 |
| Hill et al. | PRD 102, 043507 (2020), arXiv:2003.07355 | abstract; number via ACT DR6 | row 11 |
| Biskupek, Müller, Torre | arXiv:2012.12032 (Universe 7, 34 (2021) [R]) | abstract | rows 1, 7b |
| Williams, Turyshev, Boggs | PRL 93, 261101 (2004), gr-qc/0411113 | abstract | row 1b |
| Hofmann, Müller | CQG 35, 035015 (2018), DOI 10.1088/1361-6382/aa8f7a | **recalled** | row 1b |
| Zhu et al. | MNRAS 482, 3249 (2019), arXiv:1802.09206 | abstract | row 2 |
| Bonanno, Fröhlich | arXiv:1707.01866 | abstract | row 3 |
| Alvey et al. | EPJC 80, 148 (2020), arXiv:1910.10730 | abstract | row 4 |
| Lamine et al. | A&A 697, A109 (2025), arXiv:2407.15553 | abstract; Table 1 (HTML) | row 5 |
| Ballardini, Finelli, Sapone | JCAP 06 (2022) 004, arXiv:2111.09168 | abstract | row 5b |
| Lee et al. (Eöt-Wash) | PRL 124, 101101 (2020), arXiv:2002.11761 | abstract | row 6 |
| Bertotti, Iess, Tortora | Nature 425, 374 (2003), DOI 10.1038/nature01997 | abstract (publisher page) | row 7 |
| Yeh, Shelton, Olive, Fields | JCAP 10 (2022) 046, arXiv:2207.13133 | abstract | row 12 |
| Rogers, Peiris | PRL 126, 071302 (2021), arXiv:2007.12705 | abstract | row 14 |
| Hlozek, Marsh, Grin | arXiv:1708.05681 (MNRAS 2018) | abstract | row 14 |
| Fixsen et al. | ApJ 473, 576 (1996), astro-ph/9605054 | abstract | row 15 |
| Opher, Pelinson | MNRAS 362, 167 (2005), astro-ph/0409451 | abstract | row 15 |
| Klinkhamer, Volovik, "Dark matter from dark energy in q-theory" | JETP Lett. 105, 74 (2017), arXiv:1612.02326 | abstract | row 13 |
| Hu, Sawicki | arXiv:0705.1158 | recalled | §4.2 screening |
| Scalaron decay rate Γ ≈ m³/(192π M²) | Starobinsky 1980; e.g. Gorbunov-Panin arXiv:1009.2448 | recalled | §4.2 |

All downloads are in the scratchpad (`.../scratchpad/qtheory/ledger_dl/`). No paper text is in the repository.

## 7. Files

- `code/ledger_build.py`: every [C] number in §§3-4, 8 checks, and the 15-entry ledger. Run
  `python3 -I code/ledger_build.py results/`.
- `results/ledger_numbers.json`: the Candidate 1 numbers, q-mode bounds, checks and versions.
- `results/ledger.json`: the machine-readable ledger. Each entry has the observable, the construction quantity, the
  value, the CL, the source, the read level, `exposed: true`, `status_of_comparison: retrospective`, and the
  implication.
