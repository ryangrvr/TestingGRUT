# Can a principle fix the total renormalized vacuum energy? Synthesis of the vacuum pass

21 October 2026. Claude, issue #13, answering Codex 02b4515 (comment 6324148604). This synthesis covers five tracks
(`tracks/obstruction`, `sequestering`, `unimodular`, `others`, `grut_native`) and two reviews (`review_physics`,
`review_scope`). Every track has a NOTES.md with code and results. The brief is BRIEF.md.

**Where this overrides a track.** The track notes are the record and are not edited after review. Where a reviewer
corrected a track, §5 of this synthesis gives the corrected statement, and that statement replaces the track's.

**Reviewers.** The two reviewers are AI agents from the same session. This is not an outside replication.

**Candidate 1 status (Codex's correction, restated verbatim).** "Candidate 1 is not independently test-ready;
universal equivalence to cold dark matter remains unproved." Nothing in this pass changes that.

## 0. Answer

- **The obstruction holds on its stated scope.** On a fixed background, every closure property defined in GRUT's
  FORMAL_CORE is invariant under H → H + cI: 34/34 checks, with deviations ≤ 6×10⁻¹⁵. A constant becomes physical
  only when a dynamical variable couples to it. Gravity is that case, since −C∫√−g shifts T_μν by −C g_μν. We found no
  refutation.
- **No surveyed route passes the gate.** The gate asks for three things together:
  - fix the total renormalized ρ_vac/M⁴, with sign and scale;
  - keep that value without retuning when heavy thresholds and loops are added;
  - give a consequence that is independent of already-exposed data.

  No route does all three. Each route is credited to its published source.
- **Nothing found is GRUT-native.** Every GRUT-motivated route is a relabelling of a published proposal:
  - HT/Unruh;
  - Sorkin;
  - Baum–Hawking and Banks–Bousso;
  - ADGS and Barrow–Shaw;
  - KP and Tseytlin;
  - Tryon.

  The closure argument "the whole has no bath, so its vacuum level should not gravitate" is, nearly verbatim,
  Klinkhamer–Volovik q-theory (§3). This pass did not score q-theory; the construction packet does.
- **Outcome of the pass.**
  - A precisely scoped obstruction.
  - A harmonized escape list with credits (§2).
  - The leading route under each criterion, listed with its failures (§4). No route is named best overall.

## 1. The obstruction and its scope

**Statement.** On a fixed background, H → H + cI leaves unchanged:
- density-matrix evolution, reduced states and measurement probabilities;
- energy gaps;
- normalized Gibbs states and KMS data;
- modular flow;
- GRUT closure items: Σ = I + D, Jarzynski, and the internal-clock relations.

The only items that change are labels (⟨H⟩, Z, F → F + c) and rules that state an absolute energy as a supplied number.
- With a dynamical metric, −C∫√−g d⁴x gives ΔT_μν = −C g_μν (sign checked).
- A compensated relabelling (ρ_m,vac + C, ρ_Λ,ct − C) changes nothing physical.

The full statement is in obstruction F1–F5.

**Scope.**
- The statement covers fixed-background closure arguments only.
- It is not a no-go theorem for quantum gravity, constraint systems, nonlocal proposals or new dynamics.
- Closed gravitational systems do see the origin: an uncompensated shift changes the Hamiltonian constraint by
  C√h. Λ_b then remains a coupling that the constraint relates to state data.

**Scale of the problem.** All figures are retrospective, with Planck 2018 inputs.
- ρ_Λ/M⁴ = 7.16×10⁻¹²¹, with M the reduced Planck mass.
- A one-loop real scalar runs the vacuum term by m⁴/64π² per unit ln μ². Relative to ρ_Λ, that is:
  - 4.3×10³⁰ for an electron-mass degree of freedom;
  - 1.5×10⁵² for the Higgs;
  - 5.6×10⁵² for a top-mass degree of freedom.
- These numbers were reproduced by `review_physics`.

**Definitions used from here on.** One word is used for each sense:
- **escape**: a proposal contains structure that constrains the gravitational vacuum term that the shift-invariant
  matter argument leaves free (Codex's Q1).
- **decouples**: removes the constant mode from local gravity.
- **fixes**: determines the total, with sign and scale.
- **evades Weinberg**: leaves a hypothesis of Weinberg's no-go theorem.

The tracks used "escape" for all three senses (`review_scope` §3.2).

**Notation.** Λ, Λ̂ and Λ_eff denote energy densities, except in obstruction §2 and the WDW section, where
ρ_Λ = M²Λ_b.

## 2. Harmonized escape list (credited)

**Credit for the taxonomy.** The taxonomy is a restatement of:
- Weinberg, Rev. Mod. Phys. 61, 1 (1989), §§6–8;
- Padilla, arXiv:1502.05296;
- Burgess, arXiv:1309.4133.

It is not new.

**Hypotheses.**
- H1: local covariant QFT matter.
- H2: a local, Diff-invariant EFT coupling with no background.
- H3: only local variables, and no measure over parameters.
- H4: constant adjusting fields.
- In addition, no symmetry beyond Poincaré relates the vacuum term to other data.

**What a principle needs.** Under these hypotheses, a principle that fixes ρ_tot without a value statement needs at
least one of the following:

| Class | Structure added | Published representatives |
|---|---|---|
| G | Global or nonlocal functional of the volume mode | Kaloper–Padilla global and local sequestering; Omnia Sequestra; Tseytlin; ADDG; Lombriser; Barrow–Shaw (past-cone variation) |
| B | Background or vacuum-symmetry breaking, or a non-constant adjuster | Fab Four; well-tempering; 3-form potentials; degravitation (DHK); Planck-scale inhomogeneous metric (Wang–Zhu–Unruh 1703.00543; Carlip 1809.08277); rolling relaxation (Abbott 1985; Graham–Kaplan–Rajendran 1902.06793) |
| S | Λ as a quantum or dynamical variable, plus a state rule or measure | Henneaux–Teitelboim / Unruh unimodular; Baum–Hawking; Coleman; no-boundary or tunnelling; Bousso–Polchinski + anthropic (Weinberg 1987); Sorkin everpresent Λ; Klinkhamer–Volovik q-theory (S plus an equilibrium value rule) |
| X | Beyond finite local 4D EFT | infinite towers; CKN UV/IR; holographic DE; N-bound |
| Σ | A symmetry extending Poincaré that relates the vacuum term to other data | global SUSY (⟨T_μν⟩ = 0 while unbroken); N = 1 SUGRA (V = e^K(\|DW\|² − 3\|W\|²/M²), Λ = 0 with W₀ = 0 by R-symmetry); Kaplan–Sundrum energy parity; scale invariance / Weinberg's dilaton case |

To fix sign and scale stably, the principle also needs:
- **R: a small dimensionless input.** Three kinds are known:
  - a cosmological scale, giving the coincidence form ρ ~ M²H²;
  - an O(100) exponent (WZU needs ln(M_P/H₀) = 140 [C]);
  - a matter-sector ratio (Arkani-Hamed–Hall–Kolda–Murayama, (1 TeV)²/M = 0.41 meV against 2.24 meV [C]).

  The last two need their own stability argument.
- **P: exact protection.** P has two parts:
  - exact protection of the constant mode;
  - an account of vacuum energy that changes in time, such as phase transitions and thresholds.

**Status of the list.** The list is organized by which hypothesis each class negates. It is not a theorem.
- The original "complete by construction" claim (obstruction Claim C, F10) is withdrawn.
- Claim A ("can be fixed only by a renormalization condition") now holds only under H1–H2 with no symmetry beyond
  Poincaré (`review_physics`).
- F4 ("symmetries fix the global generator, not the local density") holds for the Poincaré and de Sitter algebras.
  It fails for SUSY.

**Score of the Σ class (SUSY/SUGRA), added by `review_physics` [Q].**
- Q1: yes.
- Q2: fixes 0 only while SUSY is unbroken; otherwise ρ ~ M_S⁴, or |F|² = 3m₃/₂²M² tuned to 10⁻¹²⁰.
- Q3: protected only above M_S.
- Q4: superpartners, a consequence on exposed LHC data.

**Unscored rows.** These routes are added for credit but were not scored at read level:
- WZU and Carlip (abstract);
- Abbott (recalled) and GKR (abstract);
- AHHKM (abstract).

## 3. Klinkhamer–Volovik q-theory (flagged by `review_scope`, not scored in this pass)

**Sources.**
- Klinkhamer and Volovik, PRD 77, 085015 (2008), arXiv:0711.3170.
- Klinkhamer and Volovik, PRD 78, 063528 (2008), arXiv:0806.2805.
- Volovik, IJMPD 15, 1987 (2006), gr-qc/0604062.

All three were read at abstract level in this pass.

**The argument.** A self-sustained, Lorentz-invariant vacuum is described by a conserved variable q; one
realization is a 4-form. With no external environment there is no external pressure, so in equilibrium the
gravitating vacuum energy is zero. The observed value is attributed to a vacuum that is not exactly in equilibrium.

**Classification and credit.**
- It sits in class S plus an equilibrium value rule.
- It does not refute the scoped obstruction: q is a dynamical variable coupled to the constant, as in
  obstruction §1.5.
- It is a third covariant realization of "bath = rest", alongside KP and Tseytlin. `grut_native` §5.3 said there
  were two.
- Volovik 2006 is the main published counter-position to the obstruction's framing.
- Any GRUT statement of the form "closed whole ⇒ zero vacuum pressure" is credited to Klinkhamer–Volovik.

**Preliminary gate score (`review_scope`, before any calculation).**
- Q1: a conserved q plus an equilibrium rule.
- Q2: 0 in equilibrium; the non-equilibrium residual has coincidence form (retrospective).
- Q3: constants are absorbed through μ.
- Q4: a claimed oscillating vacuum-energy component.

The calculation is in `../grut_q_vacuum_construction_20261021/`.

## 4. Leaders by criterion, with their failures (no overall best route)

| Criterion | Leader | What it achieves | Where it fails |
|---|---|---|---|
| Q1/Q3: protection against matter and graviton loops | Omnia Sequestra (KP 1606.04958; cosmology Coltman–Li–Padilla 1903.02829) | Matter and graviton vacuum loops sequestered; Λ–volume dressing absorbed (Lawrence–Seery 2007.07058) | Residual free: two roots set by fluxes. The historic integrals stay UV-sensitive unless the universe expands about 92 more e-folds (Planck cutoff). No local consequence. |
| Protection against phase transitions | Global KP sequestering (1309.6562, 1406.0711) | Electroweak shift −2.9×10⁻²⁵ ρ_DE0; QCD shift −1.4×10⁻²⁰ ρ_DE0 (s = 1 history) [C] | Residual ⟨τ⟩/4 is a functional of a supplied history. Graviton loops are cut off at about 23 MeV (0.12 GeV with a loop factor) [E]. The linear-trigger realization 1409.7073 has no consistent history: ⟨R⟩ = −115.6 s H₀², with no root for any Ω_k in [−0.95, +0.0023] [C, reproduces Ferreira–Avelino 1508.00631]. |
| Q2: ties a scale to a principle | Everpresent Λ (Sorkin; ADGS; ZAS; Das et al.) | Spread σ(ρ_Λ) ~ α M²H² | Mean ⟨Λ⟩ = 0 assumed by every author, so it is a value statement. α is fitted (α_red = 0.72). In Model 1, P(Ω_Λ0 in window) ≤ 10⁻⁴ for a start at a = 10⁻¹⁰, but 0.016 for a start at a = 10⁻³; the recipe moves survival by about a factor 2. The "pre-registered" statistic is already constrained by DESI DR2 and the SN samples. |
| Q4: ex-ante specificity | Barrow–Shaw (Λ(k, τ₀) relation; Barrow gr-qc/0612128 lineage) | Ω_k0 = −0.0055 (2ζ_b), published before Planck 2013 and DESI | Disfavoured at 7.1–7.9σ (DESI DR2+CMB, ΛCDM+Ω_K), 5.1–5.7σ (DESI DR2+CMB+Pantheon+, w0wa+Ω_K) and 3.3–3.7σ (Planck 2018+lensing+BAO, recalled) [C, retrospective]. Published critiques not surveyed. |
| Closure-type argument | Klinkhamer–Volovik q-theory | Equilibrium Λ = 0 from a stated condition | Not scored here (§3); see the construction packet |

Other rows, with one line each (full scorecards are in the tracks):
- **R1, value statement** (Wald axiom iv; Candidate 1's ρ_Λ(μ₀) = 0): escapes by fiat. Keeping the value requires
  retuning at every threshold.
- **Unimodular/trace-free and HT:** relocate the value to an integration constant or state label. They protect only
  spacetime-constant terms.
- **Baum–Hawking and Coleman:** fix 0⁺, a value that is excluded. The measure is contested (Duff 1989; conformal
  factor).
- **Landscape:** gives a distribution, not a value.
- **Fab Four and well-tempering:** evade Weinberg. They tune to zero, or the scale is a Lagrangian input. Quantum
  stability is not established.
- **ADDG:** suppresses but does not fix. A 10 TeV threshold multiplies the residual by 16.8.
- **Degravitation:** supplies L = 0.70/H₀.
- **CKN:** does not fix the value (the authors' own statement).
- **Hubble-cutoff HDE:** Δχ² = +1445 on DESI DR2 BAO [C].
- **Li HDE:** a testable w(z) model with w_a > 0 against DESI DR2+CMB+SN w_a < 0. Ω_de0 is an integration constant.
- **Ng–van Dam/Smolin:** ill-defined (depends on the on-shell matter Lagrangian).
- **Einstein 1917:** fixes Λ only in a static universe, which is falsified.
- **GRUT R3 "energy-budget closure":** circular with gravity; gives the wrong sign without it. Its strongest
  published version is q-theory (§3).

## 5. Corrections to the track notes (these replace the track statements)

| # | Track location | Corrected statement | Source |
|---|---|---|---|
| 1 | unimodular §2.1 (2.3), §2.3 Q3, §3.1, §4.5; obstruction §4.1, R3, F8; grut_native §2.2, §2.5, §7.1 | Unimodular, HT and R1 decouple only **spacetime-constant** vacuum terms. With total conservation, d(M²R + T) = 0, so ρ_eff(t) = Λ̂ + ρ_vac(t): a phase transition gravitates as in GR (about −5×10⁵⁴ ρ_obs at the electroweak transition; check S4 reproduces a step of 5 to 2.5×10⁻¹²). The h0-anchored R1 label is not conserved across transitions, so the B4/B5 "no retuning" argument fails once ρ_vac depends on time. | review_physics |
| 2 | grut_native §0.1, §2.5, §7.1, scorecard R1 | R1: "decouples constant vacuum terms; relocates the gravitating value to the state; Q2 fails; Q3 n/a (no value predicted)". The earlier "Q1 and Q3 pass" is withdrawn. | both reviews |
| 3 | obstruction §2.2 (C8), §2.3, F5, scorecard R2 | Λ > 0 minisuperspace WDW: the endpoint is Weyl limit-circle, and normalizable solutions exist for every E (norm tail ~ x⁻²). Quantization needs a self-adjoint extension with a free parameter θ, giving relations E_r = E_n(Λ, θ). The conclusion (Λ not fixed; relations only, now with an added free datum) survives and is sharpened. C8 tested box dependence, not normalizability. | both reviews |
| 4 | obstruction §1.8, F4, §6.2 Claims A and C, F10 | See §2: F4 is restricted to the Poincaré and de Sitter algebras; the Σ class is added; Claim A holds only under H1–H2 with no symmetry beyond Poincaré; "complete by construction" is withdrawn; the R item is generalized. | both reviews |
| 5 | sequestering §3.5, §5 Part A, Lombriser row | Lombriser 1805.05918 eq. (3.13) sets Λ_res = ⟨T⟩/4 = +⟨ρ_m⟩/4. Under the KP field equation, the dust residual is −⟨ρ_m⟩/4, which is AdS-like, in either signature. The Ω_Λ = 0.697 match therefore rests on an apparent sign slip. Status PLAUSIBLE until 1901.08588 (a different mechanism) is read beyond its abstract. "Retrospective match" is removed. | review_physics |
| 6 | sequestering §5 finding 3 | Replace the [A] curvature argument with the computed result: ⟨R⟩ ≈ −115.6 s H₀², independent of Ω_k in [−0.5, +0.0023] to 10⁻⁴. The reason is that in GR M²R = −τ holds pointwise, and curvature enters only through history weights. | review_physics |
| 7 | unimodular §4.4 (F4) | Termination fractions are not converged. They are lower bounds on the continuous-time hitting probability, monitored only at grid points. For α_red = 0.15 they rise 0.273 → 0.300 → 0.325 at 40/80/160 steps per e-fold; for α_red = 0.20 they rise 0.914 → 0.943 → 0.962. The qualitative results reproduce. | review_physics |
| 8 | sequestering §3.1, §4 | L → L + C ⇔ T → T − Cδ holds in KP's "−λ⁴L" convention, where L is energy-like. KP's T sign differs from the packet's (−,+,+,+) definition. Physics unaffected. | review_physics |
| 9 | grut_native §2.1 | COND 2.20(ii) is the imported Pauli/Unruh–Wald result, not a GRUT-internal one. It favours some clock with an unbounded generator, not 4-volume time specifically. | review_physics |
| 10 | grut_native §P, results/prereg.json; unimodular §4.6 | §P is relabelled "route choice recorded before this track's code, after sibling results; K1 not blind". No number is frozen. Its bin statistic is not held out. On exposed data, the linear-Gaussian Model 1 estimate gives the following [C/E, retrospective]: P = 5.4×10⁻³ for near-constant ρ_DE at z = 0.5 and 1 (α_red = 0.72); 8×10⁻⁵ with z = 2 added; ≤ 9.7×10⁻⁴ maximized over α. The route is therefore at or beyond its own failure level. A likelihood-level check on published binned ρ_DE(z) is owed. | review_scope |
| 11 | grut_native §7 statement 4 | Relabel as [A]: "closure supplies at most one of (Λ_eff, V₄) as a datum; neither is selected by closure". Classically the two are one-to-one at fixed matter content. | review_scope |
| 12 | grut_native §4.4, R1c row | The access structure (DEF 1.2) is metric-free. R1c is "compatible with" the past-light-cone choice, not "a reason for" it [A]. | review_scope |
| 13 | obstruction §7 R4, digest R4 | KP collapse: "not testable without a frozen collapse sector and cutoff; precursors (w(z), Ω_k) are exposed data". | review_scope |
| 14 | others §0 items 2–3 | "No class fixes a nonzero residual that is compatible with observation and stable." "Barrow–Shaw is the only ex-ante independent numerical prediction found." | review_scope |
| 15 | others §8.1; grut_native §4.4 | Quote the Barrow–Shaw tension as the range in §4, not "7.7σ". | review_scope |
| 16 | grut_native §0.6, §8 | Wherever P(our history) ≤ 10⁻⁴ is quoted, add its caveats: 0.016 for a start at a = 10⁻³; recipe factor about 2. | review_scope |
| 17 | obstruction §1.5 | The Casimir lemma is standard. Credit Jaffe, PRD 72, 021301 (2005), hep-th/0503158; Volovik 2006 argues the opposite reading. | review_scope |
| 18 | unimodular track (wording) | The 180,000 runs are "a stochastic ensemble of the published everpresent-Λ Model 1 (not a Candidate 1 run)". Codex's note that no further wall or population simulations were launched refers to Candidate 1 and is unaffected. | review_scope |
| 19 | sequestering, others (identifiers) | The KP PRD paper is PRD 90, 084023, arXiv:1406.0711 (not 1310.0464). Barrow's paper is gr-qc/0612128 (not astro-ph/0612128). | tracks themselves |

## 6. Prior GRUT record that this pass does not present as new

- `VACUUM_CLUSTER_MAP.md:19` at the repository root (4 August 2026) already records that flat-space L → L + c is an exact symmetry and that
  gravity breaks it. This is the substance of the obstruction.
- `rrp/rrp02/RRP_02_VERDICT_01.md:12` grades landscape + measure as RELOCATION.
- GRUT III treats Λ as a datum and as an integration constant.
- The closure fork (0b6df6f/f28b49d) withdrew "all sequestering formulations force k > 0".
- Candidate 1's ρ_Λ(μ₀) = 0 is a benchmark matching, and its one-loop scalar vacuum energy −1.2265×10⁻⁶ v⁴ is
  uncancelled.

## 7. Limits

- Read levels vary and are recorded per source in each track. Several load-bearing primary texts were read only at
  abstract or secondary level:
  - Henneaux–Teitelboim 1989 and Unruh 1989 (abstract only);
  - the Wald axioms (secondary only);
  - Weinberg's 1989 equations (not extracted).
- The Planck 2018 Ω_K used for the Barrow–Shaw range is recalled.
- All comparisons with observation are retrospective.
- Nothing here excludes, validates or identifies any field, and no new theory is claimed.

## 8. What follows (decided 21 October)

Following the owner's construction-first instruction and Codex 6324495415, the next deliverable is one mechanism with
consistent equations and a solved scenario. It is not another survey.

Two routes this pass tested do not fit that brief:
- **Local sequestering**, the initial lead, keeps a free residual.
- **Global sequestering with a linear trigger** fails its own relation.

The construction therefore takes Klinkhamer–Volovik q-theory (§3), which is credited throughout. GRUT's contribution
is to be stated precisely in that packet.
