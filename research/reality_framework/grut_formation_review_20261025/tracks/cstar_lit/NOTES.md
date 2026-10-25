# Track `cstar_lit`: literature claims and corrections in 8986d34's C* audit

Claude review, 25 October 2026. Read-only on Codex's branches. No commits.

**Inputs.** `git archive 8986d34 research/reality_framework/grut_cstar_audit_20261024/` into the scratchpad
(`review/cstar_lit/packet/`). I read CLAUDE_HANDOFF.md, REPORT.md, MATHEMATICS.md, GLOBAL_CHARGE.md,
CSTAR_HYPOTHESIS.md, LITERATURE_AUDIT.md and `code/closure_checks.py` with CLOSURE_CHECKS.json. The X-charge
radiation numbers are covered by track `xcharge`; this track does not recheck them.

**Independence.** Bibliographic identities come from my own script (`code/verify_refs.py`, arXiv export API and
Crossref). The elementary identities and a Monte Carlo reproduction of the Harlow–Usatyuk–Zhao code formulas
come from `code/identities.py`. I ran no Codex code. Codex statements are quoted [Q].

Labels: [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued, [R] recalled.
Scope: this is a literature-scope review of a research postulate. Nothing here validates or excludes C*, and
nothing here is evidence for the GRUT fields.

## 0. Verdict summary

| # | Claim (8986d34) | Verdict | Note |
|---|---|---|---|
| L1 | All nine arXiv identifiers in LITERATURE_AUDIT.md, plus the Page–Wootters DOI, resolve to the stated papers, authors and journals | CONFIRMED [C] | `results/refs.json`. 2403.13047v2 is the published JHEP version; Codex gives no journal (cosmetic) |
| L2 | Marolf–Maxfield 2002.08950: null states, α-sectors, ensemble vs fixed theory | CONFIRMED | Read §1, §2.3. Precise statement: the Ẑ[J] commute, α-states are superselection sectors, and H_BU is one-dimensional exactly when the path-integral amplitudes factorize (eqs. 2.16–2.23) [Q] |
| L3 | Usatyuk–Wang–Zhao 2402.00098: unique non-perturbative closed-universe state in 2D JT/topological models | CONFIRMED | Read §1. Add: the cosmologies are **Λ < 0** big-bang/big-crunch universes. The general higher-D argument relies on a formal continuation of divergent moments (their fn. 2) [Q] |
| L4 | Usatyuk–Zhao 2403.13047: factorization ⇒ one state in a fixed theory | CONFIRMED | Abstract only |
| L5 | Harlow–Usatyuk–Zhao 2501.02359: 1D real space; non-isometric encoding; observer assumptions | CONFIRMED | Read §1–4.1 and §6. Three stated assumptions: observer dimension e^{S_Ob}, classical pointer basis, entanglement entropy ~S_Ob. Calculations use Λ < 0; a dS realization is "hoped" [Q] |
| L6 | Harlow 2602.03835: distinguishes canonical, baby-universe and holographic spaces; calls the observer proposals speculative | CONFIRMED | Read §1, §1.1, §2.1, §3.1. A 1D space arises only in averaged/factorizing holography; canonical QG and BUFT have large spaces. α-states are well defined only in 1+1 D. The disconnected (Boltzmann-brain) contribution dominates for dS [Q] |
| L7 | Bousso hep-th/0010252: N-bound on observable entropy for Λ > 0; no global count for a fading vacuum | CONFIRMED | Read §1. The bound "becomes trivial in the limit of vanishing cosmological constant", and Banks's Λ–N relation "does not solve the cosmological constant problem except by fiat" [Q] |
| L8 | Klinkhamer–Volovik 0711.3170: q-vacuum equilibrium P = 0 ⇒ ε − qε′ = 0; extra structure | CONFIRMED | Read §I–II.B, eqs. (2)–(11). The argument is the thermodynamics of an isolated, variable-volume medium in flat space, not a GR closure result |
| L9 | Harlow–Ooguri 1810.05337: no global symmetries under AdS/CFT assumptions | CONFIRMED | Abstract |
| L10 | Kawasaki–Nakatsuka 1912.06993: A-term context only | CONFIRMED | Abstract |
| L11 | DESI 30 July 2026 Ly-α page agrees with the standard model and leaves interpretations open | CONFIRMED | Page fetched |
| L12 | Literature coverage is adequate for the C* claims | **CORRECTED** | See §2: missing are the competing observer rules, which give **different dimensions**, the published dissent, the de Sitter observer constructions, and the original finite-Hilbert-space conjecture |
| C1 | Hamiltonian constraint ≠ measured zero sum of matter + gravitational energy (Codex correction 1) | CONFIRMED, refined | §3 |
| C2 | 1D fundamental space and e^{10^122} effective states need different spaces and a dictionary (correction 2) | CONFIRMED, refined | §4. The cited dictionary (HUZ) gives dim ≈ e^{S_Ob}, not e^{S_dS} |
| C3 | Fixed dS capacity and a vacuum relaxing to zero are alternatives; closure supplies neither (correction 3) | CONFIRMED, strengthened | §5 |
| C4 | Gauss law: total gauge charge = 0 on a compact boundaryless space; no constraint on global U(1)_X | CONFIRMED [known] | §6 |
| T1 | Clause tags: 1 known+conjecture, 2 bet, 3 bet+conjecture, 4 conjecture, 5 conjecture+bet | CONFIRMED except clause 5 | Clause 5's "residue = departure from equilibrium" is Klinkhamer–Volovik's own published proposal. Tag it [conjecture, prior art], not a GRUT bet. Clause 1 under-credits the semiclassical relational-clock results (§7) |

## 1. Reference verification [C]

`python3 -I code/verify_refs.py` → `results/refs.json`. All 41 arXiv identifiers and 5 DOIs resolve. Crossref
bibliographic queries also returned the linearization-instability papers (Moncrief JMP 16, 493 (1975) and 17,
1893 (1976); Arms–Marsden–Moncrief Ann. Phys. 144, 81 (1982); Higuchi CQG 8, 1961 and 1983 (1991)). The
identities of these papers are verified. Their content is [R].

Read levels, all from PDFs retrieved into the scratchpad. HUZ: pp. 1–19, 39–41. Harlow 2026: pp. 1–8, 16–17.
Abdalla–Antonini–Iliesiu–Levine (AAIL) 2501.02632: pp. 1–6 and 66–68 (§7.3). UWZ: pp. 1–7. Marolf–Maxfield:
pp. 1–5 and 14–17. Bousso: pp. 1–5. KV: pp. 1–6. Zhao 2602.05939: pp. 1–7. CLPW 2206.10780: §2.3–2.4, read
through a tool summary of the ar5iv HTML; treat as [Q, secondary]. All other follow-ups: abstract only.
MTW is quoted through Guth's MIT 8.962 slides, lectures 23–24 (2018) [Q, secondary]. p. 463: "The spacetime must
be asymptotically flat if there is to be any possibility of defining energy and angular momentum." p. 467: local
gravitational energy-momentum "is not observable." p. 457 ("There is no such thing as 'the energy (or angular
momentum, or charge) of a closed universe'") comes from a web-search summary of the same notes and is
[Q, tertiary].

## 2. What the closed-universe results show, and what is missing from the audit

- **One state.** The one state is precise only in 1+1 D JT and topological models (MM 2020; UWZ 2024) and in a
  fixed AdS/CFT with crunching Λ < 0 cosmologies (UZ 2024). It is conditional on a well-defined gravitational path
  integral with a sum over topologies (or on factorization). McNamara–Vafa 2004.06738 gave a prior Swampland
  argument for a 1D baby-universe space, which the audit omits. Harlow 2026 shows that canonical QG and baby-universe
  field theory are **not** 1D.
- **Contested.** Antonini–Rath–Sasieta–Swingle–Vilar López 2507.10649 read the 1D space as a limit on CFT access,
  not on the universe. Liu 2509.14327 argues for a semiclassical closed universe at finite G. Engelhardt–Gesteau
  2504.14586 argue the opposite for AdS/CFT baby universes. None is cited.
- **The observer dictionary is not unique.** HUZ use a cloned observer and get dim ≈ e^{S_Ob}. AAIL 2501.02632
  require the observer worldline to connect bra to ket and get dim = d² with d ~ e^{O(1/G)}. For dS-JT they
  *suggest* O(e^{2S_0}) ("plan to explicitly compute"). Akers et al. 2503.09681 show that the two rules give
  different dimensions. Nomura–Ugajin 2505.20390 and 2602.13387 condition on observations in an unconstrained
  space without an external observer. Zhao 2602.05939 (JHEP 09 (2026) 252) solves a 1D topological dS model in which
  an inside observer has a Hilbert space of dimension Z_L (the sphere-partition-function analogue) while the
  baby-universe space is classical statistics. Only the audit's sources HUZ and Harlow are cited.

## 3. Hamiltonian constraint (C1)

On compact slices the canonical Hamiltonian is ∫(N𝓗 + N^i𝓗_i), with no surface term (ADM 1962; DeWitt 1967; compare
Regge–Teitelboim 1974 for the asymptotically flat surface term) [R/D]. Its vanishing is a gauge identity. Local
gravitational energy has no gauge-invariant density, so "matter + gravity = 0" is not a sum of measurable parts.
MTW: the total degenerates to 0 = 0 [Q, secondary]. Codex's counterexample is correct: closed dS, a = cosh(Ht)/H,
satisfies every constraint for any Λ > 0 with p = −ρ [D; C below].

**Refinement.** There *is* a precise sense of "zero total charge" for gravity. In perturbation theory about closed
backgrounds with isometries, linearization-stability constraints make the total Killing charges of all
perturbations vanish (Moncrief; Arms–Marsden–Moncrief; Higuchi) [R]. CLPW impose exactly this as Ĥ = H + H_obs and
find the observer-free algebra trivial [Q, secondary]. These are constraints on states. They are not measurements,
and they do not fix Λ. Tryon 1973 is the origin of the popular "zero-energy universe" phrasing [R].

## 4. Finite information and the 1D space (C2)

Bousso derives S_obs ≤ N = 3π/Λ for Λ > 0, resting on the covariant bound. That bound is proved only in limited
regimes (Flanagan–Marolf–Wald hep-th/9908070; Bousso–Casini–Fisher–Maldacena 1404.5635). Banks hep-th/0007146,
with Fischler, conjectured that **the fundamental** Hilbert space has dimension e^N. That is a rival to the 1D claim
if both refer to the same space. In a reconciliation, e^{S_dS} is reinterpreted as an observer/patch dimension. Two
published routes do this, in toy models only: AAIL §7.3 and Zhao 2026. CLPW give the semiclassical analogue, a
Type II₁ algebra with a maximum-entropy state. That algebra is infinite-dimensional and fixes entropy only up to a
constant, so it does not supply the number e^{10^122}.

Under HUZ's rule the effective space is the clone space, so its dimension is ≤ e^{S_Ob} exactly. Inner-product
errors are ≈ e^{−S₂(Ob)}(1 + |⟨φ*|ψ⟩|²). An e^{10^122}-state effective description in this dictionary needs
S_Ob ≳ 10^122 nats, an observer saturating the dS bound. Under AAIL the space is large, but what the observer can
learn is still ≤ log dim H_obs (their fn. 3). Codex's correction stands. Its sentence "the cited research supplies
proposals for that dictionary" is incomplete: the cited proposal yields e^{S_Ob}, and the S_dS-sized proposals
(AAIL, Zhao) are uncited and 2D/1D only.

## 5. Capacity versus relaxation (C3)

The identity ρ_Λ = 24π²M⁴/S_dS holds only for a constant Λ > 0. Bousso: the N-bound is trivial as Λ → 0, and Λ–N
fixes Λ "by fiat". A vacuum relaxing to zero therefore removes the bound that would set the capacity, so the two
cannot share one saturated relation. KV's equilibrium (0711.3170) and relaxation dynamics (0806.2805) need a
conserved q-variable and a specified action. Neither uses closure. The closed dS solution shows that GR closure
imposes no P = 0. CONFIRMED.

## 6. Gauss law (C4)

∫_Σ √h D_iE^i = 0 on a compact Σ with ∂Σ = ∅. For a Yang–Mills theory, physical states are singlets under constant
gauge transformations [D/R]. References: Pedram, AJP 78, 403 (2010), 0912.0225; 1508.06910. The contrast with C1 is
that charge density is a local gauge-invariant observable, so Q = 0 is a real neutrality statement. It covers only
gauged charges and the whole Σ; net charge inside our patch is allowed. For the proposed Z₆ repair, a residual
discrete gauge charge must also vanish mod 6 globally, unless compensated [A].

## 7. Status of the five clauses

| Clause | Precise published realization | Conjectural | GRUT bet |
|---|---|---|---|
| 1 Internal clocks | Constraint structure of canonical GR [known]. Semiclassical relational clocks: CLPW 2206.10780 (observer in dS), Chen–Penington 2406.02116 (inflaton/BH clock, no external observer), Kaplan–Marolf–Yu–Zhao 2410.00111, Page–Wootters toy | Non-perturbative 4D quantization, physical inner product | That GRUT's own matter supplies the clock |
| 2 Compact space | — | — | Yes (topology/size untested) |
| 3 Internal observers | 1D space in 2D/AdS toy models (MM, UWZ, UZ). Observer rules in JT/codes (HUZ, AAIL, Zhao) | 1D for 4D Λ > 0; which rule; contested | That GRUT records realize one rule |
| 4 Finite information | N-bound for Λ > 0 (conditional on covariant bound); CLPW Type II₁ maximum-entropy state | Finite fundamental dim e^{S_dS} (Banks/Fischler) | A capacity-selection law (absent) |
| 5 Vacuum equilibrium | KV q-theory equilibrium and relaxation (model-level) | That nature's vacuum is a q-vacuum; residue = imperfect-vacuum offset (KV's own proposal) | Link to closure; specific source law |

## 8. Computed checks [C] (`results/identities.json`)

| Check | Mine [C] | Codex [Q] | Verdict |
|---|---|---|---|
| S_dS = πc⁵/(ħGH_Λ²) = 8π²M²/H_Λ² = Bousso N = 3π/(Λl_P²), with H0 = 67.4, Ω_Λ = 0.685 | 3.3072e122 nats (all three agree to 1e-9) | 3.3063e122 | CONFIRMED. The 0.03% gap comes from Codex rounding M = 2.435e18 GeV (CODATA gives 2.4353e18) |
| ρ_Λ = 24π²M⁴/S | 2.519e-47 GeV⁴, identity residual 0 | 2.5186e-47 | CONFIRMED |
| Closed dS a = cosh(Ht)/H, k = +1 constraint and acceleration equation | residuals 3.3e-7 (finite-difference truncation) | 2.2e-16 (analytic) | CONFIRMED |
| a ∝ t²: w, horizon d_E/t | w = −2/3; d_E/t = 1.0000000002 at t = 1, 3, 10 | −2/3, 1 | CONFIRMED |
| Gauss law on lattice T³ | Σ div E / Σ\|div E\| ≤ 8.5e-18; a lone charge has a k = 0 mode of 1, so no solution | — | CONFIRMED |
| Toy constraint H = Z⊗1 − 1⊗Z | kernel dimension 2; conditional ⟨σx⟩ = cos θ | 2; cos θ | CONFIRMED |
| Rank-1 criterion Tr(M)ⁿ = Tr(Mⁿ) (UWZ/HUZ eq. 2.3) | rank 1: violation ≤ 1.1e-15; rank 2: ≥ 0.26 | — | CONFIRMED (illustration) |
| HUZ code, fluctuation of ⟨φ\|V̂†V̂\|ψ⟩ (eqs. 1.3, 4.2) | Monte Carlo matches the exact Haar formula, Var = (d/(d+2))(1+\|⟨φ*\|ψ⟩\|²)/d_obs − 2\|⟨φ\|ψ⟩\|²/(d+2) [D], within MC noise (≤ 1.7%). For d = 4096, d_obs = 16: 0.0936 vs HUZ leading 0.0938 | not in packet | HUZ formula REPRODUCED |
| HUZ effective dimension | images of 6 (12) generic states have rank exactly d_obs = 2 (4) | not in packet | dim_eff ≤ e^{S_Ob} CONFIRMED |

## 9. Limits

These are section-level reads, not replications. The CLPW and MTW content is secondary. Moncrief/Higuchi content is
[R]. I did not check the X-charge numbers (track `xcharge`). Literature after 2026-10-14 was not searched.

## Reproduction

```
python3 -I code/verify_refs.py   # network: arXiv API + Crossref -> results/refs.json
python3 -I code/identities.py    # offline, numpy only -> results/identities.json
```
