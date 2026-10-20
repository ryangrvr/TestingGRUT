# Campaign ledger entry — s05_grut_III: GRUT III / observer bath (P0–P5, P1 identification, modular KMS, GR recovery, universal coupling)

**Verified and corrected entry (adversarial re-read, 2026-10-20).** Repository `/home/user/TestingGRUT`, branch `codex/grut-closure-fork-20261019`. The working tree is at `a6be0e4` (a WIP descendant of the handoff commit `96a7e7e`); `git diff --stat 96a7e7e HEAD -- GRUT_III research/reality_framework/grut_closure_fork_20261019 research/reality_framework/grut_closed_whole_handoff_20261019` is empty, so every sector source below is byte-identical to its `96a7e7e` version. Paths are relative to the repository root unless prefixed by a branch (`branch:path:lines`, read with `git show`). Nothing was executed, modified or rerun; `GRUT_III/RESULTS.md` is the run record generated 2026-10-12 (`GRUT_III/RESULTS.md:3`). Abbreviations: `01` = `GRUT_III/01_PRINCIPLES.md`, `02` = `GRUT_III/02_DERIVATIONS.md`, …, `RTR` = `GRUT_III/RESPONSE_TO_REVIEWS.md`, `RES` = `GRUT_III/RESULTS.md`, `CSR` = `research/reality_framework/grut_closure_fork_20261019/CLOSURE_SELECTION_RESPONSE.md`, `HANDOFF` = `research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md`, `CP04` = `origin/codex/grut-causal-response-05:research/constructive_physics/attempt04/REPORT.md`, `CP05` = `…/attempt05/REPORT.md`.

## 0. Verifier's corrections to the earlier entry

The earlier entry was substantially correct: every cited line range I re-read supports the claim attached to it, and all but one quoted number match the source exactly. The following were found and are corrected below.

| # | Earlier claim | Problem | Deciding check | Fix |
|---|---|---|---|---|
| 1 | "G⁺ quadrature vs Planck rel.err ≤ 1.55e-31" | Number mismatch. 1.55e-31 is the w = 1 entry; the w = 5 entry prints 1.17e-17, so the bound is false. | `RES:32-38` (rel.err 0.0, 0.0, 1.55e-31, 0.0, **1.17e-17**, 1.07e-31, 0.0) | Worst rel.err 1.17e-17 (w = 5); all seven PASS. |
| 2 | "two closed forms agree to 1e-26" (v7) | Precision conflation. The two closed forms (position-space vs momentum-space) agree symbolically; 1.14e-26 is quadrature vs closed form at t = 0.05L. The entry followed `08:221-222`, which itself conflates the two. | `RES:250` (symbolic PASS, no number) vs `RES:258` (rel.err 1.14e-26) | Stated separately. |
| 3 | "smeared KMS to ≈ 2.3e-13" (v7 §C) | Precision conflation. ~2.3e-13 is the residual of the FDT identity N_L = ½coth·ρ_L and of ρ_L against part A; the direct KMS ratio G⁺/G⁻·e^{−βω} − 1 reaches 5.14e-11 (w = 2, L = 1). `02:134` says "verified to 1e-13". | `RES:317-385`, worst KMS-ratio line `RES:332-333` | Both numbers given. |
| 4 | "Dominated by recovered known physics (D1–D4 … IMPORTED)" | Misclassification of D2: its tier is "DERIVED; corollary of standard results", not IMPORTED. D5 is listed as a conditional constraint, but its derived content is "contained in HLZ … a verification, not a discovery". | `02:25`; `08:126-127` | D2 = DERIVED corollary; D5 derived content = recovered known physics (scope caveats retained). |
| 5 | 08 "gives CP-04 only the pattern 'CONDITIONAL ONLY' rather than its actual verdict" | Overstated. CP-04 prints two verdict labels: the header (line 3) and §6 (line 156), and §6 literally begins "CONDITIONAL ONLY as a published vacuum-constraint construction. NOT YET ESTABLISHED as a viable complete GRUT theory." So `08:87`'s "verdicts 'CONDITIONAL ONLY'" matches one of CP-04's own labels. What `08` omits is the header label. | `CP04:3`, `CP04:156`; `08:86-87` | Discrepancy restated. |
| 6 | "Modular KMS is automatic for any faithful state" | The source says "faithful normal state". | `CSR:42-43` | "faithful normal state". |
| 7 | Source range `REPORT.md:170-185` | File has 181 lines. | `wc -l` → 181 | `REPORT.md:168-181`. |
| 8 | Omitted caveats | D11's second-order formula is complete only under the full adiabatic hierarchy (ε₃ = η̇/Hη also small); the massless velocity correlator is IR-finite only for ε < 2/3; no Fourier spectrum at t₀ is computed. | `11:125-126`; `11:158-159`, `RES:668-669`; `02:186`, `11:59-60` | Added to §2/§7. |
| 9 | Omitted failures | (a) `10_SECULAR_IR_SURVEY.md` carries **two conflicting graviton verdicts**: the graviton reading report's bottom line "the literature's resummations do NOT yet agree … This is a genuine, open window, not a closed one" (`10:104`) versus the author's §4 verdict "a quantum-gravity question, not a GRUT one" (`10:138`); and the synthesis agent received only the two scalar reports intact (`10:7-9`). (b) `13` §8's adversarial review withdrew three of `13`'s own statements ("GRUT derives the initial state"; N1 originally asked an already-answered question; the DMP ν conflated with the slow-roll index) (`13:344-353`). (c) `06`'s first revision-2 text said "undetermined coefficient", corrected from Milgrom's paper (`06:18`). (d) Garbrecht–Prokopec's second-order response is non-thermal, so "the record's P1 must be understood as KMS of the first-order kernel" (`10:130`). (e) Abbott–Wise 1984 and Lucchin–Matarrese 1985 were never read (`10:110`, `10:167`). (f) Minor open v8 items and a unit-label misprint (`11:42-45`, `11:139`). | as cited | Added to §8. |
| 10 | Omitted: framework-level falsifiers | `03` lists falsifiers of the framework itself distinct from F1–F6 (Lorentz violation at horizon-thermality level; area-law failure; contextual outcome frequencies; non-polynomial pointwise dissipative kernel). | `03:40-44` | Added to §3.14. |
| 11 | Omitted: third horizon-temperature candidate | `09` records T₀ = 1/(2πR_A) = H/2π (Cai–Kim; Cai–Cao–Hu; Tian–Booth) as "likewise a geometric label"; T_det/T₀ − 1 = −ε exactly. | `09:30-33`; `RES:676-678` | Added to §3.7. |
| 12 | Unreconciled tension not flagged | `09:25` says P1 assigns no temperature inside R_A for ε ≠ 0 "an application of P1's definition, not a computation", while `13:163-167` says KPV 2021 proves a scaling-limit KMS temperature at the Hayward surface gravity for trapping horizons (cosmological case open). The sources do not say whether P1's geometric list (`01:35-39`) should be extended by that theorem. | `09:23-30`; `11:156-157`; `13:163-167` | Flagged in §7. |
| 13 | Stale values inside the sources (not entry errors; recorded so they are not quoted) | `02:177` and `RTR:189,194` still print "c = 1 to 3e-5"; the Richardson-weight fix gives c = 1 to 1e-9 and `RES` prints ≈1e-9. `12:4` says "fourteen runnable scripts"; `RES:5-22` lists sixteen (v15, v16 added in rev. 2.5). `02:5` says "Scripts: calc/v1–v11". | `11:121-123`; `RES:611-614`; `CHANGELOG.md:14-21` | Noted in §3.5/§3.17. |

Everything else in the earlier entry (P1 status, modular-KMS statement, GR-recovery and universal-coupling imports, F1–F6, the v11/v12/v13/v14/v15 numbers, CP-04/CP-05 labels, the closure-fork Obadia finding, the handoff row) was re-checked line by line and stands.

## 1. Question asked in this sector

The program's parent question — what a closed system whose every subsystem's bath is the rest of the universe forces beyond known physics — was asked here in its most precise algebraic form. GRUT III replaced GRUT I/II's supplied mechanism by an identification: "the medium an observer sees is the triple (local algebra, restricted state, modular flow); the bath is the commutant" (`GRUT_III/README.md:13-16`). The sector therefore asked: once "bath = rest of the universe" is realised as P1, do P0–P5 force any state, kernel, coupling, constant, arrow or outcome rule that GR + QFT in curved spacetime does not already contain? Sub-questions for this ledger: is P1 an identification or a law; is modular KMS a temperature for a physical clock; what is imported in GR recovery and universal coupling; what do the GRUT III predictions depend on.

The handoff's closure framing for this sector: "GRUT III's P1 identifies the bath with a commutant, equating it to the causal-complement algebra only where Haag duality holds. Local type-III algebras do not have the elementary reduced-density-matrix construction." (`HANDOFF:92-94`).

## 2. Assumptions / model class

- Tiers: POSTULATE (counted), IMPORTED (zero credit), DERIVED (computed in `calc/`), DATUM, CONJECTURE (`01:3-6`).
- Data: ħ, c, G, k_B, Λ, the particle spectrum; "this is also where the gravitational coupling enters the framework (through P2 and P3)"; "GRUT III does not derive G and says so" (`01:10-12`).
- **P0** (POSTULATE, +1): local relativistic QFT on a Lorentzian manifold, a net of von Neumann algebras A(O) with Einstein causality, covariance, spectrum condition, Hadamard states; "This single postulate bundles a great deal: Lorentzian signature, causal locality, the Hilbert-space structure, and the field content" (`01:14-21`).
- **P1** (IDENTIFICATION, 0): medium = (A(O), ω|_{A(O)}, σ^ω_t); bath = commutant A(O)′ (= A(O′) where Haag duality holds) (`01:23-27`). Local algebras are type III₁: no tensor factorisation, no reduced density matrix, "P1 supplies no influence functional and no dynamics" (`01:29-31`).
- **P2** (POSTULATE, +1; G enters): S = ηA and Clausius δQ = TδS across every local causal horizon; area law IMPORTED (Bombelli–Koul–Lee–Sorkin 1986; Srednicki 1993) with UV-divergent coefficient; "η = 1/(4Għ) is inserted"; "The Clausius relation is the one structural postulate" (`01:45-52`).
- **P3** (POSTULATE, +1; κ enters): vertex (κ/2)h_μν T^μν with κ² = 32πG; "The coupling itself is inserted; its universality is the operational equivalence principle" (`01:54-59`).
- **P4** (POSTULATE, +1): the traced-out horizon is the future horizon, β > 0; "the record's Past Hypothesis, located; it is not derived" (`01:61-64`).
- **P5** (POSTULATE, +1): non-contextuality; Gleason/Busch give Tr(ρE), IMPORTED; the HJW route is withdrawn (`01:66-74`).
- Accounting: framework postulates P0, P2, P3, P4, P5 = 5; P1 an identification; "No kernel is supplied anywhere; nothing is earned by that fact alone, because the kernels that replace it are standard physics" (`01:78-88`).
- Model class of every computed result: free conformal or massless minimally coupled scalar, Bunch–Davies / conformal vacuum, geodesic/comoving worldline, exact dS or power-law FRW with constant ε, two-point (Gaussian) level, fixed-mean-time (Wigner) spectra at ω ≠ 0 (`02:259-266`; `07:126-128`; `11:126-128` adds: conformally coupled massless scalar in flat FRW, "composite, higher-spin and smeared observables not computed"). The D11 second-order formula is complete only under the full adiabatic hierarchy, ε₃ = η̇/Hη also small (`11:125-126`); the massless velocity correlator is IR-finite only for ε < 2/3 (`RES:668-669`; `11:158-159`).

## 3. Executed evidence, with sources and exact numbers

### 3.1 P1 is an identification, not a law, and selects no state

- Labelled "IDENTIFICATION, 0; stated algebraically" (`01:23`); "P1 replaces the record's `rung1_ontology_finite_memory` and `background_time_translation_flow` by naming the same objects; it earns nothing by itself" (`01:42-43`). Ledger migration: the system/bath split and the background flow are IDENTIFICATION, new Δ = 0 (`04:8`, `04:26`).
- The author's retraction: "naming the bath as the commutant algebra supplies no dynamics, and for type III algebras supplies no factorisation either" (`00:19-21`); "the frontier question is dissolved" → "relocated" (`03:62`; `RTR:141-144`).
- Current closure-fork status (19 October 2026): "C2 identifies the appropriate observer-state object, but P1 does not currently add a state-selection condition. There is no stronger, hidden closure condition in the `0b6df6f` derivations" (`CSR:7-11`). The restriction map ω_O(A) = ω(A) "exists for every state on the encompassing algebra … Thus P1 alone cannot fix even the local one-point function. This is an elementary limitation of the definition, not a new quantum-cosmology calculation" (`CSR:28-40`). An extension of C2 "must name a condition on the local state or gluing prescription that is not already restriction, trace matching, positivity, normalization, or modular KMS. No such extension has been derived here" (`CSR:84-88`).

### 3.2 Modular KMS: automatic in the modular parameter; a temperature only where the flow is geometric

- "by Tomita–Takesaki the restricted state is KMS at β = 1 in the modular parameter automatically, so 'the medium is thermal' is a tautology; the physical content is whether the modular flow is geometric, which is exactly when a temperature is a number" (`01:30-33`). The cases, each IMPORTED: Minkowski wedge β = 2π/a (Bisognano–Wichmann); bifurcate Killing horizon β = 2π/κ (Kay–Wald 1991); dS geodesic observer T = H/2π (Gibbons–Hawking 1977); CFT double cones (Hislop–Longo); "generic regions: non-geometric, no temperature" (`01:35-39`).
- Closure fork: "for a faithful normal state the KMS property with respect to its own modular flow is automatic. It does not select the state or set a temperature in physical time. Equating that flow to a separately specified geometric observer evolution would add content; P1 has no such general condition for an inflating subregion" (`CSR:42-47`).
- Owner's charter DEF-16 (located in the owner's `vsctestinggrut` repository, quoted via `08`): KMS "relative to a flow constructed from the state itself … KMS-type credit or RECOVERY needs a time evolution certified independently" (`08:71-73`); bearing: "modular-flow thermality (P1) is self-referential and earns nothing" (`08:78-80`).

### 3.3 Exact de Sitter, conformal field (D1, D2; `v1`, `v8`)

- Worldline Wightman function g(τ) = −(H²/16π²)/sinh²(H(τ−iε)/2); G⁺(ω) = (ω/2π)/(1 − e^{−2πω/H}) (`02:9-16`). Quadrature vs Planck formula: seven PASS lines with rel.err 0.0, 0.0, 1.55e-31, 0.0, 1.17e-17, 1.07e-31, 0.0 at w = 0.1, 0.5, 1, 2, 5, −0.5, −2 (`RES:32-38`); worst 1.17e-17.
- Record correction: the record's `calc/worldline_reduction.py` "symmetrised" S(w) = 2Re∫₀^∞ equals G⁺(w); its value 0.03393 at w = 0.1H is G⁺(0.1H) = 0.034116 "up to a 0.55% ε-extrapolation error" (`02:18-23`; `RES:40-43` prints 0.034115944 and 0.545%; `RTR:112-117`).
- Stress-tensor spectrum G⁺_T(ω) = ω(ω²+H²)(ω²+4H²)(ω²+9H²)/(10080π³(1 − e^{−2πω/H})) (`02:37`): KMS ratio G⁺_T(1)/G⁺_T(−1) = 535.491655525 (rel.err 1.89e-31); G⁺_T(0) = H⁷/(560π⁴) = 1.83321111691e-5; asymptotic decay rate 4.00015702692 → 4H (`RES:58,62-63`). ρ_T is an odd polynomial; 4H "is not an autocorrelation time — the pointwise correlator diverges as t⁻⁸ and has no normalisable one" (`02:42-48`).
- Priority: the position-space kernel is Cho–Hu (2014) eqs. (29)+(30); the worldline sinh⁻⁸ form "is a two-line corollary and is not new"; G⁺_T(ω), the odd polynomial and the 4H rate "are not printed in those papers and are recorded as corollaries/placements" (`02:52-58`; `08:28-34`).

### 3.4 FDT identity and the residual instrument (`07` T1–T2; `v8`)

- "If ω is KMS at β with respect to the chosen time, then … N_f(ω) = ½coth(βω/2)ρ_f(ω) for every f. The noise is never an independent datum once the state is thermal" (`07:20-23`); trichotomy T2 (`07:34-38`). Consequence: "In the regime its own principle defines, GRUT III can only be standard" (`03:8`).
- `v8`: KMS defect of the Planckian spectrum 9.86e-32; mutation at T/2 flagged with defect 0.00344529; KMS defect of G⁺_T 2.41e-35 (`RES:401-402,413`). Template on the conformal medium: VERDICT EMPTY (`RES:473-476`). The record's rung-3 kernel: IR exponent 3 vs the medium's 1, VERDICT SUPPLIED; S5-1 Lorentzian: IR exponent 0, max |Δρ|/ρ_ref = 255, VERDICT SUPPLIED (`RES:430-431,444-445,477-484`). "C1's omega_c is set to H here by fiat" (`RES:470`). Scope: "two-point level, pointwise T00 of one conformally coupled scalar in Bunch-Davies" (`RES:417-418`).

### 3.5 Slow-roll conformal sector (D11; `v11` rev. 3)

- Theorem: for any a(t) with ε = −Ḣ/H² and η = ε̇/Hε both first-order small, W(t,t′) = W_dS(Δ; H_eff) + O(ε²), H_eff = H̄(1 − ε̄) = (d/dt)ln ȧ, "exactly KMS to first order at T_eff = (ä/ȧ)/2π (a statement at fixed x = HΔ, i.e. for ω ≫ εH; … about the fixed-mean-time Wigner spectrum" (`02:160-170`). Second order: W/W_dS(H_eff) = 1 + ε²[x²/4 − (x³/12)coth(x/2)]; "KMS is first broken at O(ε²)" (`02:172-174`).
- Numbers: D/ε → 2πω/H̄ with rel 3.76e-10, 8.2e-10, 9.1e-10, 9.38e-10 at ω/H = 0.5, 1, 2, 3 (`RES:611-614`); D₂ = 2.501629, 10.99679, 31.47289, 52.14295; D₂/(2πω/H) 0.796293 → 2.76627 (`RES:616-619`). Stale value: `02:177` still prints "c = 1 to 3e-5"; the Richardson-weight fix gives 1e-9 (`11:121-123`).
- The revision-2 headline (coefficient 1.01, "1%-of-ε frequency dependence") was "O(ε²) contamination … They were right" (`02:155-158`; `RTR:185-199`). Independent check (branch `codex/grut-causal-state-review`, commit `0aca9229`): exact power-law correlator vs W_dS(H(1−ε)) differs by 2.77e-4, 6.93e-5, 1.74e-5, 4.34e-6 at ε = 0.04, 0.02, 0.01, 0.005 (`08:149-154`).

### 3.6 The present universe (`v11` D)

- ε₀ = 0.4725, η₀ = −2.055, H_eff(t₀) = 0.5275 H₀; W_dS(H_eff)/W_exact − 1 = 0.0012174, 0.005563, 0.036387, 0.1916 at ΔH₀ = 0.25, 0.5, 1.0, 1.5; "the T = 0 (Minkowski) form is closer at every grid point" (`RES:623-634`). A Fourier spectrum at t₀ "requires a window and is not computed" (`02:186`). Auditor: "in the midpoint convention the exact kernel is KMS at no β — it has imaginary period 4π/(√Ω_Λ H₀) with a pole inside the strip" (`02:195-196`; `11:139-141`). "factors 4–18" dropped (convention-dependent by 1/(1−ε₀)² = 3.6) (`02:192-193`).

### 3.7 Two temperatures (`09`, `v12`; corrected in `13`)

- T_KH = (H/2π)(1 − ε/2), T_det = (H/2π)(1 − ε); (L_K h)_tr = −a²rḢ, "K is Killing iff Hdot = 0 iff eps = 0" (`RES:645-652`). For ε ≠ 0 "P1 assigns no temperature to the region inside R_A — an application of P1's definition, not a computation" (`09:23-25`; `11:156-157`). Horizon first laws (Cai–Kim; Akbar–Cai with the signed κ) are identities of the Einstein equations (`09:39-42`; `RES:656-663`). A third candidate, T₀ = 1/(2πR_A) = H/2π (Cai–Kim; Cai–Cao–Hu; Tian–Booth), is "likewise a geometric label"; T_det/T₀ − 1 = −ε exactly (`09:30-33`; `RES:676-678`).
- `13` correction: Kurpicz–Pinamonti–Verch 2021 prove Hadamard states become KMS w.r.t. the Kodama flow at the Hayward surface gravity in a near-horizon scaling limit for spherically symmetric outer trapping horizons; "Whether this extends to the cosmological apparent horizon is open"; T_KH and T_det are "two temperatures of two different observables, both conventional" (`13:163-167`; inline in `09:26-30`). "H₀/2π ∼ 10⁻³⁰ K, is unobservable regardless" (`09:54`).

### 3.8 Massless minimally coupled sector (D13; `v14`)

- Representation via the `v12` lemma ν = 3/2 + ε/(1−ε); first-order theorem W = (1 + 3ε/2)W_mmc^{dS}(Δ; ä/ȧ) + ε(H̄²/8π²)S₂ + c₀ + c₂Δ² + O(ε²) (`02:223-235`). Residual halving ratios → 4 only at H_eff with S₂ (worst |ratio − 4| = 0.0423); ≈ 2 with S₂ dropped, at H̄, or at Kodama–Hayward (`RES:758-782`). Exact-form KMS defect ratios 4.008, 4.004 at β_eff; 2.014, 2.007 at β_KH; 2.007, 2.004 at β̄ (Δ = 0.8) (`RES:788-796`).
- Shape change: r(0) = 2 − 2 ln 2 = 0.6137, r(H̄) = 0.46077899, r(4H̄) = 0.085242243, r(8H̄) = 0.022893644 (`RES:807-811`); D = (H̄³/4π²)[1 + c_D ε], c_D = −(1 + 2 ln 2) = −2.386, matching the stochastic c₁(σ) only at σ* = e^{1−γ_E} = 1.52621 (`RES:816-818`).
- Convention: statements are about the fixed-mean-time (Wigner) spectrum at ω ≠ 0; an instantaneous one-sided rate differs at O(ε) (auditor's conformal-sector numbers c(E) = 0.24, 0.93, 3.97 at E/H = 0.25, 0.5, 1, "not re-derived by the author"); the massless linearly coupled detector is non-thermal already at ε = 0 (`02:213-221`; `11:93-101`). Scope: general slow-roll a(t) "not established"; "no tensor computation was done"; "All of this is conventional QFT; it establishes nothing about P1 beyond consistency" (`02:259-266`; `RES:826`).

### 3.9 Zero-mode temperature route (`v13` rev. 2)

- c₁(σ) = 1 − 2γ_E − 2 ln(2σ); c₁(1) = −1.5407257; c_X(σ) = 0.4685 + ½ ln σ: = 0 at σ = 0.3918, 1/2 at 1.065, 1 at 2.895 (`RES:700-723`). "The sharp-window coefficient is scheme-dependent and cannot test horizon temperatures; rev. 1's 'no temperature matches' is withdrawn as an inference (the arithmetic at sigma = 1 stands)" (`RES:735-737`). Audit: a σ-free invariant exists (D/H³ closed form, O(ε) coefficient −(1 + 2 ln 2)), so "not a physical invariant" was too strong; its closed form "is the auditor's and is not re-derived in the repository" (`11:163-176`); Cho–Kim–Kitamoto's late-time solution has no stationary weight at O(ε) and their drift is a leading-IR approximation (`11:169-171`).

### 3.10 Prior source located by the closure fork (Obadia 2008)

- Obadia, PRD 78, 083532 (2008), eqs. (20)–(22): T_coefficient = √(H² + 2Ḣ)/2π = (H/2π)[1 − ε + O(ε²)]; "The coefficient alone is not a thermal-state theorem … D2 is withdrawn as a claimed new result and retained only as an appropriately scoped consistency check. No novelty certificate is issued" (`CSR:90-111`; `REPORT.md:71-76`). The label "D2" here is an owner-brief label whose definition is not in the repository (the D1–D4 brief is missing, `SOURCES.md:50-54`); that it denotes the first-order coefficient claim is inferred from `REPORT.md:71-76` and `CSR:97-110`.

### 3.11 GR recovery (D3, D12)

- D3 IMPORTED (Jacobson 1995, reproduced in `v2`): Clausius with T = κ/2π and S = ηA gives T_kk = (η/2π)R_kk, "κ cancels, and nothing refers to any observer's physical acceleration"; null-cone lemma + Bianchi give G_ab + Λg_ab = 8πG T_ab "with η = 1/4G inserted, Λ an integration constant" (`02:60-69`). SVD nullspace singular value 6.192e-15, spanned by η_ab (`RES:80-88`). Jacobson full text: "no condition on any acceleration appears"; "D3 is entirely contained in Jacobson 1995" (`08:129-138`). P2's state content is the Hadamard condition, with exact KMS only at bifurcate Killing horizons (Kay–Wald) and otherwise in a scaling limit (`13:156-162`; `13:333`).
- D12 (where G comes from; `v9`): per species δ(1/4G): minimal scalar 1/(48πε²), Dirac 1/(24πε²), conformal scalar 0, Maxwell −1/(12πε²); SM sums Σk₁ = +1/6 (ξ = 0), −1/2 (ξ = 1/6); SM + 3ν_R at ξ = 1/6: exactly 0; Λ_UV = √(12π) M_P = 6.14 M_P "(above M_P: one-loop dominance fails)" (`02:199-209`; `RES:544-558`). "P2 is consistent, not derived; the finite G needs the UV completion" (`02:208-209`; `RES:567`; `05:7-17`).

### 3.12 Universal coupling (P3, D5; `v4`)

- From the inserted vertex: Im Π(s) = c_eff(s)Gs², w = 1 : 3 : 12, c_SM = 283/120 = 2.358 (`02:80-89`); Dirac/scalar ratio 6.0000000000, photon-pair ratio 12.0000000000, fermion law β³(1 + 8r/3) to 2.22e-16 (`RES:129-140`). HLZ 1999 eqs. (46), (52), (53): "D5's threshold laws and weights are contained in HLZ; `v4` reproduces them … (a verification, not a discovery)"; the massive-vector polynomial "not computed" (`08:119-127`; `02:93-94`). "the variable is the timelike invariant mass √s, not a GW frequency … the Standard-Model staircase in `v4` is illustrative" (`02:94-97`; `RES:142-143,160`). FRW remainder Γ = 5.2e-162 s⁻¹ at ω = 2π × 100 Hz (`RES:161-162`).

### 3.13 Horizon transport, access geometry, gravitational sector (D4, D9, D10)

- D4 IMPORTED: η_s = c³/(16πG) = 8.03 × 10³³ kg s⁻¹ (revision 1 mislabelled Pa·s); τ_β = 1/H₀ = 4.58 × 10¹⁷ s; 4GM/c³ = 19.7 µs; T_dS(H₀) = 2.655e-30 K; "neither number is a 'memory time'" (`02:71-78`; `RES:100-112`).
- D9 (`v7`): χ_L(t) = θ(t)·t·e^{−t²/4L²}/(8π^{3/2}L³), ⟨t⟩ = √π L (`02:117-121`); for the stress tensor the two closed forms agree symbolically (`RES:250`) and quadrature vs closed form agrees to 1.14e-26 at t = 0.05L (`RES:258`); smeared KMS: FDT residual ≈ 2.3e-13, direct KMS-ratio residual up to 5.14e-11 (`RES:317-385`); "apparent memory of order L/c can always be manufactured by the access profile and is never evidence for a vacuum relaxation law" (`02:135-136`).
- D10 (`v10`, IMPORTED form): stationary part 2πi/H-periodic ("full KMS analyticity in the strip is not verified here"); Brownian covariance D = H³/4π² "an asymptotic statement, not an exact covariance at all separations"; d⟨h_ij h^ij⟩/dt = 2H³/(π²M̄_P²) = H·P_T (`02:138-151`; `RES:573-588`).

### 3.14 Predictions, falsifiers, what they depend on

- "Therefore GRUT III makes no prediction that distinguishes it from GR + QFT in curved spacetime" (`03:9`). F1–F6 "are shared with ΛCDM / GR + QFT. They refute the minimal framework; they do not single GRUT III out" (`03:13-24`): F1 w ≡ −1 with DESI DR2 Table VI 3.1σ (DESI+CMB), 2.8σ / 3.8σ / 4.2σ with Pantheon+ / Union3 / DESY5, "no combination exceeds 4.2σ", w₀ = −0.752 ± 0.057, w_a = −0.86 (DESI+CMB+DESY5), "No DR3 cosmology release exists. The live risk" (`03:19`; `08:99-104`); F2 tensor friction below ∼10⁻¹⁶⁰ s⁻¹ (`03:20`); F3 untested at precision; F4 D1–D2/D11 spectra; F5 "Not accessible"; F6 consistent (`03:21-24`).
- Falsifiers of the framework itself: Lorentz violation at the horizon-thermality level (P0/P1); failure of the area law across smooth surfaces (P2); contextual outcome frequencies (P5); a non-polynomial pointwise dissipative kernel for conformal matter in dS (D2) (`03:40-44`).
- Price of the extra-light-field escape: Cassini γ − 1 = (2.1 ± 2.3) × 10⁻⁵ ("a measurement with a 1σ error, not by itself a confidence-qualified upper bound"); Eöt-Wash λ < 38.6 µm at 95% for α = 1 (`03:26-38`).
- Dependence: P0's field content and Hadamard states; the Bunch–Davies / conformal vacuum; the imported horizon-thermality theorems; the inserted G (P2) and κ (P3); for D11/D13 the slow-roll hypothesis (ε, η first-order small), the Wigner convention and ω ≠ 0; P4 and P5 as postulates (`01:14-74`; `02:160-170,213-221,259-266`).

### 3.15 The one parameter-free state-selection candidate computed (`v15`, Sorkin–Johnston)

- UV: 16k²|β|²(a_eH)⁻² = 0.99981888, 0.99999821, 1.0 at k/(a_eH) = 100, 1000, 10000 (conformal), "power law: not Hadamard"; energy density grows like Λ² — violates P2 (`RES:842-851`). IR: R_k ∝ (k/a_eH)^0.5 (conformal), ^1.5 (massless); over the inflationary slab power × 3.4e-33..1.0e-39 and tilt +1.5 vs Planck n_s − 1 = −0.035 ± 0.004 (`RES:853-868`). Region dependence: R(η = −1) = 3.96332e-5 → 1.25331e-6 → 3.96361e-8 as a_eH goes 1 → 10 → 100 (`RES:870-873`). "A region ending today is not computed here" (`RES:877-878`; `13:213-214`).
- The DMP horizon state "comes closest to GRUT's P1" but its Hadamard proof needs Re ν < 3/2, excluding the massless minimally coupled field (`13:224-231`); "*If* P1 can be shown to select it — not yet done" (`13:283-284`); "a programme, not a result" (`13:281`).

### 3.16 CP-04 and CP-05 primary reports (located on `origin/codex/grut-causal-response-05`)

- Five attempt directories exist on that branch (attempt01–05), plus `cp04_followup/REVIEW.md` and `PROGRAM.md` (git ls-tree). The branch README head, written at CP-04, says "no rerun or CP-05" and "no automatic CP-05" (`origin/codex/grut-causal-response-05:research/README.md:1,5`), although attempt05 is present on the same branch.
- CP-04 header: "Author verdict: CONDITIONAL PUBLISHED MECHANISM RETAINED; NOT A GENERATED SPACETIME OR NEW GRUT THEORY. One selected action, no repairs. Eight author check groups pass. Independent scientific review pending" (`CP04:3`); §6: "CONDITIONAL ONLY as a published vacuum-constraint construction. NOT YET ESTABLISHED as a viable complete GRUT theory" (`CP04:156`). Not earned: "physical time/causal signature, origin of Einstein dynamics/gauge symmetry, necessity of sigma, …" (`CP04:160`).
- CP-05 header: "Verdict: WORKING LEADING SCALAR–TENSOR MODEL; NO NEW GRUT LAW OR EMPIRICAL VALIDATION. … Seven author check groups pass" (`CP05:3`); η = 2g_c²/(M² + 6g_c²), 0 ≤ η < 1/3 (`CP05:100,108`); "The bound eta<1/3 is not a theorem for all gravity or a novel GRUT law" (`CP05:125`); control M = 1, g_c = 1/2, m = 1 → Z = 5/2, μ² = 2/5, η = 1/5, slip ≈ 0.863, enhancement ≈ 14.7%, "Units and numbers are invented controls" (`CP05:141`); "Author evidence only; independent review pending. … No automatic CP-06 or Stage-3 Card 3" (`CP05:153`). GRUT III's digest truncates the CP-05 label to "…NO NEW GRUT LAW" (`08:92`) and omits CP-04's header label, giving the common pattern "CONDITIONAL ONLY" (`08:86-87`), which is CP-04's §6 wording.

### 3.17 Verification record

- Three-lens adversarial verification of `v7`–`v10` "did not run: the four verifier agents failed at launch with a usage-limit error"; "Independent of the author's agents: none yet" (`08:216-231`). Round 6: "First launch: 3 of 10 round-6 agents and 0 of 3 `v14` agents ran … relaunched … all thirteen agents have reported" (`11:9-10`). "Nothing in the arithmetic of the branch was refuted … the claims made about them had to be cut back to what was computed" (`11:180-188`). All sixteen scripts PASS (v6 INFO) on 2026-10-12 (`RES:3-22`); `12:4` still says "fourteen runnable scripts" (written at rev. 2.4; v15, v16 added in rev. 2.5, `CHANGELOG.md:14-21`).

### 3.18 The secular-infrared survey carries two graviton verdicts (`10`)

- Author's synthesis: scalars "no opening"; gravitons "a live disagreement, but not between competing resummations … a question about the definition of observables in perturbative quantum gravity on de Sitter; none of P0–P5 speaks to it" (`10:12-21`; `10:136-138`).
- The graviton reading report, appended "not seen by the synthesis agent", ends: "Bottom line for GRUT III item 3: the literature's resummations do NOT yet agree … This is a genuine, open window, not a closed one" (`10:96-104`), recording Fröb's one-loop secular effect of opposite sign to Tsamis–Woodard in constant-ε inflation (`10:101`). The author's §4 verdict re-frames this as "unresolved, and it is a quantum-gravity question, not a GRUT one" (`10:138`); `05` O5 item 3 and `07` §6 adopt the author's framing (`05:70-73,119-120`; `07:131-135`). Process caveat: "the synthesis agent received only the two scalar reports intact" (`10:7-9`). Garbrecht–Prokopec's second-order response is non-thermal, so "the record's P1 must be understood as KMS of the first-order kernel" (`10:130`).

## 4. Answer

(1) **P1 is an identification, not a law.** Zero postulate credit, no dynamics, no influence functional, no type-III factorisation (`01:23-43`; `00:19-21`), and no state-selection condition: the restriction map exists for every global state, so P1 cannot fix even a local one-point function (`CSR:7-11,28-40`).

(2) **Modular KMS is a temperature for a physical clock only where the modular flow is geometric.** For a faithful normal state KMS w.r.t. its own modular flow is automatic and "does not … set a temperature in physical time" (`CSR:42-45`); the geometric cases are each imported theorems (`01:35-39`). In FRW with ε ≠ 0 the apparent horizon is not Killing and P1's definition assigns the region no temperature (`09:21-25`; `RES:649-655`); what is computable is conventional QFT: the comoving worldline's Wigner spectrum is KMS at (ä/ȧ)/2π to O(ε) in both the conformal (`v11`) and massless (`v14`) sectors, a one-sided instantaneous rate differs at O(ε), T_KH is a separate near-horizon scaling-limit KMS temperature (KPV 2021, cosmological case open), T₀ = 1/(2πR_A) a geometric label, and the exact ΛCDM kernel today is KMS at no β (`RES:607-620,644-680,787-798`; `02:194-196,213-221`; `09:30-33`; `13:163-167`).

(3) **GR recovery imports** the area law, Jacobson's 1995 derivation and the Hadamard/Kay–Wald thermality, and **inserts** η = 1/(4Għ) via the Clausius postulate P2; Λ is an integration constant; one-loop induced gravity with SM content gives Σk₁ = +1/6 or −1/2, so G is not derived in that class (`01:45-52`; `02:60-69,199-209`; `08:129-138`; `13:156-162`). **Universal coupling inserts** κ² = 32πG (P3); its derived content Im Π(s) = c_eff(s)Gs² with weights 1 : 3 : 12 is contained in Han–Lykken–Zhang 1999 (`01:54-59`; `02:80-99`; `08:119-127`).

(4) **The GRUT III predictions depend on** P0's field content and Hadamard states, the Bunch–Davies/conformal vacuum, the imported horizon-thermality theorems, the inserted G and κ, and (D11/D13) the slow-roll hypothesis in the Wigner convention at ω ≠ 0; F1–F6 are all shared with ΛCDM / GR + QFT and "GRUT III makes no prediction that distinguishes it from GR + QFT in curved spacetime" (`03:9,13-24`); the framework-level falsifiers of `03:40-44` test P0/P1, P2, P5 and D2, not a GRUT-specific claim.

## 5. Closure verdict

**Nothing additional was forced in the stated class** (covariant QFT with Hadamard states; bath = commutant / causal complement; KMS states), and the sector contains the computations that show it. Realising "bath = rest of the universe" as P1 forces only tautological KMS in the modular parameter (`01:30-33`); by the FDT identity any extra noise at the same β is extra dissipation, i.e. a change of field content (`03:5-9`; `07:20-38`). In every regime computed — exact dS (`v8`: KMS defect 2.41e-35, residual template EMPTY), quasi-dS conformal sector through O(ε²) (`v11`: D/ε → 2πω/H̄ to ≲1e-9), quasi-dS massless sector at O(ε) (`v14`: residual halving ratios → 4 only at H_eff), horizon first laws (`v12`: identities) — "the answer is fixed by conventional theory and P0–P5 add no selection" (`05:123-125`). No state is selected (restriction exists for every state; the parameter-free SJ candidate violates P2 and the CMB; the DMP horizon state is unproven at ν = 3/2), no constant is fixed (G inserted, Λ an integration constant, Sakharov fails with SM content), and P4 (arrow) and P5 (Born) remain postulates (`05:27-40`). "before there is a geometry there is no bath — the universe has no outside" (`13:61-62`). What survives is GR + QFT with a thermodynamic reading, the Δρ residual instrument, and a priced ten-line constraint sheet for any sixth postulate (`12:36-80`). This is not a theorem that no future postulate could have content: the state door (N1) and the renormalisation-constant door (N2) are pre-registered and unexecuted (`13:357-375`), and the sector's own survey file carries an unreconciled reading-report verdict that the graviton infrared is "a genuine, open window" (`10:104`) against the author's "not a GRUT one" (`10:138`).

## 6. Classification: mixed

- **Recovered known physics**: D1, D3, D4, D7, D8, D10 (form), D12 (IMPORTED, `02:7,60,71,107,111,138,199`); D2 (DERIVED, "corollary of standard results", `02:25`); D5's weights and threshold laws (contained in HLZ, "a verification, not a discovery", `08:126-127`); D9; T1–T6 ("all imported/standard", `07:40`); the `v11`/`v14`/`v12`/`v13` computations ("conventional QFT", `02:265-266`; `12:15-25`).
- **Conditional constraints**: D11 first-order thermality (ε, η first-order small; Wigner convention; ω ≫ εH; second order under the full adiabatic hierarchy; not today) (`02:160-170`; `11:125-126`); D13 (constant ε, comoving worldline, ω ≠ 0, massless scalar, ε < 2/3 for the velocity correlator; general slow roll and gravitons not established) (`02:259-266`); D5's scope (√s variable; illustrative staircase; massive-vector polynomial not computed) (`02:93-97`); F1–F6 as falsifiers of the minimal framework only (`03:13-15`); CP-05's η < 1/3 within its kinetic family (`CP05:125`).
- **Supplied model choices**: P0's content and the data (`01:10-21`); η = 1/(4Għ) (`01:49`); κ² = 32πG (`01:56`); P4 (`01:61-64`); P5 (`01:66-74`); the Bunch–Davies/Hadamard state; the fixed-mean-time convention ("a choice, but not a free parameter", `02:242-243`).
- **Rejected mechanisms**: a₀ via Jacobson; memory time 1/(4H); bulk Planckian fluid in Pa·s; c_eff at LIGO frequencies; Born via HJW; w ≡ −1 as discriminator; "+16 → +5 compression", "GR hosted", "frontier dissolved" (`03:53-63`); v11 rev. 2 coefficient 1.01 (`02:155-158`); v13 rev. 1 inference (`RES:735-737`); the record's rung-3 and S5-1 kernels (SUPPLIED, `RES:477-484`); the sharp Sorkin–Johnston state (`RES:875-878`); the closure fork's D2 novelty, sign-only curvature falsifier and all-sequestering-positive-curvature claim (`CSR:164-167`).
- **Conjecture**: P1 selecting the DMP horizon state (`13:281-286`); a non-inserted residual in the interacting gravitational infrared (`07:131-135`; `10:104`); N1, N2 (`13:359-368`).
- **Unexecuted tasks**: O7 items 1, 2, 4 (`05:143-155`); literature check for r(ω) and −(1 + 2 ln 2) and the first-principles derivation of the `v13` 13.5 invariant (`12:6-8`; `13:374-375`; `11:176`); general slow-roll and tensor extensions of D13 (`02:259-264`); O5 item 2 (`05:116-118`); SJ region ending today (`13:213-214`); N1–N4 (`13:359-373`); independent review of CP-04/CP-05 (`CP04:3`; `CP05:153`); closure fork B Step 1 (`CSR:172-174`); reading Brandenberger–Kahn 1982, Abbott–Wise 1984, Lucchin–Matarrese 1985 (`10:110-113,167`).

## 7. Remaining freedom

- The global state: no covariant state assignment (Fewster–Verch 2012 Thm 6.13); every FRW rule buys uniqueness with an extra input; the unresolved freedom is the infrared of massless fields in an accelerating universe (`13:168-171,184-192,238-262`).
- The renormalisation constants Λ, G and curvature-squared coefficients (Hollands–Wald); no horizon principle fixing one has been proposed (`13:179-182`).
- The finite value of η = 1/4G: a bare term, BSM content or a UV completion (`05:12-17`).
- P4 and P5: "No computation proposed" (`05:35,40`).
- Field content: an added light field, priced by Cassini and Eöt-Wash (`03:26-38`).
- The interacting gravitational infrared in (quasi-)dS: ΔW(t,t′) beyond the Brownian term is "the one opening" (`07:131-135`); the sector's own reading report calls it "a genuine, open window" (`10:104`) while the author calls the dispute one about observables (`10:138`; `05:119-120`); Fröb's opposite-sign one-loop secular effect in constant-ε inflation is on record (`10:101`).
- The convention at O(ε) (Wigner vs one-sided detector rate); the D11 second-order hierarchy; extension of D13 to general a(t), gravitons and Riemann-coupled detectors (`02:213-221,259-264`; `11:125-128`).
- Whether P1's geometric list should be extended by KPV 2021's scaling-limit KMS at trapping horizons; the sources leave `09:25` ("P1 assigns no temperature … by definition") and `13:163-167` side by side without reconciling them.
- Whether Haag duality holds (bath = A(O′) requires it); the regulator/split needed for any literal subregion density matrix (`01:27`; `CSR:73-82`).
- A sixth postulate (Stage-3 Card 1) only through the ten constraints of `12` §2; "I do not have a candidate that passes all ten lines" (`12:36-80,103-105`).

## 8. Failures retained

1. Revision-1 errors retracted: surface viscosity mislabelled bulk with wrong units; "memory time" attached to a correlator with none; Jacobson "fails below cH₀"; "+16 → +5" read as derivation; "GR is hosted"; "frontier question dissolved" (`00:9-23`; `03:53-66`; `RTR:10-22`).
2. P5 from P0 + P1 via Hughston–Jozsa–Wootters withdrawn (`01:69-70`; `RTR:59-66`).
3. "w ≡ −1 is a sharp GRUT prediction; DESI DR3 decides" retracted; no DR3 release exists (`03:57`; `08:104`).
4. v11 revision-2 headline (coefficient 1.01, 1%-of-ε frequency dependence) refuted by two independent reviewers as O(ε²) contamination (`02:155-158`; `RTR:183-199`); `02:177` still carries the stale "3e-5" precision (`11:121-123`).
5. v11 D "factors 4–18" dropped (convention-dependent by 1/(1−ε₀)² = 3.6; with the εη term "4–22") (`02:192-193`; `11:61`); "deviates from any thermal form by up to 19%" corrected to a one-form time-domain comparison (`11:130-143`); unit-label misprint "T_eff = 0.084 H₀/2π" fixed (`11:139`).
6. v13 rev. 1 inference ("no horizon temperature matches") withdrawn; "not a physical invariant" narrowed to "the sharp-window coefficient is scheme-dependent" after the audit produced a σ-free invariant, whose closed form is the auditor's and "is not re-derived in the repository" (`RTR:254-261`; `RES:735-739`; `11:163-176,184-185`).
7. `12` §3 "no state-selection rule has been proposed" — wrong; several exist (`13:44-46`). `12` constraint 1 "dS KMS by invariance and positivity" corrected: needs Hadamard, α-vacua not KMS (`13:152-155`).
8. `09` "T_KH is a geometric label" corrected by KPV 2021 (`13:163-167`); `09`/`v12` misattributed T_KH to Cai–Kim and Cai–Cao–Hu and misdescribed Tian–Booth (`09:30-31`; `08:165-166,169-173`; `11:151-154`).
9. `13` §8 adversarial review: "GRUT derives the initial state" asserted a derivation nobody has done; the pre-registered N1 originally asked a question already answered by the IR divergence for ν ≥ 3/2; the horizon-state theorem's ν was conflated with the slow-roll mode index; "would describe a universe that is not ours" overstated; the interference argument stated as established (`13:344-353`).
10. Revision-1 priority claim that the conformal dS kernel "is in PNRV 2009" inaccurate (PNRV ξ = 0 only); SCH 2016 listed as a source for the dS kernel inaccurately; `08` misquoted Cho–Kim–Kitamoto ("rescaled time variable") (`08:15-18,44-46,185-186`; `11:171`).
11. `07` §5 misquote ("against 7", "−4") and T5 sign misprint corrected; `v8`'s printed conclusion and docstring overstated (`07:61-62,121-122`; `11:26-37`); minor v8 items left open (`11:42-45`).
12. `06`'s first revision-2 text said "undetermined coefficient"; corrected from Milgrom's paper (`06:18`).
13. Three-lens verification of v7–v10 did not run (usage limit); round-6 first launch 3 of 10 and 0 of 3 agents; "Independent of the author's agents: none yet" (`08:216-231`; `11:9-11`). `10`'s synthesis agent received only the two scalar reports intact (`10:7-9`); its reading report and author verdict disagree on the graviton window (`10:104` vs `10:138`).
14. Sorkin–Johnston candidate fails P2 (non-Hadamard) and the CMB in amplitude and tilt over the inflationary slab; region-dependent; region ending today not computed (`RES:875-878`).
15. Candidate conventional results (v14 r(ω), −(1 + 2 ln 2)) "not located in print; literature check owed" (`12:21-22`); Brandenberger–Kahn 1982, Abbott–Wise 1984, Lucchin–Matarrese 1985 not read (`10:110-113,167`).
16. Closure fork: "D2 is withdrawn as a claimed new result … No novelty certificate is issued" (Obadia 2008 prior); C2 adds an observer interpretation, not the missing selection condition; Step 1 unexecuted; sign-only curvature falsifier and the all-sequestering-forces-positive-curvature claim withdrawn (`CSR:90-111,164-177`).
17. CP-04 and CP-05 are author verdicts only, independent review pending; CP-05 "No automatic CP-06"; the -05 branch README head still says "no automatic CP-05" (`CP04:3`; `CP05:153`; `origin/codex/grut-causal-response-05:research/README.md:1,5`).

## 9. Must not promote

- F1 (w ≡ −1) into a GRUT discriminator; DESI DR2 significances (3.1σ; 2.8–4.2σ) into a current verdict either way (`03:15,19`).
- D11/D13 first-order thermality into a "temperature of the universe": valid only for ε, η first-order small, Wigner convention, ω ≠ 0; today ε₀ = 0.47, η₀ = −2.06 and "The exact kernel is needed"; H₀/2π ∼ 10⁻³⁰ K unobservable regardless (`02:184-197`; `09:53-55`).
- The "five framework postulates" count into a compression or achievement (`00:17-19`; `01:86-88`).
- v14's r(ω) and −(1 + 2 ln 2) into GRUT results — "candidate small conventional result", literature check owed (`12:21-22`).
- The Δρ residual template into a prediction — it is an instrument; its non-EMPTY verdicts on the record's kernels are "SUPPLIED" (`RES:472-489`).
- "KMS by construction ⇒ standard" into a theorem about all possible sixth postulates; it concerns the class P0–P5 with KMS states (`05:126-134`; `13:357`).
- The c_eff staircase table (illustrative; massive-vector polynomial not computed; √s ≠ GW frequency) (`RES:142-163`).
- τ_β = 4.58 × 10¹⁷ s, 19.7 µs, or 1/(4H) into "memory times" (`RES:109-112`).
- v15's 10⁻³³–10⁻³⁹ suppression and +1.5 tilt into a theorem against all state-selection rules (sharp SJ on the inflationary slab only) (`RES:865-878`).
- CP-05's η < 1/3 bound and control numbers (Z = 5/2, μ² = 2/5, η = 1/5, slip 0.863, 14.7%) into physical predictions — "invented controls" (`CP05:125,141`).
- D10's Brownian covariance into an exact covariance (asymptotic; full KMS analyticity not verified) (`02:145-148`).
- Historical ledger counts ("+16", "twenty campaigns") into current statements (`00:32`; `04:27-33`).
- The handoff's "P1 identifies the bath with a commutant" into a result; it is a definition (`HANDOFF:92-93`; `01:23`).
- The author's "not a GRUT one" verdict on the graviton infrared into a closed question, or the reading report's "genuine, open window" into a GRUT opening; both are in `10`, unreconciled (`10:104,138`).
- `02:177`'s "c = 1 to 3e-5" (stale; 1e-9 after the weight fix, `11:121-123`).

## 10. Handoff row discrepancies (`HANDOFF:118`)

1. No false statement in the row; all three clauses are supported by `01:23-59` and `CSR:7-11,42-45`. The items below are omissions, precision losses and one inter-source update.
2. "Modular KMS is not automatically a temperature for a physical clock" omits that two conventional temperatures of two different observables exist in FRW — T_det = (ä/ȧ)/2π (Wigner convention) and T_KH (near-horizon scaling-limit KMS for the Kodama flow, KPV 2021) — plus the geometric label T₀ = 1/(2πR_A), with the cosmological apparent-horizon extension open and the P1-definition/KPV tension unreconciled (`13:163-167`; `09:21-33`).
3. "GR recovery and universal coupling import physical inputs" understates the specifics: η = 1/(4Għ) inserted via the Clausius postulate; the Jacobson derivation itself IMPORTED; κ² = 32πG inserted; and the computed negative result Σk₁ = +1/6 (ξ = 0) or −1/2 (ξ = 1/6) for one-loop induced gravity with SM content (`01:45-59`; `02:199-209`; `RES:544-558`).
4. The row does not state the sector's headline closure result — "GRUT III makes no prediction that distinguishes it from GR + QFT in curved spacetime", F1–F6 shared with ΛCDM (`03:9,13-24`) — nor the framework-level falsifiers (`03:40-44`).
5. "Not a new state-selection law" should carry `13`'s caveats: rules do exist (correcting `12` §3); the only parameter-free one (SJ) fails P2 and the CMB; the GRUT-native DMP route is "a programme, not a result", unproven at ν = 3/2 (`13:44-46,194-236,281-286`).
6. `HANDOFF:187-188` (locate CP-04/CP-05 primaries before attaching labels or a count of five): both primaries exist on `origin/codex/grut-causal-response-05` (`attempt04/REPORT.md`, `attempt05/REPORT.md`), with five attempt directories on that branch; the program's CP numbering extends to CP-12 elsewhere (`13:5`; `14:1-3`), so "five" is a branch count. `08:84-95` truncates CP-05's label (primary: "…NO NEW GRUT LAW OR EMPIRICAL VALIDATION", `CP05:3`) and gives for CP-04 only the common pattern "CONDITIONAL ONLY" — which is CP-04's §6 wording (`CP04:156`) but omits its header label "CONDITIONAL PUBLISHED MECHANISM RETAINED; NOT A GENERATED SPACETIME OR NEW GRUT THEORY" (`CP04:3`). `SOURCES.md:50-51` ("not found in this checkout") is accurate for the working tree; the implementation is on the fetched remote. The -05 branch README head (written at CP-04) says "no automatic CP-05" (`research/README.md:1,5`) although attempt05 was later added on the same branch.
7. The closure fork's location of Obadia 2008 as prior source for the (H/2π)(1 − ε) coefficient (`CSR:92-111`) partially discharges `08:174-175` "literature check owed", `12:20` "not located in print" and `10:128` item (e) (coefficient only, not the thermality theorem or the massless-sector statements); the row does not record this inter-source update.
8. The row's source list omits `10_SECULAR_IR_SURVEY.md`, which holds the unreconciled graviton verdicts (`10:104` vs `10:138`) relevant to "what remains free".

## 11. Unread or missing

- CP-01, CP-02, CP-03 primary reports not read (only `research/README.md:1-40` and `PROGRAM.md:1-40` on the -05 branch); `cp04_followup/REVIEW.md` read at lines 1-40 only. No independent-review file for CP-04/CP-05 found in the attempt04/attempt05 directory listings (`05:58-59` says the CP-05 review was "relayed by the owner"); `cp04_followup/REVIEW.md:3` says a reported Claude check "is retained as user-provided review evidence, not authenticated independent scientific clearance".
- `GRUT_III/14_CP12_GRAVITY_ADJUDICATION.md` read at lines 1-30 only (CP-12 sector). `10_SECULAR_IR_SURVEY.md` read at 1-168 in full. `04`, `06`, `11`, `12`, `13`, `00`–`03`, `05`, `07`–`09`, `README.md`, `CHANGELOG.md`, `RTR`, `RES` read in full.
- `GRUT_III/calc/*.py` (19 files listed) not executed or read; `RESULTS.md` (generated 2026-10-12) taken as the run record.
- The owner's `vsctestinggrut` `STAGE3_CHARTER.md` (DEF-16, NR-12(b), CT-75, SCOPE-Q) is not in this repository; quoted only via `08:61-77`. The "Scope-Q / G2-12" ruling is "not found in this repository; the convergence is recorded on the owner's word" (`RTR:176-177`; `05:48-49`).
- The D1–D4 closure literature brief is missing (`SOURCES.md:50-54`); the "D2" withdrawn in `CSR:110` is an owner-brief label whose definition is not in the repository.
- Primary texts of Obadia 2008, Kurpicz–Pinamonti–Verch 2021, Fewster–Verch 2012, Jacobson 1995, Han–Lykken–Zhang 1999, DESI DR2, Garbrecht–Prokopec 2004, Cho–Kim–Kitamoto 2015 not read here; cited as GRUT III and the closure fork report them. Brandenberger–Kahn 1982, Abbott–Wise 1984, Lucchin–Matarrese 1985 were not read by the sector either (`10:110-113,167`).
- The published independent check `research/reality_framework/grut_iii_independent_review/REVIEW_V11_R2.md` on `codex/grut-causal-state-review` (commit `0aca9229`) was not opened (`08:149-154`).
- The original `GRUT_V7.md` is not at a repository path (`HANDOFF:174-177`); not needed for this sector.

## 12. Sources

- `/home/user/TestingGRUT/GRUT_III/00_LETTER.md:1-74`
- `/home/user/TestingGRUT/GRUT_III/01_PRINCIPLES.md:1-88`
- `/home/user/TestingGRUT/GRUT_III/02_DERIVATIONS.md:1-266`
- `/home/user/TestingGRUT/GRUT_III/03_PREDICTIONS_AND_FALSIFIERS.md:1-66`
- `/home/user/TestingGRUT/GRUT_III/04_LEDGER_MIGRATION.md:1-41`
- `/home/user/TestingGRUT/GRUT_III/05_ROAD_TO_COMPLETION.md:1-155`
- `/home/user/TestingGRUT/GRUT_III/06_CONJECTURE_A0.md:1-39`
- `/home/user/TestingGRUT/GRUT_III/07_ACCESS_RESOLVED_RESPONSE.md:1-135`
- `/home/user/TestingGRUT/GRUT_III/08_PRIORITY_AND_SOURCES.md:1-231`
- `/home/user/TestingGRUT/GRUT_III/09_TWO_TEMPERATURES.md:1-91`
- `/home/user/TestingGRUT/GRUT_III/10_SECULAR_IR_SURVEY.md:1-168`
- `/home/user/TestingGRUT/GRUT_III/11_VERIFICATION_ROUND.md:1-188`
- `/home/user/TestingGRUT/GRUT_III/12_BASELINE_AND_HANDOFF.md:1-113`
- `/home/user/TestingGRUT/GRUT_III/13_SYNTHESIS.md:1-375`
- `/home/user/TestingGRUT/GRUT_III/14_CP12_GRAVITY_ADJUDICATION.md:1-30`
- `/home/user/TestingGRUT/GRUT_III/README.md:1-70`
- `/home/user/TestingGRUT/GRUT_III/CHANGELOG.md:1-119`
- `/home/user/TestingGRUT/GRUT_III/RESPONSE_TO_REVIEWS.md:1-287`
- `/home/user/TestingGRUT/GRUT_III/RESULTS.md:1-932`
- `/home/user/TestingGRUT/research/reality_framework/grut_closure_fork_20261019/CLOSURE_SELECTION_RESPONSE.md:1-177`
- `/home/user/TestingGRUT/research/reality_framework/grut_closure_fork_20261019/REPORT.md:1-181`
- `/home/user/TestingGRUT/research/reality_framework/grut_closure_fork_20261019/SOURCES.md:1-54`
- `/home/user/TestingGRUT/research/reality_framework/grut_closure_fork_20261019/NEXT_QUESTION.md:1-51`
- `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md:1-251`
- `origin/codex/grut-causal-response-05:research/constructive_physics/attempt04/REPORT.md:1-168`
- `origin/codex/grut-causal-response-05:research/constructive_physics/attempt05/REPORT.md:1-153`
- `origin/codex/grut-causal-response-05:research/constructive_physics/cp04_followup/REVIEW.md:1-40`
- `origin/codex/grut-causal-response-05:research/constructive_physics/PROGRAM.md:1-40`
- `origin/codex/grut-causal-response-05:research/README.md:1-40`
- `git branch -r`, `git ls-tree -r --name-only origin/codex/grut-causal-response-05 -- research/constructive_physics`, `git diff --stat 96a7e7e HEAD -- GRUT_III research/reality_framework/...` (empty), `ls GRUT_III/calc/` at `a6be0e4`
